import os
import csv
from collections import namedtuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

from zigzag import preprocessing
from zigzag import merops_connector
from zigzag import sequences


# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1000000)

# General Seaborn aesthetics
sns.set(font='monospace')

########################################################################
# MAPS SINGLE-LETTER CATALYTIC TYPE CODES FROM MEROPS IDs TO DESCRIPTIVE
# CATALYTIC TYPE NAMES. ONLY USED FOR GENERATION OF LABELS ON PLOTS.
########################################################################
CATALYTIC_TYPES = {
    'A':'Aspartic',
    'C':'Cysteine',
    'G':'Glutamic',
    'M':'Metallo',
    'N':'Asparagine',
    'P':'Mixed',
    'S':'Serine',
    'T':'Threonine',
    'U':'Unknown',
}
########################################################################

CLUSTERMAP_ROW_COLORS = [
    'e1e1e1',               # Light grey color for 0.
    '424141',               # Dark grey color for 1.
]

########################################################################
### REFACTOR!!!
########################################################################
SUBSTITUTION_MATRIX = preprocessing.importSubstitutionMatrix(
    '../data/BLOSUM62_percentages.txt'
)
########################################################################

PatternFrequency = namedtuple(
    'PatternFrequency',
    ['absolute', 'percentage']
)

PatternSimilarityMatrix = namedtuple(
    'PatternSimilarityMatrix',
    ['matrix', 'num_clustered', 'num_unclustered']
)

def generate_position_labels(sequence_df):
    """
    Generate label string for pattern postions in y-ticks of protease
        pattern heatmap.

    Parameters:
        sequence_df -- Pandas DataFrame.

    Returns:
        position_labels -- String.
    """

    position_labels = []
    for i, position_label in enumerate(sequence_df.columns.tolist()):
        position_labels.append(position_label.center(5))
        # Add vert bar for cleavage site if even-width sequences.
        if i == ((len(sequence_df.columns) / 2) - 1):
            position_labels.append('|'.center(5))
    position_labels = ''.join(position_labels)
    position_labels += '    { Abs Fq (Pct Fq) }'
    
    return position_labels


def generate_protease_labels(protease_patterns):
    """
    Generate protease labels of format:
        "Catalytic type  -  Family  -  Protease" for all proteases with
        pre-computed patterns in MEROPS database instance.

    Parameters:
        protease_patterns -- Dict.

    Returns:
        protease_labels --
    """
    
    # Get sorted list of proteases.
    proteases = sorted(list(protease_patterns.keys()))

    # Build dictionary for hierarchical protease labeling.
    protease_groups = {}
    for protease in proteases:
        family = merops_connector.retrieve_protease_family_code(protease)
        catalytic_type = CATALYTIC_TYPES[family[0]]
        if catalytic_type in protease_groups:
            if family in protease_groups[catalytic_type]:
                protease_groups[catalytic_type][family].append(protease)
            else:
                protease_groups[catalytic_type][family] = [protease]
        else:
            protease_groups[catalytic_type] = {family: [protease]}

    # Generate sorted list of protease labels.
    protease_labels = []
    # Loop through catalytic types in alphabetical order.
    for catalytic_type in sorted(protease_groups.keys()):
        # Loop through families in alphabetical order.
        for family in sorted(protease_groups[catalytic_type].keys()):
            family_proteases = sorted(protease_groups[catalytic_type][family])
            # Generate label for each protease in family and add to list.
            protease_labels += [
                (catalytic_type.ljust(
                    len(max(protease_groups.keys(), key=len)) + 2))
                + '-  ' + family
                + '  -  ' + protease
                for protease in family_proteases
            ]
        
    return protease_labels


def generate_non_exact_protease_pattern_matrix(patterns,
                                                protease_patterns,
                                                protease_labels):
    """
    Calculate pattern similarity scores matching foreground patterns
        to MEROPS protease substrate patterns. Similarity scores are
        equal to the fraction of fixed foreground pattern positions that
        are also present MEROPS protease patterns. Foreground patterns
        must be at least as specific as protease patterns. The maximum
        similarity score is taken for each protease in cases where
        proteases have more than one corresponding pattern.

    Parameters:
        patterns -- PatternContainer instance.
        protease_patterns -- Dict.
        protease_labels -- List. Protease pattern labels generated by
                            protease_labels function.

    Returns:
        non_exact_protease_pattern_matrix -- List.
    """

    non_exact_protease_pattern_matrix = [[0] * len(protease_labels)]
    for pattern in patterns.pattern_list:
        # Score foreground patterns against protease patterns.
        pattern_protease_scores = []
        for protease_label in protease_labels:
            # Find protease name in protease label string.
            protease = protease_label[
                (len(protease_label)
                - protease_label[::-1].find('  -  ')):].strip()
            # Score foreground pattern against all protease patterns.
            protease_pattern_scores = []
            for protease_pattern in protease_patterns[protease]:
                # Convert protease pattern string to matrix.
                protease_pattern_matrix = sequences.vectorize_pattern(
                    pd.Series(protease_pattern.split(',')),
                    patterns.background,
                    add_compound_residues=False
                ).reshape(pattern.pattern_matrix.shape)
                # Get constituents from compound residues in protease pattern.
                try:
                    protease_constituent_pattern = sequences.get_pattern_constituents(
                        protease_pattern_matrix,
                        pattern.background
                    )
                except AttributeError:
                    pass

                # Score foreground pattern against protease pattern.
                protease_pattern_intersection = np.sum(
                    np.logical_and(
                        pattern.pattern_matrix,
                        protease_constituent_pattern
                    )
                )

                if (
                    protease_pattern_intersection == np.sum(
                        np.any(protease_pattern_matrix, axis=1)
                    )
                ):
                    protease_pattern_score = protease_pattern_intersection / np.sum(
                        np.any(pattern.pattern_matrix, axis=1)
                    )
                else:
                    protease_pattern_score = 0

                protease_pattern_scores.append(protease_pattern_score)
            pattern_protease_scores.append(max(protease_pattern_scores))
        non_exact_protease_pattern_matrix.append(pattern_protease_scores)

    return non_exact_protease_pattern_matrix


def generate_protease_pattern_frequency_matrix(patterns,
                                                protease_patterns,
                                                protease_labels):
    """
    Frequency matrix for pattern-protease matches.
    
    Parameters:
        patterns -- PatternContainer instance.
        protease_patterns -- Dict.
        protease_labels -- List. Protease pattern labels generated by
                            protease_labels function.

    Returns:
        frequency_matrix -- List.
    """
    frequency_matrix = [[0] * len(protease_labels)]
    for pattern in patterns.pattern_list:
        # Get pattern string.
        pattern_string = ''.join(pattern.character_pattern().tolist())
        # Frequency encode protease hits
        pattern_protease_vector = []
        for protease_label in protease_labels:
            # Find protease name in protease label string.
            protease = protease_label[
                (len(protease_label)
                - protease_label[::-1].find('  -  ')):].strip()
            if pattern_string in protease_patterns[protease]:
                pattern_protease_vector.append(pattern.subset_tensor.shape[0])
            else:
                pattern_protease_vector.append(0)
        frequency_matrix.append(pattern_protease_vector)
    
    return frequency_matrix


def generate_protease_pattern_percentage_frequency_matrix(absolute_frequency_matrix,
                                                            num_sequences):
    """
    Frequency of protease pattern matches in aligned data set in
        percentage form.

    Parameters:
        absolute_frequency_matrix -- List.
        num_sequences -- Int.
    Returns:
        percentage_frequency_matrix -- List.
    """
    
    percentage_frequency_matrix = np.divide(
        np.array(absolute_frequency_matrix, dtype=np.float64),
        num_sequences
    ).tolist()

    return percentage_frequency_matrix


def calculate_pattern_similarity_matrix(sequence_df,
                                        pattern_list,
                                        pattern_labels,
                                        subs):
    """
    Calculates mean positional similarity between each input sequence
        in aligned data set and each pattern detected in the data.

    Parameters:
        sequence_df -- Pandas DataFrame. One sequence per row.
        pattern_list -- List. Contains pattern instances.
        pattern_labels -- List.
        subs -- Pandas DataFrame. Contains pairwise
                                substitution probabilities.

    Returns:
        pattern_similarity_matrix -- Pandas DataFrame.
    """

    pattern_similarity_matrix = []
    for i, sequence in sequence_df.iterrows():
        sequence_scores = []
        for j, pattern in enumerate(pattern_list):
            position_scores = []
            for position, residue in pattern.character_pattern().iteritems():
                stripped_residue = residue.strip('[]')
                if stripped_residue != '.':
                    if stripped_residue in pattern.background.alphabet:
                        try:
                            position_scores.append(
                                subs[sequence[position]][stripped_residue]
                            )
                        except KeyError:
                            position_scores.append(0.0)
                            continue
                    else:
                        try:
                            constituent_residues = pattern.background.compound_residues[
                                stripped_residue].residues
                        except KeyError:
                            print('Residue not found:  {}'.format(stripped_residue))
                            continue
                        else:
                            if sequence[position] in constituent_residues:
                                position_scores.append(1.0)
                            else:
                                # Get all pairwise substitution probabilities
                                # between sequence position and compound residue.
                                try:
                                    constituent_substitution_scores = [
                                        subs[sequence[position]][constituent]
                                        for constituent in constituent_residues
                                    ]
                                except KeyError:
                                    position_scores.append(0.0)
                                    continue
                                else:
                                    # Select highest substitution probability.
                                    position_scores.append(
                                        max(constituent_substitution_scores)
                                    )

            sequence_scores.append(np.mean(position_scores))
        pattern_similarity_matrix.append(sequence_scores)

    pattern_similarity_matrix = pd.DataFrame(
        pattern_similarity_matrix,
        columns=pattern_labels[1:]
    )

    return pattern_similarity_matrix


def generate_pattern_labels(position_labels, patterns):
    """
    Generate pattern labels containing aligned positional residues and
        and pattern frequency in absolute and percentage forms.

    Parameters:
        position_labels -- String.
        patterns -- PatternContainer instance.
    Returns:
        pattern_labels -- List.
    """

    pattern_labels = [position_labels]
    for pattern in patterns.pattern_list:
        character_pattern = pattern.character_pattern().tolist()
        pattern_label = ''
        # Add pattern positions to label.
        for i, residue in enumerate(character_pattern):
            pattern_label += residue.center(5)
            # Add cleavage site indicator to even-width sequences.
            if i == ((len(character_pattern) / 2) - 1):
                pattern_label += '|'.center(5)
        # Add pattern frequency data to label.
        pattern_frequency = calculate_pattern_frequency(pattern)
        pattern_label += (
            '    {{ {absolute:6d} ({percentage:6.2f}) }}'.format(
                absolute=pattern_frequency.absolute,
                percentage=pattern_frequency.percentage
            )
        )
        pattern_labels.append(pattern_label)

    return pattern_labels


def calculate_pattern_frequency(pattern):
    """
    Calculates frequency of pattern in input data set.

    Parameters:
        pattern -- Pattern instance.

    Returns:
        pattern_frequency -- PatternFrequency instance.
    """

    # Get constituent pattern.
    try:
        constituent_pattern = sequences.get_pattern_constituents(
            pattern.pattern_matrix,
            pattern.background
        )
    except AttributeError:
        constituent_pattern = pattern.pattern_matrix

    # Intersect pattern with sample sequence tensor.
    sequence_tensor = pattern.pattern_container.sample.sequence_tensor
    positional_matches = np.sum(
        np.any(
            np.logical_and(
                sequence_tensor,
                constituent_pattern
            ),
            axis=2
        ).astype(np.int8),
        axis=1
    )

    # Get number of matching sequences from intersection.
    absolute_frequency = positional_matches[
        positional_matches == np.sum(np.any(constituent_pattern, axis=1).astype(np.int8))
    ].shape[0]

    # Calculate percentage frequency from absolute frequency.
    percentage_frequency = absolute_frequency / sequence_tensor.shape[0]

    pattern_frequency = PatternFrequency(
        absolute=absolute_frequency,
        percentage=percentage_frequency
    )

    return pattern_frequency


def generate_row_colors(frequency_matrix, protease_labels):
    """
    Generate row color annotation dataframe for sequence pattern
        similarity clustermap. Annotations communicate proteases
        with matching pre-computed patterns.

    Parameters:
        frequency_matrix -- List.
        protease_labels -- List. Protease pattern labels generated by
                            protease_labels function.

    Returns:
        row_colors -- Pandas DataFrame.
    """

    # Convert frequency matrix to binary matrix.
    binary_protease_pattern_matrix = np.array(
        frequency_matrix[1:], dtype=np.bool).astype(np.int8)

    # Boolean array to drop proteases with no pattern hits.
    protease_hits = np.any(binary_protease_pattern_matrix, axis=0)

    # Matrix containing protease hits only.
    binary_protease_pattern_matrix = binary_protease_pattern_matrix[:, protease_hits]

    # Drop protease labels with no hits.
    protease_labels = protease_labels[protease_hits]

    # Map row color annotations values to corresponding colors.
    annotation_values = set(binary_protease_pattern_matrix.flatten())
    lut = dict(zip(annotation_values), CLUSTERMAP_ROW_COLORS)

    # Map annotation colors binary protease hit matrix.
    row_colors = []
    for i in binary_protease_pattern_matrix.shape[1]:
        protease = pd.Series(
            binary_protease_pattern_matrix[:, i],
            name=protease_labels[i]
        ).map(lut)
        row_colors.append(protease)

    row_colors = pd.DataFrame(row_colors).transpose()
    row_colors.columns = protease_labels

    return row_colors


def generate_sequence_clustermap(title,
                                    pattern_similarity_matrix,
                                    output_path,
                                    row_colors=None,
                                    row_color_ratio=None,
                                    annotate_rows=False):
    """
    Clustermap with sequences along x-axis and patterns along y-axis.
        Optional row colors and accompanying labels to annotate
        proteases with matching patterns.

    Parameters:
        title -- String.
        patterns -- PatternContainer instance.
        pattern_similarity_matrix -- Pandas DataFrame.
        position_labels -- String.
        output_path -- String.
        row_colors -- Pandas DataFrame.
        row_color_ratio -- 
        annotate_rows -- Boolean.

    Returns:
        None
    """

    # Set up figure and axes.
    #clustermap_fig = plt.figure(figsize=(24,12))
    #clustermap_ax = clustermap_fig.add_subplot(111)
    
    if len(pattern_similarity_matrix) > 1:
        cluster_columns = True
    else:
        cluster_columns = False
    if len(pattern_similarity_matrix.columns) > 1:
        cluster_rows = True
        ytick_rotation = 90
    else:
        cluster_rows = False
        ytick_rotation = 270

    # Plot clustermap.
    clustermap = sns.clustermap(
        pattern_similarity_matrix.transpose(),
        method='complete',
        metric='cosine',
        cmap='Oranges',
        row_cluster=cluster_rows,
        col_cluster=cluster_columns,
        row_colors=row_colors,
        #row_color_ratio=row_color_ratio,
        annot=False,
        xticklabels=False,
        yticklabels=True
    )

    plt.setp(clustermap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    # Compute clustermap height and width.
    dpi = 300
    # fontsize = plt.rcParams['ytick.labelsize']
    # heatmap_height = (fontsize * len(pattern_similarity_matrix.columns)) / dpi

    ## Configure heatmap aesthetics.
    clustermap.fig.suptitle(title)
    clustermap.ax_heatmap.set(xlabel='Sequences', ylabel='Patterns')
    #loc, labels = plt.xticks()
    #clustermap.ax_heatmap.set_yticklabels(labels[::-1], rotation=ytick_rotation)
    # heatmap_position = clustermap.ax_heatmap.get_position()
    # row_dendrogram_position = clustermap.ax_row_dendrogram.get_position()
    # clustermap.ax_heatmap.set_position([
    #     heatmap_position.x0,
    #     heatmap_position.y0,
    #     heatmap_position.width,
    #     heatmap_height
    # ])
    # clustermap.ax_row_dendrogram.set_position([
    #     row_dendrogram_position.x0,
    #     row_dendrogram_position.y0,
    #     row_dendrogram_position.width,
    #     row_dendrogram_position.height * (heatmap_position.height / heatmap_height)
    # ])
    clustermap.savefig(output_path, dpi=dpi)


def generate_protease_pattern_heatmap(title,
                                        patterns,
                                        scoring_matrix,
                                        protease_labels,
                                        position_labels,
                                        output_path):
    """
    Generate and save pattern-to-protease heatmap.
    
    Parameters:
        title -- String.
        patterns -- PatternContainer instance.
        scoring_matrix -- List. Contains protease pattern frequencies.
        protease_labels -- List. Protease pattern labels generated by
                            protease_labels function.
        output_path -- String. Output path for heatmap SVG file.

    Returns:
        None
    """
    
    # Set up figure and axes.
    heatmap_fig = plt.figure(figsize=(24,12))
    heatmap_ax = heatmap_fig.add_subplot(111)
        
    # Generate pattern labels.
    pattern_labels = generate_pattern_labels(position_labels, patterns)
    
    # Plot heatmap on axes.
    sns.heatmap(scoring_matrix,
                   cmap='Oranges',
                   xticklabels=protease_labels,
                   yticklabels=pattern_labels,
                   linewidths=0.2,
                   linecolor='#a5a091',
                   ax=heatmap_ax,
                   vmin=0.0,
                   vmax=(100.0 if (np.amax(scoring_matrix).astype(np.float64) == 0)
                            else np.amax(scoring_matrix).astype(np.float64)))
    
    # Dynamic label coloring.
    protease_pattern_scores = np.array(scoring_matrix,
                                            dtype=np.float64)
    max_pattern_scores = np.amax(protease_pattern_scores, axis=1)
    max_protease_scores = np.amax(protease_pattern_scores, axis=0)
    
    # Color matching pattern labels.
    for i, max_pattern_score in enumerate(max_pattern_scores):
        if max_pattern_score == 1.0:
            heatmap_ax.get_yticklabels()[i].set_color("red")
    # Color detected protease labels.
    for i, max_protease_score in enumerate(max_protease_scores):
        if max_protease_score == 1.0:
            heatmap_ax.get_xticklabels()[i].set_color("red")

    # Configure heatmap aesthetics.
    plt.setp(heatmap_ax.xaxis.get_majorticklabels(), rotation=270)
    plt.setp(heatmap_ax.yaxis.get_majorticklabels(), rotation=0)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)


def generate_figures(sequence_df,
                        patterns,
                        analysis_title,
                        output_dir,
                        annotate_clustermap=False):
    """
    Helper function used to generate all default output figures.

    Parameters:
        sequence_df -- Pandas DataFrame.
        patterns -- PatternContainer instance.
        analysis_title -- String.
        output_dir -- String. Directory for figures.
        annotate_clustermap -- Boolean.

    Returns:
        None
    """
    
    # Create output file name prefix from analysis title.
    output_prefix = analysis_title.replace(" ", "_")

    # Create output directory.
    try:
        os.makedirs(output_dir + '/figures')
    except FileExistsError:
        pass

    # Get positions from sample data frame.
    position_labels = generate_position_labels(sequence_df)

    # Get protease patterns from MEROPS database instance.
    protease_patterns = merops_connector.retrieve_protease_patterns()
    # Generate protease labels.
    protease_labels = generate_protease_labels(protease_patterns)

    # Calculate non-exact protease pattern intersection scores.
    non_exact_scoring_matrix = generate_non_exact_protease_pattern_matrix(
        patterns,
        protease_patterns,
        protease_labels
    )

    # Set absolute frequency protease pattern heatmap title.
    protease_pattern_heatmap_title = (
        analysis_title + ' - Protease Pattern Matches (Percent Positions Matched)'
    )
    
    # Set output path for absolute frequency protease pattern heat map.
    protease_pattern_heatmap_output_path = (
        output_dir
        + '/figures/'
        + output_prefix
        + '_protease_pattern_heatmap.svg'
    )

    # Generate absolute frequency protease pattern heatmap.
    generate_protease_pattern_heatmap(
        protease_pattern_heatmap_title,
        patterns,
        non_exact_scoring_matrix,
        protease_labels,
        position_labels,
        protease_pattern_heatmap_output_path
    )

    # Compute protease frequency matrix.
    absolute_frequency_matrix = generate_protease_pattern_frequency_matrix(
        patterns,
        protease_patterns,
        protease_labels
    )

    """
    # Set absolute frequency protease pattern heatmap title.
    absolute_heatmap_title = (
        analysis_title + ' - Protease Pattern Frequency (Absolute)'
    )
    
    # Set output path for absolute frequency protease pattern heat map.
    absolute_heatmap_output_path = (
        output_dir
        + '/figures/'
        + output_prefix
        + '_absolute_heatmap.svg'
    )

    # Generate absolute frequency protease pattern heatmap.
    generate_protease_pattern_heatmap(
        absolute_heatmap_title,
        patterns,
        absolute_frequency_matrix,
        protease_labels,
        position_labels,
        absolute_heatmap_output_path
    )
    """

    # Generate percentage frequency matrix from absolute frequency matrix.
    percentage_frequency_matrix = generate_protease_pattern_percentage_frequency_matrix(
        absolute_frequency_matrix,
        len(sequence_df)
    )

    """
    # Set percentage frequency protease pattern heatmap title.
    percentage_heatmap_title = (
        analysis_title + ' - Protease Pattern Frequency (Percentage)'
    )

    # Set output path for absolute frequency protease pattern heat map.
    percentage_heatmap_output_path = (
        output_dir
        + '/figures/'
        + output_prefix
        + '_percentage_heatmap.svg'
    )

    # Generate percentage frequency protease pattern heatmap.
    generate_protease_pattern_heatmap(
        percentage_heatmap_title,
        patterns,
        percentage_frequency_matrix,
        protease_labels,
        position_labels,
        percentage_heatmap_output_path
    )
    """

    # Generate sequence-pattern similarity matrix.
    pattern_labels = [
        label[:(label.find('{') - 4)] for label in generate_pattern_labels(
            position_labels,
            patterns
        )
    ]
    pattern_similarity_matrix = calculate_pattern_similarity_matrix(
        sequence_df,
        patterns.pattern_list,
        pattern_labels,
        SUBSTITUTION_MATRIX
    )

    # If clustermap annotations is enabled, generate row colors.
    row_color_ratio = 2
    if annotate_clustermap:
        row_colors = generate_row_colors(absolute_frequency_matrix, protease_labels)
        row_color_ratio = 2
    else:
        row_colors = None

    # Set clustermap title.
    clustermap_title = (
        analysis_title
        + ' - Mean Sequence-Pattern Positional Substitution Probability'
    )

    # Set clustermap output path.
    clustermap_output_path = (
        output_dir
        + '/figures/'
        + output_prefix
        + 'sequence_clustermap.svg'
    )

    # Pattern similarity clustermap.
    if np.any(pattern_similarity_matrix.to_numpy().astype(np.bool)):
        try:
            generate_sequence_clustermap(
                clustermap_title,
                pattern_similarity_matrix,
                clustermap_output_path,
                row_colors=row_colors,
                row_color_ratio=row_color_ratio,
                annotate_rows=False
            )
        except ValueError as e:
            print('Clustermap failed:\n')
            print(e)
            pass