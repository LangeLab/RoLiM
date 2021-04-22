import os
import re
import itertools
import io
import uuid
from collections import namedtuple

import numpy as np
import pandas as pd
import scipy.stats as st
from scipy import special
import weblogo as w

from zigzag import vis
from zigzag import merops_connector
from zigzag import amino_acid_encoding_converter
from zigzag import sequences
from zigzag import merops_connector

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 100000)

ParentStats = namedtuple(
    'ParentStats',
    ['pattern_id', 'size', 'bonferroni_m', 'expected_frequency']
)


def pattern_sequence_match(pattern_matrix, background, sequence_tensor):
    """
    Matches pattern to all sequences in a sample. Returns a list of the same
        length as the sample, containing an index-matched 0 for non-matching
        sequences or 1 for matching sequences.

    Parameters:
        pattern -- Pattern instance.

    Returns:
        pattern_sequence_matches -- List.
    """
        
    # Get constituent pattern.
    try:
        constituent_pattern = sequences.get_pattern_constituents(
            pattern_matrix,
            background
        )
    except AttributeError:
        constituent_pattern = pattern_matrix

    # Intersect pattern with sample sequence tensor.
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
    pattern_sequence_matches  = np.zeros_like(positional_matches)
    pattern_sequence_matches[
        positional_matches == np.sum(np.any(constituent_pattern, axis=1).astype(np.int8))
    ] = 1

    return pattern_sequence_matches


def pattern_to_regular_expression(character_pattern):
    '''
    Converts series-style patterns to regular expressions.

    Parameters:
        patterns -- List.
    
    Returns:
        regular_expressions -- List.
    '''

    regular_expression = r"(?=("
    separation = 0
    for position in character_pattern.tolist():
        stripped_position = position.replace('[', '').replace(']', '')
        if stripped_position != '.':
            if separation != 0:
                regular_expression += (
                    r"(.){"
                    + re.escape(str(separation))+ r"}" 
                    + position
                )
                separation = 0
            else:
                regular_expression += position
        else:
            separation += 1
    regular_expression += r"))"
    
    return regular_expression


def extract_protease_substrate_patterns(background=None,
                                        percentage_frequency_cutoff=0.05,
                                        max_depth=None,
                                        p_value_cutoff=0.001,
                                        minimum_occurrences=2,
                                        fold_change_cutoff=1,
                                        multiple_testing_correction=True,
                                        positional_weighting=True,
                                        allow_compound_residue_decomposition=True,
                                        enable_compound_residues=True,
                                        position_specific=True,
                                        width=8,
                                        center=True):
    '''
    Run pattern extraction on MEROPS protease substrate data sets.
    
    Parameters:
        

    Returns:
        None
    '''

    # Replace existing protease_patterns table with empty table.
    connection = merops_connector.connect_to_merops_database()

    if enable_compound_residues:
        protease_pattern_table = 'protease_patterns'
    else:
        protease_pattern_table = 'protease_patterns_no_cr'

    with connection.cursor() as cursor:
        cursor.execute("SHOW TABLES LIKE '{}'".format(protease_pattern_table))
        if cursor.fetchone():
            cursor.execute("DROP TABLE {}".format(protease_pattern_table))
            connection.commit()
    merops_connector.create_protease_patterns_table(
        enable_compound_residues=enable_compound_residues
    )

    if enable_compound_residues:
        compound_residues = sequences.COMPOUND_RESIDUES
    else:
        compound_residues = None

    # If no background is specified, use default SwissProt background.
    if background == None:
        # Load context sequences from fasta.
        context = sequences.import_fasta(
            os.path.join('media', 'defaults', 'uniprot.fasta')
        )
        # Generate new Background instance.
        background = sequences.Background(
            context['sequence'].tolist(),
            compound_residues=compound_residues,
            position_specific=position_specific,
            width=width,
            center=center
        )

    # Get substrates for all MEROPS proteases and loop through the sets.
    pattern_containers = []
    protease_substrates = merops_connector.retrieve_substrates()
    for protease, substrates in protease_substrates.items():
        print('{}\n'.format(protease))
        # Initialize output directory name.
        if enable_compound_residues:
            protease_output_directory = os.path.join(
                'media', 'merops', 'proteases', protease.replace(" ", "_")
            )
        else:
            protease_output_directory = os.path.join(
                'media', 'merops', 'proteases_no_cr', protease.replace(" ", "_")
            )
        # Prepare substrate sets for pattern extraction.
        single_letter_substrates = amino_acid_encoding_converter.convert_encoding(
            substrates,
            1
        )
        single_letter_substrates.drop_duplicates(inplace=True)
        substrate_tensor = sequences.vectorize_sequences(
            single_letter_substrates,
            background
        )
        protease_substrate_sample = sequences.Sample(
            sequence_df=single_letter_substrates,
            sequence_tensor=substrate_tensor
        )
        
        min_occurrences = (
            percentage_frequency_cutoff
            * protease_substrate_sample.sequence_tensor.shape[0]
        )
        if min_occurrences < 2:
            min_occurrences = 2
        
        # Run pattern extraction on substrate set.
        patterns = PatternContainer(
            protease_substrate_sample,
            background,
            protease,
            protease_output_directory,
            minimum_occurrences=min_occurrences
        )

        # Post-process extracted patterns and save pattern outputs.
        patterns.prune_patterns()

        if len(patterns.pattern_list) > 0:
            # Prepare output directories.
            try:
                os.makedirs(protease_output_directory)
            except FileExistsError:
                pass
            patterns.generate_pattern_outputs()

            # Generate clustermap for protease substrate set. 
            clustermap_title = (
                protease
                + ' Substrates -'
                + ' Mean Sequence-Pattern Positional Substitution Probability'
            )
            clustermap_output_directory = protease_output_directory + '/figures'
            try:
                os.makedirs(clustermap_output_directory)
            except:
                pass
            clustermap_output_path = clustermap_output_directory + '/clustermap.svg'

            position_labels = vis.generate_position_labels(
                protease_substrate_sample.sequence_df
            )
            pattern_labels = [
                label[:(label.find('{') - 4)] for label in vis.generate_pattern_labels(
                    position_labels,
                    patterns
                )
            ]
            pattern_similarity_matrix = vis.calculate_pattern_similarity_matrix(
                protease_substrate_sample.sequence_df,
                patterns.pattern_list,
                pattern_labels,
                vis.SUBSTITUTION_MATRIX
            )

            if np.any(pattern_similarity_matrix.to_numpy().astype(np.bool)):
                try:
                    vis.generate_sequence_clustermap(
                        clustermap_title,
                        pattern_similarity_matrix,
                        clustermap_output_path
                    )
                except ValueError:
                    print('clustermap for {} failed.'.format(protease))

            clusters = [
                ','.join(pattern.character_pattern().tolist())
                for pattern in patterns.pattern_list
            ]

            substrate_patterns = {protease: clusters}
            merops_connector.insert_protease_patterns(
                substrate_patterns,
                enable_compound_residues=enable_compound_residues
            )
            pattern_containers.append(patterns)

    for pattern_container in pattern_containers:
        generate_merops_heatmap(pattern_container, position_labels)


def generate_merops_heatmap(pattern_container, position_labels):
    """
    Generate heatmap for each pattern in MEROPS database with detected
        patterns.

    Parameters:
        pattern_container -- PatternContainer instance.
        position_labels -- String.

    Returns:
        None
    """
    protease_pattern_heatmap_title = (
        pattern_container.title + ' - Protease Pattern Matches (Percent Positions Matched)'
    )
    output_prefix = re.sub(r'\W+', ' ', pattern_container.title).strip().replace(" ", "_")
    # Set output path for absolute frequency protease pattern heat map.
    protease_pattern_heatmap_output_path = (
        pattern_container.output_directory
        + '/figures/'
        + output_prefix
        + '_protease_pattern_heatmap.svg'
    )
    if pattern_container.background.compound_residues is not None:
        enable_compound_residues = True
    else:
        enable_compound_residues = False
    protease_patterns = merops_connector.retrieve_protease_patterns(
        enable_compound_residues=enable_compound_residues
    )
    protease_labels = vis.generate_protease_labels(protease_patterns)
    non_exact_scoring_matrix = vis.generate_non_exact_protease_pattern_matrix(
        pattern_container,
        protease_patterns,
        protease_labels
    )
    # Generate absolute frequency protease pattern heatmap.
    vis.generate_protease_pattern_heatmap(
        protease_pattern_heatmap_title,
        pattern_container,
        non_exact_scoring_matrix,
        protease_labels,
        position_labels,
        protease_pattern_heatmap_output_path
    )


class Pattern:
    """

    Attributes:
        subset_tensor -- 3D NumPy Array. Dtype = int8.
        pattern_matrix -- 2D NumPy Array. Dtype = int8.
        removed_positional_residues -- 2D NumPy Array. Dtype = int8.

    Methods:

    """

    def __init__(self,
                    pattern_container,
                    parent_stats,
                    subset_tensor,
                    initial_pattern,
                    removed_positional_residues,
                    current_depth,
                    max_depth,                      
                    p_value_cutoff=0.001,
                    fold_change_cutoff=1.0,
                    minimum_occurrences=20,
                    positional_weighting=True,
                    multiple_testing_correction=True,
                    allow_compound_residue_decomposition=True,
                    set_reduction=True):
        """

        Parameters:
            pattern_container -- PatternContainer instance.
            subset_tensor -- 3D NumPy Array. Dtype = int8.
            pattern_matrix -- 2D NumPy Array. Dtype = int8.
            current_depth -- Int.
            max_depth -- Int.
            p_value_cutoff -- Float.
            minimum_occurrences -- Int.
            fold_change_cutoff -- Float.
            positional_weighting -- Bool.
            multiple_testing_correction -- Bool.

        Returns:
            None
        """

        # Initialize pattern attributes.
        self.pattern_id = uuid.uuid4().hex
        self.pattern_container = pattern_container
        self.pattern_container.pattern_ids[self.pattern_id] = None
        self.background = pattern_container.background
        self.subset_tensor = subset_tensor
        self.pattern_matrix = initial_pattern.copy()
        self.parent_stats = parent_stats
        self.invalid_pattern = False

        print(' '*current_depth*4 + 'new:\n' + ' '*current_depth*4 + ''.join(self.character_pattern().tolist()))

        if allow_compound_residue_decomposition == True:
            self.removed_positional_residues = removed_positional_residues.copy()
        else:
            # Mask positions containing fixed residue
            self.removed_positional_residues = np.logical_or(
                np.any(self.pattern_matrix, axis=1)[:, np.newaxis],
                removed_positional_residues
            ).astype(np.int8)

        """
        # Remove positional compound residues corresponding to extracted
        # positional simple residues.
        self.removed_positional_residues = (np.logical_or(
            self.pattern_matrix,
            np.concatenate(
                (self.pattern_matrix[:, :len(self.background.alphabet)],
                np.amax(
                    (self.background.compound_residue_matrix.to_numpy()
                    * self.pattern_matrix[:, :len(self.background.alphabet)][:, np.newaxis]),
                    axis=2
                )),
                axis=1
            )
        ).astype(np.int8))
        """

        # Iteratively extract significantly enriched positional residues.
        if current_depth < max_depth:
            while self.subset_tensor.shape[0] >= minimum_occurrences:
                # Remove pattern from tree if
                # prior set reduction has reduced subset tensor below
                # significance threshold in context of the parent tensor
                # in which the pattern was originally detected.
                if (
                    (np.any(self.pattern_matrix))
                    and (parent_stats.expected_frequency != parent_stats.size)
                ):
                    # Calculate probability of more extreme outcome via regularized
                    # incomplete beta function.
                    pattern_p_value = special.betainc(
                        self.subset_tensor.shape[0] + 1.0,
                        parent_stats.size - self.subset_tensor.shape[0],
                        parent_stats.expected_frequency
                    )
                    print(''.join(self.character_pattern()))
                    print(parent_stats.expected_frequency)
                    print(parent_stats.size)
                    print(parent_stats.expected_frequency * parent_stats.size)
                    print(self.subset_tensor.shape[0])
                    print(pattern_p_value)
                    print(p_value_cutoff / parent_stats.bonferroni_m)
                    print(self.removed_positional_residues)
                    # Remove from pattern list if p-value above cutoff.
                    if (
                        (pattern_p_value > (p_value_cutoff / parent_stats.bonferroni_m))
                        or (
                            self.subset_tensor.shape[0]
                            <= (parent_stats.expected_frequency * parent_stats.size)
                        )
                    ):
                        """
                        if (
                            self.pattern_container.pattern_ids[parent_stats.pattern_id]
                            is not None
                        ):
                            self.pattern_container.pattern_ids[
                                parent_stats.pattern_id
                            ] = np.concatenate(
                                (
                                    self.pattern_container.pattern_ids[parent_stats.pattern_id],
                                    self.subset_tensor
                                ),
                                axis=0
                            )
                        else:
                            self.pattern_container.pattern_ids[
                                parent_stats.pattern_id
                            ] = self.subset_tensor
                        """
                        print('Invalid pattern.')
                        self.invalid_pattern = True

                # Fix positional residues present in all sequences.
                self.fix_homogeneous_positional_residues()
                print(
                    (' ' * current_depth * 4)
                    + 'homogeneous:\n'
                    + (' ' * current_depth * 4)
                    + ''.join(self.character_pattern().tolist())
                )
                
                # Calculate positional residue frequencies from subset_tensor.
                positional_residue_frequencies = self.calculate_frequencies()
                
                if self.background.position_specific:
                    positional_background_frequencies = (
                        self.calculate_position_specific_background_frequencies(
                            allow_compound_residue_decomposition
                        )
                    )
                else:
                    # Calculate positional background frequency adjustments.
                    positional_background_frequencies = (
                        self._adjust_positional_background_frequencies()
                    )
                    #positional_background_frequencies[
                    #    positional_background_frequencies > 1.0] = 1.0

                # Calculate positional residue fold change.
                positional_residue_fold_change = np.divide(
                    (positional_residue_frequencies / self.subset_tensor.shape[0]),
                    positional_background_frequencies,
                    out=np.ones_like(positional_background_frequencies),
                    where=(positional_background_frequencies != 0.0)
                )

                # Calculate binomial p values from positional residues
                # and adjusted positional residue background frequencies.
                binomial_p_values = self.calculate_binomial_p_values(
                    positional_residue_frequencies,
                    positional_background_frequencies
                )
                # If p-value is smaller than float 64 precision,
                # replace p-value with smallest float 64 decimal.
                binomial_p_values[binomial_p_values == 0.0] = np.nextafter(
                    0,
                    1,
                    dtype=np.float64
                )

                # Mask removed positional residue p-values for
                # significance test.
                removed_positional_residue_mask = np.invert(
                    self.removed_positional_residues.astype(np.bool)
                )
                binomial_p_values = (
                    binomial_p_values * removed_positional_residue_mask
                )
                binomial_p_values[binomial_p_values == 0.0] = 1.0

                # Multiple testing correction calculations.
                if multiple_testing_correction:
                    bonferroni_m = np.count_nonzero(
                        positional_residue_frequencies
                        * removed_positional_residue_mask
                    )
                    # If only ambiguous positions remain, set m=1.
                    if bonferroni_m == 0.0:
                        bonferroni_m = 1.0
                    binomial_p_values[binomial_p_values
                                        > (p_value_cutoff/bonferroni_m)] = 1.0
                
                # Fold change and minimum occurrence enforcement.
                binomial_p_values[positional_residue_fold_change
                    <= fold_change_cutoff] = 1.0
                binomial_p_values[positional_residue_frequencies
                    < minimum_occurrences] = 1.0

                # Exit loop if no significantly enriched positional
                # residues detected.
                if (binomial_p_values[(binomial_p_values > 0)
                    * (binomial_p_values < p_value_cutoff)].size == 0):
                    print('end {}'.format(''.join(self.character_pattern())))
                    break

                # Calculate positional weights.
                positional_weights = self.calculate_positional_weights(
                    positional_residue_frequencies,
                    positional_residue_fold_change,
                    positional_background_frequencies
                )

                # Calculate enrichments.
                if positional_weighting:
                    enrichment_values = self.calculate_enrichment(
                        binomial_p_values,
                        positional_weights
                    )
                else:
                    enrichment_values = (
                        -1.0 * np.log(binomial_p_values,
                                        out=np.zeros_like(binomial_p_values),
                                        where=binomial_p_values!=0.0)
                    )

                # Select most enriched positional residue and generate
                # new pattern. If multiple positional residues have
                # equal maximal enrichment, select all and generate
                # new pattern for each. 
                max_positional_residue_enrichments = self.maximum_enrichment(
                    enrichment_values
                )

                # Loop through all maximally enriched positional
                # residue index values.
                """
                if len(np.transpose(
                        np.nonzero(max_positional_residue_enrichments))) > 1:
                    print([tuple(r) for r in np.transpose(
                        np.nonzero(max_positional_residue_enrichments))])
                    print(enrichment_values)
                    input()
                """
                # Get list of maximally enriched positional residue
                # coordinates.
                enriched_positional_residues = [
                    tuple(enriched_positional_residue)
                    for enriched_positional_residue
                    in np.transpose(np.nonzero(max_positional_residue_enrichments))
                ]
                
                new_patterns = []
                prior_removed_positional_residues = self.removed_positional_residues.copy()
                for m, enriched_positional_residue in enumerate(enriched_positional_residues):
                    if m > 0:
                        print('simultaneous branch point.')
                    # Set enriched residue coordinates.
                    residue_coordinates = enriched_positional_residue
                    # Generate a new pattern for each maximally enriched
                    # positional residue and add to PatternContainer
                    # pattern list.
                    new_pattern = np.zeros_like(self.pattern_matrix)
                    new_pattern[residue_coordinates] = 1
                        
                    # Get subset of sequences matching new pattern
                    matching_sequences = self.subset_tensor[
                        tuple(
                            np.where(
                                np.any(
                                    np.any(
                                        self.subset_tensor * new_pattern,
                                        axis=1
                                    ),
                                    axis=1)
                                )
                            )
                    ]

                    # Add new positional residue to removed positional
                    # residue array.
                    new_pattern_removed_positional_residues = prior_removed_positional_residues.copy()
                    self.removed_positional_residues[residue_coordinates] = 1
                    new_pattern_removed_positional_residues[residue_coordinates] = 1

                    # Add simple residues to removed positional residue
                    # array if enriched compound residue extracted.
                    if (residue_coordinates[1]
                            >= len(self.background.alphabet)):
                        if not allow_compound_residue_decomposition:
                            new_pattern_removed_positional_residues = (
                                self._remove_compound_residue_constituents(
                                    residue_coordinates,
                                    new_pattern_removed_positional_residues
                                )
                            )
                        else:
                            self._remove_compound_residue_constituents(
                                residue_coordinates,
                                new_pattern_removed_positional_residues
                            )

                    # Add parent pattern residues to new pattern.
                    prior_positions = self.pattern_matrix.copy()
                    prior_positions[residue_coordinates[0]] = 0
                    new_pattern = np.logical_or(
                        prior_positions,
                        new_pattern
                    ).astype(np.int8)

                    try:
                        stats_bonferroni_m = bonferroni_m
                    except UnboundLocalError:
                        stats_bonferroni_m = 1.0
                    stats = ParentStats(
                        pattern_id=self.pattern_id,
                        size=self.subset_tensor.shape[0],
                        bonferroni_m=stats_bonferroni_m,
                        expected_frequency=positional_background_frequencies[residue_coordinates]
                    )

                    # Instantiate new Pattern object with new pattern,
                    # matching sequences, and current depth.
                    pattern_container.add_new_pattern(
                        Pattern(
                            pattern_container,
                            stats,
                            matching_sequences,
                            new_pattern,
                            new_pattern_removed_positional_residues,
                            current_depth + 1,
                            max_depth,
                            p_value_cutoff=p_value_cutoff,
                            fold_change_cutoff=fold_change_cutoff,
                            minimum_occurrences=minimum_occurrences,
                            positional_weighting=positional_weighting,
                            multiple_testing_correction=multiple_testing_correction,
                            allow_compound_residue_decomposition=allow_compound_residue_decomposition,
                            set_reduction=set_reduction
                        )
                    )
                    new_patterns.append(new_pattern)

                # Remove sequences matching new pattern(s) from subset.
                if set_reduction:
                    self._remove_sequences(max_positional_residue_enrichments)
                """
                # Add sequences from invalid patterns back to tensor.
                if self.pattern_container.pattern_ids[self.pattern_id] is not None:
                    self.subset_tensor = np.concatenate(
                        (
                            self.subset_tensor,
                            self.pattern_container.pattern_ids[self.pattern_id]
                        ),
                        axis=0
                    )
                """


    def calculate_position_specific_background_frequencies(self,
                                                            allow_compound_residue_decomposition):
        # Get downstream fixed residues that do not include pattern.
        if self.background.compound_residues is not None:
            constituent_pattern = sequences.get_pattern_constituents(
                self.pattern_matrix,
                self.background
            )
            if allow_compound_residue_decomposition:
                constituents_only = constituent_pattern - self.pattern_matrix
                downstream_fixed_positional_residues = (
                    self.removed_positional_residues
                    + constituents_only
                    - constituent_pattern
                )
            else:
                downstream_fixed_positional_residues = (
                    self.removed_positional_residues
                    - constituent_pattern
                )
        else:
            downstream_fixed_positional_residues = (
                self.removed_positional_residues
                - self.pattern_matrix
            )

        # Memory-intensive, high-performance version.
        if self.background.fast:
            # Restrict background to sequences matching pattern.
            background_tensor = self.background.background_tensor[
                np.logical_and(
                    self.background.background_tensor,
                    np.logical_or(
                        self.pattern_matrix,
                        np.invert(
                            self.pattern_matrix.any(axis=1)
                        )[..., np.newaxis]
                    )
                ).any(axis=(2)).all(axis=1),
                :,
                :
            ]
            # Eliminate background sequences which have been removed downstream.
            background_tensor = self.background_tensor[
                ~np.logical_and(
                    background_tensor,
                    np.logical_or(
                        downstream_fixed_positional_residues,
                        np.invert(
                            downstream_fixed_positional_residues.any(axis=1)
                        )[..., np.newaxis]
                    )
                ).any(axis=2).all(axis=1),
                :,
                :
            ]
            positional_background_frequencies = np.divide(
                np.sum(
                    background_tensor,
                    axis=0,
                    dtype=np.float64
                ),
                background_tensor.shape[0]
            )
        # Light-memory, slow version.
        else:
            # Get character representation of pattern.
            if self.background.compound_residues is not None:
                character_pattern = self.character_pattern(
                    pattern_matrix=constituent_pattern
                )
            else:
                character_pattern = self.character_pattern()
            
           #  Get subset of background matching pattern.
            background_df = self.background.background_df
            for position, residues in character_pattern.iteritems():
                if residues == '.':
                    continue
                stripped_residues = residues.replace('[', '').replace(']', '')
                background_df = background_df[
                    background_df[position].isin(list(stripped_residues))
                ]
            
            # Eliminate background sequences matching dowstream fixed residues.
            fixed_downstream_characters = self.character_pattern(
                pattern_matrix=downstream_fixed_positional_residues
            )
            for position, residues in fixed_downstream_characters.iteritems():
                if residues == '.':
                    continue
                stripped_residues = residues.replace('[', '').replace(']', '')
                background_df = background_df[
                    ~background_df[position].isin(list(stripped_residues))
                ]

            # Calculate absolute simple residue background frequencies.
            positional_background_frequencies = []
            for column in background_df.columns:
                positional_frequencies = background_df[column].value_counts()
                positional_frequencies = positional_frequencies[
                    positional_frequencies.index.isin(self.background.alphabet)
                ]
                positional_background_frequencies.append(
                    pd.Series(
                        [0] * len(self.background.alphabet),
                        index=self.background.alphabet
                    ).add(positional_frequencies).tolist()
                )
            positional_background_frequencies = np.array(
                positional_background_frequencies,
                dtype=np.float64
            )

            # Add absolute compound residue background frequencies.
            if self.background.compound_residues is not None:
                compound_residue_frequencies = np.inner(
                    positional_background_frequencies,
                    self.background.compound_residue_matrix
                )
                positional_background_frequencies = np.concatenate(
                    (positional_background_frequencies, compound_residue_frequencies),
                    axis=1
                )

            # Convert absolute frequencies to percentage frequencies.
            positional_background_frequencies = (
                positional_background_frequencies
                / len(background_df)
            )
            np.nan_to_num(positional_background_frequencies, copy=False)

        return positional_background_frequencies

    def fix_homogeneous_positional_residues(self):
        """
        Fix positions containing either a single simple or compound
            residue in pattern matrix and add to removed positional
            residue matrix.

        Parameters:
            None

        Returns:
            None
        """
        # Set num_residue_constituents array for simple residues.
        num_residue_constituents = np.ones_like(
            self.background.background_vector,
            dtype=np.int8
        )
        # Set num_residue_constituents for compound residues.
        try:
            num_compound_residue_constituents = np.sum(
                self.background.compound_residue_matrix.to_numpy(),
                axis=1,
                dtype=np.int8
            )
            num_residue_constituents[
                len(self.background.alphabet):] = num_compound_residue_constituents
        except AttributeError:
            pass

        # Set homogeneous positional residues to 1.0, others to 0.0.
        homogeneous_positional_residues = self.calculate_frequencies()
        homogeneous_positional_residues[
            homogeneous_positional_residues!=self.subset_tensor.shape[0]] = 0.0
        homogeneous_positional_residues[
            homogeneous_positional_residues==self.subset_tensor.shape[0]] = 1.0
        # Proceed if any homogeneous positional residues exist.
        if np.any(homogeneous_positional_residues):
            # Calculate most specific homogeneous residue for each position.
            specificity_ratio = homogeneous_positional_residues / num_residue_constituents
            homogeneous_positional_residues = (
                specificity_ratio == np.amax(
                    specificity_ratio,
                    axis=1,
                    keepdims=True
                )
            ).astype(np.int8) * homogeneous_positional_residues.astype(np.int8)
            homogeneous_positional_residues = (
                homogeneous_positional_residues * np.invert(
                    np.any(self.pattern_matrix, axis=1, keepdims=True)
                )
            )
            # Add homogeneous positional residues to pattern matrix.
            self.pattern_matrix = np.logical_or(
                self.pattern_matrix,
                homogeneous_positional_residues
            ).astype(np.int8)
            # Add homogeneous positional residues to removed positional
            # residues.
            self.removed_positional_residues = np.logical_or(
                self.removed_positional_residues,
                homogeneous_positional_residues
            ).astype(np.int8)

    def _adjust_positional_background_frequencies(self):
        '''
        Dynamically calculate expected residue frequencies based on
            fixed pattern residues and residues removed from non-fixed
            positions in prior iterations.

        Parameters:
            None

        Returns:
            vectorized_background_frequency_adjustment_values --
                                1D NumPy Array. Dtype = float64.
        '''

        # Get matrix of all removed simple residues.
        removed_simple_residues = self.removed_positional_residues[
            :,
            :len(self.background.alphabet)
        ]
        
        # Calculate simple residue percentage background frequencies.
        positional_background_frequencies = self._calculate_background_frequencies(
            removed_simple_residues
        )

        # Replace removed positional simple residues frequencies with 0.
        positional_background_frequencies = (
            positional_background_frequencies
            * np.invert(
                self.removed_positional_residues[
                    :,
                    :len(self.background.alphabet)
                ].astype(np.bool))
        )

        # Replace fixed position frequencies with 0.
        # fixed_positions = self.pattern_matrix[
        #     :,
        #     :len(self.background.alphabet)
        # ].copy()
        # fixed_positions[np.any(fixed_positions, axis=1)] = 1
        
        # Set all frequencies to 0.0 in fixed positions.
        fixed_positions = np.zeros_like(self.pattern_matrix)
        fixed_positions[np.any(self.pattern_matrix, axis=1)] = 1
        fixed_positions = fixed_positions[:, :len(self.background.alphabet)]
        positional_background_frequencies = (
            positional_background_frequencies
            * np.invert(fixed_positions.astype(np.bool)).astype(np.float64)
        )

        try:
            # Get compound residue matrix from background.
            compound_residue_matrix = self.background.compound_residue_matrix.to_numpy()
        except AttributeError:
            # Fixed position background frequencies should be 1.0 for
            # fixed residues and 0.0 elsewhere.
            fixed_position_background_frequencies = self.pattern_matrix.astype(np.float64)
        else:
            # Merge compound and simple residue percentage frequencies.
            positional_background_frequencies = self._compound_residue_background_frequencies(
                compound_residue_matrix,
                positional_background_frequencies
            )

            ## Get all simple residues from pattern.
            # Extend compound residues with blank.
            constituents = np.concatenate(
                (
                    compound_residue_matrix,
                    np.zeros((1, compound_residue_matrix.shape[1]), dtype=np.int8)
                ),
                axis=0
            )

            # Generate index array from pattern matrix.
            compound_residue_index = np.argmax(
                self.pattern_matrix,
                axis=1
            ) - len(self.background.alphabet)
            
            # Map simple positions to blank in compound residue matrix.
            compound_residue_index[compound_residue_index < 0] = len(compound_residue_matrix)

            # Get compound residue constituents using index array.
            positional_residue_constituents = constituents[compound_residue_index]

            # Drop removed constituents from previous decomposition.
            postional_residue_constituents = (
                positional_residue_constituents
                - np.logical_and(
                    positional_residue_constituents,
                    self.removed_positional_residues[
                        :,
                        :len(self.background.alphabet)
                    ]
                ).astype(np.int8)
            )
            # Calculate background frequencies for constituents.
            unadjusted_constituent_frequencies = (
                positional_residue_constituents
                * self.background.background_vector[:len(self.background.alphabet)]
            )
            # Calculate total positional background frequencies.
            total_positional_constituent_frequencies = np.sum(
                unadjusted_constituent_frequencies,
                axis=1
            )[:, np.newaxis]
            # Calculate weighted constituent background frequencies.
            adjusted_constituent_frequencies = np.divide(
                unadjusted_constituent_frequencies,
                total_positional_constituent_frequencies,
                out=np.zeros_like(unadjusted_constituent_frequencies),
                where=total_positional_constituent_frequencies!=0
            )
            # Pad array with zeros in compound residue positions.
            fixed_position_background_frequencies = np.zeros_like(positional_background_frequencies)
            fixed_position_background_frequencies[
                :,
                :len(self.background.alphabet)
            ] = adjusted_constituent_frequencies

        # Merge fixed positions back into background frequencies.
        positional_background_frequencies = (
            positional_background_frequencies
            + fixed_position_background_frequencies
        )

            # # Combine constiuents of fixed compounds with fixed simple.
            # pattern_constituent_matrix = np.logical_or(
            #     self.pattern_matrix[:, :len(self.background.alphabet)],
            #     positional_residue_constituents
            # )

            # # Get complement of fixed residues.
            # removed_fixed_position_residues = np.invert(
            #     pattern_constituent_matrix
            # ).astype(np.int8)

            # Calculate fixed positional background frequencies.
            # fixed_simple_residue_background_frequencies = self._calculate_background_frequencies(
            #     removed_fixed_position_residues
            # )

            # # Concatenate fixed position compound residues.
            # fixed_position_background_frequencies = self._compound_residue_background_frequencies(
            #     compound_residue_matrix,
            #     fixed_simple_residue_background_frequencies
            # )

        # # Merge adjusted fixed- and non-fixed position background frequencies.
        # positional_background_frequencies = (
        #     positional_background_frequencies
        #     + fixed_position_background_frequencies
        # )



        return positional_background_frequencies

    def _compound_residue_background_frequencies(self,
                                                    compound_residues,
                                                    background_frequencies):
        """

        Parameters:
            compound_residues -- Numpy Array.
            background_frequencies -- Numpy Array.
        
        Returns:
            positional_background_frequencies -- Numpy Array.
        """

        # Calculate adjusted compound residue background frequencies.
        positional_compound_residue_frequencies = np.sum(
            (compound_residues
            * background_frequencies[:, np.newaxis]),
            axis=2
        )

        # Merge non-fixed position compound and simple residues.
        positional_background_frequencies = np.concatenate(
            (background_frequencies,
                positional_compound_residue_frequencies),
            axis=1
        )

        return positional_background_frequencies

    def _calculate_background_frequencies(self, removed_positional_residues):
        """

        Parameters:
            removed_positional_residues -- Numpy Array
        
        Returns:
            background_frequencies -- Numpy Array
        """

        # Calculate cumulative extracted background percentage.
        positional_penalties = np.sum(
            (removed_positional_residues
                * self.background.background_vector[:len(self.background.alphabet)]),
            axis=1
        )

        background_frequencies = np.outer(
            np.reciprocal(
                (100 - positional_penalties),
                out=np.zeros_like(positional_penalties),
                where=positional_penalties < 100
            ),
            self.background.background_vector[:len(self.background.alphabet)]
        )

        return background_frequencies

    def calculate_frequencies(self):
        '''

        Parameters:
            vectorized_sample -- 3D NumPy Array. Dtype = int8.
    
        Returns:
            vectorized_effective_frequencies -- 2D NumPy Array.
                                                        Dtype = int8.
        '''

        positional_residue_frequencies = np.sum(self.subset_tensor, axis=0)

        return positional_residue_frequencies

    def calculate_binomial_p_values(self,
                                    positional_residue_frequencies,
                                    positional_background_frequencies):
        '''
        # apply positional background frequency adjustments to
        background_frequencies to generate vectorized background
        frequencies array
        # calculate binomial p-value for all positional residues
        using positional_residue_frequencies as x, n, and
        adjusted background frequency array as expected

        Parameters:
            positional_residue_frequencies -- 2D NumPy Array.
                                                        Dtype = int8.
            positional_background_frequencies --
                                    2D NumPy Array. Dtype = float64.

        Returns:
            binomial_p_values -- 2D NumPy Array. Dtype = float64.
        '''

        flattened_positional_background_frequencies = positional_background_frequencies.flatten()
        binomial_p_values = []
        
        for i, frequency in enumerate(positional_residue_frequencies.flatten()):
            if frequency == self.subset_tensor.shape[0]:
                binomial_p_values.append(np.nextafter(0, 1, dtype=np.float64))
            elif frequency == 0:
                binomial_p_values.append(1.0)
            else:
                # Calculate probability of more extreme outcome via regularized
                # incomplete beta function.
                binomial_p_value = special.betainc(
                    frequency + 1.0,
                    self.subset_tensor.shape[0] - frequency,
                    flattened_positional_background_frequencies[i]
                )
                binomial_p_values.append(binomial_p_value)

        binomial_p_values = np.array(
            binomial_p_values,
            dtype=np.float64
        ).reshape(positional_residue_frequencies.shape)

        return binomial_p_values

    def calculate_positional_weights(self,
                                        positional_residue_frequencies,
                                        positional_residue_fold_change,
                                        positional_background_frequencies):
        '''

        Parameters:
            positional_residue_frequencies -- 2D NumPy Array.
                                                Dtype = int8.

        Returns:
            positional_weights -- 1D NumPy Array.
                                    Dtype = float64.
            positional_residue_fold_change -- 2D Numpy Array.
            positional_background_frequencies -- 2D Numpy Array.
        '''

        ################################################################
        ### NON-NORMALIZED POSITIONAL WEIGHT CALCULATION
        ################################################################
        # positional_weights = (
        #     1 / np.count_nonzero(
        #         positional_residue_frequencies[:, :len(self.background.alphabet)],
        #         axis=1
        #     )
        # )
        ################################################################

        ################################################################
        ### NORMALIZED POSITIONAL WEIGHT CALCULATION
        ################################################################
        # 
        underrepresentation = np.subtract(
            1.0,
            positional_residue_fold_change[:, :len(self.background.alphabet)]
        )
        underrepresentation[underrepresentation <= 0.0] = 0.0
        positional_underrepresentation = np.sum(underrepresentation, axis=1)
        num_positional_background_residues = np.count_nonzero(
            positional_background_frequencies,
            axis=1
        )
        positional_weights = np.multiply(
            -1.0,
            np.log(
                np.subtract(
                    1.0,
                    np.divide(
                        positional_underrepresentation,
                        num_positional_background_residues
                    )
                )
            )
        )
        ################################################################

        return positional_weights

    def calculate_enrichment(self,
                                binomial_p_values,
                                positional_weights):
        '''

        Parameters:
            binomial_p_values -- 2D NumPy Array. Dtype = float64.
            positional_weights -- 1D NumPy Array. Dtype = float64.

        Returns:
            enrichment_values -- 2D NumPy Array. Dtype = float64.
        '''

        enrichment_values = np.transpose(
            np.multiply(
                positional_weights,
                -1.0 * np.log(
                    binomial_p_values.transpose(),
                    out=np.zeros_like(binomial_p_values.transpose()),
                    where=binomial_p_values.transpose()!=0
                ),
                out=np.zeros_like(binomial_p_values.transpose()),
                where=binomial_p_values.transpose()!=0
            )
        )

        return enrichment_values

    def maximum_enrichment(self, enrichment_values):
        '''

        Parameters:
            enrichment_values -- 2D NumPy Array. Dtype = float64.

        Returns:
            max_positional_residue_enrichments -- 2D NumPy Array.
                                                    Dtype = int8.
        '''

        max_positional_residue_enrichments = np.zeros_like(
            enrichment_values
        ).astype(np.int8)
        max_positional_residue_enrichments[
            enrichment_values == np.nanmax(enrichment_values)
        ] = 1

        return max_positional_residue_enrichments

    def _remove_compound_residue_constituents(self,
                                                residue_coordinates,
                                                new_pattern_removed_positional_residues):
        """
        Add simple residue constituents of compound residue to removed
            residue array.
        
        Parameters:
            residue_coordinates -- Tuple.
            background -- Background instance.

        Returns:
            removed_positional_residues -- 2D Numpy Array.
        """

        removed_positional_residues = new_pattern_removed_positional_residues.copy()
        # Get row index of compound residue in compound residue matrix.
        ind = residue_coordinates[1] - len(self.background.alphabet)
        # Get row from compound residue matrix as array.
        compound_residue_constituents = (
            self.background.compound_residue_matrix.iloc[ind]
        ).astype(np.int8)
        # Set consituents to 1 in removed positional residue array.
        self.removed_positional_residues[
            residue_coordinates[0],
            :len(self.background.alphabet)] = np.logical_or(
                self.removed_positional_residues[
                    residue_coordinates[0],
                    :len(self.background.alphabet)],
                compound_residue_constituents
        ).astype(np.int8)
        # Set constituents to 1 in new pattern removed positional residue
        # array.
        removed_positional_residues[
            residue_coordinates[0],
            :len(self.background.alphabet)] = np.logical_or(
                removed_positional_residues[
                    residue_coordinates[0],
                    :len(self.background.alphabet)],
                compound_residue_constituents
        ).astype(np.int8)

        return removed_positional_residues

    def _remove_sequences(self, max_positional_residue_enrichments):
        """
        Remove sequences containing enriched residue from subset.

        Parameters:
            max_positional_residue_enrichments --

        Returns:
            None
        """

        self.subset_tensor = self.subset_tensor[
            tuple(
                np.where(
                    np.any(
                        np.any(
                            self.subset_tensor
                            * max_positional_residue_enrichments,
                            axis=1
                        ),
                        axis=1
                    ) == False
                )
            )
        ]

    def character_pattern(self, pattern_matrix=None):
        '''

        Parameters:
            None

        Returns:
            pattern_df -- Pandas DataFrame.
        '''

        def get_positional_residue(row):
            '''

            Parameters:
                row -- Pandas Series.

            Returns:
                positional_residue -- String.
            '''

            positional_residue = '.'

            residue_codes = []
            for col in pattern_df.columns:
                if row[col] == 1:
                    residue_codes.append(col)

            if residue_codes:
                positional_residue = '[' + ''.join(residue_codes) + ']'

            return positional_residue

        if pattern_matrix is None:
            pattern_matrix = self.pattern_matrix

        pattern_df = pd.DataFrame(pattern_matrix,
                                    columns=self.background.ordered_residues)
        pattern_series = pattern_df.apply(get_positional_residue, axis=1)
        pattern_series.index = self.pattern_container.sample.sequence_df.columns

        return pattern_series

    def generate_sequence_strings(self, output_directory):
        """
        Decode each sequence in subset_tensor array to string format.
            Return all sequence strings in a list.

        Parameters:
            output_directory -- String.

        Sets:
            sequence_strings -- List. Contains sequences from
                                subset_tensor in string form.
        """
        
        def get_positional_residue(row):
            '''

            Parameters:
                row -- Pandas Series.
            Returns:
                positional_residue -- String.
            '''

            positional_residue = 'X'

            for col in sequence_df.columns:
                if row[col] == 1:
                    positional_residue = col
                    break

            return positional_residue

        sequence_strings = []
        for sequence_matrix in self.subset_tensor:
            sequence_df = pd.DataFrame(
                sequence_matrix[:, :len(self.background.alphabet)],
                columns=self.background.alphabet
            )
            sequence_series = sequence_df.apply(get_positional_residue, axis=1)
            sequence_strings.append(''.join(sequence_series.tolist()))
        self._sequence_strings = sequence_strings

        # Write sequence strings to text file.
        output_file = os.path.join(
            output_directory,
            ''.join(
                self.character_pattern()
            ).replace('[', '').replace(']', '').replace('.', 'x') + '_sequences.txt'
        )
        with open(output_file, 'w') as fout:
            for sequence in sequence_strings:
                fout.write('{}\n'.format(sequence))

    def generate_logo_map(self, output_directory):
        """
        Generate sequence logo map for set of sequences matching pattern.
        
        Parameters:
            output_directory -- String.

        Sets:
            logo_map -- String. Path to logo map file in SVG format.
        """


        fasta = '> \n' + '\n> \n'.join(self._sequence_strings)

        logo_alphabet = w.Alphabet(
            self.background.alphabet
            + ['X'] if 'X' not in self.background.alphabet
            else background.alphabet
        )

        baserules = [
            w.SymbolColor("GSTYC", "green", "polar"),
            w.SymbolColor("NQ", "purple", "neutral"),
            w.SymbolColor("KRH", "blue", "basic"),
            w.SymbolColor("DE", "red", "acidic"),
            w.SymbolColor("PAWFLIMV", "black", "hydrophobic")
        ]

        colorscheme = w.ColorScheme(baserules, alphabet = logo_alphabet)

        seqs = w.read_seq_data(io.StringIO(fasta), alphabet=logo_alphabet)
        data = w.LogoData.from_seqs(seqs)
        
        options = w.LogoOptions()
        options.color_scheme = colorscheme
        mformat = w.LogoFormat(data, options)
        
        output_file = os.path.join(
            output_directory,
            ''.join(
                self.character_pattern()
            ).replace('[', '').replace(']', '').replace('.', 'x') + '_logo_map.pdf'
        )
        with open(output_file, "wb") as f:
            f.write(w.pdf_formatter(data, mformat))

        self._logo_map = output_file

    @property
    def logo_map(self):
        return self._logo_map
    


class PatternContainer:
    '''
    Pattern tree data structure and accompanying methods.

    Attributes:
        pattern_list -- List.

    Methods:
        add_new_pattern --
    '''

    def __init__(self,
                    sample,
                    background,
                    title,
                    output_directory,
                    initial_pattern=None,
                    initial_removed_positional_residues=None,
                    max_depth=None,
                    p_value_cutoff=0.001,
                    minimum_occurrences=20,
                    fold_change_cutoff=1,
                    multiple_testing_correction=True,
                    positional_weighting=True,
                    allow_compound_residue_decomposition=True,
                    set_reduction=True):
        '''
        Initialize PatternContainer object and start recursive pattern
            extraction.

        Parameters:
            sample -- Sample instance.
            background -- Background class instance.
            initial_pattern -- 2D Numpy Array.
            max_depth -- Int.
            p_value_cutoff -- Float.
            minimum_occurrences -- Int.
            fold_change_cutoff -- Float.
            multiple_testing_correction -- Bool.
            positional_weighting -- Bool.

        Returns:
            None
        '''

        self.pattern_ids = {}
        self.pattern_list = []
        self.sample = sample
        self.background = background
        self.title = title
        self.output_directory = output_directory
        self.minimum_occurrences = minimum_occurrences

        # Initialize pattern template.
        if initial_pattern == None:
            initial_pattern = np.zeros((self.sample.sequence_tensor.shape[1],
                                    len(background.ordered_residues)),
                                    dtype=np.int8)
        if initial_removed_positional_residues == None:
            initial_removed_positional_residues = initial_pattern

        if multiple_testing_correction:
            bonferroni_m = np.count_nonzero(np.sum(sample.sequence_tensor, axis=0))
            if bonferroni_m == 0.0:
                bonferroni_m = 1.0
        else:
            bonferroni_m = 1.0
            
        parent_stats = ParentStats(
            pattern_id=None,
            size=sample.sequence_tensor.shape[0],
            bonferroni_m=bonferroni_m,
            expected_frequency=sample.sequence_tensor.shape[0]
        )

        # Begin recursive pattern tree construction.
        self.add_new_pattern(
            Pattern(
                pattern_container=self,
                parent_stats=parent_stats,
                subset_tensor=self.sample.sequence_tensor,
                initial_pattern=initial_pattern,
                removed_positional_residues=initial_removed_positional_residues,
                # Starting depth equal number of positions in initial_pattern.
                current_depth=np.count_nonzero(np.any(initial_pattern, axis=1)),
                # Max_depth either set by argument, or as len(initial_pattern)
                max_depth=(max_depth if max_depth is not None else len(initial_pattern)),
                p_value_cutoff=p_value_cutoff,
                minimum_occurrences=minimum_occurrences,
                fold_change_cutoff=fold_change_cutoff,
                positional_weighting=positional_weighting,
                multiple_testing_correction=multiple_testing_correction,
                set_reduction=set_reduction
            )
        )

    def add_new_pattern(self, new_pattern):
        '''
        Add new Pattern instance to pattern list.

        Parameters:
            new_pattern --

        Returns:
            None
        '''
        
        patterns = [
            np.array_str(pattern.pattern_matrix.flatten())
            for pattern in self.pattern_list
        ]
        
        if (
            new_pattern.pattern_matrix.any()
            and new_pattern.subset_tensor.shape[0] >= self.minimum_occurrences
            and np.array_str(new_pattern.pattern_matrix.flatten()) not in patterns
        ):
            
            self.pattern_list.append(new_pattern)

    def prune_patterns(self):
        """
        Drop patterns with fewer than required minimum sequences.

        Parameters:
            None

        Returns:
            None
        """

        drop_patterns = []
        for pattern in self.pattern_list:
            if pattern.subset_tensor.shape[0] < self.minimum_occurrences:
                drop_patterns.append(pattern.pattern_id)
            elif pattern.invalid_pattern:
                drop_patterns.append(pattern.pattern_id)

        self.pattern_list = [
            pattern for pattern in self.pattern_list
            if pattern.pattern_id not in drop_patterns
        ]

    def post_processing(self,
                        proteolysis_data=True,
                        cluster_sequences=True,
                        logo_maps=True):
        '''
        Removes patterns if their exclusive sequence subset is two
            sequences or less.

        Parameters:
            None
        
        Returns:
            None
        '''

        # Drop patterns with fewer than required minimum sequences.
        self.prune_patterns()

        # Generate outputs for each pattern.
        self.generate_pattern_outputs(logo_maps=logo_maps)

        # Generate heatmaps and clustermap.
        vis.generate_figures(
            self,
            proteolysis_data=proteolysis_data,
            cluster_sequences=cluster_sequences,
            annotate_clustermap=False
        )
        
        # Generate tabular output.
        self.save_summary_table()

        return self.sequence_summary_table

    @staticmethod
    def generate_summary_table(original_sequences, pattern_list):
        """

        Parameters:
            original_sequences -- Pandas DataFrame. Contains all rows
                                    from the supplied foreground data
                                    set in their original order, and
                                    supplemental data from the alignment
                                    and extension process if relevant.
            pattern_list -- List. Contains a Pattern object for each
                        unique pattern.
        
        Returns:
            summary_table -- Pandas DataFrame. Cotains original_sequences
                                columns, and an additional binary
                                columns for each unique pattern based on
                                sequence match.
        """

        # Match aligned sequences to detected patterns.
        pattern_rows = []
        for i, row in original_sequences.iterrows():
            if not isinstance(row['aligned_sequence'], str):
                pattern_rows.append([np.nan] * len(pattern_list))
            else:
                matching_patterns = []
                for pattern in pattern_list:
                    pattern_re = pattern_to_regular_expression(
                        pattern.character_pattern(
                            pattern_matrix=pattern.constituent_pattern
                        )
                    )
                    pattern_match = 0
                    for sequence in row['aligned_sequence'].split(';'):
                        if re.match(pattern_re, sequence):
                            pattern_match = 1
                            break
                    matching_patterns.append(pattern_match)
                pattern_rows.append(matching_patterns)

        # Generate tabular output DataFrame.
        pattern_columns = [
            ''.join(pattern.character_pattern()) for pattern in pattern_list
        ]
        summary_table = pd.DataFrame(
            pattern_rows,
            columns=pattern_columns,
            index=original_sequences.index
        )
        summary_table = pd.concat(
            [original_sequences, summary_table],
            axis=1
        )

        return summary_table

    def save_summary_table(self):
        """
        Output table conserving original input data set rows with
            addtional output columns.
        Parameters:

        Returns:

        """

        if self.sample.original_sequences is not None:
            # Match aligned sequences to detected patterns.
            summary_table = self.generate_summary_table(
                self.sample.original_sequences,
                self.pattern_list
            )

            # Save summary table to tab-delimited file.
            try:
                os.makedirs(os.path.join(self.output_directory, 'summary'))
            except FileExistsError:
                pass
            summary_table_path = os.path.join(
                self.output_directory,
                'summary',
                self.title.replace(' ', '_') + '_summary_table.txt'
            )
            summary_table.to_csv(
                summary_table_path,
                sep='\t',
                header=True,
                index=True
            )
            self.sequence_summary_table = summary_table

    def generate_pattern_outputs(self, logo_maps=True):
        """Last steps for retained patterns."""
        

        pattern_directory = os.path.join(
            self.output_directory,
            'patterns'
        )
        try:
            os.makedirs(pattern_directory)
        except FileExistsError:
            pass
        
        # Generate pattern summary table.
        pattern_labels = []
        sample_frequencies = []
        foreground_frequencies = []
        foreground_sizes = []
        enrichments = []
        for pattern in self.pattern_list:
            # Generate pattern label.
            pattern_label = ''.join(pattern.character_pattern())
            pattern_labels.append(pattern_label)
            
            # Save logo map and seqeunce strings for each pattern.  
            pattern.generate_sequence_strings(pattern_directory)
            if logo_maps:
                pattern.generate_logo_map(pattern_directory)
            
            # Calculate pattern statistics.
            sample_matches = pattern_sequence_match(
                pattern.pattern_matrix,
                self.background,
                self.sample.sequence_tensor
            )
            sample_frequency = np.sum(sample_matches)
            sample_frequencies.append(sample_frequency)
            foreground_frequencies.append(pattern.subset_tensor.shape[0])
            foreground_sizes.append(pattern.parent_stats.size)

            # Generate constituent pattern.
            try:
                pattern.constituent_pattern = sequences.get_pattern_constituents(
                    pattern.pattern_matrix,
                    self.background
                )
            except AttributeError:
                pattern.constituent_pattern = pattern.pattern_matrix

            # Calculate pattern enrichment if position-specific background.
            if self.background.position_specific:
                sample_percentage_frequency = (
                    sample_frequency / len(sample_matches)
                )
                pattern_regular_expression = pattern_to_regular_expression(
                    pattern.character_pattern(pattern_matrix=pattern.constituent_pattern)
                )
                
                background_frequency = 0
                for sequence in self.background.background_sequences:
                    background_matches = re.findall(pattern_regular_expression, sequence)
                    background_frequency += len(background_matches)
                background_percentage_frequency = (
                    background_frequency / len(self.background.background_df)
                )
                try:
                    enrichment = (
                        sample_percentage_frequency / background_percentage_frequency
                    )
                except:
                    enrichment = 0.0
            else:
                enrichment = np.nan    
            enrichments.append(enrichment)

        self.pattern_summary_table = pd.DataFrame(
            {
                'Pattern': pattern_labels,
                'Sample Frequency': sample_frequencies,
                'Foreground Frequency': foreground_frequencies,
                'Foreground Size': foreground_sizes,
                'Enrichment (Sample Frequency / Background Frequency)': enrichments  
            }
        )
        self.pattern_summary_table.to_csv(
            os.path.join(pattern_directory, f'{title}_pattern_summary_table.csv'),
            sep=',',
            header=True,
            index=True
        )