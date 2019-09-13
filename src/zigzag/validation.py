import os

import numpy as np
import pandas as pd

from zigzag.config import *
from zigzag.preprocessing import *
from zigzag.pattern_extraction import *
from zigzag.merops_connector import *
from zigzag import sequences

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1000000)


# Set up context and background frequencies.
swissprot = sequences.import_fasta('../data/uniprot/uniprot.fasta')
compound_residue_background = sequences.Background(swissprot['sequence'].tolist())
no_compound_residue_background = sequences.Background(
    swissprot['sequence'].tolist(),
    compound_residues=None
)
background_frequencies = no_compound_residue_background.background_dict
print(background_frequencies)

# Load residues and frequencies.
alphabet = [
    residue for residue in no_compound_residue_background.alphabet
    if residue != 'U'
]
background_frequencies_no_u = {
    residue: frequency for residue, frequency in background_frequencies.items()
    if residue != 'U'
}
background_frequency_no_u_rounding_adjustment = np.divide(
    (100 - np.sum(list(background_frequencies_no_u.values())))
    , len(background_frequencies_no_u)
)
for residue in background_frequencies_no_u.keys():
    background_frequencies_no_u[residue] += background_frequency_no_u_rounding_adjustment

background_vector = np.array(
    [
        fq for i, fq in enumerate(no_compound_residue_background.background_vector)
        if i != no_compound_residue_background.alphabet.index('U')
    ],
    dtype=no_compound_residue_background.background_vector.dtype
)

# Identify most/least abundant amino acids in background alphabet.
most_abundant_residue = alphabet[np.argmax(background_vector)]
least_abundant_residue = alphabet[np.argmin(background_vector)]

# Initialize expected patterns.
most_abundant_residue_pattern = [pd.Series(['.' for i in range(8)])]
most_abundant_residue_pattern[0][3] = '[' + most_abundant_residue + ']'
least_abundant_residue_pattern = [pd.Series(['.' for i in range(8)])]
least_abundant_residue_pattern[0][3] = '[' + least_abundant_residue + ']'

background_frequencies_no_l = {
    residue:((frequency / (100 - background_frequencies_no_u['L'])) * 100)
    for residue, frequency in background_frequencies_no_u.items()
    if residue != 'L'
}
background_frequency_no_l_rounding_adjustment = np.divide(
    (100 - np.sum(list(background_frequencies_no_l.values())))
    , len(background_frequencies_no_l)
)
for residue in background_frequencies_no_l.keys():
    background_frequencies_no_l[residue] += background_frequency_no_l_rounding_adjustment
print(np.sum(list(background_frequencies_no_l.values())))

background_frequencies_no_w = {
    residue:((frequency / (100 - background_frequencies_no_u['W'])) * 100)
    for residue, frequency in background_frequencies_no_u.items()
    if residue != 'W'
}
background_frequency_no_w_rounding_adjustment = np.divide(
    (100 - np.sum(list(background_frequencies_no_w.values()))),
    len(background_frequencies_no_w)
)
for residue in background_frequencies_no_w.keys():
    background_frequencies_no_w[residue] += background_frequency_no_w_rounding_adjustment
background_frequencies_no_w['A'] -= 0.00000000000001
print(np.sum(list(background_frequencies_no_w.values())))

alphabet_no_l = [residue for residue in alphabet if residue != 'L']
alphabet_no_w = [residue for residue in alphabet if residue != 'W']

background_vector_no_l = np.array(
    [
        fq for i, fq in enumerate(background_vector)
        if i != alphabet.index('L')
    ],
    dtype=background_vector.dtype
)
background_vector_no_l = np.array(
    [
        fq for i, fq in enumerate(background_vector)
        if i != alphabet.index('W')
    ],
    dtype=background_vector.dtype
)


def generate_random_sample(output_filename,
                            sample_size,
                            group_size,
                            n_patterns,
                            positions='random'):

    '''

    Parameters:
        output_filename --
        background_frequencies --
        sample_size --
        group_size --
        n_patterns --
        positions --

    Returns:

    '''

    pattern_template = pd.Series(['.' for i in range(8)])

    positional_residues = []
    patterns = []

    sequences = []

    residues = [
                residue for residue in list(background_frequencies.keys())
                if residue != 'U'
    ]

    for p in range(n_patterns):
        # Generate random position / residue pair. If pair violates uniqueness criterion for test condition, try again.
        while True:
            
            if positions in ['random','different']:
                positional_residue = (np.random.choice(7),
                                    np.random.choice(residues))
            elif positions == 'same':
                positional_residue = (3, np.random.choice(residues))
            
            if ((positional_residue not in positional_residues)
                        and not (positions == 'different' and positional_residue
                                in [i[0] for i in positional_residues])):
                positional_residues.append(positional_residue)
                pattern = pd.Series(['.' for i in range(8)])
                pattern[positional_residue[0]] = '[' + positional_residue[1] + ']'
                patterns.append(pattern)
                break

        for i in range(group_size):
            sequence = pd.Series(
                [np.random.choice(
                    no_compound_residue_background.alphabet,
                    p=[(background_frequencies[k] / 100)
                        if k in background_frequencies
                        else 0.0 for k in no_compound_residue_background.alphabet]
                )
                for j in range(8)]
            )
            sequence[positional_residue[0]] = positional_residue[1]
            sequences.append(sequence)

    for i in range(sample_size-(group_size*n_patterns)):
        sequences.append(
            pd.Series(
                [np.random.choice(
                    no_compound_residue_background.alphabet,
                    p=[(background_frequencies[k] / 100)
                        if k in background_frequencies
                        else 0.0 for k in no_compound_residue_background.alphabet]
                )
                for j in range(8)]
            )
        )

    df = pd.DataFrame(sequences)

    with open(output_filename,'w') as fout:
        for i,sequence in df.iterrows():
            fout.write(''.join(sequence.tolist()) + '\n')

    return patterns

"""
def vary_absolute_frequency(initial_size=1,
                            step_size=1,
                            num_steps=20,
                            num_iterations=100,
                            configuration='removal'):
    '''
    Determine minimum absolute frequency for detection of enriched pattern.

    Parameters:
        initial_size --
        step_size --
        num_steps --
        num_iterations --
        configuration --

    Returns:
        f1_scores --
    '''

    experiment_dir = '../results/validation/min_absolute_frequency'

    # Initialize output list.
    f1_scores = []

    for i in range(initial_size, (num_steps*step_size)+initial_size, step_size):

        for j in range(num_iterations):
            # Configure output and sample metadata.
            output_dir = experiment_dir + '/n_{}/{}'.format(i, j)
            os.makedirs(output_dir)
            
            sample_filename = output_dir + '/sample.txt'
            sample_metadata = [i, j]
            
            # Generate sample.
            expected_patterns = generate_random_sample(sample_filename,
                                                        i,
                                                        i,
                                                        1)
            # Compute results.
            results = calculate_results(sample_filename,
                                        vectorized_background_frequencies,
                                        expected_patterns,
                                        output_dir
                                        )
            # Generate result row and append to output list.
            results_row = sample_metadata + results
            f1_scores.append(results_row)

    # Generate output data frame and save to CSV file.
    columns = ['num_sequences',
                'sample',
                'tp',
                'fp',
                'fn',
                'precision',
                'recall',
                'f1_score'
                ]
    f1_scores = pd.DataFrame(f1_scores, columns=columns)
    f1_scores.to_csv(experiment_dir + '/f1_scores.csv', index=False, header=True)

    return f1_scores
"""

def min_absolute_fixed_residues(initial_size=1,
                                step_size=1,
                                num_steps=20,
                                num_iterations=100,
                                set_reduction=True):
    """
    Test minimum absolute frequency on most abundant () and least
        abundant amino acids. Test should reflect that detection floor
        is lowest for the least abundant residue, highest for the most
        abundant residue, and consequently between the two for all others.

    Parameters:
        initial_size -- Int.
        step_size -- Int.
        num_steps -- Int.
        num_iterations -- Int.
        set_reduction -- Boolean.

    Returns:
        f1_scores -- Pandas DataFrame.
    """

    base_output_dir = '../results/validation/min_absolute_frequency'
    
    f1_scores = []
    # Iterate over number of sequences.
    for n in range(1, (num_steps + 1), step_size):
        # Iterate over replicates.
        for i in range(num_iterations):
            # Initialize output directories.
            experiment_dir = base_output_dir + '/n_{}/{}'.format(n, i)
            most_abundant_output_directory = experiment_dir + '/most_abundant'
            least_abundant_output_directory = experiment_dir + '/least_abundant'
            os.makedirs(most_abundant_output_directory)
            os.makedirs(least_abundant_output_directory)

            # Initialize samples.
            most_abundant_residue_sample = []
            least_abundant_residue_sample = []
            # Generate random sequences.
            for s in range(n):
                sequence = generate_random_sequence()
                # Append sequence to respective sample w/ substituted
                # fixed residue.
                most_abundant_sequence = sequence.copy()
                most_abundant_sequence[3] = most_abundant_residue
                most_abundant_residue_sample.append(most_abundant_sequence)
                least_abundant_sequence = sequence.copy()
                least_abundant_sequence[3] = least_abundant_residue
                least_abundant_residue_sample.append(least_abundant_sequence)

            # Write samples to files
            most_abundant_residue_sample_df = pd.DataFrame(most_abundant_residue_sample)
            most_abundant_output_filename = most_abundant_output_directory + '/sample.txt'
            with open(most_abundant_output_filename,'w') as fout:
                for k, sequence in most_abundant_residue_sample_df.iterrows():
                    fout.write(''.join(sequence.tolist()) + '\n')
            
            least_abundant_residue_sample_df = pd.DataFrame(least_abundant_residue_sample)
            least_abundant_output_filename = least_abundant_output_directory + '/sample.txt'
            with open(least_abundant_output_filename,'w') as fout:
                for k, sequence in least_abundant_residue_sample_df.iterrows():
                    fout.write(''.join(sequence.tolist()) + '\n')

            # Calculate least abundant residue results.
            least_abundant_residue_results = calculate_results(
                least_abundant_output_filename,
                compound_residue_background,
                least_abundant_residue_pattern,
                least_abundant_output_directory
            )
            sample_metadata = ['Least abundant residue: {}'.format(least_abundant_residue), n, i]
            results_row = sample_metadata + least_abundant_residue_results
            f1_scores.append(results_row)
            
            # Calculate most abundant residue results.
            most_abundant_residue_results = calculate_results(
                most_abundant_output_filename,
                compound_residue_background,
                most_abundant_residue_pattern,
                most_abundant_output_directory
            )
            sample_metadata[0] = 'Most abundant residue: {}'.format(most_abundant_residue)
            results_row = sample_metadata + most_abundant_residue_results
            f1_scores.append(results_row)
    
    # Save F1 scores to CSV file.
    columns = [
        'abundance',
        'num_sequences',
        'sample',
        'tp',
        'fp',
        'fn',
        'precision',
        'recall',
        'f1_score'
    ]
    f1_scores = pd.DataFrame(f1_scores, columns=columns)
    f1_scores.to_csv(
        base_output_dir + '/f1_scores.csv',
        index=False,
        header=True
    )

    return f1_scores


def min_pct_frequency_fixed_residues(sample_sizes=[10, 100, 1000, 10000],
                                        percent_frequency_range=(10, 100),
                                        percent_step_size=10,
                                        num_iterations=100,
                                        width=8,
                                        set_reduction=True):
    """


    Parameters:
        sample_sizes -- List-like.
        percent_frequency_range -- List-like.
        percent_step_size -- Int.
        num_iterations -- Int.
        set_reduction -- Boolean.

    Returns:
        f1_scores -- Pandas DataFrame.
    """
    base_output_dir = '../results/validation/min_percentage_frequency'

    f1_scores = []
    for sample_size in sample_sizes:
        for percent in range(
            percent_frequency_range[0],
            percent_frequency_range[1] + percent_step_size,
            percent_step_size
        ):
            num_sequences = int(sample_size * (percent / 100))
            for i in range(num_iterations):
                # Initialize experiment directories.
                experiment_dir = base_output_dir + '/sample_size_{}/{}_percent/{}'.format(
                    sample_size,
                    percent,
                    i
                )
                most_abundant_output_dir = experiment_dir + '/most_abundant'
                least_abundant_output_dir = experiment_dir + '/least_abundant'
                os.makedirs(most_abundant_output_dir)
                os.makedirs(least_abundant_output_dir)

                fixed_residue_count = 0
                sequence_strings = []
                most_abundant_residue_sample = []
                least_abundant_residue_sample = []
                most_abundant_count = 0
                least_abundant_count = 0
                for j in range(sample_size):
                    while True:
                        sequence = generate_random_sequence()
                        sequence_string = ''.join(
                            [p for q, p in enumerate(''.join(sequence)) if q != 3]
                        )
                        if sequence_string not in sequence_strings:
                            most_abundant_residue_sample.append(sequence.copy())
                            least_abundant_residue_sample.append(sequence.copy())
                            if sequence[3] == most_abundant_residue:
                                most_abundant_count += 1
                            elif sequence[3] == least_abundant_residue:
                                least_abundant_count += 1
                            sequence_strings.append(sequence_string)
                            break
                
                for sequence in most_abundant_residue_sample:
                    if most_abundant_count < num_sequences:
                        if sequence[3] != most_abundant_residue:
                            sequence[3] = most_abundant_residue
                            most_abundant_count += 1
                    else:
                        break
                for sequence in least_abundant_residue_sample:
                    if least_abundant_count < num_sequences:
                        if sequence[3] != least_abundant_residue:
                            sequence[3] = least_abundant_residue
                            least_abundant_count += 1
                    else:
                        break

                # Generate sample data frames and save to text files.
                most_abundant_sample_df = pd.DataFrame(most_abundant_residue_sample)
                most_abundant_sample_filename = most_abundant_output_dir + '/sample.txt'
                with open(most_abundant_sample_filename, 'w') as fout:
                    for i,sequence in most_abundant_sample_df.iterrows():
                        fout.write(''.join(sequence.tolist()) + '\n')

                least_abundant_sample_df = pd.DataFrame(least_abundant_residue_sample)
                least_abundant_sample_filename = least_abundant_output_dir + '/sample.txt'
                with open(least_abundant_sample_filename, 'w') as fout:
                    for i,sequence in least_abundant_sample_df.iterrows():
                        fout.write(''.join(sequence.tolist()) + '\n')

                # Calculate least abundant residue results.
                least_abundant_residue_results = calculate_results(
                    least_abundant_sample_filename,
                    compound_residue_background,
                    least_abundant_residue_pattern,
                    least_abundant_output_dir
                )
                sample_metadata = [
                    'Least abundant residue: {}'.format(least_abundant_residue),
                    sample_size,
                    percent,
                    i
                ]
                results_row = sample_metadata + least_abundant_residue_results
                f1_scores.append(results_row)
                
                # Calculate most abundant residue results.
                most_abundant_residue_results = calculate_results(
                    most_abundant_sample_filename,
                    compound_residue_background,
                    most_abundant_residue_pattern,
                    most_abundant_output_dir
                )
                sample_metadata[0] = 'Most abundant residue: {}'.format(most_abundant_residue)
                results_row = sample_metadata + most_abundant_residue_results
                f1_scores.append(results_row)

    # Save F1 scores to CSV file.
    columns = [
        'abundance',
        'sample_size',
        'percentage',
        'sample',
        'tp',
        'fp',
        'fn',
        'precision',
        'recall',
        'f1_score'
    ]
    f1_scores = pd.DataFrame(f1_scores, columns=columns)
    f1_scores.to_csv(
        base_output_dir + '/f1_scores.csv',
        index=False,
        header=True
    )

    return f1_scores


def vary_n_patterns_constant_frequency(sample_size=1000,
                                        pattern_size=50,
                                        num_iterations=100,
                                        positions=['random', 'same']):
    '''
    Hold sample size and pattern size constant. Vary number of patterns.
        If sum frequency of pattern exemplar sequences is less than sample
        size, fill remainder of sample with random sequences drawn according
        to amino acid background frequencies.

    Parameters:
        sample_size -- Int.
        pattern_size -- Int.
        num_iterations -- Int.
        positions -- String.

    Returns:
        f1_scores -- Pandas Series.
    '''
    base_output_dir = '../results/validation/constant_pattern_size'
    
    # Initialize F1 scores list.
    f1_scores = []
    for p in positions:
        for n in range(1, (sample_size // pattern_size) + 1):
            for i in range(num_iterations):
                # Set up output directory.
                output_dir = (
                    base_output_dir
                    + '/{}_positions/sample_size_{}/pattern_size_{}/{}_patterns/{}'.format(
                        p,
                        sample_size,
                        pattern_size,
                        n,
                        i
                    )
                )
                os.makedirs(output_dir)
                sample_filename = output_dir + '/sample.txt'
                sample_metadata = [sample_size, p, n, i]

                # Generate pattern set and sample.
                expected_patterns = generate_random_sample(
                    sample_filename,
                    sample_size,
                    pattern_size,
                    n,
                    p
                )

                results = calculate_results(
                    sample_filename,
                    compound_residue_background,
                    expected_patterns,
                    output_dir
                )

                results_row = sample_metadata + results
                f1_scores.append(results_row)
    
    # Save F1 scores to CSV file.
    columns = [
        'sample_size',
        'positions',
        'num_patterns',
        'sample',
        'tp',
        'fp',
        'fn',
        'precision',
        'recall',
        'f1_score'
    ]
    f1_scores = pd.DataFrame(f1_scores,columns=columns)
    f1_scores.to_csv(
        base_output_dir + '/f1_scores.csv',
        index=False,
        header=True
    )

    return f1_scores


def vary_fold_difference(fixed_pattern_size=20,
                            fold_differences=(1, 10, 100, 1000),
                            num_iterations=100,
                            set_reduction=[True, False]):
    """

    Parameters:
        background -- Background instance.
        fixed_pattern_size --
        fold_differences --
        num_iterations --
    
    Returns:
        f1_scores --
    """

    # Experiment directory set-up
    base_output_dir = '../results/validation/vary_fold_difference'

    positional_residues = []
    expected_patterns = []

    # Generate patterns.
    for i in range(2):
        while True:
            positional_residue = (
                np.random.choice(7),
                np.random.choice(
                    [residue for residue in list(background_frequencies.keys())
                    if residue != 'U']
                )
            )
            if (positional_residue not in positional_residues):
                positional_residues.append(positional_residue)
                pattern = pd.Series(['.' for i in range(8)])
                pattern[positional_residue[0]] = '[' + positional_residue[1] + ']'
                expected_patterns.append(pattern)
                break

    # Initialize list for F1 scores.
    f1_scores = []
    for fold_difference in fold_differences:
        for i in range(num_iterations):
            # Set up experiment directory.
            experiment_dir = base_output_dir + '/fold_difference_{}/{}'.format(
                fold_difference,
                i
            )
            os.makedirs(experiment_dir)

            # Initialize sequence list. 
            sequences = []
            sequence_strings = []
            # Generate fixed-size pattern sequences
            for n in range(fixed_pattern_size):
                positional_residue = positional_residues[0]
                while True:
                    sequence = generate_random_sequence()
                    sequence[positional_residue[0]] = positional_residue[1]
                    sequence_string = ''.join(sequence)
                    if sequence_string not in sequence_strings:
                        sequences.append(sequence)
                        break

            # Generate variable-size pattern sequences.
            for n in range(fixed_pattern_size * fold_difference):
                positional_residue = positional_residues[1]
                while True:
                    sequence = generate_random_sequence()
                    sequence[positional_residue[0]] = positional_residue[1]
                    sequence_string = ''.join(sequence)
                    if sequence_string not in sequence_strings:
                        sequences.append(sequence)
                        break

            
            # Generate sequence data frame and save to text file.
            df = pd.DataFrame(sequences)
            sample_filename = experiment_dir + '/sample.txt'

            with open(sample_filename, 'w') as fout:
                for s, sequence in df.iterrows():
                    fout.write(''.join(sequence.tolist()) + '\n')

            for s in set_reduction:
                # output directory set-up
                output_dir = (
                    experiment_dir
                    + '/set_reduction_{}'.format('enabled' if s else 'disabled')
                )
                os.makedirs(output_dir)

                # Assign metadata for replicate.
                sample_metadata = [
                    'Set reduction enabled' if s else 'Set reduction disabled',
                    fold_difference,
                    i
                ]
                
                results = calculate_results(
                    sample_filename,
                    compound_residue_background,
                    expected_patterns,
                    output_dir,
                    alpha=0.001,
                    set_reduction=s
                )
                
                results_row = sample_metadata + results
                f1_scores.append(results_row)
    
    # F1 scores output
    columns = [
        'set_reduction',
        'fold_difference',
        'sample',
        'tp',
        'fp',
        'fn',
        'precision',
        'recall',
        'f1_score',
    ]
    f1_scores = pd.DataFrame(f1_scores, columns=columns)
    f1_scores.to_csv(
        experiment_dir + '/f1_scores.csv',
        index=False,
        header=True
    )

    return f1_scores


def false_positive_analysis(sample_sizes=[10, 100, 1000, 10000],
                            num_iterations=100,
                            width=8,
                            alpha=0.001,
                            set_reduction=True):
    """

    Parameters:
        sample_sizes -- List-like.
        num_iterations -- Int.

    Returns:
        f1_scores -- Pandas DataFrame.
    """

    # Initialize output directory.
    base_output_dir = '../results/validation/false_positive_analysis'
    
    fp = []
    for sample_size in sample_sizes:
        for i in range(num_iterations):
            # Initialize output directories.
            experiment_dir = base_output_dir + '/sample_size_{}/{}'.format(
                sample_size,
                i
            )
            swissprot_output_dir = experiment_dir + '/swissprot'
            in_silico_output_dir = experiment_dir + '/in_silico'
            os.makedirs(swissprot_output_dir)
            os.makedirs(in_silico_output_dir)

            # Generate and save Swiss-Prot sample.
            swissprot_output_filename = swissprot_output_dir + '/sample.txt'
            generate_random_swissprot_sample(
                sample_size,
                swissprot_output_filename,
                width=width
            )
            # Run pattern extraction.
            for c in range(2):
                if c == 0:
                    pattern_output_dir = (
                        swissprot_output_dir
                        + '/compound_resides_enabled'
                    )
                else:
                    pattern_output_dir = (
                        swissprot_output_dir
                        + '/compound_resides_disabled'
                    )
                os.makedirs(pattern_output_dir)

                patterns = run_pattern_extraction(
                    swissprot_output_filename,
                    compound_residue_background if c == 0 else no_compound_residue_background,
                    alpha=alpha,
                    set_reduction=set_reduction
                )
                # Save patterns.
                process_patterns(pattern_output_dir, 'output', patterns)
                # Append results.
                fp.append(
                    [
                        'Swiss-Prot',
                        ('Compound residues enabled' if c == 0
                            else 'Compound residues disabled'),
                        sample_size,
                        i,
                        len(patterns)
                    ]
                )

            # Generate and save in silico sample.
            in_silico_output_filename = in_silico_output_dir + '/sample.txt'
            generate_in_silico_sample(
                sample_size,
                in_silico_output_filename,
                width=8
            )
            # Run pattern extraction.
            for c in range(2):
                if c == 0:
                    pattern_output_dir = (
                        in_silico_output_dir
                        + '/compound_resides_enabled'
                    )
                else:
                    pattern_output_dir = (
                        in_silico_output_dir
                        + '/compound_resides_disabled'
                    )
                os.makedirs(pattern_output_dir)

                patterns = run_pattern_extraction(
                    in_silico_output_filename,
                    (compound_residue_background if c == 0
                        else no_compound_residue_background),
                    alpha=alpha,
                    set_reduction=set_reduction
                )
                # Save patterns.
                process_patterns(pattern_output_dir, 'output', patterns)
                # Append results.
                fp.append(
                    [
                        'In silico',
                        ('Compound residues enabled' if c == 0
                            else 'Compound residues disabled'),
                        sample_size,
                        i,
                        len(patterns)
                    ]
                )
    
    columns = [
        'sequence_type',
        'compound_residues',
        'sample_size',
        'sample',
        'fp',
    ]

    # FP output
    fp_df = pd.DataFrame(fp, columns=columns)
    fp_df.to_csv(
        base_output_dir + '/f1_scores.csv',
        index=False,
        header=True
    )

    return fp_df


def generate_random_sequence(width=8):
    """
    Generates random sequence using background frequencies.

    Parameters:
        width -- Int.

    Returns:
        Pandas Series object.
    """
    return pd.Series([generate_random_amino_acid() for j in range(width)])


def generate_random_amino_acid(alphabet=alphabet,
                                background_frequencies=background_frequencies_no_u):
    """Generates random amino acid using background frequencies."""
    return np.random.choice(
        alphabet,
        p=[((background_frequencies[k]) / 100) 
        if k in background_frequencies
        else 0.0 for k in alphabet]
    )

def generate_random_min_pct_residue(abundance):
    """

    Parameters:
        abundance -- String.

    Returns:
        residue -- String.
    """
    if abundance == 'most_abundant':
        residue = generate_random_amino_acid(
            alphabet=alphabet_no_l,
            background_frequencies=background_frequencies_no_l
        )
    elif abundance == 'least_abundant':
        residue = generate_random_amino_acid(
            alphabet=alphabet_no_w,
            background_frequencies=background_frequencies_no_w
        )

    return residue


def run_pattern_extraction(sample_path,
                            background,
                            alpha=0.001,
                            minimum_occurrences=2,
                            multiple_testing_correction=True,
                            positional_weighting=True,
                            mod_type='merops',
                            set_reduction=True):
    '''
    Runs pattern extraction algorithm on specified sample of sequences.
        Returns detected patterns.

    Parameters:
        sample -- String. Path to sample file in text form.
        background -- NumPy Array.
        mod_type -- String.
    Returns:
        detected_patterns -- List.
    '''

    # Select residue- or bond-centered sequences.
    if mod_type == 'phospho':
        center = False
    else:
        center = True

    # Load sample from file.
    sample = sequences.load_prealigned_file(sample_path, background, center=center)

    # Run pattern extraction.
    pattern_container = PatternContainer(
        sample,
        background,
        title=None,
        output_directory=None,
        p_value_cutoff=alpha,
        minimum_occurrences=minimum_occurrences,
        multiple_testing_correction=multiple_testing_correction,
        positional_weighting=positional_weighting,
        set_reduction=set_reduction
    )
    # Cleanup patterns and return list of character patterns.
    pattern_container.prune_patterns()
    detected_patterns = [
        pattern.character_pattern() for pattern in pattern_container.pattern_list
    ]

    return detected_patterns


def process_patterns(output_dir, pattern_type, patterns):
    '''
    Generate set of unique patterns. Write patterns to CSV file. Return
        set containing patterns.

    Parameters:
        output_dir -- String. Path to output directory for num_iteration
                                    of experiment.
        pattern_type -- String. "input" or "output".
        patterns -- List.

    Returns:
        pattern_set -- Set.
    '''

    output_filename = output_dir + '/{}_patterns.csv'. \
                        format('input' if pattern_type == 'input' else 'output')
    
    pattern_set = {''.join(pattern.tolist()) for pattern in patterns}
    pd.Series(list(pattern_set)).to_csv(output_filename, index=False)

    return pattern_set


def calculate_outcomes(expected_pattern_set,
                        detected_pattern_set):
    '''

    Parameters:
        output_dir -- String.
        expected_pattern_set -- Set.
        detected_pattern_set -- Set.

    Returns:
        outcomes -- Pandas Series.
    '''

    true_positives = len(expected_pattern_set & detected_pattern_set)
    false_positives = len(detected_pattern_set - expected_pattern_set)
    false_negatives = len(expected_pattern_set - detected_pattern_set)
    outcomes = pd.Series([true_positives, false_positives, false_negatives],
                                    index = ['TP','FP','FN'])

    return outcomes


def calculate_precision(confusion_values):
    '''

    Parameters:
        confusion_values -- Pandas Series.

    Returns:
        precision -- Int.
    '''

    # Update cumulative precision and recall lists.
    try:
        precision = confusion_values['TP'] / (confusion_values['TP']
                                                + confusion_values['FP']
                                                )
    except ZeroDivisionError:
        precision = 0
    except:
        print('Unexpected error')

    return precision


def calculate_recall(confusion_values):
    """

    Parameters:
        confusion_values -- Pandas Series.

    Returns:
        recall -- Int.
    """

    try:
        recall = confusion_values['TP'] / (confusion_values['TP']
                                            + confusion_values['FN']
                                            )
    except ZeroDivisionError:
        recall = 0
    except:
        print('Unexpected error')

    return recall


def f1_score(precision, recall):
    """
    Calculates F1 score from precision and recall.

    Parameters:
        precision -- Float.
        recall -- Float.

    Returns:
        f1_score -- Float.
    """
    
    try:
        f1_score = 2 * ((np.mean(precision) * np.mean(recall)) 
                        / (np.mean(precision) + np.mean(recall)
                        ))
    except ZeroDivisionError:
        f1_score = 0
    except:
        print('Something went wrong during this iteration.')
    else:
        return f1_score


def calculate_results(sample_filename,
                        background,
                        expected_patterns,
                        output_dir,
                        removal='removal',
                        alpha=0.001,
                        minimum_occurrences=2,
                        multiple_testing_correction=True,
                        positional_weighting=True,
                        set_reduction=True,
                        mod_type='merops'):
    """

    Parameters:
        sample_filename --
        expected_patterns --

    Returns:
        results --
    """

    
    detected_patterns = run_pattern_extraction(
        sample_filename,
        background,
        alpha=alpha,
        minimum_occurrences=minimum_occurrences,
        multiple_testing_correction=multiple_testing_correction,
        positional_weighting=positional_weighting,
        set_reduction=set_reduction
    )
    
    # Generate and save pattern sets for iteration.
    expected_pattern_set = process_patterns(
        output_dir,
        'input',
        expected_patterns
    )
    print(expected_pattern_set)
    detected_pattern_set = process_patterns(
        output_dir,
        'ouput',
        detected_patterns
    )
    print(detected_pattern_set)

    # Calculate outcomes.
    outcomes = calculate_outcomes(
        expected_pattern_set,
        detected_pattern_set
    )

    # Calculate precision and recall.
    precision = calculate_precision(outcomes)
    recall = calculate_recall(outcomes)

    results = [*outcomes, precision, recall, f1_score(precision, recall)]

    return results


def load_pattern_set(path):
    '''
    Loads patterns from CSV file containing one pattern per line to set of
        unique patterns.

    Parameters:
        path -- String. Path to CSV-formatted file containing pattern set.

    Returns:
        pattern_set -- Set. Set containing patterns from file at path.
    '''

    patterns = pd.read_csv(path, header=None)
    pattern_set = {''.join(pattern.tolist()).upper() for i, pattern in patterns.iterrows()}

    return pattern_set


def pease_test_cases(data):
    '''

    Parameters:
        data --

    Returns:
        None
    '''

    # Load MoMo w/o Harvard patterns.
    momo_patterns = load_pattern_set('../data/Pease/momo_wo_harvard.csv')
    momo_pattern_set = set(momo_patterns)
    print(momo_pattern_set)

    # Load MoMo w/ Harvard patterns.
    motif_x_patterns = load_pattern_set('../data/Pease/momo_w_harvard.csv')
    motif_x_pattern_set = set(motif_x_patterns)
    print(motif_x_pattern_set)

    # Run pattern extraction.
    detected_patterns = run_pattern_extraction('../data/Pease/pease_2.txt')
    detected_pattern_set = set(detected_patterns)
    print(detected_pattern_set)

    # Compute intersections of pattern sets.
    # Generate Venn diagram of patterns.


    return None


# Function to generate dummy samples from Swiss-Prot.
def generate_random_swissprot_sample(sample_size, output_filename, width=8):
    """Generates sample of n randomly selected sequences."""
    sample = []
    for i in range(sample_size):
        while True:
            while True:
                protein = swissprot.sample(n=1)['sequence'].to_string(index=False).strip()
                if len(protein) >= width:
                    break
            start = np.random.randint(len(protein) - width + 1)
            stop = start + width
            sequence = protein[start:stop]
            if sequence not in sample:
                sample.append(sequence)
                break
    
    df = pd.DataFrame(sample)
    with open(output_filename,'w') as fout:
        for i, sequence in df.iterrows():
            fout.write(''.join(sequence.tolist()) + '\n')


def generate_random_plasmodium_sample(sample_size, output_filename, width=9):
    """Generates n randomly selected Plasmodium falciparum sequences."""
    plasmodium_context = sequences.import_fasta('../data/Pease/Plasmodium_falciparum.EPr1.pep.all.fa')
    plasmodium_background = sequences.Background(
        plasmodium_context['sequence'].tolist(),
        compound_residues=None
    )

    sample = []
    for i in range(sample_size):
        while True:
            while True:
                protein = plasmodium_context.sample(n=1)['sequence'].to_string(index=False).strip()
                if len(protein) >= width:
                    break
            start = np.random.randint(len(protein) - width + 1)
            stop = start + width
            sequence = protein[start:stop]
            if (sequence not in sample) and (sequence[6] == 'Y'):
                sample.append(sequence)
                break
        
    df = pd.DataFrame(sample)

    print(df.columns)
    with open(output_filename,'w') as fout:
        for i, sequence in df.iterrows():
            fout.write(''.join(sequence.tolist()) + '\n')


def generate_in_silico_sample(sample_size, output_filename, width=8):
    """Generates random sample of specified width."""
    sample = []
    sequence_strings = []
    for i in range(sample_size):
        while True:
            sequence = generate_random_sequence(width=width)
            sequence_string = ''.join(sequence)
            if sequence_string not in sequence_strings:
                sample.append(sequence)
                sequence_strings.append(sequence_string)
                break

    df = pd.DataFrame(sample)
    with open(output_filename,'w') as fout:
        for i, sequence in df.iterrows():
            fout.write(''.join(sequence.tolist()) + '\n')


def calculate_min_pct_results_from_patterns():
    """Calculate min pct F1 scores and write to file using output files"""

    min_pct_dir = '../results/validation/min_percentage_frequency'

    f1_scores = []
    for root, dirs, files in os.walk(min_pct_dir):
        if 'input_patterns.csv' not in files:
            continue

        # Load input patterns.
        expected_patterns = pd.read_csv(os.path.join(root, files[0]), sep=',', header=None, index_col=False)
        expected_pattern_set = {pattern[0] for pattern in expected_patterns.iterrows()}

        # Load output patterns.
        try:
            detected_patterns = pd.read_csv(os.path.join(root, files[1]), sep=',', header=None, index_col=False)
        except:
            results = [0, 0, 1, 0, 0, 0]
        else:
            detected_pattern_set = {pattern[0] for pattern in detected_patterns.iterrows()}

            # Calculate outcomes.
            outcomes = calculate_outcomes(
                expected_pattern_set,
                detected_pattern_set
            )

            # Calculate precision and recall.
            precision = calculate_precision(outcomes)
            recall = calculate_recall(outcomes)
            results = [*outcomes, precision, recall, f1_score(precision, recall)]

        # Get sample metadata.
        root_components = root.split('/')
        if root_components[7] == 'most_abundant':
            abundance = 'Most abundant residue: L'
        else:
            abundance = 'Least abundant residue: W'
        sample_size = root_components[4][(root_components[4].rfind('_') + 1):]
        percentage = root_components[5][:root_components[5].find('_')]
        sample = root_components[6]
        sample_metadata = [abundance, sample_size, percentage, sample]

        # Store results row.
        f1_scores.append(sample_metadata + results)

    # Save F1 scores to CSV file.
    columns = [
        'abundance',
        'sample_size',
        'percentage',
        'sample',
        'tp',
        'fp',
        'fn',
        'precision',
        'recall',
        'f1_score'
    ]
    f1_scores = pd.DataFrame(f1_scores, columns=columns)
    f1_scores.to_csv(
        min_pct_dir + '/f1_scores.csv',
        index=False,
        header=True
    )


def generate_plasmodium_blank():
    """Generate set of random Plasmodium falciparum sequences."""
    pass
