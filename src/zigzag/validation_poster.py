from zigzag.validation import *

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1e6)

def vary_pattern_proportion(fixed_pattern_size=20,
                            fold_differences=(1, 10, 100, 1000),
                            num_iterations=100):
    '''

    Parameters:
        fixed_pattern_size --
        fold_differences --
        num_iterations --
    Returns:
        f1_scores --
    '''

    # Experiment directory set-up
    experiment_dir = '../results/validation/vary_pattern_proportion_1000x'

    positional_residues = []
    expected_patterns = []

    # Generate patterns.
    for i in range(2):
        while True:
            positional_residue = (np.random.choice(7),
                                    np.random.choice(list(background_frequencies.keys()))
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
            
            # output directory set-up
            output_dir = experiment_dir \
                        + '/fold_difference_{}/{}'.format(fold_difference,
                                                            i)
            os.makedirs(output_dir)

            # Initialize sequence list. 
            sequences = []

            # Generate fixed-size pattern sequences
            for i in range(fixed_pattern_size):
                positional_residue = positional_residues[0]
                sequence = pd.Series([np.random.choice(SINGLE_LETTER_CODES,
                                        p=[(background_frequencies[k] / 100) \
                                            if k in background_frequencies \
                                            else 0.0 for k in SINGLE_LETTER_CODES]) \
                                        for j in range(8)])
                sequence[positional_residue[0]] = positional_residue[1]
                sequences.append(sequence)

            # Generate variable-size pattern sequences.
            for i in range(fixed_pattern_size * fold_difference):
                positional_residue = positional_residues[1]
                sequence = pd.Series([np.random.choice(SINGLE_LETTER_CODES,
                                p=[(background_frequencies[k]/100) \
                                    if k in background_frequencies \
                                    else 0.0 for k in SINGLE_LETTER_CODES]) \
                                for j in range(8)])
                sequence[positional_residue[0]] = positional_residue[1]
                sequences.append(sequence)

            
            # Generate sequence data frame and save to text file.
            df = pd.DataFrame(sequences)
            sample_filename = output_dir + '/sample.txt'

            with open(sample_filename,'w') as fout:
                for i,sequence in df.iterrows():
                    fout.write(''.join(sequence.tolist()) + '\n')

            # Assign metadata for replicate.
            sample_metadata = [fold_difference, i]
            
            results = calculate_results(sample_filename,
                                        vectorized_background_frequencies,
                                        expected_patterns,
                                        output_dir
                                        )
            
            results_row = sample_metadata + results
            f1_scores.append(results_row)
    
    # F1 scores output
    columns = ['fold_difference',
                'sample',
                'tp',
                'fp',
                'fn',
                'precision',
                'recall',
                'f1_score']
    f1_scores = pd.DataFrame(f1_scores, columns=columns)
    f1_scores.to_csv(experiment_dir + '/f1_scores.csv',
                        index=False,
                        header=True)

    return f1_scores

def minimum_relative_frequency(sample_sizes=[10, 100, 1000, 10000],
                                percent_frequency_range=(5, 100),
                                percent_step_size=5,
                                num_iterations=100,
                                positions='random'):
    """

    Parameters:
        sample_sizes -- List-like.
        percent_frequency_range -- List-like.
        percent_step_size --  Int.
        num_iterations -- Int.
        positions -- String.

    Returns:
        f1_scores -- Pandas DataFrame.
    """

    # Set up directory for experiment.
    experiment_dir = '../results/validation/min_relative_frequency'
    
    # Initialize list for all F1 scores.
    f1_scores = []
    
    for sample_size in sample_sizes:

        # for each percent frequency condition (5%-25%, 5% increments)
        for percent_frequency in range(percent_frequency_range[0],
                                        percent_frequency_range[1]
                                        + percent_step_size,
                                        percent_step_size):

            # calculate pattern size based on percent frequency for each sample size
            pattern_size = round(sample_size * (percent_frequency / 100))

            for i in range(num_iterations):
                
                # Output configuration.
                output_dir = experiment_dir \
                            + ('/sample_size_{}'
                                '/percent_frequency_{}/{}').format(sample_size,
                                                                    percent_frequency,
                                                                    i)
                os.makedirs(output_dir)
                sample_filename = output_dir + '/sample.txt'

                
                sample_metadata = [sample_size,
                                    percent_frequency,
                                    i]

                # Generate sample and run pattern extraction.
                expected_patterns = generate_random_sample(sample_filename,
                                                            sample_size,
                                                            pattern_size,
                                                            1,
                                                            positions
                                                            )

                results = calculate_results(sample_filename,
                                            vectorized_background_frequencies,
                                            expected_patterns,
                                            output_dir
                                            )

                results_row = sample_metadata + results                
                
                f1_scores.append(results_row)

    # Generate data frame for F1 scores and save to CSV file.
    columns = ['sample_size',
                'percent_frequency',
                'sample',
                'tp',
                'fp',
                'fn',
                'precision',
                'recall',
                'f1_score']
    f1_scores = pd.DataFrame(f1_scores, columns=columns)
    f1_scores.to_csv(experiment_dir + '/f1_scores.csv',
                        index=False,
                        header=True
                        )

    return f1_scores