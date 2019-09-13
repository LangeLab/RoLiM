import itertools
import io
import uuid

from django.conf import settings

import numpy as np
import pandas as pd
import scipy.stats as st
from scipy import special
import weblogo as w

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 100000)


class Pattern:
    '''

    Attributes:
        subset_tensor -- 3D NumPy Array. Dtype = int8.
        pattern_matrix -- 2D NumPy Array. Dtype = int8.
        removed_positional_residues -- 2D NumPy Array. Dtype = int8.

    Methods:

    '''

    def __init__(self,
                    pattern_container,
                    subset_tensor,
                    initial_pattern,
                    current_depth,
                    max_depth,                      
                    p_value_cutoff=0.01,
                    fold_change_cutoff=1.0,
                    occurrences=2,
                    positional_weighting=True,
                    multiple_testing_correction=True):
        '''

        Parameters:
            pattern_container -- PatternContainer instance.
            subset_tensor -- 3D NumPy Array. Dtype = int8.
            pattern_matrix -- 2D NumPy Array. Dtype = int8.
            current_depth -- Int.
            max_depth -- Int.
            p_value_cutoff -- Float.
            occurrences -- Int.
            fold_change_cutoff -- Float.
            positional_weighting -- Bool.
            multiple_testing_correction -- Bool.

        Returns:
            None
        '''

        # Initialize pattern attributes.
        self.pattern_id = str(uuid.uuid4())
        self.background = pattern_container.background
        self.subset_tensor = subset_tensor
        self.pattern_matrix = initial_pattern.copy()
        
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
        ).astype('int8'))

        # Iteratively extract significantly enriched positional residues.
        if current_depth < max_depth:
            while self.subset_tensor.shape[0] > 1:
                # Calculate positional background frequency adjustments.
                positional_background_frequencies = (
                    self.adjust_positional_background_frequencies()
                )
                positional_background_frequencies[
                    positional_background_frequencies > 1.0] = 1.0

                # Calculate positional residue frequencies from subset_tensor.
                positional_residue_frequencies = self.calculate_frequencies()

                # Calculate positional residue fold change.
                positional_residue_fold_change = (
                    positional_residue_frequencies
                    / (positional_background_frequencies
                        * self.subset_tensor.shape[0])
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
                    dtype='float64'
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
                    if bonferroni_m == 0:
                        bonferroni_m = 1
                    binomial_p_values[binomial_p_values
                                        > (p_value_cutoff/bonferroni_m)] = 1.0
                
                # Fold change and minimum occurrence enforcement.
                binomial_p_values[positional_residue_fold_change <= fold_change_cutoff] = 1.0
                binomial_p_values[positional_residue_frequencies < occurrences] = 1.0

                # Exit loop if no significantly enriched positional
                # residues detected. Instantiate new search on reduced
                # set if ehaustive search is enabled.
                if (binomial_p_values[(binomial_p_values > 0)
                    * (binomial_p_values < p_value_cutoff)].size == 0):
                    break

                # Calculate positional weights.
                positional_weights = self.calculate_positional_weights(
                    positional_residue_frequencies
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
                extracted_subset_tensor_hashes = []
                for enriched_positional_residue in np.transpose(
                        np.nonzero(max_positional_residue_enrichments)):
                    # Set enriched residue coordinates.
                    residue_coordinates = tuple(enriched_positional_residue)
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

                    # Prevent creation of identical branches.
                    matching_sequences_hash = hash(matching_sequences.tostring())
                    if matching_sequences_hash in extracted_subset_tensor_hashes:
                        continue
                    else:
                        extracted_subset_tensor_hashes.append(matching_sequences_hash)

                    # Add new positional residue to removed positional
                    # residue array.
                    self.removed_positional_residues[residue_coordinates] = 1

                    # Add simple residues to removed positional residue
                    # array if enriched compound residue extracted.
                    if (residue_coordinates[1]
                            >= len(self.background.alphabet)):
                        self._remove_compound_residue_constituents(residue_coordinates)

                    # Add parent pattern residues to new pattern.
                    new_pattern = np.logical_or(self.pattern_matrix,
                                                new_pattern).astype('int8')

                    # Instantiate new Pattern object with new pattern,
                    # matching sequences, and current depth.
                    pattern_container.add_new_pattern(Pattern(pattern_container,
                                                            matching_sequences,
                                                            new_pattern,
                                                            current_depth + 1,
                                                            max_depth,
                                                            p_value_cutoff))

                # Remove sequences matching new pattern(s) from subset.
                self._remove_sequences(max_positional_residue_enrichments)

    def adjust_positional_background_frequencies(self):
        '''

        Parameters:
            None

        Returns:
            vectorized_background_frequency_adjustment_values --
                                1D NumPy Array. Dtype = float64.
        '''

        # Calculate cumulative extracted background percentage.
        positional_penalties = np.sum(
            (self.removed_positional_residues[:, :len(self.background.alphabet)]
                * self.background.background_vector[:len(self.background.alphabet)]),
            axis=1
        )

        # Calculate adjusted background frequencies.
        positional_background_frequencies = np.outer(
            np.reciprocal(100 - positional_penalties),
            self.background.background_vector[:len(self.background.alphabet)]
        )

        # Mask removed positional simple residues.
        positional_background_frequencies = (
            positional_background_frequencies
            * np.invert(
                self.removed_positional_residues[:, :len(self.background.alphabet)].astype(
                    np.bool))
        )

        # Broadcast positional simple residues to comound residues.
        try:
            positional_compound_residue_frequencies = np.sum(
                (self.background.compound_residue_matrix.to_numpy()
                * positional_background_frequencies[:, np.newaxis]),
                axis=2
            )
        # Pass if compound residues are not enabled.
        except NameError:
            pass
        # Merge simple and compound residue background frequencies.
        else:
            positional_background_frequencies = np.concatenate(
                (positional_background_frequencies,
                    positional_compound_residue_frequencies),
                axis=1
            )

        return positional_background_frequencies

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
        
        for i,frequency in enumerate(positional_residue_frequencies.flatten()):
            if frequency == self.subset_tensor.shape[0]:
                binomial_p_values.append(np.nextafter(0, 1, dtype='float64'))
            elif frequency == 0:
                binomial_p_values.append(1.0)
            else:
                # Calculate probability of more extreme outcome via regularized
                # incomplete beta function.
                binomial_p_value = special.betainc(frequency + 1.0,
                                        self.subset_tensor.shape[0] - frequency,
                                        flattened_positional_background_frequencies[i])
                binomial_p_values.append(binomial_p_value)

        binomial_p_values = np.array(binomial_p_values, dtype='float64').reshape(
                                            positional_residue_frequencies.shape)

        return binomial_p_values

    def calculate_positional_weights(self,
                                        positional_residue_frequencies):
        '''

        Parameters:
            positional_residue_frequencies -- 2D NumPy Array.
                                                            Dtype = int8.

        Returns:
            vectorized_positional_weighting_terms -- 1D NumPy Array.
                                                        Dtype = float64.
        '''

        vectorized_positional_weighting_terms = (
            1 / np.count_nonzero(
                positional_residue_frequencies[:, :len(self.background.alphabet)],
                axis=1
            )
        )

        return vectorized_positional_weighting_terms

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

        max_positional_residue_enrichments = np.zeros_like(enrichment_values).astype('int8')
        max_positional_residue_enrichments[enrichment_values == np.nanmax(enrichment_values)] = 1

        return max_positional_residue_enrichments

    def _remove_compound_residue_constituents(self,
                                                residue_coordinates):
        """
        Add simple residue constituents of compound residue to removed
            residue array.
        
        Parameters:
            residue_coordinates -- Tuple.
            background -- Background instance.

        Returns:
            None
        """

        # Get row index of compound residue in compound residue matrix.
        ind = residue_coordinates[1] - len(self.background.alphabet)
        # Get row from compound residue matrix as array.
        compound_residue_constituents = (
            self.background.compound_residue_matrix.iloc[ind]
        ).astype('int8')
        # Set consituents to 1 in removed positional residue array.
        self.removed_positional_residues[
            residue_coordinates[0],
            :len(self.background.alphabet)] = np.logical_or(
                self.removed_positional_residues[
                    residue_coordinates[0],
                    :len(self.background.alphabet)],
                compound_residue_constituents
        ).astype('int8')

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

    def character_pattern(self):
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

        pattern_df = pd.DataFrame(self.pattern_matrix,
                                    columns=self.background.ordered_residues)
        pattern_df = pattern_df.apply(get_positional_residue, axis=1)

        return pattern_df

    def generate_sequence_strings(self):
        """
        Decode each sequence in subset_tensor array to string format.
            Return all sequence strings in a list.

        Parameters:
            None

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
        
        output_file = output_directory + (
            '/' + self.pattern_id + '_logo_map.pdf'
        )
        
        with open(output_file, "wb") as f:
            f.write(w.pdf_formatter(data, mformat))

        self._logo_map = output_file

    @property
    def logo_map(self):
        return self._logo_map
    


class PatternContainer:
    '''

    Attributes:
        pattern_list -- List.

    Methods:
        add_new_pattern --
    '''

    def __init__(self,
                    sequence_tensor,
                    background,
                    output_directory='',
                    initial_pattern=None,
                    max_depth=None,
                    p_value_cutoff=0.01,
                    occurrences=2,
                    fold_change_cutoff=1,
                    multiple_testing_correction=True,
                    positional_weighting=True):
        '''
        Initialize PatternContainer object and start recursive pattern
            extraction.

        Parameters:
            sequence_tensor -- 3D Numpy Array.
            background -- Background class instance.
            initial_pattern -- 2D Numpy Array.
            max_depth -- Int.
            p_value_cutoff -- Float.
            occurrences -- Int.
            fold_change_cutoff -- Float.
            multiple_testing_correction -- Bool.
            positional_weighting -- Bool.

        Returns:
            None
        '''

        self.output_directory = os.path.join(output_directory)
        self.pattern_list = []
        self.sequence_tensor = sequence_tensor
        self.background = background
        
        # Initialize pattern template.
        if initial_pattern == None:
            # Null pattern of shape (sequence_tensor.shape[1],
            #                       len(background.ordered_residues))
            initial_pattern = np.zeros((sequence_tensor.shape[1],
                                    len(background.ordered_residues)),
                                    dtype='int8')
        
        # Begin recursive pattern tree construction.
        self.add_new_pattern(Pattern(
            pattern_container=self,
            subset_tensor=self.sequence_tensor,
            initial_pattern=initial_pattern,
            # Starting depth equal number of positions in initial_pattern.
            current_depth=np.count_nonzero(np.any(initial_pattern, axis=1)),
            # Max_depth either set by argument, or as len(initial_pattern)
            max_depth=(max_depth if max_depth is not None else len(initial_pattern)),
            p_value_cutoff=p_value_cutoff,
            occurrences=occurrences,
            fold_change_cutoff=fold_change_cutoff,
            positional_weighting=positional_weighting,
            multiple_testing_correction=multiple_testing_correction))

    def add_new_pattern(self, new_pattern):
        '''

        Parameters:
            new_pattern --

        Returns:
            None
        '''
        
        patterns = [np.array_str(pattern.pattern_matrix.flatten())
                    for pattern in self.pattern_list]
        
        if new_pattern.pattern_matrix.any() and \
            new_pattern.subset_tensor.shape[0] > 2 and \
            np.array_str(new_pattern.pattern_matrix.flatten()) not in patterns:
            
            self.pattern_list.append(new_pattern)

            print(len(self.pattern_list))

        return None

    def post_processing(self):
        '''
        Removes patterns if their exclusive sequence subset is two
            sequences or less.

        Parameters:
            None
        
        Returns:
            None
        '''

        drop_patterns = []
        for pattern in self.pattern_list:
            if pattern.subset_tensor.shape[0] < 2:
                drop_patterns.append(pattern)

        self.pattern_list = [i for i in self.pattern_list
                            if i not in drop_patterns]

        self._generate_pattern_outputs()
        
        return None

    def _generate_pattern_outputs(self):
        """Last steps for retained patterns."""
        for pattern in self.pattern_list:
            pattern.generate_sequence_strings()
            pattern.generate_logo_map(output_directory)