import re
import csv
import io
from collections import Counter
from collections import namedtuple

import numpy as np
import pandas as pd

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1000000)

CompoundResidue = namedtuple(
    'CompoundResidue',
    ['description', 'residues']
)

compound_residues = {
    '1': CompoundResidue(
        # Encoded as 
        description='helix breaker',
        residues=[
            'G',
            'P',
        ]
    ),
    '2': CompoundResidue(
        # Encoded as 
        description='aliphatic',
        residues=[
            'P',
            'A',
            'L',
            'V',
            'I',
        ]
    ),
    '3': CompoundResidue(
        # Encoded as 
        description='beta-branched (rigid)',
        residues=[
            'V',
            'I',
            'T',
        ]
    ),
    '4': CompoundResidue(
        # Encoded as 
        description='hydrophobic',
        residues=[
            'M',
            'F',
            'W',
            'Y',
            'P', 
            'A',
            'L',
            'V',
            'I',
            'T',
            'C',
        ]
    ),
    '5': CompoundResidue(
        # Encoded as 
        description='nucleophile in proteases',
        residues=[
            'C',
            'S',
        ]
    ),
    '6': CompoundResidue(
        # Encoded as delta.
        description='polar',
        residues=[
            'N',
            'Q',
            'S',
            'C',
            'T',
            'Y',
            'H',
            'R',
            'K',
            'D',
            'E',
        ]
    ),
    '7': CompoundResidue(
        # Encoded as 
        description='aromatic',
        residues=[
            'F',
            'W',
            'Y',
            'H',
        ]
        
    ),
    '8': CompoundResidue(
        # Encoded as 
        description='proton donor/acceptor in catalysis, Fe2+, Zn2+ binder',
        residues=[
            'C',
            'H',
        ]
    ),
    '9': CompoundResidue(
        # Encoded as 
        description='base, DNA, RNA binder',
        residues=[
            'R',
            'K',
        ]
    ),
    '10': CompoundResidue(
        # Encoded as 
        description='acid, Ca2+, Mg2+ binder',
        residues=[
            'D',
            'E',
        ]

    ),
    '11': CompoundResidue(
        # Encoded as 
        description='charged',
        residues=[
            'H',
            'R',
            'K',
            'D',
            'E',
        ]
    ),
}


def vectorize_sequences(data,
                        background,
                        empty_position_value=0):
    '''

    Parameters:
        data -- Pandas DataFrame.
        background -- Instance of Background class.
        empty_position_value -- Int. Expects value of 0 or 1. If 0,
                                matching amino acids coded as 1. If 0,
                                matching amino acids coded as 1.

    Returns:
        vectorized_sample -- 3D NumPy Array.
    '''

    vectorized_sample = np.array([vectorize_pattern(sequence, background)
                                    for i, sequence
                                    in data.iterrows()]).reshape(len(data),
                                                            len(data.columns),
                                                            len(background.ordered_residues))

    return vectorized_sample


def vectorize_pattern(pattern,
                        background,
                        empty_position_value=0):
    '''
    Takes pattern of Pandas Series type, converts single
        letter encoding to binary encoding. Returns vectorized
        pattern as Numpy array.

    Parameters:
        pattern -- Pandas Series.
        background -- Instance of Background class.
        empty_position_value -- Int. Expects value of 0 or 1. If 0,
                                matching amino acids coded as 1. If 0,
                                matching amino acids coded as 1.
    
    Returns:
        vectorized_pattern -- Numpy Array.
    '''

    vectorized_pattern = []
    
    for position in pattern:
        residue = position.strip('[]')
        
        position_vector = pd.Series(empty_position_vector(
                                len(background.alphabet), empty_position_value),
                            index=background.alphabet)

        # Encode residue with non-null value.
        if residue in background.alphabet:
            position_vector[residue] = int(not empty_position_value)
        # Code under-represented residues to -1.
        elif residue.upper() in background.alphabet:
            position_vector[residue] = -1

        positional_residues = position_vector.tolist()

        try:
            positional_residues += np.sum(position_vector
                * background.compound_residue_matrix, axis=1).tolist()
        except NameError:
            pass

        vectorized_pattern += positional_residues

    vectorized_pattern = np.array(vectorized_pattern, dtype='int8')

    return vectorized_pattern


def empty_position_vector(length, empty_position_value=0):
    '''

    Parameters:
        length -- Int. Length of alphabet.
        empty_position_value -- Int. Expects value of 0 or 1. If 0,
                                matching amino acids coded as 1. If 0,
                                matching amino acids coded as 1.

    Returns:
        empty_position -- 1D Numpy Array.
    '''

    if empty_position_value == 1:
        empty_position = np.ones(length, dtype='int8')
    elif empty_position_value == 0:
        empty_position = np.zeros(length, dtype='int8')

    return empty_position


def align_sequences(context, sequences, width):
    """
    Take unaligned sequences and context. Map unaligned sequences to
        context. Truncated prime segment to half width. Extend sequence
        by half width in non-prime direction based on context. Replace
        ambiguous residues in non-prime segment with "X". Fill missing
        positions with "-".

    Parameters:
        context -- Pandas DataFrame. Expected columns: [
                        'id' (ID from supplied FASTA file),
                        'sequence' (sequence strings),
                    ]
        sequences -- Pandas DataFrame. Expected columns: [
                        'context_id' (IDs matching context
                                        elements, sep=';'),
                        'sequence' (sequence strings),
                    ]
        width -- Int. Number of positions in aligned sequence output.

    Returns:
        aligned_sequences -- List. Aligned sequences of specified width
                                containing prime segment from sequences
                                and non-prime segment from matching
                                context.
    """

    non_prime_width = width // 2
    
    # Generate aligned sequence by mapping sequence segments to context.
    aligned_sequences = []
    for i, sequence in sequences.iterrows():
        # Select pre-mapped context elements if available.
        try:
            context_ids = sequence['context_id'].replace(' ', '').split(';')
        except KeyError:
            context_elements = context['sequence'].tolist()
        else:
            # Currently only supporting protein accesssion ID lookup for
            # Swiss-Prot. Any other IDs will fail and revert to default
            # behavior of attempting to match sequence segment against
            # all context sequences.
            swissprot_ids = []
            for context_id in context_ids:
                if context_id.startswith('sp|'):
                    swissprot_ids.append(context_id)
            if swissprot_ids:
                context_elements = context[context['swissprot_id'].isin(
                    context_ids)]['sequence'].tolist()
            else:
                context_elements = context['sequence'].tolist()
        # Regular expression matching of sequence to context.
        matches = set()
        for context_element in context_elements:
            matches |= match_segment_to_context(sequence['sequence'],
                                                    context_element,
                                                    width,
                                                    non_prime_width)
        if matches:
            # Convert matches to data frame for disambiguation.
            matches = pd.DataFrame([list(match) for match in matches])
            # Merge unique matches with 'X' in ambiguous positions.
            merged_sequence = ''.join((matches.loc[0, i] if unambiguous else 'X'
                            for i, unambiguous
                            in enumerate(matches.nunique(axis=0) == 1)))
            aligned_sequences.append(merged_sequence)
        else:
            continue

    return aligned_sequences


def match_segment_to_context(segment, context_element, width, non_prime_width):
    """
    Regular expression match of segment string to context element
        string. Return matches.

    Parameters:
        segment -- String.
        context_element -- String.
        width -- Int.
        non_prime_width -- Int.

    Returns:
        matches -- Set.
    """

    matches = set()
    # Set of unique non-prime matches in context sequence.
    for match in re.finditer(segment, context_element):
        match_index = match.start(0)
        start = (match_index - non_prime_width
                    if match_index >= non_prime_width else 0)
        end = match_index + width - non_prime_width
        non_prime = context_element[start:match_index].rjust(non_prime_width, '-')
        prime = context_element[match_index:end].ljust(width - non_prime_width, '-')
        matches.add(non_prime + prime)

    return matches


def remove_genomic_n_termini(aligned_sequences):
    """Removes any genomic N-terminal sequences from input data set."""
    pass


def generate_positions(center, num_positions):
    """
    Generate position labels for sequence dataframe columns or series
        indexes.

    Parameters:
        center -- Boolean. If true, column labels will count away from
                        center (leftward direction = non-prime,
                        rightward direction = prime). If false, column
                        labels will count from left to right beginning
                        with 1.
        num_positions -- Int. Total number of position labels.

    Returns:
        positions -- List. Contains position labels.
    """
    if center == True:
        positions = []

        p1_pos = int(num_positions / 2)

        for i in range(num_positions):
            if (i<p1_pos):
                position_label = 'p' + str(p1_pos-i)
                positions.append(position_label)
            else:
                position_label = 'p' + str(abs(p1_pos-(i+1))) + "'"
                positions.append(position_label)
    else:
        positions = ['p{}'.format(i) for i in range(1, num_positions+1)]

    return positions


def sequences_to_df(sequences, center=True):
    """
    Split sequence strings to Pandas DataFrame with one column per
        position.

    Parameters:
        sequences -- File-like object. Contains one sequence per row.
        center -- Boolean. If true, column labels will count away from
                        center (leftward direction = non-prime,
                        rightward direction = prime). If false, column
                        labels will count from left to right beginning
                        with 1.
    Returns:
        sequence_df -- Pandas DataFrame. Contains input data set split
                        to columns by position and accompanying fields
                        as specified above.
    """
    
    split_sequences = [list(sequence.rstrip().upper()) for sequence in sequences]
    cols = generate_positions(center,len(split_sequences[0]))
    sequence_df = pd.DataFrame(split_sequences, columns=cols)

    return sequence_df


def import_fasta(fasta_path, swissprot=True):
    """
    Import fasta to Pandas dataframe with ID column and sequence column.
        To convert pre-aligned sequence FASTA to input dataframe format,
        supply io.StringIO('\n'.join(<fasta_df>['sequence'].tolist()))
        to sequences_to_df.
         
    Parameters:
        fasta_path -- String. Path to fasta file.

    Returns:
        fasta_df -- Pandas DataFrame.
    """
    with open(fasta_path) as fasta:
        sequences = []
        for line in fasta:
            # Skip leading comments (not common any more).
            if line[0] == ';':
                continue
            elif line[0] == '>':
                # Add sequence and ID to sequences if not first ID.
                try:
                    sequences.append([sequence_id, sequence, swissprot_id])
                except NameError:
                    pass
                finally:
                    sequence_id = line.rstrip()[1:]
                    sequence = ''
                    swissprot_id = sequence_id[:sequence_id.find(' ')]
            else:
                # Add characters to sequence if sequence has ID.
                try:
                    sequence += line.rstrip().upper()
                except NameError:
                    continue
        sequences.append([sequence_id, sequence])

    cols = ['id', 'sequence', 'swissprot_id']
    fasta_df = pd.DataFrame(sequences, columns=cols)

    return fasta_df


def select_maxquant_evidence_experiments(evidence, experiments):
    """Select experiments from evidence file."""
    experiments = set(experiments.replace(' ', '').split(','))
    
    return evidence[evidence['Experiment'].isin(experiments)]


def select_maxquant_evidence_modifications(evidence, modifications):
    """Select peptides based on modification."""
    return evidence[evidence['Modifications'].isin(modifications)]


def subset_maxquant_evidence_by_experiment(evidence):
    """Subset evidence by filename. Return as dict of data frames."""
    experiment_names = set(evidence['Raw file'])
    experiment_dict = {experiment: evidence[evidence['Raw file']
                            == experiment].copy()
                        for experiment in experiment_names}
    
    return experiment_dict


def import_maxquant_evidence(evidence_path,
                                drop_reverse_peptides=True,
                                drop_contaminants=True):
    """
    Import and (optionally) expand peptides from MaxQuant
        'evidence.txt' output file.

    Parameters:
        evidence_path -- String. Path to evidence file.
        drop_reverse_peptides -- Boolean.
        drop_contaminants -- Boolean.

    Returns:
        evidence_df -- Pandas DataFrame. Contains rk on
    """

    # Load all experiments to Pandas Data Frame.
    evidence_df = pd.read_csv(evidence_path,
                                sep='\t',
                                low_memory=False,
                                error_bad_lines=False)
    # Drop reverse peptides and contaminants unless disabled.
    if drop_reverse_peptides:
        evidence_df = evidence_df[evidence_df['Reverse'] != '+']
    if drop_contaminants:
        evidence_df = evidence_df[evidence_df['Potential contaminant'] != '+']

    return evidence_df


def expand_maxquant_evidence_sequences(evidence_df,
                                        context,
                                        width=8):
    """
    Expands sequences from MaxQuant experiment dictionary.

    Parameters:
        eidence_dict -- Dict.
        context -- List.
        width -- Int. Default 8.
        center -- Boolean. If true, column labels will count away from
                        center (leftward direction = non-prime,
                        rightward direction = prime). If false, column
                        labels will count from left to right beginning
                        with 1.

    Returns:
        expanded_evidence_dict -- Dict.
    """

    evidence_df.drop_duplicates(inplace=True)
    prime_sequences = evidence_df.copy()[['Proteins', 'Sequence']]
    prime_sequences.columns = ['context_id', 'sequence']
    expanded_sequences = align_sequences(context, prime_sequences, width)
    sequence_df = sequences_to_df(expanded_sequences)
    
    return sequence_df


def fasta_to_sequences(fasta_df, center=False):
    """
    Convenience function to convert FASTA dataframe to dataframe
        containing one column per sequence position.

    Parameters:
        fasta_df -- Pandas DataFrame. Contains FASTA sequences and IDs.
        center -- Boolean. If true, column labels will count away from
                        center (leftward direction = non-prime,
                        rightward direction = prime). If false, column
                        labels will count from left to right beginning
                        with 1.
    Returns:
        sequence_df -- Pandas DataFrame. Contains input data set split
                        to columns by position and accompanying fields
                        as specified above.
    """
    sequence_df = sequences_to_df(fasta_df['sequence'].to_list(), center)

    return sequence_df


class Background:
    """
    Methods and data structures for background frequency data.
    
    Attributes:
        background -- List. 
        background_dict -- Dict. Percentage frequencies mapped to single
                            letter codes.
        background_vector -- 1D Numpy Array. 
                                Dtype = np.float64
                                Percentage frequencies.
    """

    def __init__(self,
                    sequences, 
                    background_dict=None,
                    compound_residues=compound_residues):
        """
        Default background mode. Initialize from list of sequences.

        Parameters:
            background -- List. Background sequence strings.

        Returns:
            None
        """

        # Load single residue background frequencies from context.
        self._background_sequences = sequences
        if self.background_sequences is not None:
            self._background_dict = self._calculate_percentage_frequencies()
        elif background_dict is not None:
            self._background_dict = background_dict
        else:
            raise Exception('No background provided.') 
        
        # Initialize single residue alphabet from context sequences.
        self._alphabet = sorted(list(self.background_dict.keys()))

        # Pre-process compound residues if specified.
        if compound_residues is not None:
            self._background_dict.update(
                self._generate_compound_residues(compound_residues)
            )

        # Generate background frequency vector.
        self._background_vector = self._vectorize_background()

    @classmethod
    def from_csv(cls, background_path):
        """
        Alternative background mode. Initialize from list of background
            frequencies.

        Parameters:
            background_csv --

        Returns:
            Background Instance
        """
        background_dict = {}
        with open(background_path) as background_file:
            reader = csv.reader(background_file,delimiter=',')
            for row in reader:
                background_dict[row[0]] = float(row[1])
        
        return cls(None, background_dict=background_dict)

    @classmethod
    def from_json(cls, background_path):
        """
        Generate background frequencies dictionary from JSON file
            containing percentage frequencies mapped to single letter
            codes.

        Parameters:
            background_path -- String. Path to background frequency CSV
                                file.

        Returns:
            Background Instance
        """
        with open(background_path) as background_file:
            background_dict = json.load(background_file)

        return cls(None, background_dict=background_dict)

    @property
    def background_sequences(self):
        """Get background sequences."""
        return self._background_sequences

    @property
    def compound_residues(self):
        return self._compound_residues
    
    @property
    def compound_residue_codes(self):
        return self._compound_residue_codes
    
    @property
    def compound_residue_frequencies(self):
        return self._compound_residue_frequencies
    
    @property
    def compound_residue_matrix(self):
        return self._compound_residue_matrix

    @property
    def ordered_residues(self):
        return self._ordered_residues
    
    
    @background_sequences.setter
    def background_sequences(self, sequences):
        """
        Set background, background dictionary, and background vector
            from supplied sequences.

        Parameters:
            sequences -- 

        Returns:
            None
        """
        self._background_sequences = sequences
        self.background_dict = self._calculate_percentage_frequencies()
        self.background_vector = self._vectorize_background()

    @property
    def background_dict(self):
        """Get background frequency dictionary."""
        return self._background_dict

    @property
    def background_vector(self):
        """Get background frequency vector."""
        return self._background_vector

    @property
    def alphabet(self):
        return self._alphabet
    
    def _calculate_percentage_frequencies(self):
        """
        Calculate background frequencies from list of sequences.

        Parameters:
            None

        Returns:
            percentage_frequencies -- Dict. Letter code frequencies.
        """
        # Use Counter to compute cumulative frequencies for each residue.
        absolute_frequencies = Counter()
        for sequence in self.background_sequences:
            absolute_frequencies += Counter(sequence)
        # Calculate percentage frequencies from absolute frequencies.
        total = sum(absolute_frequencies.values())
        percentage_frequencies = {residue: ((absolute_frequency / total) * 100)
                                    for residue, absolute_frequency
                                    in absolute_frequencies.items()}

        return percentage_frequencies

    def _generate_compound_residues(self, compound_residues):
        """
        Load compound residues from dictionary of CompoundResidue
            objects. Calculate expected compound residue frequencies
            as sum of expected constituent single residue frequencies.
            Return dictionary of expected frequencies, which is used
            to update background frequency dictionary.

        Parameters:
            compound_residues -- Dict. Contains CompoundResidue objects.
        Returns:
            compound_residue_frequencies -- Dict. Contains expected
                                        compound residue frequencies.
        """

        self._compound_residues = compound_residues
        self._compound_residue_codes = sorted(list(compound_residues.keys()))
        
        self._compound_residue_frequencies = {}
        compound_residue_vectors = []
        
        for compound_residue in self.compound_residue_codes:
            # Initialize vector for constituent single residues.
            compound_residue_vector = pd.Series(
                np.zeros(len(self.alphabet),
                            dtype='int8'),
                index=self.alphabet
            )
            # Initialize expected compound residue frequency to 0.
            compound_residue_frequency = 0
            # Iterative over constituent single residues.
            for residue in compound_residues[compound_residue].residues:
                # Set compound residue vector element.
                compound_residue_vector[residue] = 1
                # Add expected constituent frequency to total.
                compound_residue_frequency += self.background_dict[residue]
            # Add compound residue to frequencies dictionary.
            self._compound_residue_frequencies[compound_residue] = compound_residue_frequency
            # Add compound residue to single residue vector dictionary.
            compound_residue_vectors.append(compound_residue_vector)
        
        self._compound_residue_matrix = pd.DataFrame(compound_residue_vectors,
                                                        index=self.compound_residue_codes)

        return self._compound_residue_frequencies

    def _vectorize_background(self):
        """
        Convert background frequencies dictionary to Numpy array.

        Parameters:
            None

        Returns:
            background_frequency_vector -- 1D Numpy Array. 
                                            Dtype = np.float64
                                            Percentage frequencies.
        """
        
        try:
            self._ordered_residues = self.alphabet + self.compound_residue_codes
        except NameError:
            self._ordered_residues = self.alphabet

        background_frequency_vector = np.array(
            [self.background_dict[residue] for residue in self.ordered_residues],
            dtype=np.float64)
        
        return background_frequency_vector