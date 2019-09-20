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

ExtendedSequences = namedtuple(
    'ExtendedSequences',
    ['n_term', 'c_term']
)

"""
compound_residues = {
    '1': CompoundResidue(
        description='helix breaker',
        residues=[
            'G',
            'P',
        ]
    ),
    '2': CompoundResidue(
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
        description='beta-branched (rigid)',
        residues=[
            'V',
            'I',
            'T',
        ]
    ),
    '4': CompoundResidue(
        description='hydrophobic',
        residues=[
            'M',
            'F',
            'W',
            'P', 
            'A',
            'L',
            'V',
            'I',
        ]
    ),
    '5': CompoundResidue(
        description='transitional polarity',
        residues=[
            'C',
            'Y',
            'T',
        ]
    ),
    '6': CompoundResidue(
        description='polar',
        residues=[
            'N',
            'Q',
            'S',
            'H',
        ]
    ),
    '7': CompoundResidue(
        description='aromatic',
        residues=[
            'F',
            'W',
            'Y',
            'H',
        ]
        
    ),
    '8': CompoundResidue(
        description='base, DNA, RNA binder',
        residues=[
            'R',
            'K',
        ]
    ),
    '9': CompoundResidue(
        description='acid, Ca2+, Mg2+ binder',
        residues=[
            'D',
            'E',
        ]

    )
}
"""

compound_residues = {
    '1': CompoundResidue(
        description='Nonpolar, aliphatic R groups',
        residues=[
            'G',
            'A',
            'V',
            'L',
            'M',
            'I',
        ]
    ),
    '2': CompoundResidue(
        description='Aromatic R groups',
        residues=[
            'F',
            'Y',
            'W',
        ]
    ),
    '3': CompoundResidue(
        description='Polar, uncharged R groups',
        residues=[
            'S',
            'T',
            'C',
            'P',
            'N',
            'Q',
        ]
    ),
    '4': CompoundResidue(
        description='Positively charged R groups',
        residues=[
            'K',
            'R',
            'H',
        ]
    ),
    '5': CompoundResidue(
        description='Negatively charged R groups',
        residues=[
            'D',
            'E',
        ]
    )
}

Sample = namedtuple(
    'Sample',
    ['sequence_df', 'sequence_tensor']
)

SWISSPROT_ACCESSION_PATTERN = re.compile(
    "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)

PYRO_GLU = [
    'Glu->pyro-Glu',
    'Gln->pyro-Glu',
    '2 Oxidation (M),Gln->pyro-Glu',
    'Oxidation (M),Gln->pyro-Glu',
    'Oxidation (M),Glu->pyro-Glu',
]


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

    vectorized_sample = np.array([
        vectorize_pattern(sequence, background) for i, sequence in data.iterrows()
    ]).reshape(len(data), len(data.columns), len(background.ordered_residues))

    return vectorized_sample


def vectorize_pattern(pattern,
                        background,
                        empty_position_value=0,
                        add_constituents=True,
                        add_compound_residues=True):
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
        pattern_matrix -- Numpy Array.
    '''

    pattern_positions = []
    for position in pattern:
        residue = position.strip('[]').upper()
        
        position_vector = pd.Series(
            empty_position_vector(
                len(background.ordered_residues),
                empty_position_value
            ),
            index=background.ordered_residues
        )

        # Encode residue with non-null value.
        if residue in background.ordered_residues:
            position_vector[residue] = int(not empty_position_value)

        # Add positional residues to pattern.
        positional_residues = position_vector.tolist()
        pattern_positions += positional_residues

    # Convert pattern to numpy array.
    pattern_matrix = np.array(pattern_positions, dtype='int8').reshape(
        len(pattern),
        len(background.ordered_residues)
    )

    try:
        if add_constituents:
            # Add constituents to positions containing compound residue
            constituent_pattern_matrix = get_pattern_constituents(pattern_matrix, background)
        else:
            background.compound_residue_matrix
    except AttributeError:
        pass
    else:
        if add_compound_residues:
            # Add compound residues matching constituents.
            compound_residues = np.inner(
                pattern_matrix[:, :len(background.alphabet)],
                background.compound_residue_matrix
            )
            compound_residue_pattern_matrix = pattern_matrix.copy()
            compound_residue_pattern_matrix[:, len(background.alphabet):] = compound_residues

            # Combine compound residue and constituent pattern
            pattern_matrix = np.logical_or(
                constituent_pattern_matrix,
                compound_residue_pattern_matrix
            ).astype(np.int8)
        elif add_constituents:
            pattern_matrix = constituent_pattern_matrix

    return pattern_matrix


def get_pattern_constituents(pattern, background):
    """
    Return pattern featuring all original residues, including compound
        residues, plus compound residue constituents.

    Parameters:
        pattern -- 2D Numpy Array.
        background -- Background instance.

    Returns:
        constituent_pattern -- 2D Numpy Array.
    """

    try:
        # Extend compound residues with blank.
        constituents = np.concatenate(
            (
                background.compound_residue_matrix,
                np.zeros(
                    (1, background.compound_residue_matrix.shape[1]), dtype=np.int8
                )
            ),
            axis=0
        )
    except AttributeError as e:
        raise e
    else:
        # Generate index array from pattern matrix.
        compound_residue_index = np.argmax(pattern, axis=1) - len(background.alphabet)

        # Map simple positions to blank in compound residue matrix.
        compound_residue_index[
            compound_residue_index < 0] = len(background.compound_residue_matrix)
        
        # Get compound residue constituents using index array.
        positional_residue_constituents = constituents[compound_residue_index]

        # Combine constituents with simple fixed residues.
        fixed_position_residues = np.logical_or(
            pattern[:, :len(background.alphabet)],
            positional_residue_constituents)

        constituent_pattern = np.concatenate(
            (
                fixed_position_residues,
                pattern[:, len(background.alphabet):]
            ),
            axis=1
        )

        return constituent_pattern


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


def align_sequences(context, sequences, width, terminal):
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
        terminal -- String. 'n' for extension of N-term sequences,
                        'c' for extension of C-term sequences,
                        'both' for extension of intact peptides.

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
                elif SWISSPROT_ACCESSION_PATTERN.match(context_id):
                    swissprot_ids.append(context_id)
            # Get protein(s) matching peptide from context proteome.
            if swissprot_ids:
                context_elements = context[context['swissprot_id'].isin(
                    context_ids)]['sequence'].tolist()
            else:
                context_elements = context['sequence'].tolist()
        # Regular expression matching of sequence to context.
        extended_sequences = ExtendedSequences(n_term=set(), c_term=set())
        for context_element in context_elements:
            matches = match_segment_to_context(sequence['sequence'].rstrip().upper(),
                                                    context_element,
                                                    width,
                                                    non_prime_width,
                                                    terminal)
            for n_term_match in sorted(list(matches.n_term)):
                extended_sequences.n_term.add(n_term_match)
            for c_term_match in sorted(list(matches.c_term)):
                extended_sequences.c_term.add(c_term_match)

        if extended_sequences.n_term:
            aligned_sequences.append(merge_sequences(extended_sequences.n_term))
        if extended_sequences.c_term:
            aligned_sequences.append(merge_sequences(extended_sequences.c_term))

    return aligned_sequences


def merge_sequences(matches):
    """
    If multiple possible sequence extensions are found in proteome,
        combine matches to single sequence with ambiguous positions
        coded as 'X'.

    Parameters:
        matches -- Set.

    Returns:
        merged_sequence -- String.
    """

    # Convert matches to data frame for disambiguation.
    matches = pd.DataFrame([list(match) for match in sorted(list(matches))])
    # Merge unique matches with 'X' in ambiguous positions.
    merged_sequence = ''.join((matches.loc[0, i] if unambiguous else 'X'
                    for i, unambiguous
                    in enumerate(matches.nunique(axis=0) == 1)))

    return merged_sequence


def match_segment_to_context(segment,
                                context_element,
                                width,
                                non_prime_width,
                                terminal):
    """
    Regular expression match of segment string to context element
        string. Return matches.

    Parameters:
        segment -- String.
        context_element -- String.
        width -- Int.
        non_prime_width -- Int.
        terminal -- String. 'n', 'c', or both.

    Returns:
        matches -- Set.
    """

    matches = ExtendedSequences(n_term=set(), c_term=set())

    # Set of unique non-prime matches in context sequence.
    for match in re.finditer(segment, context_element):
        if terminal == 'n':
            matches.n_term.add(extend_sequence(match, context_element, width, non_prime_width, 'n'))
        elif terminal == 'c':
            matches.c_term.add(extend_sequence(match, context_element, width, non_prime_width, 'c'))
        elif terminal == 'both':
            matches.n_term.add(extend_sequence(match, context_element, width, non_prime_width, 'n'))
            matches.c_term.add(extend_sequence(match, context_element, width, non_prime_width, 'c'))

    return matches


def extend_sequence(match, context_element, width, non_prime_width, terminal):
    """
    Extend sequence in X-terminal direction.

    Parameters:
        match -- Match object.
        width -- Int.
        non_prime_width -- Int.
        terminal -- String.

    Returns:
        extended_sequence -- String.
    """
    
    if terminal == 'n':
        match_index = match.start(0)
    elif terminal == 'c':
        match_index = match.end(0)

    start = (match_index - non_prime_width
                if match_index >= non_prime_width else 0)
    end = match_index + width - non_prime_width
    non_prime = context_element[start:match_index].rjust(non_prime_width, '-')
    prime = context_element[match_index:end].ljust(width - non_prime_width, '-')
    extended_sequence = non_prime + prime

    return extended_sequence


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
    sequence_df.drop_duplicates(inplace=True)

    return sequence_df


def parse_swissprot_accession_number(identifier):
    """
    Parse Swiss-Prot accession number(s) from identifier string.

    Parameters:
        identifier -- String.

    Returns:
        swissprot_accession_number -- String.
    """

    swissprot_accession_number = ';'.join(
        [identifier[match.start():match.end()]
            for match
            in re.finditer(
                SWISSPROT_ACCESSION_PATTERN,
                identifier
            )
        ]
    )

    return swissprot_accession_number


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
                    swissprot_id = parse_swissprot_accession_number(sequence_id)
            else:
                # Add characters to sequence if sequence has ID.
                try:
                    sequence += line.rstrip().upper()
                except NameError:
                    continue
        sequences.append([sequence_id, sequence, swissprot_id])

    cols = ['id', 'sequence', 'swissprot_id']
    fasta_df = pd.DataFrame(sequences, columns=cols)

    return fasta_df


def import_peptide_list(peptide_list_file, delimiter='\t'):
    """
    Import peptide list and row-matched protein IDs from text file.

    Parameters:
        peptide_list_file -- File-like object.
        delimiter -- String.

    Returns:
        peptide_list -- Pandas DataFrame.
    """

    peptide_list = []

    # Parse peptides and accession numbers in input file.
    peptide_reader = csv.reader(peptide_list_file, delimiter=delimiter)
    for row in peptide_reader:    
        context_id = parse_swissprot_accession_number(row[1])
        if not context_id:
            continue
        sequence = row[0].upper()
        peptide_list.append([sequence, context_id])

    cols = ['sequence', 'context_id']
    peptide_list = pd.DataFrame(peptide_list, columns=cols)

    return peptide_list


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
    evidence_df = pd.read_csv(
        evidence_path,
        sep='\t',
        low_memory=False,
        error_bad_lines=False
    )
    # Drop reverse peptides and contaminants unless disabled.
    if drop_reverse_peptides:
        evidence_df = evidence_df[evidence_df['Reverse'] != '+']
    if drop_contaminants:
        evidence_df = evidence_df[evidence_df['Potential contaminant'] != '+']

    return evidence_df


def expand_maxquant_evidence_sequences(evidence_df,
                                        context,
                                        width=8,
                                        terminal='n'):
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
    for i, sequence in prime_sequences.iterrows():
        protein_id = sequence['context_id'].split(';')[0]
        swissprot_id = parse_swissprot_accession_number(protein_id)
        prime_sequences.loc[i, 'context_id'] = swissprot_id
    expanded_sequences = align_sequences(context, prime_sequences, width, terminal=terminal)
    print(expanded_sequences)
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


def peptides_to_sample(peptides, context, background, center=True, width=8, terminal='n'):
    """
    Align peptides and return data frame of positional residues.

    Parameters:
        peptides --
        context --
        width --
        terminal --

    Returns:
        sample -- Sample instance.
    """
    aligned_sequences = align_sequences(context, peptides, width, terminal)
    sequence_df = sequences_to_df(aligned_sequences)
    sequence_tensor = vectorize_sequences(sequence_df, background)
    sample = Sample(sequence_df=sequence_df, sequence_tensor=sequence_tensor)

    return sample


def load_fasta_peptides(peptide_fasta_path, context, background, center=True, width=8, terminal='n'):
    """
    Top-level helper function to load and extend peptides from FASTA
        file.
    """
    # call import_fasta
    # call alignment functions
    
    sample = peptides_to_sample(peptides, context, width, terminal)
    
    return sample


def load_peptide_list_file(peptide_list_path, context, background, center=True, width=8, terminal='n'):
    """
    Top-level helper function to load and extend peptides from text
        file.

    Parameters:
        peptide_list_path -- String.
        context -- Pandas DataFrame.
        terminal -- String.

    Returns:
        sample -- Sample instance.
    """

    # Detect delimiter from file extension (supports comma or tab).
    file_extension = peptide_list_path[-peptide_list_path[::-1].find('.'):].lower()
    if file_extension == 'csv':
        delimiter = ','
    else:
        delimiter = '\t'

    # open peptide list file
    with open(peptide_list_path, 'r') as peptide_list_file:
        peptides = import_peptide_list(peptide_list_file, delimiter)
    
    sample = peptides_to_sample(peptides, context, background, width, terminal)

    return sample


def load_peptide_list_field(peptide_list_field, context, background, center=True, width=8, terminal='n'):
    """
    Top-level helper function to load and extend peptides from text
        field.

    Parameters:

    Returns:
        sample -- Sample instance.
    """
    # create fileio object from peptide list field
    # call import_peptide_list
    peptides = import_peptide_list(peptide_list_file, delimiter)
    sample = peptides_to_sample(peptides, context, width, terminal)

    return sample


def load_prealigned_file(prealigned_file_path, background, center=True):
    """
    Top-level helper function to load pre-aligned sequences from text
        file.
    """
    
    with open(prealigned_file_path, 'r') as prealigned_file:
        sequence_df = sequences_to_df(prealigned_file, center=center)
    sequence_tensor = vectorize_sequences(sequence_df, background)
    sample = Sample(sequence_df=sequence_df, sequence_tensor=sequence_tensor)

    return sample


def load_prealigned_field(prealigned_field, background, center=True):
    """
    Top-level helper function to load pre-aligned sequences from text
        field.
    """
    # create fileio object from field
    # call sequences_to_df
    # return df
    pass


def load_prealigned_fasta(prealigned_fasta_path, background, center=True):
    """
    Top-level helper function to load pre-aligned sequences from FASTA
        file.
    """
    # call import_fasta
    # get sequences from fasta df
    # call sequences_to_df
    # return df
    pass


def load_maxquant_evidence_file(maxquant_evidence_path,
                                context,
                                background,
                                width=8,
                                eliminate_pyro_glu=True,
                                experiment_groups=None):
    """
    Top-level helper function to load and extend peptides from FASTA
        file.
    """

    # Load evidence from MaxQuant evidence.txt file.
    evidence = import_maxquant_evidence(
        evidence_path,
        drop_reverse_peptides=True,
        drop_contaminants=True
    )
    # Remove pyro-glu from peptide list if enabled.
    if eliminate_pyro_glu:
        evidence = evidence[~evidence['Modifications'].isin(PYRO_GLU)]

    # Subset experiments by group if groups supplied by user.
    if experiment_groups:
        pass

    sequence_df = expand_maxquant_evidence_sequences(peptides, context, width)
    sequence_tensor = vectorize_sequences(sequence_df, background)
    sample = Sample(sequence_df=sequence_df, sequence_tensor=sequence_tensor)

    samples = {}

    return samples


class Background:
    """
    Methods and data structures for background frequency data.
    
    Attributes:
        sequences -- List. 
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
        try:
            return self._compound_residues
        except AttributeError as e:
            raise e
    
    @property
    def compound_residue_codes(self):
        try:
            return self._compound_residue_codes
        except AttributeError as e:
            raise e

    
    @property
    def compound_residue_frequencies(self):
        try:
            return self._compound_residue_frequencies
        except AttributeError as e:
            raise e

    @property
    def compound_residue_matrix(self):
        try:
            return self._compound_residue_matrix
        except AttributeError as e:
            raise e

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
    
    def _calculate_percentage_frequencies(self, remove_ambiguous_residues=True):
        """
        Calculate background frequencies from list of sequences.

        Parameters:
            remove_ambiguous_residues = Boolean.

        Returns:
            percentage_frequencies -- Dict. Letter code frequencies.
        """
        # Use Counter to compute cumulative frequencies for each residue.
        absolute_frequencies = Counter()
        for sequence in self.background_sequences:
            absolute_frequencies += Counter(sequence)
        # Remove ambiguous residues.
        if remove_ambiguous_residues:
            try:
                del absolute_frequencies['X']
            except KeyError:
                pass
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
        except AttributeError:
            self._ordered_residues = self.alphabet

        background_frequency_vector = np.array(
            [self.background_dict[residue] for residue in self.ordered_residues],
            dtype=np.float64)
        
        return background_frequency_vector