import itertools
import pandas as pd
import pymysql.cursors
import re
import math
from zigzag.preprocessing import *


def generate_residue_product(residues,depth):
    '''
    Takes list of residues and generates a list of all possible
        permutations. 

    Parameters:
        residues -- List. Contains list of amino acids.
        depth -- Int. Specifies number of amino acids per combination.

    Returns:
        permutations -- List. Contains tuples for each combination of
                                residues with m specified by depth.
    '''

    # Generate permutation from residue list with specified number of members
    product = [i for i in itertools.product(residues,repeat=depth)]

    return product


def generate_positional_combinations(width,depth):
    '''

    Parameters:
        width -- Int.
        depth -- Int.

    Returns:
        positional_combinations -- List.
    '''

    positional_combinations = [i for i in itertools.combinations([j for j in range(width)],depth)]

    return positional_combinations


def generate_regular_expressions(residue_product,positional_combinations):
    '''
    Generates a set of regular expressions for each supplied permutation
        of amino acids. Regular expression consists of the specified
        residues separated by a variable number of spaces such that the
        set of all positional arrangements accomodated by the width of
        the experimental data set are produced.

    Parameters:
        residue_product -- List. Contains permutations of amino acids.
        positional_combinations -- List.

    Returns:
        regular_expressions -- List. Contains regular expression for
                                    each amino acid combination.
    '''

    regular_expressions = []

    # Map amino acid permutations to positional combinations.
    residue_mappings = itertools.product(residue_product,positional_combinations)

    # Loop through mappings to build regular expressions.
    for mapping in residue_mappings:
        # Initialize regular expression.
        regular_expression = r"(?=("
        # Loop through positions to construct regular expression for mappping.
        for i,position in enumerate(mapping[0]):
            # Concatenate positional residue.
            regular_expression += re.escape(position)
            # Concatenate positional separation.
            if i < (len(mapping[1]) - 1):
                separation = mapping[1][i+1] - mapping[1][i] - 1
                if separation > 0:
                    regular_expression += r"(.){" + re.escape(str(separation)) + r"}"
        # Terminate regular expression.
        regular_expression += r"))"
        # Append regular expression to regular_expressions list.
        regular_expressions.append(regular_expression)

    return list(sorted(set(regular_expressions)))


def calculate_total_residues(background,width,depth):
    '''
    Calculates total number of possible residues or permutations of
        residues in background data set. For single amino acids this is
        equivalent to the total number of residues in the background
        data set. For permutations of residues, this is based on the
        total possible number of permutations in the background data set
        for a permutation of a given width.

    Parameters:
        background -- List.
        width -- Int.
        depth -- Int.

    Returns:
        total -- Int. Total number of combinations given specified set
                            of sequences, motif length, and number of
                            elements in each combination.
    '''

    total = 0

    for sequence in background:
        num_positions = len(sequence) - width + 1
        sequence_total = ( math.factorial(width) / ( math.factorial(depth) * math.factorial(width-depth) ) ) \
                            + ( (num_positions - 1) \
                            * ( math.factorial(width-1) / ( math.factorial(depth-1) * math.factorial(width-depth) ) ) )
        
        total += sequence_total

    return total


def calculate_background_frequency(background,regular_expressions):
    '''
    Calculate absolute and percentage frequencies for specified residues
        or permutations of residues in background data set.

    Parameters:
        background -- List.
        regular_expressions -- List.

    Returns:
        background_frequencies -- Dict.
    '''

    totals = []
    # Find all occurrences of every regular expression in each seqeunce using capturing group within lookahead.
    for regular_expression in regular_expressions:
        total = 0
        for sequence in background:
            matches = re.finditer(regular_expression,sequence)
            total += len([match.group(1) for match in matches])
        totals.append(total)

    background_frequencies = pd.Series(totals,index=regular_expressions)

    return background_frequencies


def extract_unique_residues(sequences):
    '''
    Takes a list containing amino acids sequences. Returns a list of all
        unique amino acids found in the sequences.

    Parameters:
        sequences -- List. Contains amino acid sequences.
    
    Returns:
        unique -- List. Unique amino acids occurring in list of
                        sequences.
    '''

    unique = sorted(list(set(''.join(sequences))))

    return unique


def connect_to_background_frequencies_database():
    '''
    Opens a new connection to the MySQL background_frequencies database
        using PyMySQL.

    Parameters:
        None

    Returns:
        connection -- PyMySQL Connection Object. Connection to
                                                background_frequencies
                                                database.
    '''

    connection = pymysql.connect(host='localhost', \
                                    user='zigzag', \
                                    password='zigzagpass', \
                                    db='background_frequencies', \
                                    cursorclass=pymysql.cursors.DictCursor)

    return connection


def retrieve_background_frequencies(regular_expressions):
    '''
    Queries background_frequencies database for existing
        background frequency tables.

    Parameters:
        regular_expression -- List. Identified pattern used to query
                                    background_frequencies table on
                                    primary key index.

    Returns:
        background_frequency -- String.
    '''

    conn = connect_to_background_frequencies_database()
    try:
        with conn.cursor() as cursor:
            query = "SELECT * FROM background_frequencies WHERE pattern in ('{}')".format('\',\''.join(regular_expressions))
            print(query)
            cursor.execute(query)
            background_frequency = cursor.fetchall()
    finally:
        conn.close()

    return background_frequency


def drop_background_frequency_table():
    '''
    Drops existing background frequency table from
        background_frequencies database.

    Parameters:
        None

    Returns:
        None
    '''

    conn = connect_to_background_frequencies_database()

    return None


def create_background_frequency_table():
    '''
    Creates table to store calculated background frequencies in
        background_frequencies database.

    Parameters:
    
    Returns:
        None
    '''

    try:
        conn = connect_to_background_frequencies_database()
        with conn.cursor() as cursor:
            query = "CREATE TABLE background_frequencies (pattern varchar(32) NOT NULL, absolute_frequency int NOT NULL, percentage_frequency double NOT NULL, PRIMARY KEY (pattern))"
            cursor.execute(query)
        conn.commit()
    finally:
        conn.close()

    return None


def insert_background_frequency(pattern,absolute,percentage):
    '''
    Inserts calculated background frequencies into corresponding
        tables in background_frequencies database.

    Parameters:
        pattern -- String.
        absolute -- Int.
        percentage -- Float.

    Returns:
        None
    '''

    try:
        conn = connect_to_background_frequencies_database()
        with conn.cursor() as cursor:
            query = "INSERT INTO background_frequencies (pattern,absolute_frequency,percentage_frequency) VALUES ('{}', {}, {})".format(pattern,absolute,percentage)
            print(query)
            cursor.execute(query)
        conn.commit()
    finally:
        conn.close()

    return None


def generate_background_frequencies(background_path,width,depth):
    '''

    Parameters:
        background_path --
        width --
        depth --

    Returns:
        None
    '''

    background = import_fasta(background_path)
    unique_residues = extract_unique_residues(background)
    residue_product = generate_residue_product(unique_residues,depth)
    positional_combinations = generate_positional_combinations(width,depth)
    regular_expressions = generate_regular_expressions(residue_product,positional_combinations)
    total = calculate_total_residues(background,width,depth)
    background_frequencies = calculate_background_frequency(background,regular_expressions)
    for regular_expression,background_frequency in background_frequencies.iteritems():
        insert_background_frequency(regular_expression,background_frequency/total)

    return None


def extract_aligned_peptides(pattern, background):
    """
    Apply regular expression to sequences in background data set.
        Extract matching subsequences corresponding to alignment pattern
        matching context of foreground.

    Parameters:
        pattern -- Regular expression string.
        background -- List. Contains background sequence set.

    Returns:
        aligned_background_peptides -- List. Matching peptides from
                                        background data set.
    """

    aligned_background_peptides = []
    for sequence in background:
        aligned_background_peptides += [match[0] for match in re.findall(
                                            pattern, sequence)]

    return aligned_background_peptides


def ignore_background_positions(background, ignore_positions):
    """

    Remove positions from aligned background peptides.

    Parameters:
        background -- List. Contains aligned background peptides.
        ignore_positions -- List. Contains indexes of positions to
                                    remove from background peptides.
    Return:
        sliced_background -- List. Contains aligned background peptides
                                    with ignore_positions removed.
    """

    sliced_background = [''.join(residue for i, residue in enumerate(sequence)
                            if i not in ignore_positions)
                            for sequence in background]

    return sliced_background



def background_from_aligned_peptides(background,
                                        pattern,
                                        ignore_positions):
    """

    Parameters:
        background --
        pattern --
        ignore_positions --

    Returns:
       background -- 
    """



    return


def main():
    '''

    '''
    
    background_path = '../../data/UniProt/uniprot.fasta'
    background_sequences = import_fasta(background_path)
    unique_residues = extract_unique_residues(background_sequences)
    amino_acid_combinations = generate_residue_product(unique_residues,2)


if __name__ == '__main__':
    main()