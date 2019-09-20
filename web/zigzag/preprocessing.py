import math
import csv
import itertools

import numpy as np
import pandas as pd

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1000000)

"""
def import_fasta(fasta_path):
    '''


    Parameters:
        fasta_path -- String. Path to background data set FASTA file.
    Returns:
        sequences -- List. Contains seqeunces from background data set
                            FASTA file.
    '''

    with open(fasta_path) as fasta:
        sequences = []
        for i,line in enumerate(fasta):
            if i==0:
                sequence = ''
                continue
            elif line[0]=='>':
                sequences.append(sequence)
                sequence = ''
                continue
            else:
                sequence += line.rstrip()
        sequences.append(sequence)

    return sequences
"""

def generate_positions(center,num_positions):
    '''

    Parameters:
        center -- Boolean. If true, column labels will count away from
                        center (leftward direction = non-prime,
                        rightward direction = prime). If false, column
                        labels will count from left to right beginning
                        with 1.
        num_positions -- Int. Total number of position labels to
                            generate.

    Returns:
        positions -- List. Contains position labels.
    '''

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
        positions = ['p{}'.format(i) for i in range(1,num_positions+1)]

    return positions


def load_sequences(data, center):
    '''
    Generate and return a pandas dataframe containing input sequences.
        Additional columns include sequence ID, cluster ID,
        alignmentOffset (+ shifts sequence left, - shifts sequence
        right), and peptide score.
    
    Parameters:
        data -- String. Path to input data set in text file form.
        center -- Boolean. If true, column labels will count away from
                        center (leftward direction = non-prime,
                        rightward direction = prime). If false, column
                        labels will count from left to right beginning
                        with 1.

    Returns:
        sequence_df -- Pandas DataFrame. Contains input data set split
                        to columns by position and accompanying fields
                        as specified above.
    '''

    split_sequences = []

    with open(data) as sequences:
        for sequence in sequences:
            split_sequences.append(list(sequence.rstrip().upper()))

    cols = generate_positions(center,len(split_sequences[0]))

    sequence_df = pd.DataFrame(split_sequences, columns=cols)

    return sequence_df


def filter_missing_residues(data,condition):
    '''
    Filters vectorized sequences based on missing residues.

    Parameters:
        data -- 3D Numpy Array.
        condition -- String. If 'any', remove sequences containing any
                        missing residues. If 'all', remove sequences
                        missing all residues.

    Returns:
        filtered_data -- 3D Numpy Array. Contains vectorized sequences
                            filtered on missing residues and user-
                            specified filtering condition.
    '''

    if condition == 'any':
        filtered_data = data[tuple(np.where(np.any(np.all(
            data==0,axis=2),axis=1)==False))]
    elif condition == 'all':
        filtered_data = data[tuple(np.where(np.all(np.all(
            data==0,axis=2),axis=1)==False))]

    return filtered_data


def filter_unknown_residues(data,condition):
    '''
    Filters character sequences based on unknown residue codes using
        list of expected residue codes from background data set.

    Parameters:
        data -- Pandas DataFrame. Contains input data set split to
                                    columns by position and accompanying
                                    fields as specified above.
        condition -- String. If 'any', remove sequences containing any
                                missing residues. If 'all', remove
                                sequences missing all residues.

    Returns:
        filtered_data -- Pandas DataFrame. Contains input data filtered
                                            based on expected residue
                                            codes and user-specific
                                            filtering condition.
    '''


    if condition == 'any':
        filtered_data = data.iloc[
            data[data.isin(SINGLE_LETTER_CODES)].dropna(how='any').index
        ]
    if condition == 'all':
        filtered_data = data.iloc[
            data[data.isin(SINGLE_LETTER_CODES)].dropna(how='all').index
        ]

    return filtered_data


def importSubstitutionMatrix(substitutionMatrixFile):
    '''
    Import substitution matrix to pandas dataframe. Takes text file as
        argument. Returns pandas dataframe.
    
    Parameters:
        substitutionMatrixFile -- String. Path to text file containing
                                            amino acid substitution
                                            probabilities.

    Returns:
        substitutionMatrixDataFrame -- Pandas DataFrame. Symmetrical
                                            matrix containing pairwise
                                            amino acid substitution
                                            probabilities.
    '''

    with open(substitutionMatrixFile) as submat:

        matrix = []

        for i,line in enumerate(submat):
            if (i==0):
                headers = line.split()
                continue

            lineValues = line.split()
            lineValues.pop(0)

            for i,value in enumerate(lineValues):
                lineValues[i] = float(value)

            matrix.append(lineValues)
        
    substitutionMatrixDataFrame = pd.DataFrame(matrix,
                                                index=headers,
                                                columns=headers)

    return substitutionMatrixDataFrame

"""
def importBackgroundFrequencies(frequencySet):
    '''
    Generates and returns pandas dataframe containing background
        frequencies derived from external source (SWISSPROT, etc.).

    Parameters:
        frequencySet -- String. Path to text file containing amino acid
                                    background frequencies.

    Returns:
        bgFrequencies -- Dictionary. Background frequencies from file
                                        mapped to amino acids.

    '''

    cols = ['aa','frequency']

    bgFrequencies = {}

    with open(frequencySet) as bgFq:
        reader = csv.reader(bgFq,delimiter=',')

        for row in reader:
            bgFrequencies[row[0]] = float(row[1])

    return bgFrequencies
"""
"""
def generatePossibleSubstitutionGroups(substitutionProbabilities,bg):
    '''
    Generates possible substitution groups. Substitution groups
        are doublets or triplets of amino acids for which the
        average enrichment score of the combination of the amino
        acids into a single group will exceed the average enrichment
        score of the amino acids individually (cutoffs being 0.5 for
        doublets and 0.33 for triplets).

    Parameters:
        substitutionProbabilities -- Pandas DataFrame. Contains pairwise
                                        substitution probabilities for
                                        all amino acids in supplied
                                        substitution probability matrix.
        bg -- Dict. Contains background frequency of each amino acid in
                    supplied background frequency file.
    Returns:
        substitutionGroupDataFrame -- Pandas DataFrame.
    '''

    # Generate amino acid list from substitution matrix.
    singletTuples = [[i,i,1,(bg[i]/100),1,0,0,0] for i in substitutionProbabilities.index.values.tolist() if i not in ['B','J','X','Z']]
    # Generate list of possible doublets and list of index-linked probabilities.
    doublets = [''.join(sorted(list(doublet))) for doublet in itertools.combinations(sorted([i[0] for i in singletTuples]),2) if substitutionProbabilities[doublet[0]][doublet[1]] > 0.5]
    doubletTuples = [[doublet,doublet,substitutionProbabilities[doublet[0]][doublet[1]],((np.sum([bg[r] for r in doublet])*substitutionProbabilities[doublet[0]][doublet[1]])/100),1,0,0,0] for doublet in doublets]

    # Combine doublets into pairs for triplet calculation.
    doubletPairs = [doubletPair for doubletPair in itertools.combinations(doublets,2) if len(set(doubletPair[0])&set(doubletPair[1])) > 0]

    # Initialize triplet list
    triplets = []
    tripletTuples = []

    # Loop through doublet pairs to find possible triplets.
    for doubletPair in doubletPairs:
        # Get union of pair. Continue if union already in triplet list.
        union = ''.join(sorted(list(set(doubletPair[0])|set(doubletPair[1]))))
        if (union in triplets):
            continue

        # Check for intersection between doublets in pair
        intersection = set(doubletPair[0])&set(doubletPair[1])

        # Skip remaining steps if doublets in pair do not have one common member.
        if ( len(intersection)<1 ):
            continue
        # Extract list of symmetric difference between doublets for combined probability calculation.    
        symmetricDifference = list(set(doubletPair[0])^set(doubletPair[1]))
        
        # Calculate combined probability as product of substitution probability for each member of symmetric difference with intersection.
        for i in intersection:
            combinedProbability = np.prod([substitutionProbabilities[i][residue] for residue in symmetricDifference])

        # If combined probability > 0.33, add union to possible triplet list.
        if ( combinedProbability > 0.33 ):
            triplets.append(union)
            # "Directional" approach
            tripletTuples.append([union,intersection,combinedProbability,((np.sum([bg[''.join(list(intersection))]*substitutionProbabilities[''.join(list(intersection))][r] for r in symmetricDifference]) + np.sum([bg[r]*substitutionProbabilities[''.join(list(intersection))][r] for r in symmetricDifference]))/100),1,0,0,0])
            # original calculation
            # tripletTuples.append([union,combinedProbability,((np.sum([bg[r] for r in union])*combinedProbability)/100),0,0,0])
    substitutionGroups = singletTuples + doubletTuples + tripletTuples

    cols = [
        'substitutionGroup',
        'definingResidue',
        'probability',
        'backgroundFrequency',
        'binomial_p_value',
        'standardDeviation',
        'frequency',
        'enrichment'
    ]
    
    substitutionGroupDataFrame = pd.DataFrame(substitutionGroups,columns=cols)

    return substitutionGroupDataFrame
"""