import numpy as np
import scipy as sp
import pandas as pd
import math
import csv
import itertools
from collections import OrderedDict
import matplotlib.pyplot as plt
import seaborn as sns
import datetime
import os

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1e6)

def loadSequences(data):
    '''
    Generate and return a pandas dataframe containing input sequences.Additional columns include sequence ID, cluster ID,
        alignmentOffset (+ shifts sequence left, - shifts sequence right), and peptide score.
    
    Parameters:
        data -- String. Path to input data set in text file form.

    Returns:
        df -- Pandas DataFrame. Contains input data set split to columns by position and accompanying fields as specified above.
    '''

    seq = []

    with open(data) as peptides:
        for peptide in peptides:
            seq.append(list(peptide.rstrip()))

    cols = []

    numPositions = len(seq[0])

    p1Pos = int(numPositions / 2)

    for i in range(len(seq[0])):
        if (i<p1Pos):
            colName = 'p' + str(p1Pos-i)
            cols.append(colName)
        else:
            colName = 'p' + str(abs(p1Pos-(i+1))) + '\''
            cols.append(colName)

    df = pd.DataFrame(seq, columns=cols)
    df['clusterID'] = ''

    
    return df



def importSubstitutionMatrix(substitutionMatrixFile):
    '''
    Import substitution matrix to pandas dataframe. Takes text file as argument. Returns pandas dataframe.
    
    Parameters:
        substitutionMatrixFile -- String. Path to text file containing amino acid substitution probabilities.

    Returns:
        substitutionMatrixDataFrame -- Pandas DataFrame. Symmetrical matrix containing pairwise amino acid substitution probabilities.
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
        
    substitutionMatrixDataFrame = pd.DataFrame(matrix,index=headers,columns=headers)


    return substitutionMatrixDataFrame



def importBackgroundFrequencies(frequencySet):
    '''
    Generates and returns pandas dataframe containing background frequencies derived from external source (SWISSPROT, etc.).

    Parameters:
        frequencySet -- String. Path to text file containing amino acid background frequencies.

    Returns:
        bgFrequencies -- Dictionary. Background frequencies from file mapped to amino acids.

    '''

    cols = ['aa','frequency']

    bgFrequencies = {}

    with open(frequencySet) as bgFq:
        reader = csv.reader(bgFq,delimiter=',')

        for row in reader:
            bgFrequencies[row[0]] = float(row[1])


    return bgFrequencies



def heatmap(data,proteaseAnnotationFile,outputFileName,title):
    '''
    Generate heatmap of combinatorial clusters using Seaborn. Export heatmap as PNG image file.

    Parameters:
        data -- Pandas DataFrame.
        proteaseAnnotationFile -- String. Path to file containing protease annotations with index matching sequences data frame.
        outputFileName -- String. Path to data output file. Used to create path to image output file.


    Returns:
        None
    '''

    # Initialize set for all unique proteases in annotation file.
    uniqueProteases = set()
    # Initialize list for sequence annotations.
    sequenceAnnotations = []

    # Extract proteases from annotation file.
    with open(proteaseAnnotationFile) as annotations:
        reader = csv.reader(annotations,delimiter=',')

        for proteases in reader:
            uniqueProteases = uniqueProteases | set(proteases)
            sequenceAnnotations.append(proteases)

    # Generate column names for protease annotation data frame.
    cols = ['cluster'] + sorted(uniqueProteases)
    
    # generate binary protease matrix at sequence level.
    proteaseMatrix = []
    
    # Loop through data frame to extract cluster assignments.
    for i,row in data.iterrows():
        proteaseRow = []
        proteaseRow.append(row['clusterID'])

        # Loop through sorted unique proteases. Append one if in sequence annotation row, zero otherwise.
        for p in cols[1:]:
            if p in sequenceAnnotations[i]:
                proteaseRow.append(1)
            else:
                proteaseRow.append(0)

        # Add sequence annotation row to protease annotation binary matrix.
        proteaseMatrix.append(proteaseRow)

    # Create Pandas DataFrame object for protease annotation matrix.
    df = pd.DataFrame(proteaseMatrix,columns=cols)

    # Group sequence level protease matrix to cluster level protease matrix.
    groupedAnnotationMatrix = df.groupby('cluster').sum()

    # generate heatmap
    #sns.heatmap(groupedAnnotationMatrix,xticklabels=True,yticklabels=True)
    clusterMap = sns.clustermap(groupedAnnotationMatrix, method='average', metric='euclidean', row_cluster=True, col_cluster=True, row_colors=None)
    plt.setp(clusterMap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    clusterMap.fig.suptitle(title)
    plt.savefig(outputFileName,bbox_inches='tight',dpi=100)
    plt.show()

    return None