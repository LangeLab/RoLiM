import os
import re

import numpy as np
import pandas as pd

def import_fasta_to_df(fasta_path):
    """
    THIS IS THE ONE I ULTIMATELY WANT TO USE!!!

    Import fasta to Pandas dataframe with ID column and sequence column.

    Parameters:
        fasta_path -- String. Path to fasta file.

    Returns:
		fasta_df -- Pandas DataFrame.
    """

    with open(fasta_path) as fasta:
        sequences = []
        for i, line in enumerate(fasta):
            if i == 0:
                sequence = ''
                if line[0] == '>':
                    sequence_id = line.rstrip()
                else:
                    sequence_id = 'ID not found'
                continue
            elif line[0] == '>':
                sequences.append([sequence_id, sequence])
                sequence_id = line.rstrip()
                sequence = ''
                continue
            else:
                sequence += line.rstrip()
        sequences.append(sequence)

    cols = ['id', 'sequence']
    fasta_df = pd.DataFrame(sequences, columns=cols)

    return fasta_df