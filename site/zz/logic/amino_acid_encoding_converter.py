import pandas as pd

########################################################################
# Index-matched single-letter and three-letter amino acid alphabets.
########################################################################
SINGLE_LETTER_CODES = [
    'A','R','N','D','B','C','E','Q','Z','G','H','I','L','K',
    'M','F','P','S','T','W','Y','V',
]

THREE_LETTER_CODES = [
    'Ala','Arg','Asn','Asp','Asx','Cys','Glu','Gln','Glx',
    'Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser',
    'Thr','Trp','Tyr','Val',
]
########################################################################

def convert_encoding(sequences, encode):
    '''
    Converts amino acid from single-letter to three-letter encoding and
        vice versa. Takes and returns a Pandas DataFrame.

    Parameters:
        sequences -- Pandas DataFrame. Contains sequences in three
                        letter code.
        encode -- Int. Either 1 or 3 based on desired output encoding.
    Returns:
        converted -- Pandas DataFrame. Contains sequences from input
                        data from converted from three-letter to single-
                        letter encoding.
    '''

    single_to_three = {
        single:THREE_LETTER_CODES[i]
        for i,single in enumerate(SINGLE_LETTER_CODES)
    }
    three_to_single = {
        three:SINGLE_LETTER_CODES[i]
        for i,three in enumerate(THREE_LETTER_CODES)
    }


    if encode == 3:
        converted = sequences.applymap(lambda x: single_to_three[x])
    elif encode == 1:
        converted = sequences.applymap(lambda x: three_to_single[x.capitalize()])
    else:
        print('Please enter a valid output encoding format')    

    return converted