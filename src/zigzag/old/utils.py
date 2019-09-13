import pandas as pd

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1e6)

def load_f1_scores(data):
	'''

	Parameters:

	Returns:

	'''

    f1_scores = pd.read_csv(data, sep=',', header=None)
    f1_scores.fillna(0, inplace=True)
    f1_scores.columns = ['sequence', 'f1_score']
    
    return f1_scores