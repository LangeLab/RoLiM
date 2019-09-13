import numpy as np
import pandas as pd
import csv
from zigzag.merops_connector import *


# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1000000)

def import_annotations(annotation_file):
	'''

	Parameters:
		annotation_file: String.

	Returns:
		annotation_list: List.
	'''

	annotation_list = []

	with open(annotation_file) as annotation_csv:
		annotation_reader = csv.reader(annotation_csv,delimiter=',')

		for sequence_annotation in annotation_reader:
			annotation_list.append(sequence_annotation)


	return annotation_list


def unique_protease_annotations(annotation_list):
	'''

	Parameters:
		annotation_list -- List.

	Returns:
		unique_protease_annotations -- List.
	'''

	unique_protease_annotations = sorted({protease for proteases in annotation_list \
													for protease in proteases})


	return unique_protease_annotations



def generate_binary_protease_annotation_matrix(sequence_annotations,unique_proteases):
	'''

	Parameters:
		sequence_annotations -- List.
		unique_proteases

	Returns:
		binary_annotation_matrix -- Pandas DataFrame.
	'''

	binary_annotation_matrix = []

	for sequence,proteases in enumerate(sequence_annotations):
		sequence_vector = []
		for protease in unique_proteases:
			if protease in proteases:
				sequence_vector.append(1)
			else:
				sequence_vector.append(0)
		binary_annotation_matrix.append(sequence_vector)


	binary_annotation_matrix = pd.DataFrame(binary_annotation_matrix, \
											columns=unique_proteases)


	return binary_annotation_matrix



def generate_binary_family_annotation_matrix(sequence_annotations,unique_proteases):
	'''

	Parameters:

	Returns:

	'''
	families = {}
	for protease in unique_proteases:
		try:
			family = retrieve_protease_family_code(protease)
		except:
			continue
		if family in families:
			families[family].append(protease)
		else:
			families[family] = [protease]

	binary_family_annotation_matrix = []
	for sequence,proteases in enumerate(sequence_annotations):
		sequence_family_annotations = []
		for family in sorted(families.keys()):
			family_present = False
			for protease in proteases:
				if protease in families[family]:
					family_present = True
			if family_present == True:
				sequence_family_annotations.append(1)
			else:
				sequence_family_annotations.append(0)
		binary_family_annotation_matrix.append(sequence_family_annotations)

	binary_family_annotation_matrix = pd.DataFrame(binary_family_annotation_matrix, \
													columns=sorted(families.keys()))		


	return binary_family_annotation_matrix



def main():
	'''
	'''

	annotation_file = input('Enter path to annotation file (one sequence per row, comma-delimited annoations):  ')

	annotations = import_annotations(annotation_file)
	unique_proteases = unique_protease_annotations(annotations)
	binary_annotation_matrix = generate_binary_annotation_matrix(annotations,unique_proteases)



if __name__ == '__main__':
	main()