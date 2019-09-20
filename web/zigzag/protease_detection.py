import re
import itertools

import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.spatial import distance

from zigzag.calculate_background_frequencies import *
from zigzag.config import *
from zigzag.merops_connector import *
from zigzag.amino_acid_encoding_converter import *

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1000000)


"""
def generate_pattern_similarity_matrix(data,
                                        patterns,
                                        subs):
    '''

    Parameters:
        data -- Pandas DataFrame. One sequence per row.
                    One position per column.
        patterns -- List. Contains Pandas Series patterns.
        subs -- Pandas DataFrame. Contains pairwise substitution probabilities.

    Returns:

    '''
    
    pattern_similarity_matrix = []
    for i,sequence in data.iterrows():
        sequence_scores = []
        for pattern in patterns:
            position_scores = []
            for position,residue in pattern.iteritems():
                if residue != '.':
                    position_scores.append(subs[
                        sequence[position]][residue.strip('[]')])
            sequence_scores.append(np.mean(position_scores))
        pattern_similarity_matrix.append(sequence_scores)

    return pd.DataFrame(pattern_similarity_matrix)
"""

def pattern_to_regular_expression(patterns):
    '''
    Converts series-style patterns to regular expressions.

    Parameters:
        patterns -- List.
    
    Returns:
        regular_expressions -- List.
    '''

    regular_expressions = []

    for pattern in patterns:
        pattern_string = ''.join(pattern).strip('.')
        regular_expression = r"(?=("
        separation = 0
        for position in pattern_string:
            if position != '.':
                if separation != 0:
                    regular_expression += (r"(.){"
                        + re.escape(str(separation))+ r"}" 
                        + re.escape(position))
                    separation = 0
                else:
                    regular_expression += re.escape(position)
            else:
                separation += 1
        regular_expression += r"))"
        regular_expressions.append(regular_expression)
    
    return regular_expressions

"""
def scan_protease_substrates(proteases,
                                patterns,
                                pattern_proportion_cutoff,
                                pattern_p_value_cutoff):
    '''

    Parameters:
        proteases -- Dict.
        patterns -- List.
    Returns:
        protease_patterns -- Dict.

    '''

    protease_hits = {}

    # Scan protease substrates with sample patterns.
    for protease,substrates in proteases.items():
        substrate_array = np.array(substrates)
        for pattern in patterns:
            print('Testing:')
            print(protease)
            print(pattern)
            vectorized_pattern = vectorize_pattern(pattern,1)
            pattern_frequency = np.bincount(
                np.dot(substrate_array,vectorized_pattern),minlength=9)[8]
            if ((pattern_frequency
                    >= (pattern_proportion_cutoff*substrate_array.shape[0]))
                    and (pattern_frequency > 1)):
                print('HIT:  Frequency = %d, Percentage Frequency = %f'
                    % (pattern_frequency,
                        (pattern_frequency/substrate_array.shape[0])))
                
                if protease in protease_hits:
                    protease_hits[protease]['patterns'].append(
                        ','.join(pattern.tolist()))
                else:
                    protease_hits[protease] = {
                        'patterns':[','.join(pattern.tolist())],
                        'perfect_match':False
                    }
    
    # Flag perfect pattern intersections.
    protease_patterns = retrieve_protease_patterns()
    for protease,hits in protease_hits.items():
        try:
            if len(hits['patterns']) == len(protease_patterns[protease]):
                hits['perfect_match'] = True
        except KeyError:
            hits['perfect_match'] = 'No protease patterns found'

    return protease_hits
"""

"""
def pattern_intersection(sample_patterns, protease_patterns):
    '''

    Parameters:
        sample_patterns -- List. Contains sample patterns as
                                    Pandas Series objects.
        protease_patterns --  Dict. Contains protease substrate
                                    patterns as Pandas Series objects.

    Returns
        protease_hits --

    '''

    sample_pattern_set = {','.join((str(position)
                            for position in pattern))
                            for pattern in sample_patterns}
    protease_hits = {}
    for protease,patterns in protease_patterns.items():
        pattern_intersection = sample_pattern_set & set(patterns)
        if pattern_intersection:
            protease_hits[protease] = {
                'patterns': list(pattern_intersection)
            }

    return protease_hits
"""

"""
def vectorized_pattern_intersection(sample_patterns,
                                    protease_patterns,
                                    alphabet):
    '''

    Parameters:
        sample_patterns --
        protease_patterns --
        alphabet -- List. Ordered alphabet of single-letter codes.

    Returns:
        protease_hits --
    '''

    sample_pattern_vectors = [
        vectorize_pattern(pattern,0).reshape(
            len(POSITIONS),len(alphabet))
            for pattern in sample_patterns
    ]

    protease_hits = {}
    for protease,patterns in protease_patterns.items():
        protease_series_patterns = [
            pd.Series(pattern.strip('[]').split(','),
                index=sample_patterns[0].index) for pattern in patterns
        ]
        protease_pattern_vectors = [
            vectorize_pattern(pattern,0).reshape(len(POSITIONS),
            len(alphabet)) for pattern in protease_series_patterns
        ]
        for i,sample_pattern_vector in enumerate(sample_pattern_vectors):
            for protease_pattern_vector in protease_pattern_vectors:
                vector_product = protease_pattern_vector*sample_pattern_vector
                mask_negative = vector_product < 0
                vector_product[mask_negative] = 0
                if np.array_equal(protease_pattern_vector,sample_pattern_vector):
                #if np.array_equal(np.any(protease_pattern_vector,axis=1),np.any(vector_product,axis=1)):
                    if protease in protease_hits:
                        if (','.join(sample_patterns[i])
                                in protease_hits[protease]['patterns']):
                            continue
                        else:
                            protease_hits[protease]['patterns'].append(
                                ','.join(sample_patterns[i]))
                    else:
                        protease_hits[protease] = {
                            'patterns':[','.join(sample_patterns[i])]
                        }

    '''
    vectorized_sample_patterns = [vectorize_pattern(pattern,0) for pattern in sample_patterns]
    protease_hits = {}
    for protease,patterns in protease_patterns.items():
        protease_substrates = np.array(retrieve_vectorized_substrates(names=[protease])[protease],dtype='int8')
        split_patterns = [pd.Series(pattern.strip('[]').split(',')) for pattern in patterns]
        for pattern in split_patterns:
            sub_patterns = [pd.Series(i,index=sample_patterns[0].index) for i in itertools.product(*pattern)]
            for sub_pattern in sub_patterns:
                print(sub_pattern)
                vectorized_sub_pattern = vectorize_pattern(sub_pattern,1)
                if np.bincount(np.dot(protease_substrates,vectorized_sub_pattern),minlength=9)[8] > 0.33*len(protease_substrates):
                    for i,vectorized_sample_pattern in enumerate(vectorized_sample_patterns):
                        print(sample_patterns[i])
                        if np.dot(vectorized_sub_pattern,vectorized_sample_pattern) == sub_pattern[sub_pattern!='.'].size:
                            if protease in protease_hits:
                                protease_hits[protease]['patterns'].append(','.join(sub_pattern.tolist()))
                            else:
                                protease_hits[protease] = {'patterns':[','.join(sub_pattern.tolist())]}

    '''

    return protease_hits
"""

"""
def calculate_pattern_specific_positional_substitution_frequencies(protease,
                                                                    pattern,
                                                                    alphabet):
    '''

    Parameters:
        protease --
        pattern --
        alphabet -- List. Ordered alphabet of single-letter codes.
    Returns:
        pattern_specific_positional_substitution_frequencies --
    '''

    pattern = pd.Series([i.strip('[]') for i in pattern.split(',')])
    vectorized_pattern = vectorize_pattern(pattern,0).reshape(8,
        len(alphabet))
    protease_substrates = retrieve_substrates(names=[protease])
    substrates = convert_encoding(protease_substrates[protease],1)
    pattern.index = substrates.columns

    # Subset of substrates complementary to pattern.
    pattern_complement_index = []

    for i,substrate in substrates.mask(substrates!=pattern).fillna('.').iterrows():
        if not substrate.equals(pattern):
            pattern_complement_index.append(i)
    
    print(pattern_complement_index)
    pattern_complement = substrates.loc[pattern_complement_index]

    # Calculate position-specific substitution frequencies for protease substrate pattern.
    vectorized_pattern_complement_frequencies = np.sum(vectorize_sequences(
        pattern_complement),axis=0)
    pattern_specific_positional_substitution_frequencies = (np.divide(
        vectorized_pattern_complement_frequencies,len(pattern_complement),
        out=np.zeros_like(vectorized_pattern_complement_frequencies,
            dtype='float64'),
        where=vectorized_pattern_complement_frequencies!=0) + vectorized_pattern)
        * np.any(vectorized_pattern,axis=1).reshape(8,1)

    print(pattern_specific_positional_substitution_frequencies)

    return pattern_specific_positional_substitution_frequencies
"""

"""
def pattern_protease_distances(sample_patterns,
                                protease_patterns,
                                alphabet):
    '''

    Parameters:
        sample_patterns --
        protease_patterns --
        alphabet -- List. Ordered alphabet of single-letter codes.

    Returns:
        sample_pattern_protease_distances --

    '''

    sample_pattern_vectors = [
        vectorize_pattern(pattern,0).reshape(
            8,
            len(alphabet)) for pattern in sample_patterns
    ]
    
    sample_pattern_protease_distances = {}
    for i,sample_pattern_vector in enumerate(sample_pattern_vectors):
        pattern_protease_distances = {}
        for protease,patterns in protease_patterns.items():
            protease_series_patterns = [
                pd.Series(pattern.strip('[]').split(','),
                    index=sample_patterns[0].index) for pattern in patterns
            ]
            
            protease_pattern_vectors = [
                vectorize_pattern(pattern,0).reshape(8,len(alphabet))
                for pattern in protease_series_patterns
            ]
            
            pattern_distances = []
            for protease_pattern_vector in protease_pattern_vectors:
                pattern_distances.append(distance.euclidean(
                    sample_pattern_vector.flatten(),
                    protease_pattern_vector.flatten()))
            pattern_protease_distances.update({protease:np.mean(pattern_distances)})
        
        sample_pattern_protease_distances.update({
            ','.join(sample_patterns[i]):pattern_protease_distances})
            
    return sample_pattern_protease_distances
"""

"""
def maximum_likelihood_family(protease_distances):
    '''

    Parameters:
        protease_distances --

    Returns:
        pattern_families --

    '''

    pattern_families = {pattern:[] for pattern in protease_distances.keys()}
    for pattern,proteases in protease_distances.items():
        family_distances = {}
        for protease,distance in proteases.items():
            family = retrieve_protease_family_code(protease)
            if family in family_distances:
                family_distances[family].append(distance)
            else:
                family_distances[family] = [distance]
        min_mean_distance = 100
        for family,distances in family_distances.items():
            print(pattern)
            print(family + ':  ' + ', '.join([str(i) for i in distances]))
            mean_family_distance = np.mean(distances)
            if mean_family_distance < min_mean_distance:
                pattern_families[pattern] = [family]
                min_mean_distance = mean_family_distance
            elif mean_family_distance == min_mean_distance:
                pattern_families[pattern].append(family)

    return pattern_families
"""