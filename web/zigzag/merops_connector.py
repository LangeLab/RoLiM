import re
import itertools

import pymysql.cursors
import pandas as pd
import numpy as np

from zigzag.config import *

# Ensure that Pandas always loads full sequence.
pd.set_option('max_colwidth', 1000000)

SITES = [
    "Site_P4",
    "Site_P3",
    "Site_P2",
    "Site_P1",
    "Site_P1prime",
    "Site_P2prime",
    "Site_P3prime",
    "Site_P4prime",
]

POSITIONS = [
    "p4",
    "p3",
    "p2",
    "p1",
    "p1'",
    "p2'",
    "p3'",
    "p4'",
]


def connect_to_merops_database():
    '''

    Parameters:
        None

    Returns:
        merops_connection -- 
    '''

    # Connect to MEROPS database.
    merops_connection = pymysql.connect(host='localhost',
                                    user='zigzag',
                                    password='zigzagpass',
                                    db='merops',
                                    cursorclass=pymysql.cursors.DictCursor)

    return merops_connection


def retrieve_substrates(names=[],ids=[],organisms=['Homo sapiens']):
    '''
    Connects to MEROPS, retrieves protease substrate sequences,
        and returns a Pandas DataFrame containing the sequences
        for each protease in the Substrate_search table.
    
    Parameters:
        names -- List. Default empty. Optional list of protease
                        names used to filter results.
        ids -- List. Default empty. Optional list of substrate IDs
                        used to filter results.
        organisms -- List. Default ['Homo sapiens']. Optional list
                            of organisms used to filter results.
                            Homo sapiens must be overridden with empty
                            list or alternative organisms.

    Returns:
        protease_substrates -- Dict. Maps protease name to Pandas 
                                    DataFrame for each protease with
                                    associated substrutes in MEROPS.
    '''

    sites = ["Site_P4",
                "Site_P3",
                "Site_P2",
                "Site_P1",
                "Site_P1prime",
                "Site_P2prime",
                "Site_P3prime",
                "Site_P4prime"]

    # Connect to MEROPS database.
    connection = connect_to_merops_database()

    # Get all substrates from MEROPS Substrate_search table.
    try:
        with connection.cursor() as cursor:
            conditional_string = ""
            if names != []:
                conditional_string += (
                    " AND Substrate_search.Protease IN ('"
                    + "','".join(names)
                    + "')"
                )
            if ids != []:
                conditional_string += (
                    " AND Substrate_search.cleavageID IN ('"
                    + "','".join(ids)
                    + "')"
                )
            if organisms != []:
                conditional_string += (
                    " AND Substrate_search.organism IN ('"
                    + "','".join(organisms)
                    + "')"
                ) 
            
            for site in sites:
                conditional_string += (
                    " AND %s" % site
                    + " IN ('"
                    + "' ,'".join(THREE_LETTER_CODES)
                    + "')"
                )
            
            query = (
                "SELECT Substrate_search.Protease, {0} FROM Substrate_search"
                + " INNER JOIN (SELECT code FROM activity_status WHERE organism"
                + " = 'human' AND status = 'active') as human_active"
                + " ON Substrate_search.code = human_active.code WHERE cleavage_type"
                + " NOT IN ('synthetic','theoretical'){1}"
            ).format(", Substrate_search.".join(sites), conditional_string)
            cursor.execute(query)
            all_substrates = cursor.fetchall()
    
    finally:
        connection.close()
    
    protease_substrates = partition_proteases(all_substrates,sites,POSITIONS)

    return protease_substrates  


def retrieve_vectorized_substrates(names=[]):
    '''
    Retrieves vectorized substrates from vectorized_substrates
        table in merops MySQL database.

    Parameters:
        names -- List. Default empty. Optional list of protease
                        names used to filter results.

    Returns:
        vectorized_substrates -- Dict.
    '''

    conditional_string = ""
    if names != []:
        conditional_string += " WHERE Protease in ('%s')" % ("', '".join(names))

    connection = connect_to_merops_database()
    try:
        with connection.cursor() as cursor:
            query = "SELECT * FROM vectorized_substrates%s" % (conditional_string)
            cursor.execute(query)
            all_substrates = cursor.fetchall()
    finally:
        connection.close()

    fields = []
    for site in SITES:
        for residue in SINGLE_LETTER_CODES:
            fields.append(site+"_"+residue)

    positions = []
    for position in POSITIONS:
        for residue in SINGLE_LETTER_CODES:
            positions.append(position+"_"+residue)

    vectorized_substrates = partition_proteases(all_substrates,fields,positions)

    return vectorized_substrates


def retrieve_protease_patterns(names=[], enable_compound_residues=True):
    '''
    Retrieve protease substrate patterns from
        protease_patterns table in merops database.
        Query optionally filtered by protease name.

    Parameters:
        names -- List. Default empty.
    Returns:
        protease_patterns -- Dict,
    '''

    connection = connect_to_merops_database()

    if enable_compound_residues:
        protease_pattern_table = 'protease_patterns'
    else:
        protease_pattern_table = 'protease_patterns_no_cr'

    if names != []:
        conditional_string = " WHERE Protease in ('%s')" % ("', '".join(names))
    else:
        conditional_string = ""

    protease_patterns = {}

    try:
        with connection.cursor() as cursor:
            query = "SELECT Protease, pattern FROM {0}{1}".format(
                protease_pattern_table,
                conditional_string
            )
            cursor.execute(query)
            while True:
                row = cursor.fetchone()
                if row:
                    if row['Protease'] in protease_patterns:
                        protease_patterns[row['Protease']].append(row['pattern'])
                    else:
                        protease_patterns[row['Protease']] = [row['pattern']]
                else:
                    break
    finally:
        connection.close()

    return protease_patterns


def partition_proteases(rows,fields,columns):
    '''

    Parameters:
        rows -- 
        fields --
        columns -- 
    Returns:
        partitioned_results -- 
    '''

    # Partition result set by protease.
    proteases = {}
    for row in rows:
        df_row = []
        for field in fields:
            df_row.append(row[field])
        protease = row['Protease']
        if protease in proteases.keys():
            proteases[protease].append(df_row)
        else:
            proteases[protease] = [df_row]

    # Generate Pandas DataFrame for each protease and append
    # to list of protease substrate DataFrames.
    partitioned_results = {}    
    for protease,df_rows in proteases.items():
        partitioned_results[protease] = pd.DataFrame(df_rows,columns=columns)

    return partitioned_results


def vectorize_substrates(substrates):
    '''
    Pulls substrates from merops database Substrate_search
        table, converts to single-letter amino acid encoding,
        inserts substrates into vectorized_substrates table.
    
    Parameters:
        substrates --

    Returns:
        None
    '''

    vector_template = pd.Series(
        np.zeros(len(THREE_LETTER_CODES), dtype=np.int8),
        index=THREE_LETTER_CODES
    )

    for protease,sequences in substrates.items():
        for i,sequence in sequences.iterrows():
            skip_sequence = False
            row = [protease.strip("'\"")]
            for residue in sequence:
                positional_vector = vector_template.copy()
                if residue in positional_vector.index:
                    positional_vector[residue] = 1
                else:
                    skip_sequence = True
                    break
                row += positional_vector.tolist()
            if skip_sequence:
                print('Skipping sequence')
                continue
            row = [str(i) for i in row]
            print(row)
            insert_vectorized_substrate(row)


def insert_vectorized_substrate(substrate):
    '''

    Parameters:
        substrate --
    
    Returns:
        None
    '''

    columns = ['Protease']
    for site in SITES:
        for residue in SINGLE_LETTER_CODES:
            columns.append(site+"_"+residue)

    connection = connect_to_merops_database()

    try:
        with connection.cursor() as cursor:
            query = (
                "INSERT INTO vectorized_substrates ({0}) VALUES ('{1}}',{2})"
            ).format(", ".join(columns), substrate[0], ", ".join(substrate[1:]))
            cursor.execute(query)
        connection.commit()
    finally:
        connection.close()


def create_vectorized_substrates_table():
    '''
    Creates table for vectorized substrates containing
        a positional column for each residue in backgroun
        data.

    Parameters:
        None

    Returns:
        None
    '''

    # Generate column name and data type string.
    columns = [
        "id int NOT NULL AUTO_INCREMENT",
        "Protease varchar(100) NOT NULL",
    ]
    for site in SITES:
        for residue in SINGLE_LETTER_CODES:
            columns.append(site + "_" + residue + " int NOT NULL")
    columns.append("PRIMARY KEY (id)")

    # Open merops database connection and create vectorized_substrates table.
    connection = connect_to_merops_database()
    try:
        with connection.cursor() as cursor:
            query = (
                "CREATE TABLE vectorized_substrates ({}) "
            ).format(", ".join(columns))
            cursor.execute(query)
        connection.commit()
    finally:
        connection.close() 
    print(query)


def retrieve_protease_family_code(protease):
    '''

    Parameters:
        protease -- String.

    Returns:
        protease_family_code -- String.
    '''

    connection = connect_to_merops_database()

    try:
        with connection.cursor() as cursor:
            query = (
                "SELECT DISTINCT code FROM Substrate_search"
                + " WHERE Protease = '{}'"
            ).format(protease)
            cursor.execute(query)
            protease_code = cursor.fetchone()["code"]
    finally:
        connection.close()

    return protease_code[:protease_code.find('.')]


def create_protease_patterns_table(enable_compound_residues=True):
    '''
    Creates protease_patterns table in merops database.

    Parameters:
        None
    
    Returns:
        None
    '''

    connection = connect_to_merops_database()

    try:
        with connection.cursor() as cursor:
            if enable_compound_residues:
                protease_pattern_table = 'protease_patterns'
            else:
                protease_pattern_table = 'protease_patterns_no_cr'
            query = (
                "CREATE TABLE {} (id int NOT NULL AUTO_INCREMENT,"
                + " Protease varchar(100) NOT NULL, pattern varchar(64) NOT NULL,"
                + " PRIMARY KEY (id))"
            ).format(protease_pattern_table)
            cursor.execute(query)
        connection.commit()
    finally:
        connection.close()


def insert_protease_patterns(protease_patterns, enable_compound_residues=True):
    '''

    Parameters:
        protease_patterns -- Dict.

    Returns:
        None
    '''

    for protease,patterns in protease_patterns.items():
        for pattern in patterns:
            connection = connect_to_merops_database()
            try:
                with connection.cursor() as cursor:
                    if enable_compound_residues:
                        protease_pattern_table = 'protease_patterns'
                    else:
                        protease_pattern_table = 'protease_patterns_no_cr'
                    query = (
                        "INSERT INTO {0} (Protease, pattern)"
                        + " VALUES ('{1}', '{2}')"
                    ).format(protease_pattern_table, protease, pattern)
                    cursor.execute(query)
                connection.commit()
            finally:
                connection.close()


def swap_protease_patterns(target):
    '''

    Parameters:

    Returns:

    '''

    connection = connect_to_merops_database()

    try:
        with connection.cursor() as cursor:
            cursor.execute("DROP TABLE protease_patterns")
            cursor.execute("CREATE TABLE protease_patterns LIKE protease_patterns_06")
            cursor.execute("INSERT INTO protease_patterns SELECT * FROM %s" % target)
        connection.commit()
    finally:
        connection.close()


def main():
    '''

    '''
    
    protease_substrates = retrieve_substrates()
    print(protease_substrates)


if __name__ == '__main__':
    main()