def calculate_enrichments(data,
                            bg,
                            subs,
                            substitution_groups_enabled,
                            positions,
                            frequency_p_value_cutoff,
                            positional_p_value_cutoff):
    '''
    Generate substitution groups for each position using precomputed
        list of approved groups. Approved groups based on combined
        substitution probability exceeding 1/(number of goups).
    
    Parameters:
        data -- Pandas DataFrame. Contains sequences and metadata.
        bg --
        subs -- Pandas DataFrame. Contains possible substitution
                    groups. Serves as template DataFrame for enrichment
                    calculations.
        positions -- List. Contains column headers from sequences
                            data frame.
        p_value_cutoff --

    Returns:
        groups -- Dictionary. Contains Pandas DataFrames as items
                        mapped to positions as keys, each containing
                        substitution groups, effective frequencies,
                        and enrichment scores for a given position.
    '''

    # Initialize dictionary to contain scored groups.
    topPositionalGroups = {}
    # Loop through positions to calculate effective frequencies
    # and enrichment scores.
    for position in positions:
        # Calculate positional amino acid frequencies and prepare
        # list of positionally-present substitution groups.
        groups = subs.copy()
        frequencies = data.groupby(position).size().reset_index(name='count')
        if substitution_groups_enabled:
            positionalDoublets = [''.join(sorted(list(i))) for i in itertools.combinations(frequencies[position],2)]
            positionalTriplets = [''.join(sorted(list(i))) for i in itertools.combinations(frequencies[position],3)]
            positionalSinglets = [i for i in frequencies[position]]

            # Filter on possible doublets and triplets, doublets on
            # triplets, and singlets on doublets.
            positionalDoublets = [i for i in positionalDoublets if (i in [j for k in groups[groups['substitutionGroup'].str.len()==2].values.tolist() for j in k])]
            positionalTriplets = [i for i in positionalTriplets if (i in [j for k in groups[groups['substitutionGroup'].str.len()==3].values.tolist() for j in k])]
            doubletResidues = [j for i in positionalDoublets for j in i]
            tripletResidues = [j for i in positionalTriplets for j in i]
            positionalSinglets = [i for i in positionalSinglets if i not in doubletResidues]
            positionalDoublets = [i for i in positionalDoublets if (i[0] not in tripletResidues and i[1] not in tripletResidues)]

            # Combine triplets, doublets, and singlets into a single
            # list representing all singlets and substitution groups
            # present in position.
            positionalGroups = positionalTriplets + positionalDoublets + positionalSinglets

            # Filter possible substitution group data frame on substitution groups groups present in position.
            groups = groups[groups['substitutionGroup'].isin(positionalGroups)]
        else:
            groups = groups[groups['substitutionGroup'].isin(frequencies[position])]

        # Calculate effective group frequencies and enrichment scores.
        for i,residue in frequencies.iterrows():
            groups.loc[groups['substitutionGroup'].str.contains(residue[position]),'frequency'] += residue['count']
        
        groups = groups[groups['frequency']!=0]
        
        if (len(groups)>0):
            ## Calculate positional weighting term.
            positionalWeighting = np.mean(groups['frequency'])/len(data)
            
            for i,group in groups.iterrows():
                # binomial significance test (method used be motif-x)
                groups.loc[i,'binomial_p_value'] = st.binom_test(x=group['frequency'],n=len(data),p=group['backgroundFrequency'])

                if ( groups.loc[i,'binomial_p_value'] < frequency_p_value_cutoff ) and ( group['frequency'] > 1 ):
                    groups.loc[i,'enrichment'] = (group['frequency'] / (group['backgroundFrequency'] * len(data))) * positionalWeighting #/ groups.loc[i,'binomial_p_value']
                else:
                    groups.loc[i,'enrichment'] = 0

            if groups.loc[groups['enrichment'].idxmax(),'enrichment'] != 0:
                topPositionalGroups[position] = groups.loc[groups['enrichment'].idxmax()]

    print('****************TOP POSITIONAL GROUPS**************************')
    print(topPositionalGroups)
    if len(topPositionalGroups)>0:
        topPositionalGroups = pd.DataFrame.from_dict(topPositionalGroups,orient='index')
        topGroup = topPositionalGroups.loc[topPositionalGroups['enrichment'].idxmax(),]
        matchingSequences = data[topGroup.name][data[topGroup.name].isin([i for i in topGroup['substitutionGroup']])].index.values.tolist()
        mostEnriched = pd.DataFrame(topGroup).transpose()
        mostEnriched['matchingSequences'] = pd.Series([matchingSequences],index=mostEnriched.index)
        print("*******************************")
        print(mostEnriched)
        print("*******************************")

        return mostEnriched
    else:
        return pd.DataFrame()


def assign_sequences(data,
                        bg,
                        subs,
                        substitution_groups_enabled,
                        positions,
                        maxDepth,
                        currentDepth,
                        patternTemplate,
                        frequency_p_value_cutoff,
                        positional_p_value_cutoff):
    '''
    Recursively loops through unassigned sequences in input data
    set. Generates specificity patterns and assigns matching
    sequences to patterns.

    Parameters:
        data -- Pandas DataFrame. Contains sequences and metadata.
        bg --
        subs -- Pandas DataFrame. Contains possible substitution groups.
                                    Serves as template DataFrame for
                                    enrichment
            calculations.
        positions -- List. Contains column headers from sequences data
                            frame.
        maxDepth --
        currentDepth --
        patternTemplate --
        frequency_p_value_cutoff --
        positional_p_value_cutoff --

    Returns:
        newClusters --
        newAssignments --
        clusterPaterns --

    '''

    print('new level' + str(currentDepth))
    newAssignments = []
    kill = False
    clusterPatterns = []
    print(patternTemplate)

    while (len(newAssignments) < len(data)):
        unassigned = data.loc[~(data.index.isin(newAssignments))]
        parentEnrichment = calculate_enrichments(unassigned,bg,subs,substitution_groups_enabled,positions,frequency_p_value_cutoff,positional_p_value_cutoff)
        if not parentEnrichment.empty:
            parentPattern = patternTemplate.copy()
            parentPattern[parentEnrichment.index] = '[' + parentEnrichment['substitutionGroup'] + ']'
            print(parentEnrichment['substitutionGroup'])
            clusterPatterns.append(parentPattern)
            clusterPositions = [i for i in positions if i not in parentEnrichment.index.values.tolist()]
            print(clusterPositions)
            clusterDepth = currentDepth+1
            if (clusterDepth<maxDepth):
                childEnrichment,childSequences,childPatterns = assign_sequences(data[data.index.isin([j for i in parentEnrichment['matchingSequences'] for j in i]) & ~(data.index.isin(newAssignments))],bg,subs,substitution_groups_enabled,clusterPositions,maxDepth,clusterDepth,parentPattern,frequency_p_value_cutoff,positional_p_value_cutoff)
                if not childEnrichment.empty:
                    print('appending child')
                    parentEnrichment = pd.concat([parentEnrichment,childEnrichment])
                    clusterPatterns += childPatterns
            try:
                newClusters = pd.concat([newClusters,parentEnrichment])
            except:
                newClusters = parentEnrichment.copy()
            newAssignments = sorted(list(set(newAssignments)|set([j for i in parentEnrichment['matchingSequences'] for j in i])))
            print('new assignments: ' + str(newAssignments))
            print(newClusters)
        else:
            break


    try:
        return newClusters,newAssignments,clusterPatterns
    except:
        return pd.DataFrame(),[],[]


def score_sequences(data,patterns,subs,groups,mean_substitution_scores):
    '''

    Parameters:
        data --
        patterns --
        subs --
        groups --
        mean_substitution_scores --

    Returns:
        scoringMatrix --
    '''

    scoringMatrix = []

    ### Loop through sequences, building cluster similarity score vector for each sequence.
    for i,sequence in data.iterrows():
        sequenceRow = []
        ## Loop through clusters to calculate similarity score for sequence.
        for pattern in patterns:
            #fits_cluster = True
            clusterScores = []
            ## Loop through positions in cluster to calculate substitution probability for sequence AA and positional substitution group.
            
            for position,residues in pattern.iteritems():
                ## Skip position if no substitution group present for position in cluster.
                if (residues != '.'):
                    '''
                    if (sequence[position] not in residues):
                        fits_cluster = False
                        break
                    '''
                    positionScores = []
                    ## Assign score of 1 if sequence AA matches substitution group. Continue to next position in cluster.
                    if (sequence[position] in residues):
                        clusterScores.append(1)
                        continue
                    ## Loop through positional substitution group if sequence AA does not match. Calculate mean substitution probability.
                    positional_residues = [j for i in groups.loc[groups['substitutionGroup']==residues[1:-1],'definingResidue'] for j in i]
                    for residue in positional_residues:
                        probability = subs[sequence[position]][residue] #- mean_substitution_scores[sequence[position]]
                        '''
                        if probability >= 0.5:
                            positionScores.append(1)
                        else:
                            positionScores.append(0)
                        '''
                        positionScores.append(probability)

                    positionScore = np.mean(positionScores)#/float(groups.loc[groups['substitutionGroup']==residues[1:-1],'backgroundFrequency'])
                    clusterScores.append(positionScore)
                    
            # Append mean of positional scores to sequence cluster score vector.
            #sequenceRow.append(1 if np.mean(clusterScores) == 1 else 0.000001)
            sequenceRow.append(np.mean(clusterScores))
            '''
            if fits_cluster:
                sequenceRow.append(1)
            else:
                sequenceRow.append(0)
        # Append sequence vector to scoring matrix.
        scoringMatrix.append(sequenceRow)
    '''

    return scoringMatrix