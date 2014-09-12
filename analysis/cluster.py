
# some functions to compute connectivity etc.

def similarity(matrix, rnas, known):
    ''' computes similarity between new rna and known rna
     input: similarity matrix: list of dicts, set of known rna(indexes)
     output: dict (rna -> similarity value)
    '''
    result = {}
    for rna, similar in zip(rnas, matrix):
        if rna in known:
            continue
        score = sum([sim for r, sim in similar.iteritems() if r in known])
        result[rna] = score
        
    return result


def totalclustering(matrix, rnas):
    ''' sum of similarity to neighbors in similarity matrix'''
    result = {}
    for rna, similar in zip(rnas, matrix):
        result[rna] = sum(similar.itervalues())
    
    return result
        
#     result = {rna: sum(similar.itervalues()) for rna, similar in zip(rnas, matrix)}
        
    