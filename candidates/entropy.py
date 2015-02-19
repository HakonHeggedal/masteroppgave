'''
Created on 9. okt. 2014

@author: hakon
'''
import math

def entropy(candidates):
    print "\n\tcalculating simple entropy"
    
    for candidate in candidates:
        
        hairpin_entropy = _entropy(candidate.hairpin, 2)
        structure_entropy = _entropy(candidate.hairpin_fold, 4)
        
#         print hairpin_entropy
#         print structure_entropy
        candidate.set_entropy(hairpin_entropy, structure_entropy)


def _entropy(sequence, part_size = 2):
    ''' calculate the degree of variation between all parts of the sequence 
    '''
    parts = {}
    chars = set(sequence)
    
    for x in xrange(len(sequence) - part_size):
        part = sequence[x:x+part_size]
        parts[part] = 1 if part not in parts else parts[part]+1
    
    entropy_val = 0
    for value in parts.itervalues():
        entropy_val += value * math.log(value, 2)
        
    return entropy_val
    