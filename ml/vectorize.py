'''
Created on 4. nov. 2014

@author: hakon
'''
import numpy

def candidates_to_array(candidates):
    
    all_candidates = []
    
    for candidate in candidates:
        
        a = len(candidate.mapped_sequences)
        b = candidate.hairpin_energy
        c = candidate.entropy_nucleotides
        d = candidate.entropy_structure
#         e = candidate.heterogenity_5
#         f = candidate.heterogenity_3
        g = candidate.quality
        
        
        features = [a,b,c,d,g]
        
        all_candidates.append(features)
        
    return numpy.array(all_candidates)
