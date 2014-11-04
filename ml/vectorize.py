'''
Created on 4. nov. 2014

@author: hakon
'''
import numpy

def candidates_to_array(candidates):
    
    candidate_array = numpy.array()
    
    for candidate in candidates:
        
        a = len(candidate.mapped_sequences)
        b = candidate.hairpin_energy
        c = candidate.entropy_nucleotides
        d = candidate.entropy_structure
        e = candidate.heterogenity_5
        f = candidate.heterogenity_3
        g = candidate.quality
        
        
        
        
        
        
        features = numpy.array(a,b,c,d,e,f,g)
        numpy.concatenate(candidate_array, features)
        
    return candidate_array