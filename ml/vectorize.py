'''
Created on 4. nov. 2014

@author: hakon
'''
import numpy

def candidates_to_array(candidates):
    
    all_candidates = []
    
    for candidate in candidates:
        
        
        aa = 0.0
        for el in candidate.mapped_sequences:
            name = el.data[1]
            aa += float(name.split("-")[1])
        
        a = len(candidate.mapped_sequences)  # not scaled

        b = candidate.hairpin_energy
        c = candidate.entropy_nucleotides
        d = candidate.entropy_structure
        e = candidate.heterogenity_5_begin
        f = candidate.heterogenity_5_end
        g = candidate.heterogenity_3_begin
        h = candidate.heterogenity_3_end
        i = candidate.quality
        
        
        features = [a,aa,b,c,d,e,f,g,h,i]
        
        all_candidates.append(features)
        
    return numpy.array(all_candidates)
