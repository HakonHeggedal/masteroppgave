'''
Created on 4. nov. 2014

@author: hakon
'''
import numpy

def candidates_to_array(candidates):
    
    all_candidates = []
    
    pos_to_feat = []
    pos_to_feat.append("mapped sequences")
    pos_to_feat.append("hairpin energy")
    pos_to_feat.append("hairpin energy")

    
    for candidate in candidates:
        
        aa = 0.0
        for el in candidate.mapped_sequences:
            name = el.data[1]
            aa += float(name.split("-")[1])
        
        a = len(candidate.mapped_sequences)  # not scaled
        b = candidate.hairpin_energy
        bb = candidate.hairpin_energy_10
        c = candidate.entropy_nucleotides
        d = candidate.entropy_structure
        e = candidate.heterogenity_5_begin
        f = candidate.heterogenity_5_end
        g = candidate.heterogenity_3_begin
        h = candidate.heterogenity_3_end
        i = candidate.quality
        j = candidate.bindings_max_10
        k = candidate.overhang_level_outer_10 
        l = candidate.overhang_outer_10 
        m = candidate.overhang_level_inner_10 
        n = candidate.overhang_inner_10
        o = candidate.small_subs
        p = candidate.small_subs_5p
        q = candidate.small_subs_3p
        
        
        features = [a,aa,b,bb,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q]
        
        all_candidates.append(features)
        
    return numpy.array(all_candidates)





#     pos_to_feat[1] = "hairpin energy"
#     pos_to_feat[2] = "hairpin energy"
#     pos_to_feat[3] = "hairpin energy"
#     pos_to_feat[4] = "hairpin energy"
#     pos_to_feat[5] = "hairpin energy"
#     pos_to_feat[6] = "hairpin energy"
#     pos_to_feat[7] = "hairpin energy"
#     pos_to_feat[8] = "hairpin energy"
#     pos_to_feat[9] = "hairpin energy"
#     pos_to_feat[10] = "hairpin energy"
#     pos_to_feat[11] = "hairpin energy"
#     pos_to_feat[12] = "hairpin energy"
#     pos_to_feat[13] = "hairpin energy"
#     pos_to_feat[14] = "hairpin energy"
#     pos_to_feat[15] = "hairpin energy"
#     pos_to_feat[16] = "hairpin energy"