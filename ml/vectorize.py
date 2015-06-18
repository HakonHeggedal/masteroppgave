'''
Created on 4. nov. 2014

@author: hakon
'''
import numpy
import math

FEATURE_NAMES = []
FEATURE_NAMES.append("mapped sequences")
# FEATURE_NAMES.append("mapped sequences rpm count")
FEATURE_NAMES.append("hairpin energy")
# FEATURE_NAMES.append("hairpin + 10nt energy")
FEATURE_NAMES.append("entropy nucleotides")
FEATURE_NAMES.append("entropy structure")
FEATURE_NAMES.append("heterogenity 5p start")
FEATURE_NAMES.append("heterogenity 5p end")
FEATURE_NAMES.append("heterogenity 3p start")
FEATURE_NAMES.append("heterogenity 3p end")
FEATURE_NAMES.append("candidate quality")
# FEATURE_NAMES.append("nr of bindings hairpin+10nt")
FEATURE_NAMES.append("overhang level outer")
FEATURE_NAMES.append("overhang outer")
FEATURE_NAMES.append("overhang level inner")
FEATURE_NAMES.append("overhang inner")
# FEATURE_NAMES.append("bulge factor")
FEATURE_NAMES.append("short seq align score")
FEATURE_NAMES.append("short vs long seqs")
# FEATURE_NAMES.append("small sequences aligned 5p")
# FEATURE_NAMES.append("small sequences aligned 3p")
FEATURE_NAMES.append("loop size")
FEATURE_NAMES.append("folds 5p")
FEATURE_NAMES.append("folds 3p")
FEATURE_NAMES.append("folds before")
FEATURE_NAMES.append("folds after")
FEATURE_NAMES.append("has hairpin structure")



def one_null_bool(bool_value):
    return 1 if bool_value else 0

def candidates_to_array(candidates):
    
    all_candidates = []
    
    for candidate in candidates:
        
        a = 0.0
        for el in candidate.mapped_sequences:
            name = el.data[1]
            a += float(name.split("-")[1])
            
        a = math.log(a + 1.0)
        
#         a = len(candidate.mapped_sequences) # nr of mapped sequences
#         a = math.log(a) if a > 0 else a        
#         b = candidate.hairpin_energy

        b = candidate.hairpin_energy_10
        
        c = candidate.entropy_nucleotides
        d = candidate.entropy_structure
        
        e = candidate.heterogenity_5_begin
        f = candidate.heterogenity_5_end
        g = candidate.heterogenity_3_begin
        h = candidate.heterogenity_3_end
        
        i = candidate.quality
        
#         j = candidate.bindings_max_10 # max number of bindings in one direction in the hp-struct
        k = candidate.overhang_level_outer_10 
        l = candidate.overhang_outer_10 
        m = candidate.overhang_level_inner_10 
        n = candidate.overhang_inner_10

#         p = one_null_bool(candidate.has_short_seqs_5p)
#         pp = one_null_bool(candidate.has_short_seqs_3p)
        
        p = candidate.short_seq_align
        q = math.log(candidate.ratio_short_long)
        
        s = candidate.loop_size # noise?
        t = candidate.folds_5p
        u = candidate.folds_3p
        v = candidate.folds_before
        w = candidate.folds_after
        
#         x = one_null_bool(candidate.has_hairpin_struct_5p)
#         y = one_null_bool(candidate.has_hairpin_struct_3p)
        z = one_null_bool(candidate.has_hairpin_struct)
        
        
#         features = [a,aa,b,bb,c,d,e,f,g,h,i,j,k,l,m,n,o,q,r]
#         features = [a,aa,b,bb,c,d,e,f,g,h,i,j,k,l,m,n,o,p,pp,q,s,t,u,v,w,x,y,z]
#         features = [aa,bb,c,d,e,f,g,h,i,k,l,m,n,p,pp,q,s,t,u,v,w,x,y,z]
#         features = [a,b,c,d,e,f,g,h,i,k,l,m,n,p,q,s,t,u,v,w,x,y,z]
        features = [a,b,c,d,e,f,g,h,i,k,l,m,n,p,q,s,t,u,v,w,z]
        
        all_candidates.append(features)
        
        assert len(features) == len(FEATURE_NAMES), (len(features), len(FEATURE_NAMES))
        
    return numpy.array(all_candidates)







