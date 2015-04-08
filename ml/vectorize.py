'''
Created on 4. nov. 2014

@author: hakon
'''
import numpy
import math

# feature_names = []
# feature_names.append("log mapped sequences")
# feature_names.append("mapped sequences rpm count")
# feature_names.append("hairpin energy")
# feature_names.append("hairpin + 10nt energy")
# feature_names.append("entropy nucleotides")
# feature_names.append("entropy structure")
# feature_names.append("heterogenity 5p start")
# feature_names.append("heterogenity 5p end")
# feature_names.append("heterogenity 3p start")
# feature_names.append("heterogenity 3p end")
# feature_names.append("candidate quality")
# feature_names.append("nr of bindings hairpin+10nt")
# feature_names.append("overhang level outer")
# feature_names.append("overhang outer")
# feature_names.append("overhang level inner")
# feature_names.append("overhang inner")
# feature_names.append("bulge factor")
# feature_names.append("small sequences aligned")
# feature_names.append("small sequences aligned 5p")
# feature_names.append("small sequences aligned 3p")

def one_null_bool(bool_value):
    return 1 if bool_value else 0

def candidates_to_array(candidates):
    
    all_candidates = []
    
    for candidate in candidates:
        
        aa = 0.0
        for el in candidate.mapped_sequences:
            name = el.data[1]
            aa += float(name.split("-")[1])
        
        a = len(candidate.mapped_sequences)
        a = math.log(a) if a > 0 else a
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
        
        o = candidate.bulge_factor # skip 123
        
#         p = candidate.small_subs
#         q = candidate.small_subs_5p
#         r = candidate.small_subs_3p

        p = candidate.has_short_seqs_5p
        pp = candidate.has_short_seqs_3p
        
        q = candidate.ratio_short_long_5p
        qq = candidate.short_seq_5p_stdev
        qqq = candidate.short_seq_5p_offset
        
        r = candidate.ratio_short_long_3p
        rr = candidate.short_seq_3p_stdev
        rrr = candidate.short_seq_3p_offset
        
        s = candidate.loop_size # noise?
        t = candidate.folds_5p
        u = candidate.folds_3p
        v = candidate.folds_before
        w = candidate.folds_after
        
        x = one_null_bool(candidate.has_hairpin_struct_5p)
        y = one_null_bool(candidate.has_hairpin_struct_3p)
        z = one_null_bool(candidate.has_hairpin_struct)
        
        
#         r = candidate.junction_pos_5
#         s = candidate.junction_pos_3
        
        
#         features = [a,aa,b,bb,c,d,e,f,g,h,i,j,k,l,m,n,o,q,r]
        features = [a,aa,b,bb,c,d,e,f,g,h,i,j,k,l,m,n,o,p,pp,q,qq,qqq,r,rr,rrr,s,t,u,v,w,x,y,z]
        
        all_candidates.append(features)
        
    return numpy.array(all_candidates)







