'''
Created on 25. mar. 2015

@author: hakon
'''



def hairpin_stats(candidates, mirna, hc_mirna):
    
    for candidate in candidates:
        
        entropy_dict = candidate.bitpair_probabilities
        fold = candidate.hairpin_fold_40
        
        begin_5p = candidate.pos_5p_begin - candidate.hairpin_start + 40
        end_5p = candidate.pos_5p_end - candidate.hairpin_start + 40
        begin_3p = candidate.pos_3p_begin - candidate.hairpin_start + 40
        end_3p = candidate.pos_3p_end - candidate.hairpin_start + 40
        match_5p = -100
        match_3p = -100
        
        
        if begin_5p > -1:
            offset_5pb = _align_distance(begin_5p,entropy_dict)
            est_3pe = begin_5p + offset_5pb
            
        if end_5p > -1:
            offset_5pe = _align_distance(end_5p,entropy_dict)
            est_3pb = end_5p + offset_5pe
            
        if begin_3p > -1:
            offset_3pb = _align_distance(begin_3p,entropy_dict)
            est_5pe = begin_3p + offset_3pb

        if end_3p > -1:        
            offset_3pe = _align_distance(end_3p,entropy_dict)
            est_5pb = end_3p + offset_3pe
            
        #TODO: fix stuff

        
        print begin_5p,match_5p
        print end_3p, match_3p
        print " "*40 + "-"* (end_5p - begin_5p)
        print fold
        print "-" * match_5p
        print
        print "+" * end_3p
        print fold
        print "+" * match_3p
        print candidate.mapped_sequences
        
        assert 0


def _match_pos(pos, entropy_dict):
    
    if pos in entropy_dict:

        v = list(entropy_dict[pos].values())
        if max(v) < 0.1: return -1
        p = list(entropy_dict[pos].keys())
        m = p[v.index(max(v))]
        print "\t", pos, m, entropy_dict[pos]
        return m
    
    return -1
    



def _align_distance(position, entropy_dict, seachlen=5):
    
    area = range(position-seachlen,position+seachlen+1)
    
    print
    print position
    print area
    print entropy_dict

    
    match_area = [_match_pos(x, entropy_dict) for x in area]
    print match_area
    
    
    
    match_offset = [ p + m - position*2 for p,m in zip(area, match_area) if m != -1]
    print "\toffset", match_offset
    
    assert len(match_offset)
    
    counts = [ (match_offset.count(p), p)  for p in set(match_offset)]
    
    print counts

    
    count, best = max(counts)
    print count, best
    
    return best




def _folds(fold_seq, outer, inner):
    
    
    if inner > outer: # 5p part
        fold_sign = "("
        unfold_sign = ")"
        start = outer
        end = inner
    else:
        fold_sign = ")"
        unfold_sign = "("
        start = inner
        end = outer
        
    fold_in = fold_seq.count(fold_sign, start, end)
    fold_out = fold_seq.count(unfold_sign, start, end)
    
    return fold_in-fold_out, fold_in, fold_out


    





    