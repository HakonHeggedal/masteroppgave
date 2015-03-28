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
        
        print fold
        
        assert candidate.has_5p or candidate.has_3p
        
        if candidate.has_5p:
            assert candidate.pos_5p_begin != -1
            assert candidate.pos_5p_end != -1
            offset_5pb = _align_distance(begin_5p,entropy_dict)
            est_3pe = begin_5p + offset_5pb

            offset_5pe = _align_distance(end_5p,entropy_dict)
            est_3pb = end_5p + offset_5pe
            
        if candidate.has_3p:
            assert candidate.pos_3p_begin != -1
            assert candidate.pos_3p_end != -1
            offset_3pb = _align_distance(begin_3p,entropy_dict)
            est_5pe = begin_3p + offset_3pb
     
            offset_3pe = _align_distance(end_3p,entropy_dict)
            est_5pb = end_3p + offset_3pe
            
        
        b5 = begin_5p if candidate.has_5p else est_5pb
        e5 = end_5p if candidate.has_5p else est_5pe
        b3 = begin_3p if candidate.has_3p else est_3pb
        e3 = end_3p if candidate.has_3p else est_3pe
        
        folds_5p, folds_in_5p, folds_out_5p = _folds(fold, b5, e5)
        
        folds_3p, folds_in_3p, folds_out_3p = _folds(fold, e3, b3)
        
        folds_before, folds_before_in, folds_before_out = _folds(fold, b5-15, b5)
        folds_after, folds_after_in, folds_after_out = _folds(fold, e3+15, e3)
        
        
        loop_size = b3 - e5
        
        candidate.loop_size = loop_size
        candidate.folds_5p = folds_5p
        candidate.folds_3p = folds_3p
        candidate.folds_before = folds_before
        candidate.folds_after = folds_after
        
        print
        print fold
        print candidate.has_5p or candidate.has_3p
        print b5, e5, b3, e3
        print folds_5p, folds_in_5p, folds_out_5p
        print folds_3p, folds_in_3p, folds_out_3p
        print folds_before, folds_before_in, folds_before_out
        print folds_after, folds_after_in, folds_after_out
        print loop_size
        print "----------------"
        
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
    
    # TODO: overhang part 
    
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




def _folds(fold, outer, inner):
    
    
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
    
    
    fold_seq = fold[start:end+1] #TODO: not sure, but probably including end
        
    fold_in = fold_seq.count(fold_sign)
    fold_out = fold_seq.count(unfold_sign)
    
    print start, end, fold_seq, fold_in, fold_out, fold_sign
    
    return fold_in-fold_out, fold_in, fold_out


    

a = "123"[1:0]
print a.count("1")


    