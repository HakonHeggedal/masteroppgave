'''
Created on 25. mar. 2015

@author: hakon
'''

search_outside = 0
search_inside = 4

def hairpin_stats(candidates, mirna, hc_mirna):
    
    for candidate in candidates:
        
        entropy_dict = candidate.bitpair_probabilities
        fold = candidate.hairpin_fold_40
        
        begin_5p = candidate.pos_5p_begin - candidate.hairpin_start + 40
        end_5p = candidate.pos_5p_end - candidate.hairpin_start + 40
        begin_3p = candidate.pos_3p_begin - candidate.hairpin_start + 40
        end_3p = candidate.pos_3p_end - candidate.hairpin_start + 40
        
#         print "\n-----------------------"
#         
#         hp = begin_3p - end_5p
#         print fold
#         print " "*begin_5p+fold[begin_5p:end_5p] + " "*hp+ fold[begin_3p:end_3p]
        
        if not candidate.has_5p and not candidate.has_3p:
            # no mature seqs -> not possible to find hairpin
            print candidate.pos_5p_begin, candidate.pos_5p_end
            print candidate.pos_3p_begin, candidate.pos_3p_end
            print candidate.hairpin_start, candidate.hairpin_end
            print candidate.mapped_sequences
            print candidate.candidate_type, candidate.miRNAid, candidate.miRNAid in hc_mirna
            print begin_5p, end_5p, begin_3p, end_3p
            continue
        
        if candidate.has_5p:
            assert candidate.pos_5p_begin != -1
            assert candidate.pos_5p_end != -1
            
            assert begin_5p > -5
            assert end_5p > 0
            
#             print "5p fold\t",fold[begin_5p:end_5p]
            
            offset_5pb, dist_5pb = _align_distance(begin_5p, entropy_dict, search_outside, search_inside)
            est_3pe = begin_5p + offset_5pb

            offset_5pe, dist_5pe = _align_distance(end_5p, entropy_dict, search_inside, search_outside)
            est_3pb = end_5p + offset_5pe
            
            folds_5p, folds_in_5p, folds_out_5p = _folds(fold, begin_5p, end_5p)
            folds_before, folds_before_in, folds_before_out = _folds(fold, begin_5p-15, begin_5p)
            
            candidate.folds_5p = folds_5p
            candidate.folds_before = folds_before
            
            
            
            # is the estimated 3p seq of apx. same length as the 5p seq ? 
            if offset_5pb != -1000 or offset_5pe != -1000:
                if est_3pe > est_3pb > end_5p:
                    mature_len = end_5p - begin_5p
                    est_len = est_3pe - est_3pb
                    if abs(mature_len - est_len) < 10:
                        candidate.has_hairpin_struct_5p = True
#                     print "mature vs est", mature_len, est_len, abs(mature_len - est_len), candidate.has_hairpin_struct_5p
            
            
        if candidate.has_3p:
            assert candidate.pos_3p_begin != -1
            assert candidate.pos_3p_end != -1
            
            assert end_3p > begin_3p > 0
            
            
#             print "3p fold\t",fold[begin_3p:end_3p]
            
            
            offset_3pb, dist_3pb = _align_distance(begin_3p, entropy_dict, search_outside, search_inside)
            est_5pe = begin_3p + offset_3pb
     
            offset_3pe, dist_3pe = _align_distance(end_3p, entropy_dict, search_inside, search_outside)
            est_5pb = end_3p + offset_3pe
            
            folds_3p, folds_in_3p, folds_out_3p = _folds(fold, end_3p, begin_3p)
            folds_after, folds_after_in, folds_after_out = _folds(fold, end_3p+15, end_3p)
            
            candidate.folds_3p = folds_3p
            candidate.folds_after = folds_after
            
            
            if offset_3pb != -1000 and offset_3pe != -1000:
                if est_5pb < est_5pe < begin_3p:
                    mature_len = end_3p - begin_3p
                    est_len = est_5pe - est_5pb
                    if abs(mature_len - est_len) < 10:
                        candidate.has_hairpin_struct_3p = True

#                     print "mature vs est", mature_len, est_len, abs(mature_len - est_len), candidate.has_hairpin_struct_3p
                    
        

        if (not candidate.has_hairpin_struct_3p) and (not candidate.has_hairpin_struct_5p):
            assert not candidate.has_hairpin_struct_5p
            assert not candidate.has_hairpin_struct_3p
            print "no hairpin struct"
            continue
        
        
        if candidate.has_5p and candidate.has_3p:
            
            assert candidate.has_hairpin_struct_3p or candidate.has_hairpin_struct_5p, "hp struct"
            #===================================================================
            # calculate overhang between 5p and 3p mature seqs. 
            # When 5p and 3p does not align with each other,
            #    the fold with the highest sum of binding probabilities is chosen.
            # 
            # makes an estimate for fold probability: sum(binding)/ sum(all bindings)
            #===================================================================
            hp = begin_3p - end_5p
            print "\nOverhang"
            print fold
            print " "*begin_5p+fold[begin_5p:end_5p] + " "*hp+ fold[begin_3p:end_3p]
            print end_5p - begin_5p, end_3p - begin_3p
#             print sorted(candidate.mapped_sequences)
#             print 
              
#             assert candidate.has_hairpin_struct_5p
#             assert candidate.has_hairpin_struct_3p
              
            overhang_outer_5p = est_5pb - begin_5p 
            overhang_outer_3p = est_3pe - end_3p
              
            overhang_inner_5p = est_5pe - end_5p
            overhang_inner_3p = est_3pb - begin_3p

#             sum_outer = sum(dist_5pb.values()) + sum(dist_3pe.values())
#             sum_inner = sum(dist_5pe.values()) + sum(dist_3pb.values())
            
#             print "sum prob. outer:", sum_outer
#             print "sum prob. inner:", sum_inner
            
            # set outer overhang 
            if overhang_outer_5p == overhang_outer_3p:
#                 sum_pos = dist_5pb[offset_5pb] + dist_3pe[offset_3pe]
#                 confidence_outer = sum_pos / sum_outer
                  
                candidate.overhang_outer = overhang_outer_5p
#                 candidate.overhang_outer_conf = confidence_outer
            else:
                offset_diff = overhang_outer_3p - overhang_outer_5p
                  
                # offset predicted by the 5p side
                pos1 = dist_5pb[offset_5pb] if dist_5pb else 0.0
  
                pos1_offset = offset_3pe + offset_diff
                
                pos1_extra = 0
                if dist_3pe:
                    pos1_extra = dist_3pe[pos1_offset] if pos1_offset in dist_3pe else 0.0
                  
                pos1_sum = pos1 + pos1_extra
                  
                # offset predicted by the 3p side
                pos2 = dist_3pe[offset_3pe] if dist_3pe else 0.0
                  
                pos2_offset = offset_5pb - offset_diff
                pos2_extra = 0
                if dist_5pb:
                    pos2_extra = dist_5pb[pos2_offset] if pos2_offset in dist_5pb else 0.0
                pos2_sum = pos2 + pos2_extra
                  
                print "\tmultiple overhang positions outer"
                print "\t\t", pos1_sum, pos1, pos1_extra, pos1_offset, overhang_outer_5p
                print "\t\t", pos2_sum, pos2, pos2_extra, pos2_offset, overhang_outer_3p

                  
                # 5p or 3p fold most probable ?
                if pos1_sum > pos2_sum:
                    candidate.overhang_outer = overhang_outer_5p
                    # 5p fold is most probable
#                     print "\t5p outer"
#                     confidence_outer = pos1_sum / sum_outer
#                     candidate.overhang_outer_conf = confidence_outer
                else:
                    candidate.overhang_outer = overhang_outer_3p
#                     print "\t3p outer"
#                     confidence_outer = pos2_sum / sum_outer
#                     candidate.overhang_outer_conf = confidence_outer
                
#                 assert abs(candidate.overhang_outer) < 10 or abs(offset_diff) < 5
                
                


            # inner overhang
            if overhang_inner_5p == overhang_inner_3p:
#                 sum_pos = dist_5pe[offset_5pe] + dist_3pb[offset_3pb]
#                 confidence_inner = sum_pos / sum_inner
                  
                candidate.overhang_inner = overhang_inner_5p
#                 candidate.overhang_inner_conf = confidence_inner

            else:
                offset_diff = overhang_inner_3p - overhang_inner_5p
                  
                # 5p prediction
                pos1 = dist_5pe[offset_5pe] if dist_5pe else 0
                  
                pos1_offset = offset_3pb + offset_diff
                pos1_extra = 0
                if dist_3pb:
                    pos1_extra = dist_3pb[pos1_offset] if pos1_offset in dist_3pb else 0.0
                pos1_sum = pos1 + pos1_extra
                  
                  
                # 3p prediction
                pos2 = dist_3pb[offset_3pb] if dist_3pb else 0
                  
                pos2_offset = offset_5pe - offset_diff
                pos2_extra = 0
                if dist_5pe:
                    pos2_extra = dist_5pe[pos2_offset] if pos2_offset in dist_5pe else 0.0
                pos2_sum = pos2 + pos2_extra
                  
                  
                print "\tmultiple overhang positions inner"
                print "\t\t", pos1_sum, pos1, pos1_extra, pos1_offset, overhang_inner_5p
                print "\t\t", pos2_sum, pos2, pos2_extra, pos2_offset, overhang_inner_3p
                
                assert pos2_sum or pos1_sum
                
                if pos1_sum > pos2_sum:
                    candidate.overhang_inner = overhang_inner_5p
#                     print "5p inner"
#                     confidence_inner = pos1_sum / sum_inner
#                     candidate.overhang_inner_conf = confidence_inner
                else:
                    candidate.overhang_inner = overhang_inner_3p
#                     print "3p inner"
#                     confidence_inner = pos2_sum / sum_inner
#                     candidate.overhang_inner_conf = confidence_inner

#                 assert abs(candidate.overhang_inner) < 10 or abs(offset_diff) < 5
            
            
            
            # overhang direction depends on strand ... hotfix:
            if candidate.chromosome_direction == "-":
                candidate.overhang_inner = - candidate.overhang_inner
                candidate.overhang_outer = - candidate.overhang_outer

            
            
            if candidate.overhang_inner < -15:
                candidate.overhang_inner = - 15
            elif candidate.overhang_inner > 15:
                candidate.overhang_inner = 15
           
           
            if candidate.overhang_outer < -15:
                candidate.overhang_outer = - 15
            elif candidate.overhang_outer > 15:
                candidate.overhang_outer = 15
                
            

            
            print "outer:"
            print candidate.overhang_outer, candidate.overhang_outer_10,
            print overhang_outer_5p, overhang_outer_3p
            print "inner:"
            print candidate.overhang_inner, candidate.overhang_inner_10,
            print overhang_inner_5p, overhang_inner_3p
            print (est_5pb, begin_5p), (est_5pe, end_5p),
            print (est_3pb, begin_3p), (est_3pe, end_3p)
        
            assert -15 <= candidate.overhang_outer <= 15, candidate.overhang_outer
            assert -15 <= candidate.overhang_inner <= 15, candidate.overhang_inner
        
        assert candidate.has_hairpin_struct_5p or candidate.has_hairpin_struct_3p
#         assert candidate.has_5p == candidate.has_hairpin_struct_5p
#         assert candidate.has_3p == candidate.has_hairpin_struct_3p
        candidate.has_hairpin_struct = True
        
        b5 = begin_5p
        e5 = end_5p
        b3 = begin_3p
        e3 = end_3p
        
        # choose positions for loop size based on best information available
        if not candidate.has_hairpin_struct_5p:
            if candidate.has_hairpin_struct_3p:
                b5 = begin_5p if abs(est_5pb - begin_5p) < 10 else est_5pb
                e5 = end_5p if abs(est_5pe - end_5p) < 10 else est_5pe
            else:
                assert False, "wrong hairpin struct"
            
        if not candidate.has_hairpin_struct_3p:
            if candidate.has_hairpin_struct_5p:
                b3 = begin_3p if abs(est_3pb - begin_3p) < 10 else est_3pb
                e3 = end_3p if abs(est_3pe - end_3p) < 10 else est_3pe
            else:
                assert False, "wrong hairpin struct"
                
        

        
#         print begin_5p, end_5p, begin_3p, end_3p
#         print b5, e5, b3, e3
        
        assert b5 < e5 <= b3 < e3
        
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
        
        
#         print "stats:"
#         print b5, e5, b3, e3
#         print folds_before, folds_5p, loop_size, folds_3p, folds_after





def _match_pos(pos, align_pos, entropy_dict):
    
    distance = abs(pos - align_pos)
    relevance = 1.0 - (0.2 * distance)
    
    if pos in entropy_dict:

        scores = list(entropy_dict[pos].values())
        if max(scores) < 0.1:
#             print "\t", pos, max(entropy_dict[pos]), entropy_dict[pos]
            return -1, 0
        
        best_score = max(scores)
        positions = list(entropy_dict[pos].keys())
        best_match_pos = positions[scores.index(best_score)]
        

        score = best_score * relevance
        
#         print "\t", pos, best_match_pos, score, best_score #, entropy_dict[pos]
        return best_match_pos, score
    
    return -1, 0


    

def _align_distance(align_pos, entropy_dict, seach_before, search_after):
    
    
    # area to align against other strand
    area = range(align_pos-seach_before,align_pos+search_after+1)
    
    assert len(area) < 20, area[:20]
#     print align_pos


    # mapped positions
    match_area = [_match_pos(pos, align_pos, entropy_dict) for pos in area]
    
#     print zip(area, match_area)
    
    # sum scores for relative distances to the other strand
    offset_scores = {}
    value_sum = 0
    
    for pos, (mapped_pos, mapped_val) in zip(area, match_area):

        
        if mapped_pos != -1:
            
            offset = pos + mapped_pos - align_pos*2
            
            if offset in offset_scores:
                offset_scores[offset] += mapped_val
            else:
                offset_scores[offset] = mapped_val
            value_sum += mapped_val
    
    if not len(offset_scores):
        return -1000, None
    
    best_offset = max(offset_scores, key=offset_scores.get)
    
#     print "\t", best_offset, offset_scores[best_offset], offset_scores
    
    return best_offset, offset_scores

    





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
    
#     print start, end, fold_seq, fold_in, fold_out, fold_sign
    
    return fold_in-fold_out, fold_in, fold_out


    


# 
# def overhang(pos_5p, pos_3p, entropy_dict, search_before, search_after):
#      
#     pos_5p_area = range(pos_5p-search_before, pos_5p+search_after+1)
#     pos_3p_area = range(pos_3p-search_after, pos_3p+search_before+1)
#     
#     match_5p = [_match_pos_limit(pos_5p, align_pos, limit, limit_pos, entropy_dict):]
#     
#     
    

# def _match_pos_limit(pos, align_pos, limit, limit_pos, entropy_dict):
#     
#     distance = abs(pos - align_pos)
#     relevance = 1.0 - (0.2 * distance)
#     
#     def _scale(p):
#         dist = abs(p - limit_pos)
#         scaled = 1 - (0.2 *  dist)
#         scaled = 0.0 if scaled < 0.0 else scaled
#         return scaled
#     
#     if pos in entropy_dict:
#         
#         scores = [(v * _scale(k), k, v) for k,v in entropy_dict[pos] if k in limit]
#         
#         if not len(scores):
#             return -1, 0
# 
#         score, best_match_pos, score_unscaled = max(scores)
#         
#         if score_unscaled < 0.1:
# #             print "\t", pos, max(entropy_dict[pos]), entropy_dict[pos]
#             return -1, 0
# 
#         
#         
#         score = max(scores) * relevance
#         
# #         print "\t", pos, best_match_pos, max(scores) #, entropy_dict[pos]
#         return best_match_pos, score
#     
#     return -1, 0

#         b5 = begin_5p if candidate.has_hairpin_struct_5p else est_5pb
#         e5 = end_5p if candidate.has_hairpin_struct_5p else est_5pe
#         b3 = begin_3p if candidate.has_hairpin_struct_3p else est_3pb
#         e3 = end_3p if candidate.has_hairpin_struct_3p else est_3pe

#         if candidate.has_5p and candidate.has_3p:
#             #===================================================================
#             # calculate overhang between 5p and 3p mature seqs. 
#             # When 5p and 3p does not align with each other,
#             #    the fold with the highest sum of binding probabilities is chosen.
#             # 
#             # makes an estimate for fold probability: sum(binding)/ sum(all bindings)
#             #===================================================================
#              
#             print "\nOverhang"
#             print candidate.pos_5p_begin
#             print sorted(candidate.mapped_sequences)
#             print 
#              
#             assert candidate.has_hairpin_struct_5p == candidate.has_hairpin_struct_3p == True
#              
#             overhang_outer_5p = est_5pb - begin_5p 
#             overhang_outer_3p = est_3pe - end_3p  
#              
#             sum_outer = sum(dist_5pb.values()) + sum(dist_3pe.values())
#             print "sum prob. outer:", sum_outer
#              
#             # set outer overhang 
#             if overhang_outer_5p == overhang_outer_3p:
#                 sum_pos = dist_5pb[offset_5pb] + dist_3pe[offset_3pe]
#                 confidence_outer = sum_pos / sum_outer
#                  
#                 candidate.overhang_outer = overhang_outer_5p
#                 candidate.overhang_outer_conf = confidence_outer
#                  
#             else:
#                  
#                 offset_diff = overhang_outer_3p - overhang_outer_5p
#                  
#                 # offset predicted by the 5p side
#                 pos1 = dist_5pb[offset_5pb]
#  
#                 pos1_offset = offset_3pe + offset_diff
#                 pos1_extra = dist_3pe[pos1_offset] if pos1_offset in dist_3pe else 0.0
#                  
#                 pos1_sum = pos1 + pos1_extra
#                  
#                 # offset predicted by the 3p side
#                 pos2 = dist_3pe[offset_3pe]
#                  
#                 pos2_offset = offset_5pb - offset_diff
#                 pos2_extra = dist_5pb[pos2_offset] if pos2_offset in dist_5pb else 0.0
#                 pos2_sum = pos2 + pos2_extra
#                  
#                 print "multiple overhang positions outer"
#                 print pos1_sum, pos1, pos1_extra, pos1_offset, overhang_outer_5p
#                 print pos2_sum, pos2, pos2_extra, pos2_offset, overhang_outer_3p
#                  
#                 # 5p or 3p fold most probable ?
#                 if pos1_sum > pos2_sum:
#                     # 5p fold is most probable
#                     print "\t5p outer"
#                     confidence_outer = pos1_sum / sum_outer
#                     candidate.overhang_outer = overhang_outer_5p
#                     candidate.overhang_outer_conf = confidence_outer
#                 else:
#                     print "\t3p outer"
#                     confidence_outer = pos2_sum / sum_outer
#                     candidate.overhang_outer = overhang_outer_3p
#                     candidate.overhang_outer_conf = confidence_outer
#              
#                 print "outer overhang", candidate.overhang_outer
#                 print  "confidence", candidate.overhang_outer_conf
#                 print
#                 assert 0
#              
#             overhang_inner_5p = est_5pe - end_5p 
#             overhang_inner_3p = est_3pb - begin_3p
#              
#             sum_inner = sum(dist_5pe.values()) + sum(dist_3pb.values())
#             print "sum prob. inner:", sum_inner
#              
#              
#             # inner overhang
#             if overhang_inner_5p == overhang_inner_3p:
#                 sum_pos = dist_5pe[offset_5pe] + dist_3pb[offset_3pb]
#                 confidence_inner = sum_pos / sum_inner
#                  
#                 candidate.overhang_inner = overhang_inner_5p
#                 candidate.overhang_inner_conf = confidence_inner
#             else:
#                  
#                 offset_diff = overhang_inner_3p - overhang_inner_5p
#                  
#                 # 5p prediction
#                 pos1 = dist_5pe[offset_5pe]
#                  
#                 pos1_offset = offset_3pb + offset_diff
#                 pos1_extra = dist_3pb[pos1_offset] if pos1_offset in dist_3pb else 0.0
#                 pos1_sum = pos1 + pos1_extra
#                  
#                  
#                 # 3p prediction
#                 pos2 = dist_3pb[offset_3pb]
#                  
#                 pos2_offset = offset_5pe - offset_diff
#                 pos2_extra = dist_5pe[pos2_offset] if pos2_offset in dist_5pe else 0.0
#                 pos2_sum = pos2 + pos2_extra
#                  
#                  
#                 print "multiple overhang positions inner"
#                 print pos1_sum, pos1, pos1_extra, pos1_offset, overhang_inner_5p
#                 print pos2_sum, pos2, pos2_extra, pos2_offset, overhang_inner_3p
#                  
#                 if pos1_sum > pos2_sum:
#                     print "5p inner"
#                     confidence_inner = pos1_sum / sum_inner
#                     candidate.overhang_inner = overhang_inner_5p
#                     candidate.overhang_inner_conf = confidence_inner
#                 else:
#                     print "3p inner"
#                     confidence_inner = pos2_sum / sum_inner
#                     candidate.overhang_inner = overhang_inner_3p
#                     candidate.overhang_inner_conf = confidence_inner
#              
#                 print  
#                 print "inner overhang", candidate.overhang_inner
#                 print "confidence", candidate.overhang_inner_conf
#                 print "5p", dist_5pe
#                 print "3p", dist_3pb
#                 print overhang_inner_5p, end_5p, est_5pe, offset_3pb, dist_5pe[offset_5pe]
#                 print overhang_inner_3p, begin_3p, est_3pb, offset_5pe, dist_3pb[offset_3pb]
#                  
#                 print 
#                 assert 0
#             
# #             print overhang_outer_5p, begin_5p, est_5pb, offset_3pe, dist_5pb[offset_5pb]
# #             print overhang_outer_3p, end_3p, est_3pe, offset_5pb, dist_3pe[offset_3pe]
# #             print overhang_inner_5p, end_5p, est_5pe, offset_3pb, dist_5pe[offset_5pe]
# #             print overhang_inner_3p, begin_3p, est_3pb, offset_5pe, dist_3pb[offset_3pb]
# #             print "\t outer",  candidate.overhang_outer
# #             print "\t inner",  candidate.overhang_inner
#                         
# #             assert overhang_outer_5p == overhang_outer_3p
# #             assert overhang_inner_5p == overhang_inner_3p