'''
Created on 19. nov. 2014

@author: hakon
'''
# import SuffixTree

from SuffixTree import SubstringDict
import math


def align_small_seqs(candidates, small_seqs, small_seqs_copies):
    
    print "\n\tSmall seq analysis"
    print "\tcandidates:", len(candidates), candidates[0].hairpin
    print "\tsmall seqs:", len(small_seqs), len(small_seqs_copies)


    max_frequent_seq = 0
    c_used = set()
    s_used = set()

    # adding all canididate seqs to a substring dict
    find_candidates = SubstringDict()

    for i, candidate in enumerate(candidates):
        find_candidates[candidate.hairpin_padded_40] = i 

    
    # for each small seq. find the set of matching candidates
    for j, seq in enumerate(small_seqs):
        
        
        candidate_set = find_candidates[seq]
        
        if len(candidate_set) > max_frequent_seq:
            max_frequent_seq = len(candidate_set)
            
        # add the seq. to that candidate
        for cnr in candidate_set:
            copies = small_seqs_copies[j]
            
            if not candidates[cnr].small_subs:
                candidates[cnr].small_subs = {}
                
            candidates[cnr].small_subs[seq] = copies
            
            c_used.add(cnr)
            s_used.add(j)
        
            
    print "\tmost frequent small hit", max_frequent_seq, "times"
    print "\tsmall seq hits:", len(s_used)
    print "\tcandidates with small seqs:",len(c_used)


def small_seq_stats(candidates):
    
    
    has_small_seqs = [c for c in candidates if c.small_subs]
    print len(has_small_seqs)
    
    for c in has_small_seqs:
        print
        hairpin = c.hairpin_padded_40
        total_subseqs = sum(c.small_subs.values())
        
        start_5 = c.pos_5p_begin - c.hairpin_start + 40 if c.pos_5p_begin != -1 else 40
        end_3 = c.pos_3p_end - c.hairpin_start + 40 if c.pos_3p_end != -1 else len(hairpin) - 40
        
        end_5 = c.pos_5p_end - c.hairpin_start + 40 if c.pos_5p_end != -1 else start_5 + 24
        start_3 = c.pos_3p_begin - c.hairpin_start + 40 if c.pos_3p_begin != -1 else end_3 - 24
        
        if end_5 >= start_3:
            end_5 = (start_5 + end_3) / 2
            start_3 = end_5 + 1
        
        
        
        area_5p_short = {pos: 0.0 for pos in range(start_5, end_5)}
        area_3p_short = {pos: 0.0 for pos in range(start_3, end_3)}
        
        
        area_5p_long = set(range(start_5, end_5))
        area_5p_long_sum = 0.0
        area_3p_long = set(range(start_3, end_3))
        area_3p_long_sum = 0.0
        

#         area_5p_long = {pos: 0.0 for pos in range(start_5, end_5)}
#         area_3p_long = {pos: 0.0 for pos in range(start_3, end_3)}

        

        for seq, copies in c.small_subs.iteritems():
            start_seq = hairpin.find(seq) + 40
            
            
#             end_seq = start_seq + len(seq)
#             if start_seq > end_3:
#                 starts_after.append((start_seq, end_seq))
#                             
#             elif end_seq < start_5:
#                 ends_before.append((start_seq, end_seq))
#             print start_seq, end_seq, copies, seq
            
            if start_seq in area_5p_short:
                area_5p_short[start_seq] += copies
            elif start_seq in area_3p_short:
                area_3p_short[start_seq] += copies
            else:
                "seq not in any mature", (start_5, end_5), (start_3, end_3)
#                 assert 0, ((start_5, end_5), (start_3, end_3))
        
#         print "ends before:", start_5, sorted(ends_before) 
#         print "starts after:", end_3, sorted(starts_after) 
#         print len(c.mapped_sequences)
        
        if len(c.mapped_sequences):
            for i in c.mapped_sequences:
                 
                copies = float(i.data[1].split("-")[1])
                seq_start = i.begin-c.hairpin_start + 40
                 
                if seq_start in area_5p_long:
                    area_5p_long_sum += copies
                elif seq_start in area_3p_long:
                    area_3p_long_sum += copies
                else:
                    print seq_start,
                     
                
        area_5p_short_sum = sum(area_5p_short.values())
        ratio_short_long_5p = area_5p_short_sum * 1.0 / area_5p_long_sum if area_5p_long_sum else 0.0
        
        area_3p_short_sum = sum(area_3p_short.values())
        ratio_short_long_3p = area_3p_short_sum * 1.0 / area_3p_long_sum if area_3p_long_sum else 0.0
        
        
        
        print "\n"
        print (start_5, end_5), (start_3, end_3)
        print min(area_5p_long), max(area_5p_long),
        print min(area_3p_long), max(area_3p_long)
        print (c.pos_5p_begin - c.hairpin_start + 40, c.pos_5p_end - c.hairpin_start + 40), (c.pos_3p_begin - c.hairpin_start + 40, c.pos_3p_end - c.hairpin_start + 40)
        print total_subseqs
        print area_5p_short_sum, "\t", area_3p_short_sum
        print area_5p_long_sum, "\t", area_3p_long_sum, c.mapped_sequences
        print ratio_short_long_5p, "\t", ratio_short_long_3p
        
        c.ratio_short_long_5p = ratio_short_long_5p
        c.ratio_short_long_3p = ratio_short_long_3p
        
        print c.ratio_short_long_5p
        print c.ratio_short_long_3p

#     print len([c for c in candidates if c.ratio_short_long_5p])
    has_ratio = [math.log(c.ratio_short_long_3p) for c in candidates if c.ratio_short_long_3p]
    print len(has_ratio)
    print sum(has_ratio) / len(has_ratio)


    
    
    
    
    