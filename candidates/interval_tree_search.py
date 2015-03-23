'''
Created on 1. okt. 2014

@author: hakon
'''

from intervaltree.bio import GenomeIntervalTree
from candidates import structure
from candidates.interval_tree_misc import best_interval, filter_intervals





MAX_CANDIDATE_LEN = 100
# MIN_CANDIDATE_LEN = 41
MIN_HAIRPIN_LOOP = 10
MAX_HAIRPIN_LOOP = 20
# MIN_MATURE_SEQ = 16
# MAX_MATURE_SEQ = 30 # old?
MAX_MATURE_SEQ = 30
SEARCH_LEN = 100
                

def find_candidates_2(sequence_hits):
    global derpy_errors
    ''' finds microRNA candidates from bowtie data (using interval trees)
    
        sequence_hits -- an iterable of lists on bowtie output format:
  0          1    2                                 3           4                           5                           6
['1-15830', '-', 'gi|224589818|ref|NC_000006.11|', '72113295', 'AGCTTCCAGTCGAGGATGTTTACA', 'IIIIIIIIIIIIIIIIIIIIIIII', '0']
        returns a list of candidates, and the interval tree with all sequences
    '''
    print "nr of sequence hits: ", len(sequence_hits)
    sequence_tree = GenomeIntervalTree()
    candidate_tree = GenomeIntervalTree() # only candidates here
    candidate_list = []
    all_mapped_sequences = set()
    seq_to_candidates = {}
    
    f = 0
    a = 0
    
    print "adding all intervals to the tree"
    for prop in sequence_hits:
#         print prop
        seq_name = prop[0]
        strand_dir = prop[1] # forward: + backward: -
        genome_nr = prop[2].split("|")[3] # which genome (and version)
        genome_offset = int(prop[3]) # offset into the genome, 0-indexed
        dna_sequence = prop[4] # the dna_sequence matching this position.
        sequence_info = [strand_dir, seq_name, dna_sequence]
        
        sequence_tree.addi(genome_nr, genome_offset, genome_offset + len(dna_sequence), sequence_info)
    
    print "\tall intervals added to the tree" 
    
    interval_sum = 0.0

    # test all intervals to find candidates
    for tree in sequence_tree:
        print tree
        
        for interval in sorted(sequence_tree[tree]):
            if interval in all_mapped_sequences:
#                 print "fast skip\t", interval
                continue
            
            start_interval = interval.begin
            end_interval = start_interval + SEARCH_LEN


#             find a peak in this interval
            candidate_intervals = sequence_tree[tree][start_interval:end_interval]
            if not candidate_intervals:
                continue
            
#             filter by direction
            candidate_intervals = [s for s in candidate_intervals if s.data[0] == interval.data[0]]
            if len(candidate_intervals) <= 1:
                continue
            
#             search for more sequences close to this one
            max_end_interval = max(candidate_intervals, key=lambda x:x.end)
            max_end = max_end_interval.end
            candidate_intervals = set(candidate_intervals)
            
            while max_end + MAX_HAIRPIN_LOOP + MAX_MATURE_SEQ > end_interval: # extend search area
                
                new_seqs = sequence_tree[tree][end_interval:end_interval+SEARCH_LEN]
                new_seqs = [s for s in new_seqs if s.data[0] == interval.data[0]]
                
                if len(new_seqs) == 0:
                    break
                
                candidate_intervals.update(new_seqs)
                end_interval += SEARCH_LEN
                max_end_interval = max(new_seqs, key=lambda x:x.end)
                max_end = max_end_interval.end 
            
            
            
            all_mapped_sequences.update(candidate_intervals) # do not need to use these next iteration
            
#             if len(candidate_intervals) <= 50:
#                 continue
            
            candidate_intervals = sorted(candidate_intervals) # first interval is picked if several equal.
#             print len(candidate_intervals)
            
            i = 0
            while candidate_intervals:
                i += 1
                if i > 1: break
                
                # finding the best interval (highest peaks):
                start_peak, start_peak_val, end_peak, end_peak_val = best_interval(candidate_intervals, start_interval)

                
                # no intervals at all
                if start_peak_val == -1 or end_peak_val == -1 or start_peak == -1 or end_peak == -1:
                    break
                

                # finding interval close before best
                start_before_limit = max(-3,start_peak - MAX_HAIRPIN_LOOP - MAX_MATURE_SEQ )
                stop_before_limit = max(-3, start_peak - MIN_HAIRPIN_LOOP)
                
                five_intervals = filter_intervals(candidate_intervals, start_interval, start_before_limit, stop_before_limit)
                start_before, start_before_val, stop_before, stop_before_val = best_interval(five_intervals, start_interval)
                
                
                # interval close after best
                start_after_limit = end_peak + MIN_HAIRPIN_LOOP
                stop_after_limit = end_peak + MAX_HAIRPIN_LOOP + MAX_MATURE_SEQ
                
                three_intervals = filter_intervals(candidate_intervals, start_interval, start_after_limit, stop_after_limit)
                start_after, start_after_val, stop_after, stop_after_val = best_interval(three_intervals, start_interval)
                
     
                no_peak_after = start_after == -1 or stop_after == -1
                no_peak_after = no_peak_after or start_after_val == -1 or stop_after_val == -1
                
                no_peak_before = start_before == -1 or stop_before == -1
                no_peak_before = no_peak_before or start_before_val == -1 or stop_before_val == -1
                
                 
                p1 = p2 = p3 = p4 = -1
                a += 1
                
                if no_peak_after and no_peak_before:
                    f += 1
                    break

                # best peak is 5p
                elif no_peak_before or start_after_val + stop_after_val > start_before_val + stop_before_val:
                    
                    begin_5p = start_peak + start_interval # best peak
                    end_5p = end_peak + start_interval
                    begin_3p = start_after + start_interval # peak after is second best
                    end_3p = stop_after + start_interval
                    
                    p1 = start_peak_val # peak value, used for testing
                    p2 = end_peak_val
                    p3 = start_after_val
                    p4 = stop_after_val
                    
                    assert begin_5p < end_5p < begin_3p < end_3p
                    
                    assert start_after_val != -1
                    assert stop_after_val != -1
                    assert start_after_val + stop_after_val > start_before_val + stop_before_val
#                     assert stop_after_val > stop_before_val
                
                # best peak is 3p
                else:
                    begin_5p = start_before + start_interval # peak before
                    end_5p = stop_before + start_interval
                    begin_3p = start_peak + start_interval # best peak
                    end_3p = end_peak + start_interval
                    
                    p1 = start_before_val
                    p2 = stop_before_val
                    p3 = start_peak_val
                    p4 = end_peak_val
                    assert begin_5p < end_5p < begin_3p < end_3p
                    
                    
    #             else:
    #                 print "???"
    #                 print no_peak_after, no_peak_before
    #                 print start_before, start_before_val, stop_before, stop_before_val
    #                 print start_peak, start_peak_val, end_peak, end_peak_val
    #                 print start_after, start_after_val, stop_after, stop_after_val
    #                 assert False
    #             elif start_3p_val + end_3p_val > start_5p_val + end_5p_val:
    #                 print "best: 3p"
    #             else:   # start_5p_val + end_5p_val > start_3p_val + end_3p_val:
    #                 print "best: 5p"
    
                strand_dir = interval.data[0]
                chromosome = tree
                
    #             if end_3p - begin_5p > MAX_CANDIDATE_LEN:
    #             print
    #             print end_3p-begin_5p
    #             print begin_5p-start_interval, end_5p-start_interval, begin_3p-start_interval, end_3p-start_interval,
    #             print (start_peak, end_peak),
    #             print start_after_limit, stop_after_limit
                    
                assert end_3p - begin_5p <= MAX_CANDIDATE_LEN
                
#                 qwe = len(candidate_intervals) > 50

#                 if  i > 1:
#                     print "before....", i
# #                     for c in sorted(candidate_intervals):
# #                         print c
#                     print begin_5p, end_3p
#                     print len(candidate_intervals)
                
                # close intervals are the intervals overlapping the candidate
                close_intervals = set()
                for c in candidate_intervals:
                    if begin_5p < c.begin < end_3p or begin_5p < c.end < end_3p:
                        close_intervals.add(c)
                
#                 close_intervals = set(c for c in candidate_intervals if c.begin <= end_3p or c.end >= begin_5p)
                candidate_intervals = set(candidate_intervals)
                candidate_intervals -= close_intervals
                candidate_intervals = list(candidate_intervals)
                
#                 if  i > 1:
#                     
#                     print len(candidate_intervals)
#                     print len(close_intervals)
# #                     for c in sorted(candidate_intervals):
# #                         print c
#                     print len(candidate_intervals)
#                     print i
#                     print
                    
#                     if derp:
#                         assert False
#                     derp = True
#                     assert False
                    
#                 print "adding new canididate from", begin_5p ," sequence hits:", len(close_intervals)
                
                hairpin_start = begin_5p
                hairpin_end = end_3p

                candidate = structure.Candidate(chromosome,
                                                 strand_dir,
                                                 hairpin_start, # used as gene offset. sometimes 5p mature seq. is missing...
                                                 hairpin_end,
                                                 begin_5p,
                                                 end_5p,
                                                 begin_3p,
                                                 end_3p,
                                                 close_intervals)
                
                candidate.candidate_type = structure.TYPE_CANDIDATE
                
                candidate.peak_5b = p1
                candidate.peak_5e = p2
                candidate.peak_3b = p3
                candidate.peak_3e = p4
                
                assert candidate.pos_5p_begin < candidate.pos_5p_end < candidate.pos_3p_begin < candidate.pos_3p_end
                      
                for candidate_interval in close_intervals:
                    name = candidate_interval.data[1]
                    if name not in seq_to_candidates:
                        number_id = int(name.split("-")[0])
                        duplicates = float(name.split("-")[1])
                        interval_sum += duplicates
                        s = structure.Sequence(number_id, duplicates, candidate_interval.data[2])
                        s.add_candidate(candidate)
                        seq_to_candidates[name] = s
                    else:
                        seq_to_candidates[name].add_candidate(candidate)
                     
                candidate_tree[tree][begin_5p:end_3p] = candidate
                candidate_list.append(candidate)
     
                if len(all_mapped_sequences) == 0:
                    all_mapped_sequences = candidate.all_mapped_sequences
    

    print "find candidates 2.0"
    print "candidates:", a-f
    print "tests:", a, (a-f) * 1.0/ a
    print "fail:", f, f * 1.0 / a
    print "sum interval in candidates:", interval_sum
    
    return candidate_tree, sequence_tree, candidate_list, seq_to_candidates



#===============================================================================
# 
# def align_miRNAs(mirna_hits, hairpinID_to_mature, candidate_tree, candidate_list, sequence_tree,
#                  seq_to_candidates, miRNA_species, miRNA_high_conf):
# 
#     candidate_to_miRNAid = {}
# 
#     noseq_set = set()
# 
#     candidate_already = 0
#     
#     unique_mirnas = set()
# 
#     candidate_count = 0
#     miRNA_with_candidates = set()
#     has_seqs = []
#     noseqs = 0
#     
#     
#     for mirna_loki in mirna_hits:
#         
#         miRNAid = ">" + mirna_loki[0]
#         unique_mirnas.add(miRNAid)
# 
#         strand_dir = mirna_loki[1]
#         chromosome = mirna_loki[2].split("|")[3]
#         genome_offset = int(mirna_loki[3])
#         hairpin = mirna_loki[4]
# 
#         mature_pos = hairpinID_to_mature[miRNAid]
#         
#         mature_len = max(mature_pos[3]-mature_pos[2], mature_pos[1]-mature_pos[0])
#         if mature_len < 12:
#             mature_len = 22
#         
#         begin_5p = end_5p = begin_3p = end_3p = genome_offset
#         
#         begin_5p += mature_pos[0] if mature_pos[0] != -1 else 0
#         end_5p +=  mature_pos[1] if mature_pos[1] != -1 else  mature_len
#         begin_3p +=  mature_pos[2] if mature_pos[2] != -1 else  len(hairpin) - mature_len
#         end_3p +=  mature_pos[3] if mature_pos[3] != -1 else  len(hairpin)
#         
#         is_candidate = False
#         
# #         print str(123) + " " + miRNAid
# #         print mature_pos
# #         print mature_pos[2], mature_pos[2] == -1, len(hairpin) - mature_len, mature_len, len(hairpin)
# #         print begin_5p, end_5p, begin_3p, end_3p
# #         assert begin_5p < end_5p <= begin_3p < end_3p
# 
#         
# #         print mature_pos, mirna_loki
#         
#         tree = candidate_tree[chromosome]
#         if not tree:
# #             print "no:", chromosome
#             continue
#         candidates = tree[begin_5p:end_3p]
# 
#         
# #         if len(candidates) > 1:
# #             print
# #             for c in candidates:
# #                 print c.data.chromosome_direction
# #                 print "\t", c
# #                 print "\t", c.data.mapped_sequences
# #                 print "\t", c.data.hairpin
# # #                 print "\t", c.data.hairpin_fold_10
# #                 print "\t", miRNAid, miRNAid in miRNA_high_conf
# #                 for ses in c.data.mapped_sequences:
# #                     print "\t\t", ses
# #             print candidates
# 
#         if candidates:
#             miRNA_with_candidates.add(miRNAid)
#             candidate_already += 1
#             
#             for candidate in candidates:
#                 
#                 if candidate.data.chromosome_direction != strand_dir:
#                     continue
# 
# 
#                 candidate_count += 1
#                 
#                 len_5 = candidate.data.pos_5p_end - candidate.data.pos_5p_begin
#                 shift_5 = abs(candidate.data.pos_5p_begin - begin_5p)
#                 shift_5 += abs(candidate.data.pos_5p_end - end_5p)
#                 
#                 len_3 = candidate.data.pos_3p_end - candidate.data.pos_3p_begin
#                 shift_3 = abs(candidate.data.pos_3p_begin - begin_3p)
#                 shift_3 += abs(candidate.data.pos_3p_end - end_3p)
#                 
#                 hashval = candidate.data.chromosome + strand_dir + str(candidate.data.pos_5p_begin)
#                 
#                 wrong_shift_end = abs(candidate.data.pos_5p_begin - begin_3p)
#                 wrong_shift_end += abs(candidate.data.pos_5p_end - end_3p)
#                 
#                 wrong_shift_middle = abs(candidate.data.pos_3p_begin - begin_5p)
#                 wrong_shift_middle += abs(candidate.data.pos_3p_end - end_5p)
#                 
#                 
# 
#                 # miRNA should overlap with candidate
#                 if shift_5 < len_5*2 or shift_3 < len_3*2 or shift_5+shift_3 < (len_3+len_5)*2:
#                     candidate_to_miRNAid[hashval] = miRNAid
#                     is_candidate = True
#                     break
# 
# 
#         if is_candidate:
#             continue
#         
# #         print 
# #         assert candidate.pos_5p_begin < candidate.pos_5p_end < candidate.pos_3p_begin < candidate.pos_3p_end
#         
# #         no candidates aligns the "miRNA"
#         tree = sequence_tree[chromosome]
#         sequences = tree[begin_5p:end_3p]
#         if sequences:
#             best_start_pos, _, best_end_pos, _ = _best_interval(sequences, begin_5p)
#             
# #             print 123
# #             print best_start_pos, best_end_pos, begin_5p, end_5p, begin_3p, end_3p
#             
#             avgpos = (best_start_pos + best_end_pos) / 2.0
#             
#             offset = begin_5p
#             
#             if avgpos < (end_3p - begin_5p) / 2.0 :
#                 # peak is 5p
#                 begin_5p = best_start_pos + offset
#                 end_5p = best_end_pos + offset
#                 if end_5p > begin_3p:
#                     begin_3p = end_5p
#             else:
#                 # peak is 3p
#                 begin_3p = best_start_pos + offset
#                 end_3p = best_end_pos + offset
#                 if begin_3p < end_5p:
#                     end_5p = begin_3p
#                 
# #             print best_start_pos, best_end_pos, begin_5p, end_5p, begin_3p, end_3p
#             
#             has_seqs.append(miRNAid)
# #             print best_start_pos, best_start, best_end_pos, best_end
# #             if miRNAid in miRNA_high_conf:
# #                 print "mirna positions:", begin_5p, end_3p
# #                 for se in sequences:
# #                     print "\t", se.begin-begin_5p, se.end-begin_5p, se.data
# #                 print
# #             pass
#             
#         else:
# #             no sequences at all
#             noseqs += 1
#             noseq_set.add(miRNAid)
#             pass
#         
#         assert begin_5p < end_5p <= begin_3p < end_3p
#         
#         candidate = structure.Candidate(chromosome,
#                          strand_dir,
#                          begin_5p,
#                          end_5p,
#                          begin_3p,
#                          end_3p,
#                          sequences)
#         
# 
#         for candidate_interval in sequences:
#             name = candidate_interval.data[1]
#             if name not in sequences:
#                 number_id = int(name.split("-")[0])
#                 duplicates = float(name.split("-")[1])
# 
#                 s = structure.Sequence(number_id, duplicates, candidate_interval.data[2])
#                 s.add_candidate(candidate)
#                 seq_to_candidates[name] = s
#             else:
#                 seq_to_candidates[name].add_candidate(candidate)
#         
#         if end_3p - begin_5p > 200:
#             print "\t200+ length", begin_5p, end_5p, begin_3p, end_3p
#             assert False
#         
#         assert candidate.pos_5p_begin < candidate.pos_5p_end <= candidate.pos_3p_begin < candidate.pos_3p_end
#         
#         candidate_list.append(candidate)
#         hashval = chromosome + strand_dir + str(begin_5p)
#         
#         candidate_to_miRNAid[hashval] = miRNAid
#     
#     
#     print
#     print "nr of miRNA bowtie hits:\t", len(mirna_hits)
#     print "unique miRNAs (after bowtie):\t", len(unique_mirnas)
#     print "miRNA aligns with candidate:\t", candidate_already
#     print "Unique mirna aligning with candidate:\t", len(miRNA_with_candidates), len(miRNA_with_candidates) * 1.0 / len(unique_mirnas)
#     print
#     print "set of candidates with 1+ seq aligning:", len(set(candidate_to_miRNAid.iterkeys())), len(list(candidate_to_miRNAid.iterkeys()))
# 
#     has_seqs = set(has_seqs)
#     print "miRNA only aligning sequences:\t\t", len(set(has_seqs)), len(has_seqs) * 1.0 / len(unique_mirnas)
#     
#     has_seqs.update(miRNA_with_candidates)
#     print "miRNA with candidate or sequences:\t", len(has_seqs), len(has_seqs) * 1.0 / len(unique_mirnas)
#     print "no sequences aligning at all:\t", noseqs, len(noseq_set)
#     print "no seqs vs high confidence:\t", len(noseq_set.intersection(miRNA_high_conf) ), len(miRNA_high_conf)
#     
# #     assert False
#     return candidate_to_miRNAid
#===============================================================================


