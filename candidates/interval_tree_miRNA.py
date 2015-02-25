'''
Created on 24. feb. 2015

@author: hakon
'''

from candidates import structure, interval_tree_search


def align_miRNAs(mirna_hits, hairpinID_to_mature, candidate_tree, candidate_list, sequence_tree,
                 seq_to_candidates, miRNA_species, miRNA_high_conf):

    candidate_to_miRNAid = {}
    noseq_set = set()
    unique_mirnas = set()
    candidate_already = 0

    candidate_count = 0
    miRNA_with_candidates = set()
    has_seqs = []
    noseqs = 0
    
    
    for dead_loki in mirna_hits:
        
        miRNAid = ">" + dead_loki[0]
        unique_mirnas.add(miRNAid)

        strand_dir = dead_loki[1]
        chromosome = dead_loki[2].split("|")[3]
        genome_offset = int(dead_loki[3])
        hairpin = dead_loki[4]

        mature_pos = hairpinID_to_mature[miRNAid]
        
        mature_len = max(mature_pos[3]-mature_pos[2], mature_pos[1]-mature_pos[0])
        if mature_len < 12:
            mature_len = 22
        
        begin_5p = end_5p = begin_3p = end_3p = genome_offset
        
        begin_5p += mature_pos[0] if mature_pos[0] != -1 else 0
        end_5p +=  mature_pos[1] if mature_pos[1] != -1 else  mature_len
        begin_3p +=  mature_pos[2] if mature_pos[2] != -1 else  len(hairpin) - mature_len
        end_3p +=  mature_pos[3] if mature_pos[3] != -1 else  len(hairpin)
        
        is_candidate = False
        
#         print str(123) + " " + miRNAid
#         print mature_pos
#         print mature_pos[2], mature_pos[2] == -1, len(hairpin) - mature_len, mature_len, len(hairpin)
#         print begin_5p, end_5p, begin_3p, end_3p
#         assert begin_5p < end_5p <= begin_3p < end_3p

        
#         print mature_pos, dead_loki
        
        tree = candidate_tree[chromosome]
        if not tree:
#             print "no:", chromosome
            continue # TODO : SKIP everything? ALLWAYS?
        candidates = tree[begin_5p:end_3p]

        
#         if len(candidates) > 1:
#             print
#             for c in candidates:
#                 print c.data.chromosome_direction
#                 print "\t", c
#                 print "\t", c.data.mapped_sequences
#                 print "\t", c.data.hairpin
# #                 print "\t", c.data.hairpin_fold_10
#                 print "\t", miRNAid, miRNAid in miRNA_high_conf
#                 for ses in c.data.mapped_sequences:
#                     print "\t\t", ses
#             print candidates

        if candidates:
            miRNA_with_candidates.add(miRNAid)
            candidate_already += 1
            
            for candidate in candidates:
                
                if candidate.data.chromosome_direction != strand_dir:
                    continue

                candidate_count += 1
                
                len_5 = candidate.data.pos_5p_end - candidate.data.pos_5p_begin
                shift_5 = abs(candidate.data.pos_5p_begin - begin_5p)
                shift_5 += abs(candidate.data.pos_5p_end - end_5p)
                
                len_3 = candidate.data.pos_3p_end - candidate.data.pos_3p_begin
                shift_3 = abs(candidate.data.pos_3p_begin - begin_3p)
                shift_3 += abs(candidate.data.pos_3p_end - end_3p)
                
                hashval = candidate.data.chromosome + strand_dir + str(candidate.data.pos_5p_begin)
                
                wrong_shift_end = abs(candidate.data.pos_5p_begin - begin_3p)
                wrong_shift_end += abs(candidate.data.pos_5p_end - end_3p)
                
                wrong_shift_middle = abs(candidate.data.pos_3p_begin - begin_5p)
                wrong_shift_middle += abs(candidate.data.pos_3p_end - end_5p)
                
                

                # miRNA should overlap with candidate
                if shift_5 < len_5*2 or shift_3 < len_3*2 or shift_5+shift_3 < (len_3+len_5)*2:
                    candidate_to_miRNAid[hashval] = miRNAid
                    is_candidate = True
                    break


        if is_candidate:
            continue
        
#         print 
#         assert candidate.pos_5p_begin < candidate.pos_5p_end < candidate.pos_3p_begin < candidate.pos_3p_end
        
#         no candidates aligns the "miRNA"
        tree = sequence_tree[chromosome]
        sequences = tree[begin_5p:end_3p]
        if sequences:
            best_start_pos, _, best_end_pos, _ = interval_tree_search._best_interval(sequences, begin_5p)
            
#             print 123
#             print best_start_pos, best_end_pos, begin_5p, end_5p, begin_3p, end_3p
            
            avgpos = (best_start_pos + best_end_pos) / 2.0
            
            offset = begin_5p
            
            if avgpos < (end_3p - begin_5p) / 2.0 :
                # peak is 5p
                begin_5p = best_start_pos + offset
                end_5p = best_end_pos + offset
                if end_5p > begin_3p:
                    begin_3p = end_5p
            else:
                # peak is 3p
                begin_3p = best_start_pos + offset
                end_3p = best_end_pos + offset
                if begin_3p < end_5p:
                    end_5p = begin_3p
                
#             print best_start_pos, best_end_pos, begin_5p, end_5p, begin_3p, end_3p
            
            has_seqs.append(miRNAid)
#             print best_start_pos, best_start, best_end_pos, best_end
#             if miRNAid in miRNA_high_conf:
#                 print "mirna positions:", begin_5p, end_3p
#                 for se in sequences:
#                     print "\t", se.begin-begin_5p, se.end-begin_5p, se.data
#                 print
#             pass
            
        else:
#             no sequences at all
            noseqs += 1
            noseq_set.add(miRNAid)
            pass
        
        assert begin_5p < end_5p <= begin_3p < end_3p
        
        candidate = structure.Candidate(chromosome,
                         strand_dir,
                         begin_5p,
                         end_5p,
                         begin_3p,
                         end_3p,
                         sequences)
        

        for candidate_interval in sequences:
            name = candidate_interval.data[1]
            if name not in sequences:
                number_id = int(name.split("-")[0])
                duplicates = float(name.split("-")[1])

                s = structure.Sequence(number_id, duplicates, candidate_interval.data[2])
                s.add_candidate(candidate)
                seq_to_candidates[name] = s
            else:
                seq_to_candidates[name].add_candidate(candidate)
        
        if end_3p - begin_5p > 200:
            print "\t200+ length", begin_5p, end_5p, begin_3p, end_3p
            assert False
        
        assert candidate.pos_5p_begin < candidate.pos_5p_end <= candidate.pos_3p_begin < candidate.pos_3p_end
        
        candidate_list.append(candidate)
        hashval = chromosome + strand_dir + str(begin_5p)
        
        candidate_to_miRNAid[hashval] = miRNAid
    
    
    print
    print "nr of miRNA bowtie hits:\t", len(mirna_hits)
    print "unique miRNAs (after bowtie):\t", len(unique_mirnas)
    print "miRNA aligns with candidate:\t", candidate_already
    print "Unique mirna aligning with candidate:\t", len(miRNA_with_candidates), len(miRNA_with_candidates) * 1.0 / len(unique_mirnas)
    print
    print "set of candidates with 1+ seq aligning:", len(set(candidate_to_miRNAid.iterkeys())), len(list(candidate_to_miRNAid.iterkeys()))

    has_seqs = set(has_seqs)
    print "miRNA only aligning sequences:\t\t", len(set(has_seqs)), len(has_seqs) * 1.0 / len(unique_mirnas)
    
    has_seqs.update(miRNA_with_candidates)
    print "miRNA with candidate or sequences:\t", len(has_seqs), len(has_seqs) * 1.0 / len(unique_mirnas)
    print "no sequences aligning at all:\t", noseqs, len(noseq_set)
    print "no seqs vs high confidence:\t", len(noseq_set.intersection(miRNA_high_conf) ), len(miRNA_high_conf)
    
#     assert False
    return candidate_to_miRNAid

