'''
Created on 24. feb. 2015

@author: hakon
'''

from candidates import structure, interval_tree_search

complimented = {"A":"T", "T":"A", "G":"C", "C":"G"}
def _reverse_compliment(sequence):
    return "".join(complimented[x] for x in sequence[::-1])

def align_miRNAs(mirna_hits, hairpinID_to_mature, hpID_to_mseqs, candidate_tree, candidate_list, sequence_tree,
                 seq_to_candidates, miRNA_species, miRNA_high_conf):

    candidate_to_miRNAid = {}
    noseq_set = set()
    unique_mirnas = set()
    candidate_already = 0

    candidate_count = 0
    miRNA_with_candidates = set()
    has_seqs = []
    noseqs = 0
    
    
    for loki in mirna_hits:
        
        miRNAid = ">" + loki[0]
        unique_mirnas.add(miRNAid)

        strand_dir = loki[1]
        chromosome = loki[2].split("|")[3]
        genome_offset = int(loki[3])
        hairpin = loki[4]

        _OLDmature_pos = hairpinID_to_mature[miRNAid]
        mature_seqs = hpID_to_mseqs[miRNAid] if miRNAid in hpID_to_mseqs else []
        

#         print
#         print strand_dir, miRNAid
#         for s in mature_seqs:
#             if strand_dir == "-":
#                 print _reverse_compliment(s) in hairpin,
#                 assert _reverse_compliment(s) in hairpin
#             else:
#                 print s in hairpin,
#                 assert s in hairpin 
#         
#         print mature_seqs
        
        
        mature_pos = []
        for seq in mature_seqs:
            if strand_dir == "-":
                seq = _reverse_compliment(seq)
            pos = hairpin.find(seq)
            if pos > -1:
                mature_pos.append((pos,pos+len(seq)))
            else:
                print "mature seq not mapping:", seq, mature_seqs
            assert seq in hairpin
        
        
        if mature_pos:
            mature_pos = sorted(mature_pos)
            
        if len(mature_pos) > 1:
            if mature_pos[-2][1] > mature_pos[-1][0]:
#                 print
#                 print mature_pos
#                 print "remove last of overlapping mature seqs",
#                 print mature_pos[-2],mature_pos[-2][1], mature_pos[-1], mature_pos[-1][0],
#                 print mature_pos[-2][1] > mature_pos[-1][0]
                
                mature_pos.pop(-1)
#                 print _OLDmature_pos, mature_pos
#                 assert False
                
            elif mature_pos[0][1] > mature_pos[1][0]:
#                 print
#                 print mature_pos
#                 print "remove last of overlapping mature seqs, part 2",
#                 print mature_pos[0],mature_pos[0][1], mature_pos[1], mature_pos[1][0],
#                 print mature_pos[0][1] > mature_pos[1][0]
                
                mature_pos.pop(1)
#                 print _OLDmature_pos, mature_pos
#                 assert False

#         print _OLDmature_pos, mature_pos
        assert len(mature_pos) <= 2
        

        mature_len = max( map(len, mature_seqs)) if mature_seqs else 10
        if mature_len < 12:
            mature_len = 20
        
        begin_5p = genome_offset + 0
        end_5p =  genome_offset + mature_len
        begin_3p =  genome_offset + len(hairpin) - mature_len
        end_3p = genome_offset + len(hairpin)
        
        if len(mature_pos) == 2:
            begin_5p = genome_offset + mature_pos[0][0]
            end_5p =  genome_offset + mature_pos[0][1]
            begin_3p =  genome_offset + mature_pos[1][0]
            end_3p = genome_offset + mature_pos[1][1]
            
        elif len(mature_pos) == 1:
            avg_val = (mature_pos[0][0] + mature_pos[0][1] ) / 2.0
            if avg_val < len(hairpin) / 2.0:
                begin_5p = genome_offset + mature_pos[0][0]
                end_5p = genome_offset + mature_pos[0][1]
            else:
                begin_3p = genome_offset + mature_pos[0][0]
                end_3p = genome_offset + mature_pos[0][1]
            
        
#         print begin_5p, end_5p, begin_3p, end_3p
        assert begin_5p >= genome_offset
        assert begin_5p < end_5p <= begin_3p < end_3p

        
#         begin_5p = end_5p = begin_3p = end_3p = genome_offset
#         
#         begin_5p += mature_pos[0] if mature_pos[0] != -1 else 0
#         end_5p +=  mature_pos[1] if mature_pos[1] != -1 else  mature_len
#         begin_3p +=  mature_pos[2] if mature_pos[2] != -1 else  len(hairpin) - mature_len
#         end_3p +=  mature_pos[3] if mature_pos[3] != -1 else  len(hairpin)

        
        is_candidate = False
        
        tree = candidate_tree[chromosome]
        if not tree:
#             print "no:", chromosome
            continue # TODO : SKIP everything? ALLWAYS?
        
        candidates = tree[begin_5p:end_3p]

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

#         assert candidate.pos_5p_begin < candidate.pos_5p_end < candidate.pos_3p_begin < candidate.pos_3p_end
        
#         no candidates aligns the "miRNA"
        tree = sequence_tree[chromosome]
        sequences = tree[begin_5p:end_3p]
        sequences = [s for s in sequences if s.data[0] == strand_dir]
#         sequences = [s for s in sequences if s.begin >= genome_offset and s.end <= genome_offset+len(hairpin)]
        sequences = set(sequences)
        
        if sequences:
            best_start_pos, _, best_end_pos, _ = interval_tree_search._best_interval(sequences, begin_5p)
            
#             print
#             print sequences
#             print "seqs in hp: ", strand_dir, genome_offset, genome_offset+len(hairpin), hairpin
#             for s in sorted(sequences):
#                 print s.begin, s.end, s.data[2], s.data[2] in hairpin, s.data[0]
            
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
    
     
    
#     print "nr of miRNA bowtie hits:\t", len(mirna_hits)
#     print "unique miRNAs (after bowtie):\t", len(unique_mirnas)
#     print "miRNA aligns with candidate:\t", candidate_already
#     print "Unique mirna aligning with candidate:\t", len(miRNA_with_candidates), len(miRNA_with_candidates) * 1.0 / len(unique_mirnas)
#     print
#     print "set of candidates with 1+ seq aligning:", len(set(candidate_to_miRNAid.iterkeys())), len(list(candidate_to_miRNAid.iterkeys()))

    has_seqs = set(has_seqs)
    
    hc_len = len( miRNA_high_conf)
    sec_cand = has_seqs | miRNA_with_candidates
    
    print
    print "miRNA aligning sequences:\t\t", len(has_seqs) * 1.0 / len(unique_mirnas)
    print "miRNA aligning candidates:\t\t", len(miRNA_with_candidates) * 1.0 / len(unique_mirnas)    
    print "miRNA aligning seq/candidate:\t\t", len(sec_cand) * 1.0 / len(unique_mirnas)
    print
    print "HIGH CONFIDENCE -- sequences \t\t", len(miRNA_high_conf & has_seqs) * 1.0 / hc_len
    print "HIGH CONFIDENCE -- candidates \t\t", len(miRNA_high_conf & miRNA_with_candidates) * 1.0 / hc_len
    print "HIGH CONFIDENCE -- seq/candidates \t", len(miRNA_high_conf & sec_cand) * 1.0 / hc_len
    
#     print "no sequences aligning at all:\t",  len(noseq_set)
#     print "no seqs vs high confidence:\t", len(noseq_set.intersection(miRNA_high_conf) ), len(miRNA_high_conf)
    
#     assert False
    return candidate_to_miRNAid






