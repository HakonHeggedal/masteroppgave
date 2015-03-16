'''
Created on 24. feb. 2015

@author: hakon
'''
import interval_tree_search
import structure

complimented = {"A":"T", "T":"A", "G":"C", "C":"G"}
def _reverse_compliment(sequence):
    return "".join(complimented[x] for x in sequence[::-1])


def align_dead_miRNAs(mirna_hits, _,  id_to_mature, candidate_tree, candidate_list, sequence_tree,
                 seq_to_candidates):
    
    print "\naligning dead miRNA"
    candidate_to_dead = {}
    
    unique_mirnas = set()
    candidated = set()
    seqd = set()
    
    both_matures = 0
    
    for dead_loki in mirna_hits:
        
        miRNAid = dead_loki[0]
        if miRNAid in unique_mirnas:
            continue # use first entry only
        
#         print "\n\t new dead...",
        unique_mirnas.add(miRNAid)

        strand_dir = dead_loki[1]
        chromosome = dead_loki[2].split("|")[3]
        genome_offset = int(dead_loki[3])
        
        hairpin = dead_loki[4]
        
        is_candidate = False
        
        begin_5p =  -1
        end_5p = -1
        begin_3p = -1
        end_3p = -1
        
#         begin_5p =  genome_offset
#         end_5p = genome_offset
#         begin_3p = genome_offset
#         end_3p = genome_offset
        
#         print
#         print miRNAid, begin_5p, begin_5p + len(hairpin)
        
        # put mature seq into 5p or 3p
        if miRNAid in id_to_mature:
            mature_seq = id_to_mature[miRNAid]
            if strand_dir == "-":
                mature_seq = _reverse_compliment(mature_seq)
            
#             print "\t", strand_dir
#             print "\t", mature_seq in hairpin
#             print "\t", mature_seq
#             print "\t", hairpin
            

            begin_mature = hairpin.find(mature_seq)
            end_mature = begin_mature + len(mature_seq)
            
            avg_val = (begin_mature + end_mature ) / 2.0
            
            if avg_val < len(hairpin) / 2.0:
                begin_5p = genome_offset + begin_mature
                end_5p = genome_offset + end_mature
            else:
                begin_3p = genome_offset + begin_mature
                end_3p =genome_offset + end_mature
            
#             print "\t", begin_mature, end_mature, len(hairpin)
#             print "\t", begin_5p, end_5p, begin_3p, end_3p
        
#         print begin_5p, end_5p, begin_3p, end_3p, "\t", len(hairpin)
        
        
            
        tree = candidate_tree[chromosome]
        if tree:
            
            candidates = tree[genome_offset:genome_offset+len(hairpin)]
            
#             print len(candidates), candidates
            
            for candidate in candidates:
                
                if candidate.data.chromosome_direction != strand_dir:
                    continue
                
                hashval = candidate.data.chromosome + strand_dir + str(candidate.data.pos_5p_begin)
                
                shift_start = abs(genome_offset - candidate.data.pos_5p_begin)
                shift_end = abs( (genome_offset+len(hairpin)) - candidate.data.pos_3p_end)
                
#                 print shift_start, shift_end, len(hairpin)
                
                if shift_start + shift_end < len(hairpin) / 2:
                    candidate_to_dead[hashval] = miRNAid
                    is_candidate = True
                    candidated.add(miRNAid)
#                     print "-- candidate"
                    break
    
        else:
            print "no", tree
        
        is_both_matures = begin_5p != -1 and end_5p != -1 and begin_3p != -1 and end_3p != -1
        
        if (not is_candidate) and (not is_both_matures):
        
            tree = sequence_tree[chromosome]
            
            if not tree:
                continue
            
            sequences = tree[genome_offset:genome_offset+len(hairpin)]
            sequences = [s for s in sequences if s.data[0] == strand_dir]
#             sequences = [s for s in sequences if s.begin >= genome_offset and s.end <= genome_offset+len(hairpin)]
            sequences = set(sequences)
            if sequences:


                best_start_pos, _, best_end_pos, _ = interval_tree_search._best_interval(sequences, genome_offset)
                
                avgpos = (best_start_pos + best_end_pos) / 2.0
                halfsize = len(hairpin) / 2.0 
                
#                 offset = begin_5p
#                 print "--"
#                 print avgpos, halfsize, halfsize_derp
#                 print begin_5p, end_5p, begin_3p, end_3p
                
                if avgpos <= halfsize:
                    # peak is 5p
#                     print "<-"
                    begin_5p = genome_offset + best_start_pos
                    end_5p = genome_offset + best_end_pos
#                     if end_5p > begin_3p:
#                         begin_3p = end_5p
                else:
                    # peak is 3p
#                     print "->"
                    begin_3p = genome_offset + best_start_pos
                    end_3p = genome_offset + best_end_pos
#                     if begin_3p < end_5p:
#                         end_5p = begin_3p
            

            
        
            candidate = structure.Candidate(chromosome,
                     strand_dir,
                     genome_offset,
                     genome_offset+len(hairpin),
                     begin_5p,
                     end_5p,
                     begin_3p,
                     end_3p,
                     sequences)
            
            if sequences:
                seqd.add(miRNAid)
            
            
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

                    
            assert candidate.pos_5p_begin < candidate.pos_5p_end or candidate.pos_5p_begin == -1
            assert candidate.pos_3p_begin < candidate.pos_3p_end or candidate.pos_3p_begin == -1
            
            if begin_5p != -1 and end_3p != -1:
                both_matures += 1
#             if candidate.pos_5p_end >= candidate.pos_3p_begin + 10:
#                 print candidate.pos_5p_end, candidate.pos_3p_begin
#                 assert candidate.pos_5p_end < candidate.pos_3p_begin + 10
                
    
#             if end_3p - begin_5p > 200:
#                 print "\t200+ length", begin_5p, end_5p, begin_3p, end_3p
#                 assert False
            
    #         assert candidate.pos_5p_begin < candidate.pos_5p_end <= candidate.pos_3p_begin < candidate.pos_3p_end
            
            candidate_list.append(candidate)
            hashval = chromosome + strand_dir + str(genome_offset)
            
            candidate_to_dead[hashval] = miRNAid
#             print "123"
        

    print len(unique_mirnas)
    
    print "\ndead stats:"
    print "has both seqs:", both_matures, (both_matures + len(candidated))*1.0 / len(unique_mirnas) 
    print "aligning with candidates:\t", len(candidated), len(candidated)*1.0 / len(unique_mirnas) 
    print "aligning with seqs:\t", len(seqd), len(seqd)*1.0 / len(unique_mirnas) 
    print "aligning with either:\t", len(candidated | seqd), len(candidated | seqd)*1.0 / len(unique_mirnas) 
    
#     assert False
    return candidate_to_dead

#  a = "AAGCT"
# print _reverse_compliment(a)
# def yes():
#     while True:
#         yield "yes"
#          
# def p(x): print x
# 
# 
# map(p, yes())





















#             print begin_5p, end_5p, begin_3p, end_3p
            
#             if begin_5p == end_5p == begin_3p == end_3p:
#                 # using random values as mature length
#                 random_len = 20 if len(hairpin) > 40 else len(hairpin) / 2
#                 
#                 end_5p += random_len
#                 begin_3p += len(hairpin) - random_len
#                 end_3p += len(hairpin)
#             
#             elif begin_5p == end_5p:
#                 
#                 mature_len = end_3p - begin_3p
#                 mature_len = mature_len if 30 < mature_len > 10 else 20
#                 
#                 end_5p += mature_len            
# 
#             elif begin_3p == end_3p:
#                 
#                 mature_len = end_5p - begin_5p
#                 mature_len = mature_len if 30 < mature_len > 10 else 20
#                 
#                 begin_3p += len(hairpin) - mature_len 
#                 end_3p += len(hairpin)
                
            
#             print begin_5p, end_5p, begin_3p, end_3p
