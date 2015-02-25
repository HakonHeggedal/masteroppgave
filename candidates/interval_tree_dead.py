'''
Created on 24. feb. 2015

@author: hakon
'''

complimented = {"A":"T", "T":"A", "G":"C", "C":"G"}
def _reverse_compliment(sequence):
    return "".join(complimented[x] for x in sequence[::-1])


def align_dead_miRNAs(mirna_hits, _,  id_to_mature, candidate_tree, candidate_list, sequence_tree,
                 seq_to_candidates):
    

    unique_mirnas = set()
    
    for dead_loki in mirna_hits:
        
        miRNAid = dead_loki[0]
        if miRNAid in unique_mirnas:
            continue # use first entry only
         
        unique_mirnas.add(miRNAid)

        strand_dir = dead_loki[1]
        chromosome = dead_loki[2].split("|")[3]
        genome_offset = int(dead_loki[3])
        
        hairpin = dead_loki[4]
        
        is_candidate = False
        
        begin_5p =  genome_offset
        end_5p = genome_offset
        begin_3p = genome_offset
        end_3p = genome_offset
        
        print
        print miRNAid, begin_5p, begin_5p + len(hairpin)
        
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
                begin_5p += begin_mature
                end_5p += end_mature
            else:
                begin_3p += begin_mature
                end_3p += end_mature
            
            print "\t", begin_mature, end_mature, len(hairpin)
            print "\t", begin_5p, end_5p, begin_3p, end_3p
        
        
            
            
        tree = candidate_tree[chromosome]
        if tree:
            
            candidates = tree[genome_offset:genome_offset+len(hairpin)]
            
            print len(candidates), candidates
            
            for candidate in candidates:
                
                if candidate.data.chromosome_direction != strand_dir:
                    continue
                
                shift_start = abs(genome_offset - candidate.data.pos_5p_begin)
                shift_end = abs( (genome_offset+len(hairpin)) - candidate.data.pos_3p_end)
                hashval = candidate.data.chromosome + strand_dir + str(candidate.data.pos_5p_begin)
                
                print shift_start, shift_end, len(hairpin)
                
                if shift_start + shift_end < len(hairpin) / 2:
                    candidate_to_miRNAid[hashval] = miRNAid
                    is_candidate = True
                    break
    
        else:
            print "no", tree
            
        if not is_candidate:
            pass
        


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