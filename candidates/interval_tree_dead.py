'''
Created on 24. feb. 2015

@author: hakon
'''

complimented = {"A":"T", "T":"A", "G":"C", "C":"G"}
def _reverse_compliment(sequence):
    return "".join(complimented[x] for x in sequence[::-1])


def align_dead_miRNAs(mirna_hits, _,  id_to_mature, candidate_tree, candidate_list, sequence_tree,
                 seq_to_candidates):
    
    
#     reverse_hp = {}
    
#     def rev((k,v)):
#         reverse_hp[v] = k
    
    def _print(el):
        print el
    
#     map(rev, id_to_hairpin.iteritems())
#     map(_print, reverse_hp.iteritems())
#     assert False
#     unique_mirnas = set()
    
    
    for dead_loki in mirna_hits:
        
        miRNAid = dead_loki[0]
#         if miRNAid in unique_mirnas:
#             continue # use first position only
#         
#         unique_mirnas.add(miRNAid)

        strand_dir = dead_loki[1]
        chromosome = dead_loki[2].split("|")[3]
        genome_offset = int(dead_loki[3])
        
        hairpin = dead_loki[4]
        
        
        
        
#         true_hairpin = id_to_hairpin[miRNAid]
        
        
        begin_5p =  genome_offset
        end_5p = genome_offset
        begin_3p = genome_offset
        end_3p = genome_offset
        
        print
        print miRNAid
        
        # put find hairpin pos
        if miRNAid in id_to_mature:
            mature_seq = id_to_mature[miRNAid]
            if strand_dir == "-":
                mature_seq = _reverse_compliment(mature_seq)
            
#             print "\t", hairpin in reverse_hp
            print "\t", strand_dir
            
            print "\t", mature_seq in hairpin
            
            print "\t", mature_seq
            print "\t", hairpin

            
#             print "\t" , hairpin
#             print "\t" mature_seq, hairpin
#             assert mature_seq in hairpin
            
            begin_mature = hairpin.find(mature_seq)
            end_mature = begin_mature + len(mature_seq)
            
            print "\t", begin_mature, end_mature


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