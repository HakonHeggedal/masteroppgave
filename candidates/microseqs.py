'''
Created on 19. nov. 2014

@author: hakon
'''
import SuffixTree

from SuffixTree import SubstringDict

def align_small_seqs(candidates, small_seqs, small_seqs_copies):
    
    print "SMALL seq analysis"
    print "candidates:", len(candidates), candidates[0].hairpin
    print "small seqs:", len(small_seqs), len(small_seqs_copies)
    print small_seqs[0]
    print small_seqs_copies[0]
    print "123"
    
    c_used = set()

#     candidate_suffixes = SuffixTree.SuffixTree()
    find_candidates = SubstringDict()

    for i, candidate in enumerate(candidates):
        print
        print "hairpin", i,  candidate.hairpin
        print "padded?", candidate.hairpin_padded
        assert candidate.hairpin_padded.find(candidate.hairpin) == candidate.padding_size
        
#         candidate_suffixes.add(candidate.hairpin, i)
#         print "hp", candidate.hairpin
        find_candidates[candidate.hairpin_padded] = i 
    
    for j, seq in enumerate(small_seqs):
        
        candidate_set = find_candidates[seq]
#         print seq
#         print len(candidates)
#         if candidate_set: print candidate_set, seq
        for cnr in candidate_set:
            index = candidates[cnr].hairpin_padded.find(seq)
            copies = small_seqs_copies[j]
            candidates[cnr].set_small_seq(index, copies)
            
            c_used.add(cnr)
#             print "\tnumber?", cnr, 

            p3 = 40+candidate.pos_3p_begin-candidate.pos_5p_begin
            print
            print index, 40, p3,
            print 40 == index or p3 == index,
            print "\t", small_seqs_copies[j]
#             print candidates[cnr].hairpin
            

            
    print "candidates with small seqs:",len(c_used)
#             pos[]
            
#             if pos[candidate] # find all candidates, and update positions
        
#         match_length, suffix_node, endpos = candidate_suffixes.match(seq)
#         
#         if len(seq) is match_length:
#             #TODO find all candidates 
            
    assert False
    