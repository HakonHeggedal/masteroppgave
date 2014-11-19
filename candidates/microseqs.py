'''
Created on 19. nov. 2014

@author: hakon
'''
import SuffixTree

from SuffixTree import SubstringDict

def align_small_seqs(candidates, small_seqs, small_seqs_copies):
    
    pos = [{} for _ in xrange(len(candidates))]
#     candidate_suffixes = SuffixTree.SuffixTree()
    find_candidates = SubstringDict()

    for i, candidate in enumerate(candidates):
#         candidate_suffixes.add(candidate.hairpin, i)
        find_candidates[candidate.hairpin] = i 
    
    for j, seq in enumerate(small_seqs):
        
        candidates = find_candidates[seq]
        
        for candidate in candidates:
            print "number?", candidate 
            if pos[candidate] # find all candidates, and update positions
        
#         match_length, suffix_node, endpos = candidate_suffixes.match(seq)
#         
#         if len(seq) is match_length:
#             #TODO find all candidates 
            
        

    