'''
Created on 19. okt. 2014

@author: hakon
'''



def candidate_quality(candidates, seq_to_candidates):

    for candidate in candidates:
        total_hits = 0
       
        for interval in candidate.mapped_sequences:
            seq_name = interval.data[1]
            total_hits += len(seq_to_candidates[seq_name])
        
        quality = len(candidate.mapped_sequences)*1.0 / total_hits
        candidate.set_quality(quality)
        