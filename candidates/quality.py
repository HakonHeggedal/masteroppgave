'''
Created on 19. okt. 2014

@author: hakon
'''

def candidate_quality(candidates, seq_to_candidates):
    ''' sequences matching both this candidate and others /
        sequences matching this candidate
    '''
    
    for candidate in candidates:
        total_hits = 0
        candidate_hits = 0
       
        for interval in candidate.mapped_sequences:
            seq_name = interval.data[1]
            duplicates = int(seq_name.split("-")[1])
            candidate_hits += duplicates 
            total_hits += len(seq_to_candidates[seq_name].candidates) * duplicates
            
        quality = candidate_hits*1.0 / total_hits
        candidate.set_quality(quality)
        
        print quality, candidate_hits, total_hits