'''
Created on 19. okt. 2014

@author: hakon
'''

def candidate_quality(candidates, seq_to_candidates):
    ''' sequences matching both this candidate and others /
        sequences matching this candidate
    '''
    print "candidate qual"
    for candidate in candidates:
        total_hits = 0
        candidate_hits = 0
       
        for interval in candidate.mapped_sequences:
            seq_name = interval.data[1]
            duplicates = float(seq_name.split("-")[1])
            candidate_hits += duplicates
            print "\t", seq_name, len(candidate.mapped_sequences), seq_name in seq_to_candidates
            total_hits += len(seq_to_candidates[seq_name].candidates) * duplicates
        
        
        quality = candidate_hits*1.0 / total_hits if total_hits > 0 else 0
        candidate.set_quality(quality)
        
        print quality, candidate_hits, total_hits