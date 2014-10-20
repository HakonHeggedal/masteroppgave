'''
Created on 8. okt. 2014

@author: hakon
'''

from _collections import defaultdict



def tailing_au(candidates, sequences, align_len=8):
    ''' counts the number of a/u-tailing sequences
        compared to all mapping sequences
    '''
    
    tailing = [0]*len(candidates)   
    last_parts = defaultdict([])
    
    for nr, candidate in enumerate(candidates):
        # store parts at the end of 5' and 3'
        candidate_seq = candidate.hairpin_padded
        
        start = candidate.padding_size + candidate.pos_3_end - candidate.pos_3_begin - 1
        
        for i in xrange(start, start+align_len, 1):
            seq = candidate_seq[i:i+align_len]
            last_parts[seq].append(nr)
        
        start = candidate.padding_size + candidate.pos_5_end - candidate.pos_5_begin - 1
        
        for i in xrange(start, start+align_len, 1):
            seq = candidate_seq[i:i+align_len]
            last_parts[seq].append(nr)
            

    for sequence in sequences:
        # find sequences that map to the ends of candidates
        
        seq = sequence.nucleotides
        count = sequence.duplicates
        
        assert seq[-align_len:] not in last_parts
#             print  "already here??????"
        
        seq = _remove_trailing(seq, "AUau")
        seq = seq[-align_len:]
        
        if seq in last_parts:
            print "tailing:", seq, sequence.nucleotides
            for cnr in last_parts[seq]:
                tailing[cnr] += count
                
            
            
            
        
    for candidate, tailing in zip(candidates, tailing):
        candidate.set_tailing(tailing)









































# def three_prime_au(candidates, sequence_frequency):
#     
#     threes_to_aucount = {}
#     s = 0
#     
#     for candidate in candidates:    
#         three_seq = candidate.seq_prime3
#         three_seq = _remove_trailing(three_seq, "AUau")
#         threes_to_aucount[three_seq] = 0
#         
#         
#     for seq, count in sequence_frequency:
#         seq_no_au = _remove_trailing(seq, "AUau")
#         if seq_no_au in threes_to_aucount:
#             threes_to_aucount[seq_no_au] += count
# 
#     
#     for candidate in candidates:
#         three_seq = candidate.seq_prime3
#         three_seq = _remove_trailing(three_seq, "AUau")
#         
#         # set tailing value for all candidates
#         candidate.tailing_au = threes_to_aucount[three_seq]
#         print threes_to_aucount[three_seq]
#         s += threes_to_aucount[three_seq]
#         
#     print s
#             
#     
#     
# 
def _remove_trailing(sequence, chars):
    chars = set(chars)
    cut = len(sequence)
     
    for x in xrange(len(sequence)-1, 0, -1):
        if sequence[x] not in chars:
            break
        cut = x

    return sequence[:cut]




# t = "QWEWRRTYUUUU"
# t = "QWEWRRTAYAUAUT"
# print _remove_trailing_AU(t)

        
    
            