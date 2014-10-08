'''
Created on 8. okt. 2014

@author: hakon
'''


def three_prime_au(candidates, seq_to_frequency):
    
    threes_to_aucount = {}
    
    for candidate in candidates:    
        three_seq = candidate.seq_prime3
        three_no_au = _remove_trailing(three_seq, "AUau")
        threes_to_aucount[three_no_au] = 0
        
        
    for seq, count in seq_to_frequency.iteritems():
        seq_no_au = _remove_trailing(seq, "AUau")
        threes_to_aucount[seq_no_au] += count
    
    for candidate in candidates:
        three_seq = candidate.seq_prime3
        three_no_au = _remove_trailing(three_seq, "AUau")
        candidate.tailing_au = threes_to_aucount[three_no_au]
        print threes_to_aucount[three_no_au]
        
        
    
            
    
    

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

        
    
            