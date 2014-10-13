'''
Created on 8. okt. 2014

@author: hakon
'''


def three_prime_au(candidates, sequence_frequency):
    
    threes_to_aucount = {}
    s = 0
    
    for candidate in candidates:    
        three_seq = candidate.seq_prime3
        three_seq = _remove_trailing(three_seq, "AUau")
        threes_to_aucount[three_seq] = 0
        
        
    for seq, count in sequence_frequency:
        seq_no_au = _remove_trailing(seq, "AUau")
        if seq_no_au in threes_to_aucount:
            threes_to_aucount[seq_no_au] += count

    
    for candidate in candidates:
        three_seq = candidate.seq_prime3
        three_seq = _remove_trailing(three_seq, "AUau")
        
        # set tailing value for all candidates
        candidate.tailing_au = threes_to_aucount[three_seq]
        print threes_to_aucount[three_seq]
        s += threes_to_aucount[three_seq]
        
    print s
            
    
    

def _remove_trailing(sequence, chars):
    chars = set(chars)
    cut = len(sequence)
    
    for x in xrange(len(sequence)-1, 0, -1):
        if sequence[x] not in chars:
            break
        cut = x
    if cut is not len(sequence):
        print sequence
        print sequence[:cut]
    return sequence[:cut]




# t = "QWEWRRTYUUUU"
# t = "QWEWRRTAYAUAUT"
# print _remove_trailing_AU(t)

        
    
            