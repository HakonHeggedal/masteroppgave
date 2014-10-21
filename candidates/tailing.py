'''
Created on 8. okt. 2014

@author: hakon
'''

from _collections import defaultdict



def tailing_au(candidates, sequences, align_len=8):
    ''' counts the number of a/u-tailing sequences
        compared to all mapping sequences
    '''
    print "\ttailing------------------------------"
    full_hits = [0]*len(candidates)
    tailing = [0]*len(candidates)   
    last_parts = defaultdict(list)
    
    for nr, candidate in enumerate(candidates):
        # store parts at the end of 5' and 3'
        candidate_seq = candidate.hairpin_padded
        candidate_len = candidate.pos_3_end - candidate.pos_3_begin
        start = candidate.padding_size + candidate_len - 1
#         print "starting at", start, candidate.padding_size, candidate_len
        
        for i in xrange(start, start+4, 1):
            seq = candidate_seq[i:i+align_len]
            last_parts[seq].append(nr)
        
        start = candidate.padding_size + candidate.pos_5_end - candidate.pos_5_begin - 1
        
        for i in xrange(start, start+4, 1):
            seq = candidate_seq[i:i+align_len]
            last_parts[seq].append(nr)
            

    for sequence in sequences:
        # find sequences that map to the ends of candidates
        
        seq = sequence.nucleotides
        count = sequence.duplicates
        
        if len(seq) < 18:
            continue
        
        seq_end = seq[-align_len:]
        has_matches = False
        if seq_end in last_parts:
            # all non-tailing hits
            
            for c_nr in last_parts[seq_end]:

                if seq in candidates[c_nr].hairpin_padded:
                    full_hits[c_nr] += count
                    has_matches = True
                
#             if not has_matches:
#                 pass
#                 print "stupid", seq
            
        
        tailfree, is_tail =  _remove_trailing(seq, "AUau")
        
        if not is_tail or has_matches:
            continue
        
        tailfree_end = tailfree[-align_len:]
        
        if tailfree_end in last_parts:
            # tailing AU here :)
#             print "got tail", seq, tailfree_end
            
            for c_nr in last_parts[tailfree_end]:
                
                
#                 print "\t", candidates[c_nr].hairpin
                if tailfree in candidates[c_nr].hairpin:
                    tailing[c_nr] += count
                    
                    five = candidates[c_nr].hairpin[:candidate_len]
                    print "yes", count, seq, tailfree, five, candidates[c_nr].hairpin[-10:]
                    full_hits[c_nr] += count

                    
                    # sequence matches this candidate
            
            

            
        
    for candidate, tails, hits in zip(candidates, tailing, full_hits):
        tailing_val = tails * 1.0 / (hits+1)
        if tails > 0:
            print "yeah", hits, tails
            print tailing_val
        candidate.set_tailing(tailing_val)



    print "\ttailing finished^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"


def _remove_trailing(sequence, chars):
    chars = set(chars)
    cut = len(sequence)
    iscut = False
    
    for x in xrange(len(sequence)-1, 0, -1):
        if sequence[x] not in chars:
            break
        cut = x
        iscut = True

    return sequence[:cut], iscut


































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





# t = "QWEWRRTYUUUU"
# t = "QWEWRRTAYAUAUT"
# print _remove_trailing_AU(t)

        
    
            