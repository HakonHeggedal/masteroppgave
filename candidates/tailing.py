'''
Created on 8. okt. 2014

@author: hakon
'''

from _collections import defaultdict



def tailing_au(candidates, sequences, align_len=8):
    ''' counts the number of a/u-tailing sequences
        compared to all mapping sequences
    '''
    print "\n\tcalculating tailing for the candidates"
    full_hits = [0]*len(candidates)
    tailing = [0]*len(candidates)   
    last_parts = defaultdict(list)
    
    for nr, candidate in enumerate(candidates):
        # store parts at the end of 5' and 3'
        candidate_seq = candidate.hairpin
        candidate_len = candidate.pos_3p_end - candidate.pos_3p_begin
        start = candidate.padding_size + candidate_len - 1
        
        for i in xrange(start, start+4, 1):
            seq = candidate_seq[i:i+align_len]
            last_parts[seq].append(nr)
        
        start = candidate.padding_size + candidate.pos_5p_end - candidate.pos_5p_begin - 1
        
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

                if seq in candidates[c_nr].hairpin:
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
                    
#                     five = candidates[c_nr].hairpin[:candidate_len]
#                     print "yes", count, seq, tailfree, five, candidates[c_nr].hairpin[-10:]
                    full_hits[c_nr] += count

                    
                    # sequence matches this candidate
            
            

            
    tailcount = 0
    for i, (candidate, tails, hits) in enumerate(zip(candidates, tailing, full_hits)):
        tailing_val = tails * 1.0 / (hits+1)
        if tails > 0:
            tailcount += 1
#             print i, "yeah\t", hits, tails, tailing_val

        candidate.set_tailing(tailing_val)


    
    print "\tnr of canididates with tails:", tailcount, tailcount*1.0/len(candidates)


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






# t = "QWEWRRTYUUUU"
# t = "QWEWRRTAYAUAUT"
# print _remove_trailing_AU(t)

        
    
            