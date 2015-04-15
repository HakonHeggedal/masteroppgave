'''
Created on 8. okt. 2014

@author: hakon
'''

from _collections import defaultdict
import math
import pickle

from candidates.interval_tree_misc import reverse_compliment





def tailing_au_simple(candidates=None, sequences=None):
    
#     pickle.dump(candidates, open("candidates_test_only.p", "wb"))
#     pickle.dump(sequences, open("sequences_test_only.p", "wb"))
#     
#     assert False
# 
#     candidates = candidates[:50]
#     pickle.dump(candidates, open("candidates_test_only2.p", "wb"))
#     assert 0
    
#     print "loading data ..."
#     candidates = pickle.load( open("candidates_test_only2.p", "rb"))
#     sequences = pickle.load( open("sequences_test_only.p", "rb"))
#     print "... loaded"

    
    tail_au = []
    pre_au = []
    tail_a_list = []
#     pre_a = []
    
    min_tailfree = 15
    
    end_chars = {"A":0.0, "C":0.0, "G":0.0, "T":0.0}
    
    for s in sequences:
        n = seq = s.nucleotides[-1]
        end_chars[n] += 1
        
    print end_chars
#     assert 0
    
    
    for sequence in sequences:
        
        seq = sequence.nucleotides
        count = sequence.duplicates if sequence.duplicates <= 1 else math.log(sequence.duplicates)
        
        tailfree, is_tail_end =  _remove_trailing(seq)
        au_start, is_au_start = _remove_trailing_start(seq)
        
        if is_tail_end and len(tailfree) > min_tailfree:
            assert len(tailfree) != len(seq)
            tail_au.append((tailfree, seq, count))
            
        if is_au_start and len(au_start) > min_tailfree:
            assert len(au_start) != len(seq)
            pre_au.append((au_start, seq, count))
            
            
        tail_a, is_tail_a = _remove_trailing(seq, set("Aa"))
        
        if is_tail_a and tail_a > min_tailfree:
            tail_a_list.append((tail_a, seq, count))
            
    
    print "long tails", len([1 for (au, seq, count) in tail_au if len(seq) - len(au) > 1])
    print "long pre", len([1 for (au, seq, count) in pre_au if len(seq) - len(au) > 1])
    print "found all tails and starts", len(tail_au), len(pre_au), "only a ends:", len(tail_a_list)
    
    
#     assert 0
            
    for c in candidates:
        
        hairpin = c.hairpin_padded_40[20:-20]
        
        if c.chromosome_direction == "-":
            hairpin = reverse_compliment(hairpin)
        
        
        full_hits_tail = 0.0
        full_hits_pre = 0.0
        full_hits_bowtie = 0.0
        
        au_tails = 0.0
        au_starts = 0.0
        hit_tailing = []
        hit_pre = []
        
        for notail, seq, count in tail_au:
            if len(notail) < 15:
                continue
            
            if notail in hairpin:
                
                
                if seq not in hairpin:
                    pos = hairpin.find(notail)
                    assert pos >= 0, (notail in hairpin, pos, len(seq), notail, seq[len(notail):], hairpin[pos:len(seq)+5])
#                     print "\t", notail in hairpin, pos, len(seq), notail, seq[len(notail):], hairpin[pos:len(seq)+5]
                    au_tails += count
                    hit_tailing.append((pos, pos+len(notail)))
                else:
                    full_hits_tail += count
                        
        
        for no_au_start, seq, count in pre_au:
             
            if no_au_start in hairpin:
                if seq not in hairpin:
                    pos = hairpin.find(no_au_start)
                    assert pos >= 0, (no_au_start in hairpin, pos, len(seq), no_au_start, seq[len(no_au_start):], hairpin[pos:len(seq)+5])
                    
                    au_starts += count
                    hit_pre.append( (pos, pos+len(no_au_start)) )
                else:
                    full_hits_pre += count
                    
        
#         full_hits_bowtie = sum([float(s.data[1].split("-")[1])] for s in c.mapped_sequences)
        
        full_hits_missed = 0
        for sequence in sequences:
            
            seq = sequence.nucleotides
            count = sequence.duplicates if sequence.duplicates <= 1 else math.log(sequence.duplicates)
            
            if seq in hairpin:
                full_hits_missed += count
            
            
        
        
        
        full_hits_bowtie = 0
        for sequence in c.mapped_sequences:
            count = float(sequence.data[1].split("-")[1])
            full_hits_bowtie += math.log(count) if count > 1 else count
            
                    
        print c.chromosome_direction, len(hairpin)
        print au_tails, full_hits_tail, sorted(hit_tailing)
        print au_starts, full_hits_pre, sorted(hit_pre)
        print full_hits_bowtie, full_hits_missed
        print
        
        
        total_hits_tail = full_hits_bowtie + full_hits_missed + au_tails
        total_hits_tail = math.log(total_hits_tail) if total_hits_tail > 1 else total_hits_tail
        
        au_tails_score = math.log(au_tails) if au_tails > 1 else au_tails
        
        score_tailing = (au_tails_score + 1) / (total_hits_tail + 1)
        
        
        total_hits_tail = full_hits_bowtie + full_hits_missed + au_starts
        total_hits_tail = math.log(total_hits_tail) if total_hits_tail > 1 else total_hits_tail
        
        
        au_leading_score = math.log(au_starts) if au_starts > 1 else au_starts
        
        score_leading = (au_leading_score + 1) / (total_hits_tail + 1)
        
        
        c.tailing_au = score_tailing
        c.leading_au = score_leading
        
#         c.set_tailing(score)
        
        
    
        
    
#     for candidate in candidates:
#         
#         hairpin = candidate.hairpin_padded_40[20:-20]
# #         hairpin_padded = candidate.hairpin_padded_40
#         
# #         print "\n---------"
# #         print hairpin
# 
#         full_hits = 0.0
#         tail_end_hits = 0.0
#         tail_start_hits = 0.0
#         
#         for sequence in sequences:
#             # find sequences that map to the ends of candidates
#             
#             seq = sequence.nucleotides
#             count = sequence.duplicates
#             
#             count = math.log(count) if count > 1 else count
#             
# 
#             
#             if seq in hairpin:
#                 full_hits += count
#             else:
#                 
#                 tailfree, is_tail_end =  _remove_trailing(seq, "AUau")
#                 
#                 if is_tail_end and len(tailfree) > 8:
#                     if tailfree in hairpin:
#                         tail_end_hits += count
#                 else:
#                     pass
#                     
#                     tail_free_start, is_au_start = _remove_trailing_start(seq)
#                     if is_au_start and len(tail_free_start) > 8:
#                         if tail_free_start in hairpin:
#                             tail_start_hits += count
#                             print "\tstart:", count, tail_free_start[:8], seq[:8], seq in hairpin, hairpin
#                 
#         
# 
#         print full_hits, 
#         print tail_end_hits, 
#         print tail_start_hits
#         
#         assert not tail_start_hits
                
                





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
        
        count = math.log(count) if count > 1 else count
        
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
            
        
        tailfree, is_tail =  _remove_trailing(seq, set("AUau"))
        
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
        
        assert tails >= 0, tails
        assert hits >= 0, hits
        tailing_val = (tails+1.0) / (hits+1.0)
        if tails > 0:
            tailcount += 1
#             print i, "yeah\t", hits, tails, tailing_val

        candidate.set_tailing(tailing_val)



    
    print "\tnr of canididates with tails:", tailcount, tailcount*1.0/len(candidates)


def _remove_trailing(sequence, chars=set("Aa")):
#     chars = set(chars)
    cut = len(sequence)
    iscut = False
    
    for x in xrange(len(sequence)-1, 0, -1):
        if sequence[x] not in chars:
            break
        cut = x
        iscut = True

    return sequence[:cut], iscut


def _remove_trailing_start(sequence, chars=set("AUTaut")):
    
    i = 0
    for i, c in enumerate(sequence):
        if c not in chars:
            break
    
    return sequence[i:], i != 0




# tailing_au_simple()

# 
# print _remove_trailing_start("AUBCDE")
# print _remove_trailing_start("BAUBCDE")

# print _remove_trailing("GCAU")

    


# t = "QWEWRRTYUUUU"
# t = "QWEWRRTAYAUAUT"
# print _remove_trailing_AU(t)

        
    
            