'''
Created on 19. nov. 2014

@author: hakon
'''


from SuffixTree import SubstringDict
import math
import numpy
from colorama.ansi import Fore

def align_small_seqs(candidates, small_seqs, small_seqs_copies):
    
    print "\n\tSmall seq analysis"
    print "\tcandidates:", len(candidates), candidates[0].hairpin
    print "\tsmall seqs:", len(small_seqs), len(small_seqs_copies)


    max_frequent_seq = 0
    c_used = set()
    s_used = set()

    # adding all canididate seqs to a substring dict
    find_candidates = SubstringDict()

    for i, candidate in enumerate(candidates):
        find_candidates[candidate.hairpin_padded_40] = i

    
    # for each small seq. find the set of matching candidates
    for j, seq in enumerate(small_seqs):
        
        
        candidate_set = find_candidates[seq]
        
        if len(candidate_set) > max_frequent_seq:
            max_frequent_seq = len(candidate_set)
            
        # add the seq. to that candidate
        for cnr in candidate_set:
            copies = small_seqs_copies[j]
            
            if not candidates[cnr].small_subs:
                candidates[cnr].small_subs = {}
                
            candidates[cnr].small_subs[seq] = copies
            
            c_used.add(cnr)
            s_used.add(j)
        
    print "\tmost frequent small hit", max_frequent_seq, "times"
    print "\tsmall seq hits:", len(s_used)
    print "\tcandidates with small seqs:",len(c_used)


def small_seq_stats(candidates):
    '''
    compares the amount of small sequences vs long sequences in 3p and 5p positions
    uses log values of values larger than 1 RPM in each sum
    adds 1.0 to each sum to avoid division against 0-scores
    stores the log score, as some sums are still very large
    '''
    
    
    def log_over_one(val):
        return math.log(val) if val > 1.0 else val
    
    def log_sum(list_to_sum):
        log_list = [log_over_one(x) for x in list_to_sum]
        list_sum = sum(log_list)
        return list_sum
    
    
    has_small_seqs = [c for c in candidates if c.small_subs]
    print "small seqs:", len(has_small_seqs)
    
    for c in has_small_seqs:

        hairpin = c.hairpin_padded_40
        total_subseqs = sum(c.small_subs.values())
        
        start_5 = c.pos_5p_begin - c.hairpin_start + 40 if c.has_5p else 40
        end_5 = c.pos_5p_end - c.hairpin_start + 40 if c.has_5p else start_5 + 24
        
        end_3 = c.pos_3p_end - c.hairpin_start + 40 if c.has_3p else len(hairpin) - 40
        start_3 = c.pos_3p_begin - c.hairpin_start + 40 if c.has_3p else end_3 - 24
        
        if end_5 >= start_3:
            end_5 = (start_5 + end_3) / 2
            start_3 = end_5 + 1
        
        
        area_5p_short = {pos: 0.0 for pos in range(start_5, end_5)}
        area_3p_short = {pos: 0.0 for pos in range(start_3, end_3)}
        
        
        area_5p_1115 = {pos: 0.0 for pos in range(start_5, end_5)}
        area_3p_1115 = {pos: 0.0 for pos in range(start_3, end_3)}
        
        area_5p_long = set(range(start_5, end_5))
        area_5p_long_sum = 0.0
        area_3p_long = set(range(start_3, end_3))
        area_3p_long_sum = 0.0
        
        print "\n---------------"
        print "---------------"
        print "---------------"
        print "---------------"
        print "direction:", c.chromosome_direction
        print c.pos_5p_begin - c.hairpin_start, c.pos_3p_end - c.hairpin_start
        print hairpin
        print " "*start_5 + "5"*(end_5 - start_5) + " "*(start_3 - end_5) + "3"*(end_3 - start_3)

        
        
        for seq in  c.mirBase_matures:
            
            start_seq = hairpin.find(seq)

            print " " * (start_seq) + seq, "\t", (start_seq)

        for seq, copies in c.small_subs.iteritems():
            start_seq = hairpin.find(seq)
            
            if 10 < len(seq) < 16:
                print " " * (start_seq) + seq, "\t", (start_seq, copies)
#             print " " * (start_seq) + seq, "\t", (start_seq, copies)
            
#             end_seq = start_seq + len(seq)
#             if start_seq > end_3:
#                 starts_after.append((start_seq, end_seq))
#                             
#             elif end_seq < start_5:
#                 ends_before.append((start_seq, end_seq))
#             print start_seq, end_seq, copies, seq
            
            if start_seq in area_5p_short:
                area_5p_short[start_seq] += copies
                
                if 10 < len(seq) < 16:
                    area_5p_1115[start_seq] += copies
                    
            elif start_seq in area_3p_short:
                area_3p_short[start_seq] += copies
                
                if 10 < len(seq) < 16:
                    area_3p_1115[start_seq] += copies
                    
            else:
                "seq not in any mature", (start_5, end_5), (start_3, end_3)
#                 assert 0, ((start_5, end_5), (start_3, end_3))
        
#         print "ends before:", start_5, sorted(ends_before) 
#         print "starts after:", end_3, sorted(starts_after) 
#         print len(c.mapped_sequences)

        print "+" * len(hairpin)
        print hairpin
        if len(c.mapped_sequences):
            for i in c.mapped_sequences:
                 
                copies = float(i.data[1].split("-")[1])
                seq_start = i.begin-c.hairpin_start + 40
                
                print " " * seq_start + i.data[2], "\t", (seq_start, copies)
                 
                if seq_start in area_5p_long:
                    area_5p_long_sum += log_over_one(copies)
                elif seq_start in area_3p_long:
                    area_3p_long_sum += log_over_one(copies)
                else:
                    pass
#                     print seq_start,
            
            print
                     
        area_5p_short_sum = log_sum(area_5p_short.values()) + 1.0
        area_3p_short_sum = log_sum(area_3p_short.values()) + 1.0
        area_5p_long_sum += 1.0
        area_3p_long_sum += 1.0

#         has_small_seqs_5p = area_5p_short_sum != 1.0
#         has_small_seqs_3p = area_3p_short_sum != 1.0
        
        ratio_short_long_5p = area_5p_short_sum / area_5p_long_sum if area_5p_long_sum else 0.0
        ratio_short_long_3p = area_3p_short_sum / area_3p_long_sum if area_3p_long_sum else 0.0
        
        ratio_short_long_5p = math.log(ratio_short_long_5p)
        ratio_short_long_3p = math.log(ratio_short_long_3p)
        
        
        
        # only using length 11 to 15 for alignment stats
        # weighted average and weighted variance
        
        if sum(area_5p_1115.values()):
            short_5p_positions = area_5p_1115.keys()
            short_5p_weigths = map(log_over_one, area_5p_1115.values())
            
            avg_5p = numpy.average(short_5p_positions, weights=short_5p_weigths)
            
            offset_5p_squared = [ (p - avg_5p)**2 for p in short_5p_positions]
            biased_variance_5p = numpy.average(offset_5p_squared, weights=short_5p_weigths)
            st_dev_5p = math.sqrt(biased_variance_5p)
            
            
            c.short_seq_5p_stdev = st_dev_5p
            
            if c.has_5p:
                offset = abs(avg_5p - start_5)
                c.short_seq_5p_offset = offset
        
        if sum(area_3p_1115.values()):
            short_3p_positions = area_3p_1115.keys()
            short_3p_weigths = map(log_over_one, area_3p_1115.values())
            
            avg_3p = numpy.average(short_3p_positions, weights=short_3p_weigths)
            
            offset_3p_squared = [(p - avg_3p)**2 for p in short_3p_positions]
            biased_variance_3p = numpy.average(offset_3p_squared, weights=short_3p_weigths)
            st_dev_3p = math.sqrt(biased_variance_3p)
            
            c.short_seq_3p_stdev = st_dev_3p
            
            if c.has_3p:
                offset = abs(avg_3p -  start_3)
                c.short_seq_3p_offset = offset

#         if has_small_seqs_5p:
#             short_5p_positions = area_5p_short.keys()
#             short_5p_weigths = map(log_over_one, area_5p_short.values())
#             
#             avg_5p = numpy.average(short_5p_positions, weights=short_5p_weigths)
#             
#             offset_5p_squared = [ (p - avg_5p)**2 for p in short_5p_positions]
#             biased_variance_5p = numpy.average(offset_5p_squared, weights=short_5p_weigths)
#             st_dev_5p = math.sqrt(biased_variance_5p)
#             
#             
#             c.short_seq_5p_stdev = st_dev_5p
#             
#             if c.has_5p:
#                 offset = abs(avg_5p - start_5)
#                 c.short_seq_5p_offset = offset
#         
#         if has_small_seqs_3p:
#             short_3p_positions = area_3p_short.keys()
#             short_3p_weigths = map(log_over_one, area_3p_short.values())
#             
#             avg_3p = numpy.average(short_3p_positions, weights=short_3p_weigths)
#             
#             offset_3p_squared = [(p - avg_3p)**2 for p in short_3p_positions]
#             biased_variance_3p = numpy.average(offset_3p_squared, weights=short_3p_weigths)
#             st_dev_3p = math.sqrt(biased_variance_3p)
#             
#             c.short_seq_3p_stdev = st_dev_3p
#             
#             if c.has_3p:
#                 offset = abs(avg_3p -  start_3)
#                 c.short_seq_3p_offset = offset
        
        
        c.ratio_short_long_5p = ratio_short_long_5p
        c.ratio_short_long_3p = ratio_short_long_3p
        

#         print min(area_5p_long), max(area_5p_long),
#         print min(area_3p_long), max(area_3p_long)
#         print (c.pos_5p_begin - c.hairpin_start + 40, c.pos_5p_end - c.hairpin_start + 40), (c.pos_3p_begin - c.hairpin_start + 40, c.pos_3p_end - c.hairpin_start + 40)
#         print total_subseqs
#         print area_5p_short_sum, "\t", area_3p_short_sum
#         print area_5p_long_sum, "\t", area_3p_long_sum, c.mapped_sequences
#         print ratio_short_long_5p, "\t", ratio_short_long_3p
#         print c.ratio_short_long_5p
#         print c.ratio_short_long_3p

        
        if sum(area_5p_1115.values()):
            print "\n"
            print (start_5, end_5), (start_3, end_3)
            print c.short_seq_5p_stdev
            print "has 5p?"
            print c.short_seq_5p_offset, avg_5p, start_5, c.has_5p
            print [(k,v) for k,v in area_5p_short.items() if v > 0.0]
             
        if sum(area_3p_1115.values()):
            print "\n"
            print (start_5, end_5), (start_3, end_3)
            print c.short_seq_5p_stdev
            print "has 3p?"
            print c.short_seq_3p_offset, avg_3p, start_3, c.has_3p
            print c.short_seq_3p_offset
            print [(k,v) for k,v in area_3p_short.items() if v > 0.0]
            

#     print len([c for c in candidates if c.ratio_short_long_5p])
#     has_ratio = [math.log(c.ratio_short_long_3p) for c in candidates if c.ratio_short_long_3p]
#     print len(has_ratio)
#     print sum(has_ratio) / len(has_ratio)


# import random
# teste = {random.random():random.randint(1,100) for _ in range(1000)}
# 
# print [(k,v) for k,v in teste.items()]
# print zip(teste.keys(), teste.values())


    
    