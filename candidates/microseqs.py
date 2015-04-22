'''
Created on 19. nov. 2014

@author: hakon
'''


from SuffixTree import SubstringDict
import math
import numpy
from matplotlib import pyplot



def length_distribution(sequences, copies):
    from matplotlib import pyplot
    
    lengths = {}
    lenghts_unique = {}
    for s, c in zip(sequences, copies):
        l = len(s)
        if l in lengths:
            lengths[l] += c
        else:
            lengths[l] = c
            
        if l in lenghts_unique:
            lenghts_unique[l] += 1
        else:
            lenghts_unique[l] = 1
    
    res = sorted(lengths.items())
    lengths, rpm_sum = zip(*res)
    print res
    print lengths, rpm_sum
    
    pyplot.title("Small sequence length distribution")
    pyplot.plot(lengths, rpm_sum)
    pyplot.xlabel("sequence lenght")
    pyplot.ylabel("RPM")
    pyplot.grid(True)
    pyplot.show()


    res = sorted(lenghts_unique.items())
    lengths, uniques = zip(*res)
    print res
    print lengths, uniques
    
    pyplot.title("Short sequence length distribution")
    pyplot.plot(lengths, uniques)
    pyplot.grid(True)
    pyplot.xlabel("sequence lenght")
    pyplot.ylabel("unique sequences")
    pyplot.yticks([5000*x for x in range(0,8)])
    pyplot.show()
    
    relative_copies = [r*1.0 / u for r, u in zip(rpm_sum, uniques)]
    
    pyplot.title("RPM vs unique, short sequences")
    pyplot.plot(lengths, relative_copies)
    pyplot.grid(True)
    pyplot.xlabel("sequence lenght")
    pyplot.ylabel("unique sequences")
#     pyplot.yticks([5000*x for x in range(0,8)])
    pyplot.show()


    relative_copies = [u*1.0 / r for r, u in zip(rpm_sum, uniques)]
    
    pyplot.title("RPM vs unique, short sequences")
    pyplot.plot(lengths, relative_copies)
    pyplot.grid(True)
    pyplot.xlabel("sequence lenght")
    pyplot.ylabel("unique sequences")
#     pyplot.yticks([5000*x for x in range(0,8)])
    pyplot.show()

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
        medium_short = []
        
#         area_5p_1115 = {pos: 0.0 for pos in range(start_5, end_5)}
#         area_3p_1115 = {pos: 0.0 for pos in range(start_3, end_3)}
        
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

        
        if c.mirBase_matures:
            for seq in c.mirBase_matures:
                
                start_seq = hairpin.find(seq)
    
                print " " * (start_seq) + seq, "\t", (start_seq)

        for seq, copies in c.small_subs.iteritems():
            start_seq = hairpin.find(seq)
            end_seq = start_seq + len(seq)
            
            if 5 < len(seq) < 18:
                print " " * (start_seq) + seq, "\t", (start_seq, copies)
                medium_short.append( (len(seq), start_seq, end_seq, copies) )
            
            
            if start_seq in area_5p_short:
                area_5p_short[start_seq] += copies
                
#                 if 10 < len(seq) < 16:
#                     area_5p_1115[start_seq] += copies
                    
            elif start_seq in area_3p_short:
                area_3p_short[start_seq] += copies
                
#                 if 10 < len(seq) < 16:
#                     area_3p_1115[start_seq] += copies
                    
            else:
                "seq not in any mature", (start_5, end_5), (start_3, end_3)

        
        

        print "+" * len(hairpin)
        print hairpin
        
        unscaled_long = 0
        
        
        if len(c.mapped_sequences):
            for i in c.mapped_sequences:
                
                
                copies = float(i.data[1].split("-")[1])
                seq_start = i.begin-c.hairpin_start + 40
                
                print " " * seq_start + i.data[2], "\t", (seq_start, copies)
                 
                unscaled_long += copies
                
                if seq_start in area_5p_long:
                    area_5p_long_sum += log_over_one(copies)
                elif seq_start in area_3p_long:
                    area_3p_long_sum += log_over_one(copies)
                else:
                    pass
#                     print seq_start,
            
            print
            
        unscaled_short = sum(c.small_subs.values()) + 1.0
        unscaled_long += 1.0
        unscaled_short_fraction = unscaled_short / unscaled_long
                     
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
        
        
        starts = (start_5, start_3)
        ends = (end_5, end_3)
#         
#         small_sum = 0.0
#         align_sum = 0.0
#         
#         unscaled_sum = 0.0
#         unscaled_align_sum = 0.0
        
        

        

        distances = []
        # finding distances to start / end of mature sequences
        for (length, start, end, copies) in medium_short:
        
            start_dist = mindist(starts, start)
            end_dist = mindist(ends, end)
            min_dist = min(start_dist, end_dist)
            
            weight = math.log(copies+1.0)
            
            dist_weight_log = weight / (min_dist + 1.0)            
            dist_weigth = copies / (min_dist + 1.0)
            
            distances.append( (length, dist_weigth, dist_weight_log, copies) )
                    
        
        
        dist_11_13 = sum( w for l, w, _lw, c in distances if 9 <= l <= 13 )
        dist_11_14 = sum( w for l, w, _lw, c in distances if 9 <= l <= 14 )
        dist_11_15 = sum( w for l, w, _lw, c in distances if 9 <= l <= 15 )
        dist_11_16 = sum( w for l, w, _lw, c in distances if 9 <= l <= 16 )
        dist_11_17 = sum( w for l, w, _lw, c in distances if 9 <= l <= 17 )
        dist_11_18 = sum( w for l, w, _lw, c in distances if 9 <= l <= 18 )
        
        small_11_13 = sum( c for l, w, _lw, c in distances if 9 <= l <= 13 )
        small_11_14 = sum( c for l, w, _lw, c in distances if 9 <= l <= 14 )
        small_11_15 = sum( c for l, w, _lw, c in distances if 9 <= l <= 15 )
        small_11_16 = sum( c for l, w, _lw, c in distances if 9 <= l <= 16 )
        small_11_17 = sum( c for l, w, _lw, c in distances if 9 <= l <= 17 )
        small_11_18 = sum( c for l, w, _lw, c in distances if 9 <= l <= 18 )

        score_11_13 = dist_11_13 / small_11_13 if small_11_13 else 0
        score_11_14 = dist_11_14 / small_11_14 if small_11_14 else 0
        score_11_15 = dist_11_15 / small_11_15 if small_11_15 else 0
        score_11_16 = dist_11_16 / small_11_16 if small_11_16 else 0
        score_11_17 = dist_11_17 / small_11_17 if small_11_17 else 0
        score_11_18 = dist_11_18 / small_11_18 if small_11_18 else 0

        
        c.short_seq_align_9_13 = score_11_13
        c.short_seq_align_9_14 = score_11_14
        c.short_seq_align_9_15 = score_11_15
        c.short_seq_align_9_16 = score_11_16
        c.short_seq_align_9_17 = score_11_17
        c.short_seq_align_9_18 = score_11_18


#         dist_9_13 = sum( w for l, w, _lw, c in distances if 9 <= l <= 13 )
#         dist_9_14 = sum( w for l, w, _lw, c in distances if 9 <= l <= 14 )
#         dist_9_15 = sum( w for l, w, _lw, c in distances if 9 <= l <= 15 )
#         dist_9_16 = sum( w for l, w, _lw, c in distances if 9 <= l <= 16 )
#         dist_9_17 = sum( w for l, w, _lw, c in distances if 9 <= l <= 17 )
#         dist_9_18 = sum( w for l, w, _lw, c in distances if 9 <= l <= 18 )
#         
#         
#         small_9_13 = sum( c for l, w, _lw, c in distances if 9 <= l <= 13 )
#         small_9_14 = sum( c for l, w, _lw, c in distances if 9 <= l <= 14 )
#         small_9_15 = sum( c for l, w, _lw, c in distances if 9 <= l <= 15 )
#         small_9_16 = sum( c for l, w, _lw, c in distances if 9 <= l <= 16 )
#         small_9_17 = sum( c for l, w, _lw, c in distances if 9 <= l <= 17 )
#         small_9_18 = sum( c for l, w, _lw, c in distances if 9 <= l <= 18 )
# 
#         score_9_13 = dist_9_13 / small_9_13 if small_9_13 else 0
#         score_9_14 = dist_9_14 / small_9_14 if small_9_14 else 0
#         score_9_15 = dist_9_15 / small_9_15 if small_9_15 else 0
#         score_9_16 = dist_9_16 / small_9_16 if small_9_16 else 0
#         score_9_17 = dist_9_17 / small_9_17 if small_9_17 else 0
#         score_9_18 = dist_9_18 / small_9_18 if small_9_18 else 0
# 
#         
#         c.short_seq_align_9_13 = score_9_13
#         c.short_seq_align_9_14 = score_9_14
#         c.short_seq_align_9_15 = score_9_15
#         c.short_seq_align_9_16 = score_9_16
#         c.short_seq_align_9_17 = score_9_17
#         c.short_seq_align_9_18 = score_9_18
#             
#             
#         
#         align_score = align_sum / small_sum if small_sum else 0
#         
#         c.short_seq_5p_stdev = unscaled_align_sum / unscaled_sum  if unscaled_sum else -1
        

#         c.ratio_short_long_5p = ratio_short_long_5p
#         c.ratio_short_long_3p = ratio_short_long_3p

        
        c.ratio_short_long = ratio_short_long_5p + ratio_short_long_3p  
#         c.short_seq_align = align_score
        
        if unscaled_short > 1.0:
            c.ratio_short_long_5p = unscaled_short_fraction

#         print min(area_5p_long), max(area_5p_long),
#         print min(area_3p_long), max(area_3p_long)
#         print (c.pos_5p_begin - c.hairpin_start + 40, c.pos_5p_end - c.hairpin_start + 40), (c.pos_3p_begin - c.hairpin_start + 40, c.pos_3p_end - c.hairpin_start + 40)
#         print total_subseqs
#         print area_5p_short_sum, "\t", area_3p_short_sum
#         print area_5p_long_sum, "\t", area_3p_long_sum, c.mapped_sequences
#         print ratio_short_long_5p, "\t", ratio_short_long_3p
#         print c.ratio_short_long_5p
#         print c.ratio_short_long_3p

        
        if sum(area_5p_short.values()):
            print "\n"
            print (start_5, end_5), (start_3, end_3)
            print c.short_seq_5p_stdev
            print "has 5p?"
            print c.short_seq_5p_offset, start_5, c.has_5p
            print [(k,v) for k,v in area_5p_short.items() if v > 0.0]
             
        if sum(area_3p_short.values()):
            print "\n"
            print (start_5, end_5), (start_3, end_3)
            print c.short_seq_5p_stdev
            print "has 3p?"
            print c.short_seq_3p_offset, start_3, c.has_3p
            print c.short_seq_3p_offset
            print [(k,v) for k,v in area_3p_short.items() if v > 0.0]
            




def mindist(goalpos, pos):
    
    distances = [abs(p - pos) for p in goalpos]
    min_dist = min(distances)
    
    return min_dist
    








#     print len([c for c in candidates if c.ratio_short_long_5p])
#     has_ratio = [math.log(c.ratio_short_long_3p) for c in candidates if c.ratio_short_long_3p]
#     print len(has_ratio)
#     print sum(has_ratio) / len(has_ratio)



        # only using length 11 to 15 for alignment stats
        # weighted average and weighted variance
        
#         if sum(area_5p_1115.values()):
#             short_5p_positions = area_5p_1115.keys()
#             short_5p_weigths = map(log_over_one, area_5p_1115.values())
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
#         if sum(area_3p_1115.values()):
#             short_3p_positions = area_3p_1115.keys()
#             short_3p_weigths = map(log_over_one, area_3p_1115.values())
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








# 
# 
# res =  sorted(test.items())
# 
# a, b =  zip(*res)
# 
# pyplot.plot(a,b)
# pyplot.show()


