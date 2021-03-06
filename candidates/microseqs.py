'''
Created on 19. nov. 2014

@author: hakon
'''


from SuffixTree import SubstringDict
import math
import numpy
from matplotlib import pyplot
from candidates.interval_tree_misc import reverse_compliment



def length_distribution(sequences, copies):

    
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
    name = "Short sequence length distribution RPM"
    save_name = "short_rpm_dist"
    pyplot.title(name)
    pyplot.plot(lengths, rpm_sum)
    pyplot.xlabel("sequence lenght")
    pyplot.ylabel("RPM")
    pyplot.grid(True)
    pyplot.show()
    pyplot.savefig("figures/" + save_name + ".png")
    
    
    res = sorted(lenghts_unique.items())
    lengths, uniques = zip(*res)
    print res
    print lengths, uniques
    
    name = "Short unique sequence length distribution"
    save_name = "short_unique_dist"
    pyplot.title(name)
    pyplot.plot(lengths, uniques)
    pyplot.grid(True)
    pyplot.xlabel("sequence lenght")
    pyplot.ylabel("unique sequences")
    pyplot.yticks([5000*x for x in range(0,8)])
    pyplot.show()
    pyplot.savefig("figures/" + name + ".png")
    
    relative_copies = [r*1.0 / u for r, u in zip(rpm_sum, uniques)]
    name = "RPM vs unique, short sequences"
    pyplot.title(name)
    pyplot.plot(lengths, relative_copies)
    pyplot.grid(True)
    pyplot.xlabel("sequence lenght")
    pyplot.ylabel("unique sequences")
#     pyplot.yticks([5000*x for x in range(0,8)])
    pyplot.show()
    pyplot.savefig("figures/" + name + ".png")


    relative_copies = [u*1.0 / r for r, u in zip(rpm_sum, uniques)]
    name = "unique vs RPM, short sequences"
    pyplot.title(name)
    pyplot.plot(lengths, relative_copies)
    pyplot.grid(True)
    pyplot.xlabel("sequence lenght")
    pyplot.ylabel("unique sequences")
#     pyplot.yticks([5000*x for x in range(0,8)])
    pyplot.show()
    pyplot.savefig("figures/" + name + ".png")

def align_small_seqs(candidates, small_seqs, small_seqs_copies):
    
    print "\n\tSmall seq analysis"
    print "\tcandidates:", len(candidates), candidates[0].hairpin
    print "\tsmall seqs:", len(small_seqs), len(small_seqs_copies)


    max_frequent_seq = 0
    c_used = set()
    s_used = set()

    # adding all canididate seqs to a substring dict
    # reversing the hairpin if it is on negative strand
    
    find_candidates = SubstringDict()

    for i, candidate in enumerate(candidates):
        
        hairpin_part = candidate.hairpin_padded_40[20:-20]
        if candidate.chromosome_direction == "-":
            hairpin_part = reverse_compliment(hairpin_part)
        find_candidates[hairpin_part] = i

    
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


def small_seq_stats(candidates, sc=0, na=0, c_to_m=0):
    '''
    compares the amount of small sequences vs long sequences in 3p and 5p positions
    uses log values of values larger than 1 RPM in each sum
    adds 1.0 to each sum to avoid division against 0-scores
    stores the log score, as some sums are still very large
    '''
    
#     na = [x for li in na for x in li] # unpack folds
#     sc = zip(*sc)
#     print len(na), na[0], len(set(na))
#     print len(sc), sc[0]
#     
#     names_to_scores = {n:s for n,s in zip(na,sc)}
#     
#     print na[0]
#     print sc[0]
#     print names_to_scores[na[0]]
#     
#     
#     print len(names_to_scores)
    
    
#     def get_name_scores(c):
#         hashval = c.chromosome+c.chromosome_direction+str(c.hairpin_start)
#         
#         mi_name = c_to_m[hashval]
#         scores = names_to_scores[mi_name] if mi_name in names_to_scores else ""
#         return mi_name, scores
    
#     def _is_miRNA(c):
#         hashval = c.chromosome+c.chromosome_direction+str(c.hairpin_start)
#         return hashval in candidate_to_miRNA
    
    def log_over_one(val):
        return math.log(val) if val > 1.0 else val
     
    def log_sum(list_to_sum):
        log_list = [log_over_one(x) for x in list_to_sum]
        list_sum = sum(log_list)
        return list_sum
    
    print "candidates:", len(candidates)
    has_small_seqs = [c for c in candidates if c.small_subs]
    print "small seqs:", len(has_small_seqs)
    
    padding = 20
    
    for c in has_small_seqs:    
        
        hairpin = c.hairpin_padded_40[padding:-padding]
        folding = c.hairpin_fold_40[padding:-padding]
        if c.chromosome_direction == "-":
            hairpin = reverse_compliment(hairpin)
        
        begin_5 = c.pos_5p_begin - c.hairpin_start + padding if c.has_5p else c.estimate_5pb
        end_5 = c.pos_5p_end - c.hairpin_start + padding if c.has_5p else c.estimate_5pe
        
        begin_3 = c.pos_3p_begin - c.hairpin_start + padding if c.has_3p else c.estimate_3pb
        end_3 = c.pos_3p_end - c.hairpin_start + padding if c.has_3p else c.estimate_3pe
        
        
        # swap positions if hairpin is reversed
        has_5p = c.has_5p
        has_3p = c.has_3p
        
        if c.chromosome_direction == "-":
            b5 = len(hairpin) - end_3
            e5 = len(hairpin) - begin_3
            b3 = len(hairpin) - end_5
            e3 = len(hairpin) - begin_5
            
            begin_5 = b5 
            end_5 = e5 
            begin_3 = b3
            end_3 = e3
            
            #swap 3p 5p positions
            has_5p, has_3p = has_3p, has_5p
        
        # filter out bad positions
        if not (begin_5 < end_5 < begin_3 < end_3):
            if has_5p:
                if begin_3 < end_5 or end_3 < end_5:
                    begin_3 = -1
                    end_3 = -1
            elif has_3p:
                if end_5 > begin_3 or begin_5 > begin_3:
                    begin_5 = -1
                    end_5 = -1
            
            if begin_5 > end_5:
                begin_5 = -1
                end_5 = -1
            if begin_3 > end_3:
                begin_3 = -1
                end_3 = -1

#         assert begin_5 < end_5 < begin_3 < end_3, (begin_5, end_5, begin_3, end_3)
        
        mature_pos = [pos for pos in [begin_5, end_5, begin_3, end_3] if 0 <= pos <= len(hairpin)]     
        
        medium_short = []
        

         
        print "\n---------------"
        print "---------------"
        print [begin_5, end_5, begin_3, end_3]
        print mature_pos
        print c.has_5p
        print c.has_3p
        print "direction:", c.chromosome_direction
        print c.pos_5p_begin - c.hairpin_start, c.pos_3p_end - c.hairpin_start
#         print get_name_scores(c)
        print folding
        print hairpin
        print " "*begin_5 + "5"*(end_5 - begin_5) + " "*(begin_3 - end_5) + "3"*(end_3 - begin_3)
 
         
        if c.mirBase_matures:
            for seq in c.mirBase_matures:
                
                if c.chromosome_direction == "-":
                    seq = reverse_compliment(seq)
                
                start_seq = hairpin.find(seq)
                print " " * (start_seq) + seq, "\t", (start_seq)

        for seq, copies in c.small_subs.iteritems():
            start_seq = hairpin.find(seq)
            end_seq = start_seq + len(seq)
            
            print " " * (start_seq) + seq, "\t", (start_seq, copies)
            medium_short.append( (len(seq), start_seq, end_seq, copies) )
        
        
#         print "+" * len(hairpin)
#         print hairpin
        
        unscaled_long = 0
        scaled_long = 0
        
        if len(c.mapped_sequences):
            for i in c.mapped_sequences:
                
                
                copies = float(i.data[1].split("-")[1])

                seq = i.data[2]
                if c.chromosome_direction == "-":
                    seq = reverse_compliment(seq)

                seq_start = hairpin.find(seq)
                
                print " " * seq_start + seq, "\t", ( seq_start, copies)
                 
                unscaled_long += copies
                scaled_long += log_over_one(copies)
                
#===============================================================================
# #                 if seq_start in area_5p_long:
# #                     area_5p_long_sum += log_over_one(copies)
# #                 elif seq_start in area_3p_long:
# #                     area_3p_long_sum += log_over_one(copies)
#===============================================================================



            
        unscaled_short = sum(c.small_subs.values()) + 1.0
        unscaled_long += 1.0
        unscaled_short_fraction = unscaled_short / unscaled_long
        
        
        scaled_short = log_sum(c.small_subs.values()) + 1.0
        scaled_long += 1
        scaled_short_fraction = scaled_short / scaled_long
        
        
        c.ratio_short_long = unscaled_short_fraction
        c.ratio_short_long_logval = scaled_short_fraction
        
        
        def nondecreasing(l):
            for i, el in enumerate(l[1:]):
        
                if l[i] > el:
                    return False
            return True
        
        assert nondecreasing(mature_pos), (mature_pos, len(hairpin))


        

        

        distances = []
        # finding distances to start / end of mature sequences
        for (length, start, end, copies) in medium_short:
        
            start_dist = min_dist(mature_pos, start)
            end_dist = min_dist(mature_pos, end)
            smallest_dist = min(start_dist, end_dist)
            
            weight = math.log(copies+1.0)
            
            dist_weight_log = weight / (smallest_dist + 1.0)            
            dist_weigth = copies / (smallest_dist + 1.0)
            
            distances.append( (length, dist_weigth, dist_weight_log, copies) )
                    
        
        
        def score_distance(minlen, maxlen):
            
            distance_weighted = sum( w for l, w, _lw, c in distances if minlen <= l <= maxlen )
            all_weights = sum( c for l, w, _lw, c in distances if minlen <= l <= maxlen )
            score = distance_weighted / all_weights if all_weights else 0
            
            return score
            
        
        
        c.short_seq_align_10_13 = score_distance(10,13)
        c.short_seq_align_10_14 = score_distance(10,14)
        c.short_seq_align_10_15 = score_distance(10,15)
        c.short_seq_align_10_16 = score_distance(10,16)
        c.short_seq_align_10_17 = score_distance(10,17)
        c.short_seq_align_10_18 = score_distance(10,18)
        
        

#         c.short_seq_align_8_17 = score_distance(8,8)
#         c.short_seq_align_9_17 = score_distance(9,9)
#         c.short_seq_align_10_17 = score_distance(10,10)
#         c.short_seq_align_11_17 = score_distance(11,11)
#         c.short_seq_align_12_17 = score_distance(12,12)
#         c.short_seq_align_13_17 = score_distance(13,13)
#         c.short_seq_align_14_17 = score_distance(14,14)
#         c.short_seq_align_15_17 = score_distance(15,15)
#         c.short_seq_align_16_17 = score_distance(16,16)
#         c.short_seq_align_17_17 = score_distance(17,17)
        
        c.short_seq_align_8_17 = score_distance(8,17)
        c.short_seq_align_9_17 = score_distance(9,17)
        c.short_seq_align_10_17 = score_distance(10,17)
        c.short_seq_align_11_17 = score_distance(11,17)
        c.short_seq_align_12_17 = score_distance(12,17)
        c.short_seq_align_13_17 = score_distance(13,17)
        c.short_seq_align_14_17 = score_distance(14,17)
        c.short_seq_align_15_17 = score_distance(15,17)
        c.short_seq_align_16_17 = score_distance(16,17)
        c.short_seq_align_17_17 = score_distance(17,17)
        
        c.short_seq_align = score_distance(13,17)
        




def min_dist(goalpos, pos):
    
    distances = [abs(p - pos) for p in goalpos]
    mindist = min(distances)
    
    return mindist








# 
# print nondecreasing([1,2,3,3])
# 
# print nondecreasing([3,2,1])

# class derpy:
#     def __init__(self):
#         self.derp = 0
# 
# a = derpy()
# a.newstuff = 123
# 
# print a.newstuff  



