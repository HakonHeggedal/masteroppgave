'''
Created on 10. okt. 2014

@author: hakon
'''
import itertools



def heterogenity(candidates, offset=20):
    
    print "\n\tcalculating variation in 5p and 3p start positions"
    
    for candidate in candidates:
        sequences = candidate.mapped_sequences

        start = candidate.hairpin_start
        
        # relative positions
        five_b = candidate.pos_5p_begin - start
        five_e = candidate.pos_5p_end - start
        three_b = candidate.pos_3p_begin - start
        three_e = candidate.pos_3p_end - start
        
        # all start and end positions tuples
        start_positions = sorted([(s.begin-start, float(s.data[1].split("-")[1]) ) for s in sequences])
        end_positions = sorted([(s.end-start, float(s.data[1].split("-")[1]))  for s in sequences])
        
        def sumgroups(group):
            return (group[0], sum([x[1] for x in group[1]]) )
        
        # merge tuples on the same position
        grouped_starts = list(map(sumgroups, itertools.groupby(start_positions, lambda x: x[0])))
        grouped_ends = list(map(sumgroups, itertools.groupby(end_positions, lambda x: x[0])))
        
        
        # find all positions and values in same area (max 20 away) to 5p and 3p start/end
        close_to_5begin = [g for g in grouped_starts if abs(g[0]-five_b) <= offset ]
        close_to_5end = [g for g in grouped_ends if abs(g[0]-five_e) <= offset ]
        close_to_3begin = [g for g in grouped_starts if abs(g[0]-three_b) <= offset ]
        close_to_3end = [g for g in grouped_ends if abs(g[0]-three_e) <= offset ]
        
        
        sum_peaks_5b = sum(g[1] for g in close_to_5begin)
        sum_peaks_5e = sum(g[1] for g in close_to_5end)
        sum_peaks_3b = sum(g[1] for g in close_to_3begin)
        sum_peaks_3e = sum(g[1] for g in close_to_3end)
        
        five_b_peak = -1
        five_e_peak = -1
        three_b_peak = -1
        three_e_peak = -1
        
        for pos, val in close_to_5begin:
            if pos == five_b:
                five_b_peak = val
                break
        
        for pos, val in close_to_5end:
            if pos == five_e:
                five_e_peak = val
                break
            
        for pos, val in close_to_3begin:
            if pos == three_b:
                three_b_peak = val
                break
            
        for pos, val in close_to_3end:
            if pos == three_e:
                three_e_peak = val
                break

        
        #compute heterogenity
        het_5b_values = [abs(pos-five_b)*val for (pos,val) in close_to_5begin]
        het_5b_sum = sum(het_5b_values) + 1
        het_5b = het_5b_sum / (sum_peaks_5b+1) if sum_peaks_5b else -1
        
        het_5e_values =  [abs(pos-five_e)*val for (pos,val) in close_to_5end]
        het_5e_sum = sum(het_5e_values) + 1
        het_5e = het_5e_sum / (sum_peaks_5e+1) if sum_peaks_5e else -1

        het_3b_values = [abs(pos-three_b)*val for (pos,val) in close_to_3begin]
        het_3b_sum = sum(het_3b_values) + 1
        het_3b = het_3b_sum / (sum_peaks_3b+1) if sum_peaks_3b else -1
             
        het_3e_values = [abs(pos-three_e)*val for (pos,val) in close_to_3end]
        het_3e_sum = sum(het_3e_values) + 1
        het_3e =  het_3e_sum / (sum_peaks_3e+1) if sum_peaks_3e else -1
        

        if five_b_peak != candidate.peak_5b or five_e_peak != candidate.peak_5e or three_b_peak != candidate.peak_3b or three_e_peak != candidate.peak_3e:
            print
            print five_b, five_e, three_b, three_e
            print candidate.peak_5b, candidate.peak_5e, candidate.peak_3b, candidate.peak_3e
            print five_b_peak, five_e_peak, three_b_peak, three_e_peak
            print "candidate pos / peaks:", candidate.chromosome
            print "all intervals pos + value"
            print sorted([(s.begin-start, float(s.data[1].split("-")[1]), s.end-start ) for s in sequences])
            print "starts - grouped"
            print grouped_starts
            print "ends: - grouped"
            print grouped_ends
            print "close to pos, heterogenity"
            print five_b, close_to_5begin
            print five_e, close_to_5end
            print three_b, close_to_3begin
            print three_e, close_to_3end
            print "heterogenity position values"
            print het_5b_values
            print het_5e_values
            print het_3b_values
            print het_3e_values
            print "results??"
            print (het_5b_sum, het_5e_sum), (het_3b_sum, het_3e_sum)
            print (sum_peaks_5b, sum_peaks_5e, sum_peaks_3b, sum_peaks_3e)
            print (five_b_peak, five_e_peak), (three_b_peak, three_e_peak)
            print (het_5b, het_5e), (het_3b, het_3e)
        
        
        
        candidate.set_heterogenity(het_5b, het_5e, het_3b, het_3e)
        #candidate.set_reads(reads_5b, reads_5e, reads_3b, reads_3e)


    #assert False
        

'''     
        peak_size = end - start + offset + offset + 1
        begins = [0] * peak_size
        ends = [0] * peak_size
        
        h_5p_begin = -1
        h_5p_end = -1
        h_3p_begin = -1
        h_3p_end = -1
        
        
        # making arrays of start and end positions 
        count_sum = 0
        for sequence in sequences:
            count = float(sequence.data[1].split("-")[1])
            begin_pos = sequence.begin - start + offset
            end_pos = sequence.end - start + offset
            
            if begin_pos < 0 or end_pos >= len(ends):
#                 print "out of range", begin_pos, end_pos
                continue
            
            begins[begin_pos] += count
            ends[end_pos] += count
            count_sum += count
        
#         print begins
#         print ends
#         print
        
        if five_end != -1 and five_begin != -1:
            
            
            # 5p start
            het = 0
            reads_5b = 0
            i = 0
            
            s = begins[:2*offset + 1] # feil...
            
            for i, count in enumerate(s):
                off = abs(offset - i)
                het += off * count
                reads_5b += count
     
            h_5p_begin = begins[offset] / (reads_5b + 1.0)
            
#             if begins[offset] != max(s):
            print 
            print candidate.candidate_type, "five start"
            print (five_begin - start), "- -"
            print begins[offset], max(s), candidate.peak_5b
            print s
            print "het:", h_5p_begin
#             assert begins[offset] == max(s)
    
    
            # 5p end
            het = 0
            reads_5e = 0
            startpos = five_end - start
            
            s = ends[startpos:startpos+2*offset+1]
            
            for i, count in enumerate(s):
                off = abs(offset - i)
                het += off * count
                reads_5e += count
            
            h_5p_end = ends[startpos+offset] / (reads_5e + 1.0)
            

#             if ends[startpos+offset] != max(s):
            print
            print candidate.candidate_type, "five end"
            print (five_begin - start), five_end, "- -"
            print ends[startpos+offset], max(s), candidate.peak_5e
            print s
            print "het:", h_5p_end
#             assert ends[startpos+offset] == max(s)


        if three_begin != -1 and three_end != -1:
            # three begin
            het = 0
            reads_3b = 0
            startpos = three_begin - start 
            s = begins[startpos:startpos+2*offset+1]
    
            for i, count in enumerate(s):
                off = abs(offset - i)
                het += off * count
                reads_3b += count
    
            h_3p_begin = begins[startpos+offset] / (reads_3b + 1.0)

            #if begins[startpos+offset] != max(s):
            print
            print candidate.candidate_type, "three start"
            print "- -", (three_begin - start), (three_end, end)
            print begins[startpos+offset], max(s), candidate.peak_3b
            print s
            print "het:", h_3p_begin
            assert begins[startpos+offset] == max(s)

            #three end
            het = 0
            reads_3e = 0
            s = ends[-(2*offset + 1):] # feil...

            for i, count in enumerate(s):
                off = abs(offset - i)
                het += off * count
                reads_3e += count

            h_3p_end = ends[-offset-1] / (reads_3e + 1)

            #if ends[-offset-1] != max(s):
            print
            print candidate.candidate_type, "three end"
            print "- -", (three_end - end)
            print ends[-offset-1], max(s), candidate.peak_3e
            print s
            print "het:", h_3p_end
            assert ends[-offset-1] == max(s)
'''
        
#         print h_5p_begin, h_5p_end, h_3p_begin, h_3p_end
#         print reads_5b, reads_5e, reads_3b, reads_3e
#         print
        


def _old_heterogenity(candidates, offset=5):
    
    print "\n\tcalculating variation in 5p and 3p start positions"
    
    for candidate in candidates:
        sequences = candidate.mapped_sequences
        

        
        
        five_begin = candidate.pos_5p_begin
        five_end = candidate.pos_5p_end
        three_begin = candidate.pos_3p_begin
        three_end = candidate.pos_3p_end
        start = candidate.hairpin_start
        end = candidate.hairpin_end
        
#         start_positions = sorted([(s.begin-start, float(s.data[1].split("-")[1]) ) for s in sequences])
#         end_positions = sorted([(s.end-start, float(s.data[1].split("-")[1]))  for s in sequences])
#         
#         grouped_starts = list(map(sumgroups, itertools.groupby(start_positions, lambda x: x[0])))
#         grouped_ends = list(map(sumgroups, itertools.groupby(end_positions, lambda x: x[0])))
#         
#         print
#         print "starts:"
# #         print start_positions
#         print grouped_starts
#         print five_begin-start, three_begin-start
#         print "ends:"
# #         print end_positions
#         print grouped_ends
#         print five_end-start, three_end-start
#         print "-",
        
#         assert -1 < five_begin < five_end <= three_begin < three_end
        
        peak_size = end - start + offset + offset + 1
        begins = [0] * peak_size
        ends = [0] * peak_size
        
        h_5p_begin = -1
        h_5p_end = -1
        h_3p_begin = -1
        h_3p_end = -1
        
        
        # making arrays of start and end positions 
        count_sum = 0
        for sequence in sequences:
            count = float(sequence.data[1].split("-")[1])
            begin_pos = sequence.begin - start + offset
            end_pos = sequence.end - start + offset
            
            if begin_pos < 0 or end_pos >= len(ends):
#                 print "out of range", begin_pos, end_pos
                continue
            
            begins[begin_pos] += count
            ends[end_pos] += count
            count_sum += count
        
#         print begins
#         print ends
#         print
        
        if five_end != -1 and five_begin != -1:
            
            
            # 5p start
            het = 0
            reads_5b = 0
            i = 0
            
            s = begins[:2*offset + 1] # feil...
            
            for i, count in enumerate(s):
                off = abs(offset - i)
                het += off * count
                reads_5b += count
     
            h_5p_begin = begins[offset] / (reads_5b + 1.0)
            
            if begins[offset] != max(s):
                print 
                print candidate.candidate_type, "five start"
                print (five_begin, start), "- -"
                print begins[offset], max(s), candidate.peak_5b, s
                print begins
                print
#             assert begins[offset] == max(s)
    
    
            # 5p end
            het = 0
            reads_5e = 0
            startpos = five_end - start
            
            s = ends[startpos:startpos+2*offset+1]
            
            for i, count in enumerate(s):
                off = abs(offset - i)
                het += off * count
                reads_5e += count
            
            h_5p_end = ends[startpos+offset] / (reads_5e + 1.0)
            

            if ends[startpos+offset] != max(s):
                print
                print candidate.candidate_type, "five end"
                print (five_begin, start), five_end, "- -"
                print ends[startpos+offset], max(s), candidate.peak_5e, s
                print ends
                print
#             assert ends[startpos+offset] == max(s)


        if three_begin != -1 and three_end != -1:
            # three begin
            het = 0
            reads_3b = 0
            startpos = three_begin - start 
            s = begins[startpos:startpos+2*offset+1]
    
            for i, count in enumerate(s):
                off = abs(offset - i)
                het += off * count
                reads_3b += count
    
            h_3p_begin = begins[startpos+offset] / (reads_3b + 1.0)
            
            if begins[startpos+offset] != max(s):
                print
                print candidate.candidate_type, "three start"
                print "- -", (three_begin, start), (three_end, end)
                print begins[startpos+offset], max(s), candidate.peak_3b, s
                print begins
                print    
#             assert begins[startpos+offset] == max(s)
            
            # three end
            het = 0
            reads_3e = 0
            s = ends[-(2*offset + 1):] # feil...
            
            for i, count in enumerate(s):
                off = abs(offset - i)
                het += off * count
                reads_3e += count
            
            h_3p_end = ends[-offset-1] / (reads_3e + 1)
            
            if ends[-offset-1] != max(s):
                print
                print candidate.candidate_type, "three end"
                print "- -", (three_end, end) 
                print ends[-offset-1], max(s), candidate.peak_3e, s
                print ends
                print    
#             assert ends[-offset-1] == max(s)

        
#         print h_5p_begin, h_5p_end, h_3p_begin, h_3p_end
#         print reads_5b, reads_5e, reads_3b, reads_3e
#         print
        
        candidate.set_heterogenity(h_5p_begin, h_5p_end, h_3p_begin, h_3p_end)
        candidate.set_reads(reads_5b, reads_5e, reads_3b, reads_3e)


    assert False






test = [(-1, 5.56952055662), (0, 0.45652146621), (0, 0.45652146621), (0, 0.684782199315), (0, 0.755811053283),
        (0, 1.35776392484), (0, 2.31687990961), (0, 4.92326175656), (0, 9.3051189643), (0, 10.2477547709),
        (0, 13.1643265211), (0, 19.017983184), (0, 153.365633154), (1, 1.54427848927), (1, 4.95403198636),
        (26, 2.5365677152), (27, 1.11810256504), (27, 2.91814436947), (27, 66.3157908263), (28, 0.580339097938),
        (28, 1.39733138991), (29, 0.860401566619), (30, 1.69769221507), (30, 10.7545960107), (31, 0.502594518553),
        (31, 36.3874619966), (31, 41.250641651), (31, 61.7909565118), (31, 1125.03076517), (32, 43.7672911117),
        (32, 63.5136033375), (32, 137.627220276), (32, 1682.1063132), (33, 1.09422817087), (33, 4.78043610531),
        (33, 6.98928242747), (33, 45.9799562024), (33, 164.073065902), (34, 1.11413743851), (34, 5.14641924648),
        (34, 7.12364906516), (34, 48.9696523441), (35, 0.742758292339), (35, 1.85953914853), (35, 29.7243569883),
        (35, 76.3024667924), (36, 2.4060755665), (36, 3.76696964832), (36, 9.54552783057), (36, 38.723242964),
        (37, 1.95216025115), (37, 4.33308466534), (37, 13.0962479581), (38, 6.0919556485), (38, 6.6004246249),
        (38, 8.37959055423), (38, 30.7372608456), (39, 2.65283220981), (39, 3.90543975898), (39, 5.91568379626),
        (39, 25.8738412711), (40, 3.33437597547), (40, 6.05180425974), (40, 27.7477669639), (41, 3.91273801896),
        (41, 19.5219993697), (41, 27.6361305592), (41, 68.2266443687), (42, 3.11740663755), (42, 3.64297834827),
        (42, 20.9962993883), (43, 2.70601039643), (43, 23.0832015591), (43, 60.9988827057), (43, 216.18246544), 
        (44, 0.530624319308), (44, 1.49856934562), (44, 13.6843553261), (44, 53.2569365633), (45, 19.0953427992), 
        (45, 19.3847259274), (45, 26.3594344715), (45, 47.0693102033), (46, 6.06832114118), (46, 28.7297918683), 
        (46, 83.0357990988), (46, 108.091866931), (46, 149.302703609), (47, 1.95893055461), (47, 3.93475700397), 
        (47, 5.43879620575), (48, 5.88220576432), (48, 31.0785074767), (48, 43.1823103515), (48, 193.402280143), 
        (49, 0.580339097938), (49, 55.7390494386), (49, 81.0155826823), (49, 181.797969305), (49, 637.206306722), 
        (50, 13.8648705842), (50, 38.7514767215), (50, 41.8318136308), (50, 196.572111443), (51, 1.11413743851), 
        (51, 3.13827444737), (51, 18.4769515531), (52, 10.9853514031), (52, 50.1902875508), (52, 134.383608404), 
        (52, 313.533401227), (53, 3.17111542775), (53, 16.1525070188), (53, 29.3313544923), (53, 67.8576397159), 
        (54, 3.92621482146), (54, 7.89972188463), (54, 29.3841150397), (54, 119.445706792), (55, 0.742758292339), 
        (55, 1.60681988227), (55, 3.04315900117), (55, 25.5026212067), (56, 28.2100247961), (56, 66.105406128), 
        (57, 1.38083847274), (57, 28.1446549798), (58, 5.75901049227)]


print 1 / 0 if 0 else 123

# print list(itertools.groupby(test, lambda x: x[0]))



 
# print list(map(sumgroups, itertools.groupby(test, lambda x: x[0])))

# testing only, delete pls
# intervals = {1:"10", 2:"11", 3:"12", 11:"13", 4:"14", 1:"15"} # error
# 
# best_start_pos = 1
# 
# 
# ends = filter(lambda (k,v): k == best_start_pos, intervals.items())
# 
# print ends



