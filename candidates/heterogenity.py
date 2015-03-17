'''
Created on 10. okt. 2014

@author: hakon
'''

def heterogenity(candidates, offset=5):
    
    print "\n\tcalculating variation in 5p and 3p start positions"
    
    for candidate in candidates:
        sequences = candidate.mapped_sequences
        
        five_begin = candidate.pos_5p_begin
        five_end = candidate.pos_5p_end
        three_begin = candidate.pos_3p_begin
        three_end = candidate.pos_3p_end
        start = candidate.hairpin_start
        end = candidate.hairpin_end
        
        print
        
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

