'''
Created on 10. okt. 2014

@author: hakon
'''












def heterogenity(candidates, offset=5):
    
    for candidate in candidates:
        sequences = candidate.mapped_sequences
        
        five_begin = candidate.pos_5p_begin
        five_end = candidate.pos_5p_end
        three_begin = candidate.pos_3p_begin
        three_end = candidate.pos_3p_end
        
        assert -1 < five_begin < five_end < three_begin < three_end
        
        
        begins = [0] * (three_end-five_begin+offset+offset+1)
        ends = [0] * (three_end-five_begin+offset+offset+1)
        
        count_sum = 0
        for sequence in sequences:
            count = float(sequence.data[1].split("-")[1])
            begin_pos = sequence.begin - five_begin + offset
            end_pos = sequence.end - five_begin + offset
            
            if begin_pos < 0 or end_pos >= len(ends):
                print "out of range", begin_pos, end_pos
#                 assert False
                break
            
            begins[begin_pos] += count
            ends[end_pos] += count
            count_sum += count
            
        print begins
        print ends

        print "results:"
        het = 0
        reads_5b = 0
#         0 to 10
        i = 0
        
        s = begins[:2*offset + 1]
        print len(s), s, # five start
        for i, count in enumerate(s):
            off = abs(offset - i)
            het += off * count
            reads_5b += count
#         print reads_5b, begins[offset], begins[offset] / reads_5b, i
        
        h_5p_begin = begins[offset] / (reads_5b + 1.0)

        het = 0
        reads_5e = 0
        startpos = five_end - five_begin # five end
        s = ends[startpos:startpos+2*offset+1]
        print len(s), s,
        for i, count in enumerate(s):
            off = abs(offset - i)
            het += off * count
            reads_5e += count
#         print reads_5e, ends[startpos+offset], ends[startpos+offset] / reads_5e, i
        
        h_5p_end = ends[startpos+offset] / (reads_5e + 1.0)

        het = 0
        reads_3b = 0
        startpos = three_begin - five_begin # three begin 
        s = begins[startpos:startpos+2*offset+1]
        print len(s), s,
        for i, count in enumerate(s):
            off = abs(offset - i)
            het += off * count
            reads_3b += count
#         print reads_3b, begins[startpos+offset], begins[startpos+offset] / reads_3b, i
        h_3p_begin = begins[startpos+offset] / (reads_3b + 1.0)

        het = 0
        reads_3e = 0
        s = ends[-(2*offset + 1):]
        print len(s), s,
        for i, count in enumerate(s):
            off = abs(offset - i)
            het += off * count
            reads_3e += count
#         print reads_3e, ends[-offset-1], ends[-offset-1] / reads_3e, i
        
        h_3p_end = ends[-offset-1] / (reads_3e + 1)
#         for i, count in enumerate(ends[]) 
        
        print h_5p_begin, h_5p_end, h_3p_begin, h_3p_end
        print reads_5b, reads_5e, reads_3b, reads_3e
        print
        
        candidate.set_heterogenity(h_5p_begin, h_5p_end, h_3p_begin, h_3p_end)
        candidate.set_reads(reads_5b, reads_5e, reads_3b, reads_3e)


def frequency_counting(candidates, freq_range=5):
    ''' 
    compute quality of 5' and 3' positions
    '''
    
    
    for candidate in candidates:
        
        sequences = candidate.mapped_sequences
        count_all = 0
        misses = 0
        #TODO: lists to store position hits
        
        five_begin = candidate.pos_5_begin
        five_end = candidate.pos_5_end
        three_begin = candidate.pos_3_begin
        three_end = candidate.pos_3_end
        
        begin_5s = [0]*((2*freq_range) + 1)
        end_5s = [0]*((2*freq_range) + 1)
        begin_3s = [0]*((2*freq_range) + 1)
        end_3s = [0]*((2*freq_range) + 1)
        
        for seq in sequences:
            #TODO: find positon for start / end
            # compare to 5' and 3'
            # add to position lists
            count = int(seq.data[1].split("-")[1])
            count_all += count
            
            if (five_begin-freq_range <= seq.begin <= five_begin+freq_range and
                five_end -freq_range <= seq.end <= five_end+freq_range):
#                 print five_begin-freq_range, seq.begin, five_begin-freq_range <= seq.begin
#                 print seq.begin, five_begin+freq_range, seq.begin >= five_begin+freq_range
                pos = freq_range + seq.begin - five_begin
                begin_5s[pos] += count
                
                pos = freq_range + seq.end - five_end
                end_5s[pos] += count

            
            elif (three_begin-freq_range <= seq.begin and
                  seq.begin <= three_begin+freq_range and
                  three_end-freq_range <= seq.end and
                  seq.end <= three_end+freq_range):
                
                pos = freq_range + seq.begin - three_begin
                begin_3s[pos] += count
                
                pos = freq_range + seq.end - three_end
                end_3s[pos] += count

            else:
                misses += count
                
                
            
        print begin_5s
        print end_5s
        print begin_3s
        print end_3s
        print misses
        print count_all
        print misses * 1.0 / count_all
        print
            
        

