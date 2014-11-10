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
        
        
        begins = [0] * (three_end-five_begin+offset+offset)
        ends = [0] * (three_end-five_begin+offset+offset)
        
        count_sum = 0
        for sequence in sequences:
            count = int(sequence.data[1].split("-")[1])
            begin_pos = sequence.begin - five_begin + offset
            end_pos = sequence.end - five_begin + offset
            
            if begin_pos < 0 or end_pos >= len(ends):
                print "out of range", begin_pos, end_pos
                break
            
            begins[begin_pos] += count
            ends[end_pos] += count
            count_sum += count
            
        print begins
        print ends
        print
        
        het = 0
        reads = 0
        for i, count in enumerate(begins[:2*offset + 1]):
            off = abs(offset - i)
            het += off * count
            reads += count
        
#         for i, count in enumerate(ends[]) 
        
        
        
        


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
            
        

