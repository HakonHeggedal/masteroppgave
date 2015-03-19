'''
Created on 18. mar. 2015

@author: hakon
'''

import itertools



complimented = {"A":"T", "T":"A", "G":"C", "C":"G"}
def reverse_compliment(sequence):
    ''' returns the reverse compliment of a DNA-sequence '''
    global complimented
    return "".join(complimented[x] for x in sequence[::-1])


def filter_intervals(intervals, offset, startpos, endpos):
    if endpos <= 0: return []
    intervals = [i for i in intervals if i.begin-offset >= startpos and i.end-offset <= endpos]
    return intervals


def _sumgroups(group):
    return (group[0], sum([x[1] for x in group[1]]) )
    

def best_interval(intervals, offset, max_end=100): #TODO: max limit to ends, i.e has 0-40 (but only 0-20 allowed)
    
    if not intervals:
        return -1, -1, -1, -1
    
#     start_pos = -1
#     start_max_val = -1
    end_pos = -1
    end_max_val = -1    
    
    # tuples of (startpos, val)
    start_positions = sorted([(s.begin-offset, float(s.data[1].split("-")[1]), s.data[0] ) for s in intervals])
    # grouped tuples for every startpos
    grouped_starts = list(map(_sumgroups, itertools.groupby(start_positions, lambda x: x[0])))

        
#     for pos, val in grouped_starts:
#         if val > start_max_val:
#             start_max_val = val
#             start_pos = pos
            
    start_pos, start_max_val = max(grouped_starts, key=lambda x: x[1])
    
#     mpos,maxv = max(grouped_starts, key=lambda x: x[1])
#     if maxv != start_max_val:
#         print start_max_val, maxv
#     assert start_max_val == maxv
    
    MAX_MATURE_SEQ = 44

    # only positions overlapping sequences from start peak are viable
    end_range_pos = [s.end-offset for s in intervals if abs(s.begin-offset - start_pos) < 2]
    end_test = [s.end-offset for s in intervals if s.begin-offset == start_pos]
    end_range = range(min(end_range_pos)-3, max(end_range_pos)+4 )
    end_range = filter(lambda x: x-start_pos <=  MAX_MATURE_SEQ, end_range)
    end_range = set(end_range)
    
    # make tuples (endpos, sum_value)
    end_positions = sorted([(s.end-offset, float(s.data[1].split("-")[1])) for s in intervals if s.end-offset in end_range])
    grouped_ends = list(map(_sumgroups, itertools.groupby(end_positions, lambda x: x[0])))
    
    
#     _e_p = sorted([(s.end-offset, float(s.data[1].split("-")[1])) for s in intervals])
#     _a_e = list(map(_sumgroups, itertools.groupby(_e_p, lambda x: x[0])))

    if not grouped_ends:
        return -1, -1, -1, -1
    
    for pos, val in grouped_ends:
        if val > end_max_val:
            end_max_val = val
            end_pos = pos
    
    
    
    
    s, sv, e, ev = _old_best_interval(intervals, offset)
    
    derp = [(i.begin-offset, i.end-offset) for i in intervals]
    
#     if abs(start_max_val - 0.45652146621) < 0.01 and abs(end_max_val - 0.45652146621) < 0.01:
    
    if abs( end_max_val - ev) > 0.1:
        print
        print intervals[0].data
        print "all"
        print derp
        print "starts:"
    #     print start_positions
        print grouped_starts
        print "ends:"
        print grouped_ends
#         print "all ends:"
#         print _a_e
        print (start_pos, start_max_val), (end_pos, end_max_val)
        print (s, sv), (e, ev)
        print end_range, end_range_pos, end_test
    
    assert start_pos == s
#     assert start_max_val - 0.001 < sv
#     assert sv < start_max_val + 0.001
    
#     assert end_pos == e
#     assert end_max_val - 0.001 < ev 
#     assert ev < end_max_val + 0.001
    
    
#     return s, sv, e, ev
    return start_pos, start_max_val, end_pos, end_max_val
    

    
    

def _old_best_interval(intervals, offset):
    '''
        finding best interval region given a list of intervals
         returns the interval region with highest total frequency in start and end positions
        the end position must be close to the end of any interval starting in the best start position 
     '''
    
    if len(intervals) == 0:
        return -1, -1, -1, -1
    
    start_positions = {} # startpos positions -> frequency
    end_positions = {} # end pos --> frequency
    
    best_start = -1
    best_start_pos = -1
    best_end = -1
    best_end_pos = -1


    # find best start position
    for interval in intervals:
        
        startpos = interval.begin - offset
        _name = interval.data[1]
        freq = float(_name.split("-")[1])
        
        if startpos in start_positions:
            start_positions[startpos] += freq
        else:
            start_positions[startpos] = freq
            
        if start_positions[startpos] > best_start:
            best_start = start_positions[startpos]
            best_start_pos = startpos
    
    assert best_start == max(start_positions.values()) # test that best start is found
    
    # possible end positions set to the range of +- 3 from ends of intervals starting in best start pos
    start_end_pos = [ i.end - offset for i in intervals if i.begin - offset == best_start_pos ]
    legal_end_pos = set( range(min(start_end_pos)-3, max(start_end_pos)+3) )     
    
    
    # find best end position
    for interval in intervals:
        
        endpos = interval.end - offset
        if endpos not in legal_end_pos:
            continue
        
        _name = interval.data[1]
        freq = float(_name.split("-")[1])
        
        if endpos not in end_positions:
            end_positions[endpos] = freq
        else:
            end_positions[endpos] += freq
            
        if end_positions[endpos] > best_end:
            best_end = end_positions[endpos]
            best_end_pos = endpos
            

    best_end2 = -1
    best_end_pos2 = -1
    for endpos, freq in end_positions.iteritems():
        if freq > best_end:
            best_end2 = freq
            best_end_pos2 = endpos
    
    if max(end_positions.values()) > best_end:

        print max(end_positions.values()),  best_end, best_start_pos, best_end2, best_end_pos2
        print sorted(end_positions.items())
        print sorted(legal_end_pos)
        print 
        assert max(end_positions.values()) == best_end
    
    return best_start_pos, best_start, best_end_pos, best_end


# @DeprecationWarning
# def _old_best_interval(intervals, offset):
#     # random end position gives error
#     
#     if len(intervals) == 0:
#         return -1, -1, -1, -1
#     start_positions = {} # startpos positions -> frequency
#     best_start = -1
#     best_start_pos = -1
#     end_pos_midpoint = -1
#     best_end = -1
#     best_end_pos = -1
#     global derpy_errors
# 
# 
#     for interval in intervals:
#         
#         startpos = interval.begin - offset
#         endpos = interval.end - offset
#         
#         name = interval.data[1]
#         freq = float(name.split("-")[1])
#         
#         if startpos in start_positions:
#             start_positions[startpos] += freq
#         else:
#             start_positions[startpos] = freq
#             
#         if start_positions[startpos] > best_start:
#             best_start = start_positions[startpos]
#             best_start_pos = startpos
#             end_pos_midpoint = endpos
# #             find 3p endpos if possible
# 
#     assert best_start == max(start_positions.values())
#     
#     ends = {}
#     for interval in intervals:
#         
#         startpos = interval.begin - offset
#         endpos = interval.end - offset
#         
#         if endpos+5 < end_pos_midpoint or endpos-5 >= end_pos_midpoint:
#             continue # out of bounds
#         
#         name = interval.data[1]
#         freq = float(name.split("-")[1])
#         
#         if endpos not in ends:
#             ends[endpos] = freq
#         else:
#             ends[endpos] += freq
#             
#     for endpos, freq in ends.iteritems():
#         if freq > best_end:
#             best_end = freq
#             best_end_pos = endpos     
#         
#     return best_start_pos, best_start, best_end_pos, best_end

