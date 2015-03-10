'''
Created on 13. okt. 2014

@author: hakon
'''


def get_alignment(candidates):
    
    
    for i, candidate in enumerate(candidates):
        fold_10 = candidate.hairpin_fold_10
#         fold_40 = candidate.hairpin_fold_40
        
#         print "overhang alignment:",
#         print "\t", fold_10
#         print "\t", fold_40
#         
#         # overhang: level and overhang both padded and not
#         print
#         print len(fold_10),
#         print candidate.pos_5p_end - candidate.pos_5p_begin,
#         print candidate.pos_3p_end - candidate.pos_3p_begin,
#         print candidate.pos_3p_end - candidate.pos_5p_begin + 20,
#         print i, candidate.pos_5p_begin, candidate.pos_5p_end,  candidate.pos_3p_begin,  candidate.pos_3p_end
        
#         assert candidate.pos_5p_begin < candidate.pos_5p_end <= candidate.pos_3p_begin < candidate.pos_3p_end
        
#         if len(fold_10) != 20 + candidate.pos_3p_end - candidate.pos_5p_begin:
#             print candidate.hairpin
#             print candidate.hairpin_padded_40
#             print len(candidate.hairpin_padded_40)
#             print len(candidate.hairpin)
#         assert len(fold_10) == 20 + candidate.pos_3p_end - candidate.pos_5p_begin
        
        folds_10 = max_fold(fold_10)
#         folds_40 = max_fold(fold_40)
        
        start_5 = 10 # start point for hairpin with 10 nt padding
        end_5 = candidate.pos_5p_end - candidate.pos_5p_begin + 10
        start_3 = candidate.pos_3p_begin - candidate.pos_5p_begin + 10
        end_3 = candidate.pos_3p_end - candidate.pos_5p_begin + 10
        
        folds_out, off_out = find_overhang(fold_10, start_5, end_3)
        folds_in, off_in = find_overhang(fold_10, end_5, start_3)
        
        start_5 = 40 # 40 nt padding
        end_5 = candidate.pos_5p_end - candidate.pos_5p_begin + 40
        start_3 = candidate.pos_3p_begin - candidate.pos_5p_begin + 40
        end_3 = candidate.pos_3p_end - candidate.pos_5p_begin + 40
        
#         p_folds_out, p_off_out = find_overhang(fold_40, start_5, end_3)
#         p_folds_in, p_off_in = find_overhang(fold_40, end_5, start_3)
        
#         print "\tmaks fold:", folds_10, folds_40
#         print "\t10:", folds_out, off_out, folds_in, off_in
#         print "\t10:", p_folds_out, p_off_out, p_folds_in, p_off_in
        
        candidate.set_alignment_10(folds_10, folds_out, off_out, folds_in, off_in)
#         candidate.set_alignment_40(folds_40, p_folds_out, p_off_out, p_folds_in, p_off_in)



#     assert False

def max_fold(fold):
    
    max_fold = 0
    level = 0
    for sign in fold:
        if sign == "(":
            level += 1
            max_fold = max(max_fold, level)
        elif sign == ")":
            level -= 1
            
    return max_fold

def find_overhang(fold, five_start, three_end):
    ''' calculate alignment of 5' and 3' candidates '''
    
    #TODO vienna needs 10nt on each side of candidate
#     for candidate in candidates:
    
    five_folds = 0
    five_distance = -1
    
    three_folds = 0


    for nr, sign in enumerate(fold):
        
        if sign == "(":
            five_folds += 1
        elif sign == ")":
            five_folds -= 1

        if nr >= five_start:
            if sign == "(":
                five_distance = nr - five_start
                break
  
  
    three_end_pos = len(fold) - three_end
    align_pos = -1

    for nr, sign in enumerate(fold[::-1]):
        if sign == "(":
            three_folds -= 1
        elif sign == ")":
            three_folds += 1
            if three_folds == five_folds:
                align_pos = nr
                break

    
#     print five_folds, five_distance
#     print three_folds, "offset:", three_end_pos - align_pos
    
    return three_end_pos - align_pos, five_folds
    
    
       
    
#     
# f1 = ".(((((((((..((((((((.((.................)).))))))))..)))))))))..."
# start = 22
# end = len(f1)-10
#  
# print max_fold(f1)
# print find_overhang(f1, start, end)






