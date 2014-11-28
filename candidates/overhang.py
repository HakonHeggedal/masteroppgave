'''
Created on 13. okt. 2014

@author: hakon
'''


def get_overhangs(candidates):
    
    for candidate in candidates:
        fold = candidate.hairpin_fold
        fold_padded = None
        
        # overhang: level and overhang both padded and not
        
        hairpin_fold = max_fold(fold)
        padded_fold = max_fold(fold_padded)


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

def find_overhang(fold, startpos, endpos):
    ''' calculate alignment of 5' and 3' candidates '''
    
    #TODO vienna needs 10nt on each side of candidate
#     for candidate in candidates:
        
#     hairpin = candidate.hairpin
#     fold = candidate.hairpin_fold
    
    five_start = 4
    five_end = 20
    three_start = 30
    three_end = len(fold)-6
    
    
    five_folds = 0
    five_distance = -1
    
    three_folds = 0
    three_distance = -1

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

    
    print five_folds, five_distance
    print three_folds, "offset:", three_end_pos - align_pos
    
       
    
    
f1 = ".(((((((((..((((((((.((.................)).))))))))..)))))))))..."

find_overhang(f1)