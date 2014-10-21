'''
Created on 13. okt. 2014

@author: hakon
'''


def find_overhang(candidates):
    ''' calculate alignment of 5' and 3' candidates '''
    
    #TODO vienna needs 10nt on each side of candidate
    for candidate in candidates:
        
        hairpin = candidate.hairpin
        fold = candidate.hairpin_fold
        
        five_start = 0
        five_end = 20
        three_start = 30
        three_end = len(fold)
        
        
        five_folds = 0
        
        three_folds = 0
        
        level_found = False
        for nr, sign in enumerate(fold):
            
            
            if nr == five_start:
                level_found = True
            
            if sign == "-":
                pass
            elif sign == "(":
                five_folds += 1
            elif sign == ")":
                five_folds -= 1
            