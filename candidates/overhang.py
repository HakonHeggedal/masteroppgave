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
        for nr, sign in enumerate(fold[::-1]):
            if sign == "(":
                three_folds -= 1
            elif sign == ")":
                three_folds += 1
                
            if nr >= three_end_pos:
                if sign == ")":
                    three_distance = nr - three_end_pos
                    break        
        
        print five_folds, five_distance
        print three_folds, three_distance
            
            
            
    
    
    
f1 = ".(((((((((..((((((((.((.................)).))))))))..)))))))))..."