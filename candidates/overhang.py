'''
Created on 13. okt. 2014

@author: hakon
'''


def find_overhang(candidates):
    
    for candidate in candidates:
        
        fold = candidate.hairpin_fold
        
        five_start = 0
        
        three_end = len(fold)