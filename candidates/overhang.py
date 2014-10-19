'''
Created on 13. okt. 2014

@author: hakon
'''


def find_overhang(candidates):
    ''' calculate alignment of 5' and 3' candidates '''
    for candidate in candidates:
        
        hairpin = candidate.hairpin
        fold = candidate.hairpin_fold
        
        five_start = 0
        
        three_end = len(fold)