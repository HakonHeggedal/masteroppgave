'''
Created on 25. mar. 2015

@author: hakon
'''



def hairpin_stats(candidates, mirna, hc_mirna):
    pass



def match_pos(pos, entropy_dict):
    
    if pos in entropy_dict:
        return max(entropy_dict[pos].values())
    
    return -1
    



def find_matching_pos(candidate):
    
    
    entropy_dict = candidate.bitpair_entropy_dict
    precursor = candidate.hairpin_fold_40
    
    area = range(-5,6)
    
    
    match_area = map(match_pos, area)
    
    match_positions = [x-i for i,x in enumerate(match_area)]
    
    best = max(set(match_positions), key=match_positions.count)