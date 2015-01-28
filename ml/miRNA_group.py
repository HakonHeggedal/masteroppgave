'''
Created on 28. jan. 2015

@author: hakon
'''
import random

def split_candidates(candidates, mirna_dict, mirna_groups, folds=5, shuffle=False):
    
    training_lists = [[] for _ in xrange(folds)]
    test_lists = [[] for _ in xrange(folds)]
    
    if shuffle:
        random.shuffle(candidates)
    
    mis = []
    c_count = 0
    
    # split candidate list into only candidate and microRNAs
    for c in candidates:
        hashval = c.chromosome + str(c.pos_5p_begin)
        if hashval in mirna_dict:
            mis.append(c)
        else:
            test_lists[c_count].append(c)
            c_count += 1
    
    
    ungrouped_mirna = []
    group_to_mis = {}
    max_mi_size = 0
    ungroup_fill_size = 0
    
    for m in mis:
        hashval = m.chromosome + str(m.pos_5p_begin)
        mi_name = mirna_dict[hashval]
        
        if mi_name in mirna_groups:
            group = mirna_groups[mi_name]
            
            if group in group_to_mis:
                group_to_mis[group].append(m)
            else:
                group_to_mis[group] = [m]
                
    
    mi_set_nr = 0
    for group, mi_list in group_to_mis.iteritems():
        size = 0
        while size < max_mi_size:
            
            pass
            #TODO: legg til miRNA fra grupper i folds
            
            
        else:
            mi_set_nr += 1
        
        if mi_set_nr > folds:
            break
    
    # TODO: legg in ugrupperte i de foldene med minst
    for m in ungrouped_mirna:
        pass
            
            
    
    # order mirnas by group -> mirnas 
#     for key,miRNA_name in mis.iteritems():
        
        
        
        
        pass
            
    