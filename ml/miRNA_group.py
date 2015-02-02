'''
Created on 28. jan. 2015

@author: hakon
'''
import random

def split_candidates(candidates, mirna_dict, mirna_groups, folds=5, shuffle=False):
    
    print
    print "splitting data into", folds, "folds"
    
    training_lists = [[] for _ in xrange(folds)]
    test_lists = [[] for _ in xrange(folds)]
    
    if shuffle:
        random.shuffle(candidates)
    
    mis = []
    c_count = 0
    
    # split candidate list into only candidate and microRNAs
    for c in candidates:
        hashval = c.chromosome + c.chromosome_direction + str(c.pos_5p_begin)
        if hashval in mirna_dict:
            mis.append(c)
        else:
            test_lists[c_count].append(c)
            c_count = (c_count + 1) % folds
    
    
    ungrouped_mirna = []
    group_to_mis = {}
    
    # make miRNA groups: group -> [mirnas]
    for m in mis:
        hashval = m.chromosome + m.chromosome_direction + str(m.pos_5p_begin)
        mi_name = mirna_dict[hashval]
        
        if mi_name in mirna_groups:
            group = mirna_groups[mi_name]
            
            if group in group_to_mis:
                group_to_mis[group].append(m)
            else:
                group_to_mis[group] = [m]
        else:
            ungrouped_mirna.append(m)
                
    
    avg_size = len(mis) / folds
    nr = 0
    # add groups of miRNAS to different sets
    for group, mi_list in group_to_mis.iteritems():
        
        while len(training_lists[nr]) > avg_size:
            nr = (nr + 1) % folds
            
        training_lists[nr] += mi_list
        nr = (nr + 1) % folds


    
    nr = 0
    avg_size += 1
    for m in ungrouped_mirna:
        
        while len(training_lists[nr]) > avg_size:
            nr = (nr + 1) % folds
            
        training_lists[nr].append(m)
        nr = (nr + 1) % folds
                
    print "\nstats:"
    print "folds:", folds
    print "candidates:", len(candidates)
    print "mirna:", len(mis), len(mirna_dict)
    print
    print "mirna, sizes training:"
    for x in training_lists:
        print len(x),
    print
    print "candidates, sizes test:"
    for x in test_lists:
        print len(x),
    print
    print "ungrouped:", len(ungrouped_mirna)
    
    
    return training_lists, test_lists




    