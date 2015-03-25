'''
Created on 11. feb. 2015

@author: hakon
'''



def create_folds(candidates, mirna_dict, dead_dict, high_conf, mirna_groups, folds):
    
    print
    print "splitting data into", folds, "folds"
    
    hc_mirna_lists = [[] for _ in xrange(folds)]
    other_mirna_lists = [[] for _ in xrange(folds)]
    candidate_lists = [[] for _ in xrange(folds)]
    dead_lists = [[] for _ in xrange(folds)]
    
    mi_hq = []

    c_count = 0
    d_count = 0
    
    # split candidate list into only candidate and microRNAs
    for c in candidates:
        hashval = c.chromosome + c.chromosome_direction + str(c.hairpin_start)
        
        if hashval in dead_dict:
            dead_lists[d_count].append(c)
            d_count = (d_count + 1) % folds
        
        elif hashval in mirna_dict:
            mi = mirna_dict[hashval]
            if mi in high_conf:
                mi_hq.append(c)
            else:
                other_mirna_lists[c_count].append(c)
        else:
            candidate_lists[c_count].append(c)
            c_count = (c_count + 1) % folds
    
    
    ungrouped_mirna = []
    group_to_mis = {}
    
    # make miRNA groups: group -> [mirnas]
    for m in mi_hq:
        hashval = m.chromosome + m.chromosome_direction + str(m.hairpin_start)
        mi_name = mirna_dict[hashval]
        
        if mi_name in mirna_groups:
            group = mirna_groups[mi_name]
            
            if group in group_to_mis:
                group_to_mis[group].append(m)
            else:
                group_to_mis[group] = [m]
        else:
            ungrouped_mirna.append(m)
                
    
    avg_size = len(mi_hq) / folds
    nr = 0
    # add groups of miRNAS to different sets
    for group, mi_list in group_to_mis.iteritems():
        
        while len(hc_mirna_lists[nr]) > avg_size:
            nr = (nr + 1) % folds
            
        hc_mirna_lists[nr] += mi_list
        nr = (nr + 1) % folds


    print "not grouped:", len(ungrouped_mirna)
    nr = 0
    avg_size += 1
    # add mirna without groups
    for m in ungrouped_mirna:
        
        while len(hc_mirna_lists[nr]) > avg_size:
            nr = (nr + 1) % folds
            
        hc_mirna_lists[nr].append(m)
        nr = (nr + 1) % folds

    
    training_data = []
    annotations = []
    
    # combine hc-candidates + candidates + dead, and create annotations
    for hc_list, c_list, d_list in zip(hc_mirna_lists, candidate_lists, dead_lists):
        
        training_list = c_list + hc_list + d_list
        training_data.append(training_list)
        
        annotation = [0]*len(c_list) + [1]*len(hc_list) + [0]*len(d_list)
        annotations.append(annotation)
#         print len(training_list), len(annotation), len(hc_list), sum(annotation)
    
    
    print sum( map(len, dead_lists)  ), "dead"
    print sum( map(len, candidate_lists)  ), "candidates"
    print sum( map(len, hc_mirna_lists)  ), "high confidence"
    print sum( map(len, other_mirna_lists)  ), " low confidence"
    
#     assert False
    return training_data, annotations, other_mirna_lists




# lili = [ [1,2,3], [1,2,4,2], [1,2,5] ]
# 
# print sum( map(len, lili))




    
#     print "\nstats:"
#     print "folds:", folds
#     print "candidates:", len(candidates)
#     print "mirna:", len(mi_hq), len(mirna_dict)
#     print
#     print "mirna, sizes training:"
#     for x in hc_mirna_lists:
#         print len(x),
#     print
#     print "candidates, sizes test:"
#     for x in candidate_lists:
#         print len(x),
#     print
#     print "ungrouped:", len(ungrouped_mirna)

