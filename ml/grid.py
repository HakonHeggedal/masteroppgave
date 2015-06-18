'''
Created on 19. feb. 2015

@author: hakon
'''






def filter_miRNAs(data, annotations):
    ''' return only miRNAs and their annotations (1s) '''
    data = [c for c,a in zip(data, annotations) if a ]
    annotations = [1] * len(data)
    return data, annotations


def filter_candidates(data, annotations):
    ''' returns only canididates, and their annotations (0s) '''
    not_candidate = annotations.index(1)
    
    data_r = data[:not_candidate]
    annotations_r = annotations[:not_candidate]
    return data_r, annotations_r



def filter_dead(data, annotations):
    ''' returns only dead miRNAs, and their annotations (0s)'''
    not_candidate = annotations.index(1)
    not_miRNA = annotations.index(0, not_candidate)
    
    data_r = data[not_miRNA:]
    annotations_r = annotations[not_miRNA:]
    return data_r, annotations_r










# used with pool.map for grid search svm kernel parameters
def one_grid(param):
    
    import itertools
    from sklearn import svm
#     print param
    (all_data, all_annotations, test_fold_nr, filter_test) = param

    
#     gamma_values = [4.0 ** i for i in xrange(-4,3)] # too large differences
#     c_values = [4.0 ** i for i in xrange(-3,4)]
    
    gamma_values = [2.0 ** i for i in xrange(-18,-1)]
    c_values = [2.0 ** i for i in xrange(-1,10)]

    
    test = all_data[test_fold_nr]
    test_annotations = all_annotations[test_fold_nr]
    
#     print test
    
    if filter_test:
        test, test_annotations = filter_test(test, test_annotations)
    
#     print len(test_annotations), test_annotations
    
    test_hc, test_hc_an = filter_miRNAs(test, test_annotations)
    test_cand, test_cand_an = filter_candidates(test, test_annotations)
    test_dead, test_dead_an = filter_dead(test, test_annotations)
    
#     print len(test_hc)
#     print len(test_cand)
#     print len(test_dead)
#     print
#     
#     print "fold", test_fold_nr
#     print len(c_values), c_values, gamma_values
#     print test[0][0]
    
    p1 = list(all_data[:test_fold_nr])
    p2 = list(all_data[test_fold_nr+1:])
    train_lists = p1 + p2
    
    train_annotations_lists = all_annotations[:test_fold_nr] + all_annotations[test_fold_nr+1:]
    
    
    train = list(itertools.chain.from_iterable(train_lists))
    train_annotations = list(itertools.chain.from_iterable(train_annotations_lists))
    
    
    best_score = 0.0
    best_c = -12345
    best_gamma = -12345

    
    results = [ [ 0 for _ in xrange(len(gamma_values))] for _ in xrange(len(c_values))]
    for j, gamma_val in enumerate(gamma_values):
        for k, c_val in enumerate(c_values):
            
            
            learner = None
            learner = svm.SVC(C=c_val,gamma=gamma_val, cache_size=500)
            learner.fit(train, train_annotations)
#             score = learner.score(test, test_annotations)
            score_hc = learner.score(test_hc, test_hc_an)
            score_cand = learner.score(test_cand, test_cand_an)
            score_dead = learner.score(test_dead, test_dead_an)
            
#             score = score_hc + score_cand + score_dead
            score = score_hc * score_cand * score_dead
#             print (score_hc, score_cand, score_dead), "\t", ("gamma:",gamma_val,"C:", c_val)
            
            
            results[k][j] = score
            
            if score > best_score:
                best_c = c_val
                best_gamma = gamma_val
                best_score = score
#     print (best_score, best_c, best_gamma)
#     print
    print best_score, best_c, best_gamma
    
    assert best_c != -12345, "c has no value, default -12345"
    
    return (best_score, best_c, best_gamma, results)
    




