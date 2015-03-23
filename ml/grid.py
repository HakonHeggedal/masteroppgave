'''
Created on 19. feb. 2015

@author: hakon
'''

# used with pool.map for grid search svm kernel parameters
def one_grid(param):
    
    import itertools
    from sklearn import svm
#     print param
    (all_data, all_annotations,test_fold_nr) = param

    
    gamma_values = [4.0 ** i for i in xrange(-4,3)]
    c_values = [4.0 ** i for i in xrange(-3,4)]

    
    test = all_data[test_fold_nr]
    test_annotations = all_annotations[test_fold_nr]
    
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
            score = learner.score(test, test_annotations)
            
            results[k][j] = score
            
            if score > best_score:
                best_c = c_val
                best_gamma = gamma_val
                best_score = score
    
    return (best_score, best_c, best_gamma, results)
    

    
#     
# c_values = [10.0 ** i for i in xrange(-3,4)]
# gamma_values = [10.0 ** i for i in xrange(-3,4)]
# 
# x = 4
# a = range(10)
# print a
# print a[:x]+ a[x+1:]
# 
# 
# derp = [ [0]*len(c_values) for _ in range(len(gamma_values))]
# print derp
# derp[0][0] = 1
# print derp






