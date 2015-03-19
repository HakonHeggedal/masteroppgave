'''
Created on 19. feb. 2015

@author: hakon
'''





# used by pool.map
def one_grid(param):
    
    import itertools
    from sklearn import svm
    
    (all_data, all_annotations,test_fold_nr) = param
    
    c_values = [10.0 ** i for i in xrange(-3,4)]
    gamma_values = [10.0 ** i for i in xrange(-3,4)]
    
    
    test = all_data[test_fold_nr]
    test_annotations = all_annotations[test_fold_nr]
    
    train_lists = all_data[:test_fold_nr] + all_data[test_fold_nr+1:]
    train_annotations_lists = all_annotations[:test_fold_nr] + all_annotations[test_fold_nr+1:]
    
    train = list(itertools.chain.from_iterable(train_lists))
    train_annotations = list(itertools.chain.from_iterable(train_annotations_lists))
    
    
    best_score = 0.0
    best_c = -1
    best_gamma = -1
    
    learner = None
    
    results = [ [0]*len(c_values) for _ in range(len(gamma_values))]
    
    for j, gamma_val in enumerate(gamma_values):
        for k, c_val in enumerate(c_values):
            learner = svm.SVC(C=c_val,gamma=gamma_val, cache_size=500)
            learner.fit(train, train_annotations)
            score = learner.score(test, test_annotations)
            
            results[j][k] = score
            
            if score > best_score:
                best_c = c_val
                best_gamma = gamma_val
                
                
                
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






