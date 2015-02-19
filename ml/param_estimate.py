'''
Created on 19. feb. 2015

@author: hakon
'''



# used by pool.map
def param_estimate_fold(param):
    import numpy
    import itertools
    from sklearn import svm
    
    
    (all_data, all_annotations, test_fold_nr) = param
    
    training = numpy.delete(all_data, test_fold_nr, 0)
    training_annotation = numpy.delete(all_annotations, test_fold_nr, 0)
    
    
#         print len(training_annotation), training_annotation
    training = numpy.vstack(training)
    training_annotation = list(itertools.chain.from_iterable(training_annotation) )
    

    testing = all_data[test_fold_nr]
    testing_annotations = all_annotations[test_fold_nr]
    
#         print len(testing), testing
#         print
#         print len(training), training
#         assert False
    features = len(testing[0])
    scores = [0]*features
    
    learner = svm.SVC(probability=True, cache_size=500)
    learner.fit(training, training_annotation)
    base_score = learner.score(testing, testing_annotations)
    

    for i in xrange(features):
        train_removed = numpy.delete(training, i, 1)
        test_removed = numpy.delete(testing, i, 1)
        learner.fit(train_removed, training_annotation)
        removed_score = learner.score(test_removed, testing_annotations)
        scores[i] = base_score - removed_score
        
    return scores