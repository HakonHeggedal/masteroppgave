'''
Created on 21. mar. 2015

@author: hakon
'''
import numpy
import copy
from multiprocessing import Pool
# from multiprocessing import current_process
# from multiprocessing import cpu_count
from ml.grid import one_grid
import pickle
from ml import vectorize
from sklearn import preprocessing

print "starting"

annotations = pickle.load( open("save_an.p", "rb"))
print "loaded annotations"



# annotated_data = pickle.load( open("save_da.p", "rb"))

train = pickle.load( open("save_train.p", "rb"))
print "loaded data"

train_annotations = annotations[1:]

# train = annotated_data[1:]
# train = map(vectorize.candidates_to_array, train)
# train = map(preprocessing.scale, train)
# pickle.dump(train, open("save_train.p", "wb"))



print len(train)


folds = range(len(train))

copy_train = [numpy.copy(train) for _ in folds]
copy_train_anno = [copy.deepcopy(train_annotations) for _ in folds]



threads = len(train)
pool = Pool(threads)

zipzap = zip(copy_train, copy_train_anno, folds)
print len(zipzap)
print 


res = map(one_grid, zipzap)

print res

best_scores, best_cs, best_gammas, results = zip(*res)

print best_scores
print best_cs
print best_gammas
print results




