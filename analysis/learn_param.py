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
from matplotlib import pyplot as plot
import time
import math
from sklearn import svm
import itertools






def filter_miRNAs(data, annotations):
    ''' return only miRNAs and their annotations (1s) '''
    data = [c for c,a in zip(data, annotations) if a ]
    annotations = [1] * len(data)
    return data, annotations

def filter_candidates(data, annotations):
    ''' returns only canididates, and their annotations (0s) '''
    not_candidate = annotations.index(1)
    
    data = data[:not_candidate]
    annotations = annotations[:not_candidate]
    return data, annotations

def filter_dead(data, annotations):
    ''' returns only dead miRNAs, and their annotations (0s)'''
    not_candidate = annotations.index(1)
    not_miRNA = annotations.index(0, not_candidate)
    
    data = data[not_miRNA:]
    annotations = annotations[not_miRNA:]
    return data, annotations

# from ml import vectorize
# from sklearn import preprocessing
# import math
start_time = time.clock()


print "starting"

annotations = pickle.load( open("save_an.p", "rb"))
print "loaded annotations ", time.clock() - start_time

data = pickle.load( open("save_scaled_data.p", "rb"))

train_annotations = annotations[1:]
train = data[1:]

print "loaded data", len(train), len(train_annotations), time.clock() - start_time


folds = range(len(train))

copy_train = [numpy.copy(train) for _ in folds]
copy_train_anno = [copy.deepcopy(train_annotations) for _ in folds]


filters = [filter_miRNAs for _ in folds]
# filters = [filter_candidates for _ in folds]
# filters = [filter_dead for _ in folds]
# filters = [None for _ in folds]

threads = len(train)
pool = Pool(threads)

zipzap = zip(copy_train, copy_train_anno, folds, filters)

"ready to run the folds", time.clock() - start_time

res = pool.map(one_grid, zipzap)
# res = map(one_grid, zipzap)



best_scores, best_cs, best_gammas, results = zip(*res)

def log_avg(l): 
    '''   average of log values, then scaled back '''
    log_list = (math.log(x, 2.0) for x in l)
    log_val = sum(log_list) / len(l)
    avg_val = 2.0 ** log_val 
    return avg_val

def avg(l):
    return sum(l) / len(l)


print 
print "Results", time.clock() - start_time
print "score?\t", avg(best_scores), log_avg(best_scores), best_scores
print "c-val\t", avg(best_cs), log_avg(best_cs), best_cs
print "gamma\t", avg(best_gammas), log_avg(best_gammas), best_gammas

best_c = log_avg(best_cs)
best_gamma = log_avg(best_gammas)

np_res = numpy.array(results)
array_res = numpy.mean( np_res , axis=0 )


# gamma_values = [10.0 ** i for i in xrange(-4,5)]
# c_values = [10.0 ** i for i in xrange(-4,5)]
# gamma_values = [10.0 ** i for i in xrange(-4,5)]
# c_values = [10.0 ** i for i in xrange(-4,5)]


gamma_values = [4.0 ** i for i in xrange(-4,3)]
c_values = [4.0 ** i for i in xrange(-3,4)]


gamma_name = "_gamma_" + str(gamma_values[0]) + "-" + str(gamma_values[-1])
c_name = "C_" + str(c_values[0]) + "-" + str(c_values[-1])
f_name = "" if not filters[0] else filters[0].__name__
plot_name = "plot/Grid_avg_" + f_name + c_name + gamma_name + ".png"


row_labels = gamma_values
column_labels = c_values


#  create 
fig, ax = plot.subplots()
heatmap = ax.pcolor(array_res)
 
# put the major ticks at the middle of each cell
ax.set_xticks(numpy.arange(array_res.shape[0])+0.5, minor=False)
ax.set_yticks(numpy.arange(array_res.shape[1])+0.5, minor=False)
 
ax.set_xticklabels(row_labels, minor=False,  rotation="vertical")
ax.set_yticklabels(column_labels, minor=False)

plot.xlabel("gamma")
plot.ylabel("C")
plot.show()
plot.savefig(plot_name)

print "saved plot", time.clock() - start_time


train = list(itertools.chain.from_iterable(train))
train_annotations = list(itertools.chain.from_iterable(train_annotations))



#===============================================================================
#  first candidates, then HC mirnas, then dead miRNA: 
#  annotations is like [000000111111111000]
#===============================================================================
print "\n final testing"
test_all = data[0]
test_all_annotations = annotations[0]

test_miRNA = [c for c,a in zip(test_all, test_all_annotations) if a ]
test_miRNA_annotations = [1] * len(test_miRNA)


first_miRNA = test_all_annotations.index(1)
test_candidates = test_all[:first_miRNA]
test_candidates_annotations = test_all_annotations[:first_miRNA]

first_dead = test_all_annotations.index(0,first_miRNA)
test_dead = test_all[first_dead:]
test_dead_annotations = test_all_annotations[first_dead:]



print first_miRNA, first_dead, len(test_all_annotations)
print len(test_candidates), len(test_miRNA), len(test_dead)
print len(test_candidates) + len(test_miRNA) + len(test_dead)

learner = svm.SVC(C=best_c,gamma=best_gamma, cache_size=500, probability=True)
learner.fit(train, train_annotations)

score_all = learner.score(test_all, test_all_annotations)
score_miRNA = learner.score(test_miRNA, test_miRNA_annotations)
score_candidates = learner.score(test_candidates, test_candidates_annotations)
score_dead = learner.score(test_dead, test_dead_annotations)

print "final scores:"
print "\t total score:\t\t", score_all
print "\t HC miRNA score:\t", score_miRNA
print "\t candidate score:\t", score_candidates
print "\t dead score:\t\t", score_dead





print "finished", time.clock() - start_time


