'''
Created on 21. mar. 2015

@author: hakon
'''

import numpy
import copy
from multiprocessing import Pool
# from multiprocessing import current_process
# from multiprocessing import cpu_count
from ml.grid import one_grid, filter_candidates, filter_miRNAs
import pickle
from matplotlib import pyplot as plot
import time
import math
from sklearn import svm, metrics
import itertools
from numpy import unravel_index
# import matplotlib
from matplotlib import cm




def roc_plot(learner, test_data, test_data_annotations, name=""):
    probs = learner.predict_proba(test_data)
    fpr, tpr, _thresholds = metrics.roc_curve(test_data_annotations, probs[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    
    print 
    print "ROC PLOT:"
#     print fpr
#     print
#     print tpr
    print "area under curve:", roc_auc
    
    title = "ROC plot " + name
    
#     f = name.replace(" ", "_")
#     f = f.replace(",", "")
    filename = "results/ROC_candidates.png"

    plot.plot(fpr, tpr)
    plot.title(title)
    plot.xticks([0.1*x for x in range(0, 11)])
    plot.xlabel("False positive rate")
    plot.yticks([0.1*x for x in range(0, 11)])
    plot.ylabel("True positive rate")
    plot.ylim(ymin=0,ymax=1)
    plot.savefig(filename)
    plot.show()





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

def remove_candidates(data, annotations, keep=0):
    ''' removes canididates, and their annotations (0s) '''
    not_candidate = annotations.index(1)
    
    data = data[not_candidate-keep:]
    annotations = annotations[not_candidate-keep:]
    return data, annotations


def heatplot(array_res, c_values, gamma_values, filters):
    
    
    gamma_name = "_gamma_" + str(gamma_values[0]) + "-" + str(gamma_values[-1])
    c_name = "C_" + str(c_values[0]) + "-" + str(c_values[-1])
    f_name = "" if not filters[0] else filters[0].__name__
    plot_name = "plot/Grid_avg_" + f_name + c_name + gamma_name + ".png"
    
    
    row_labels = gamma_values
    column_labels = c_values
    
    
    #  create gridsearch heatmap
    fig, ax = plot.subplots()
    _heatmap = ax.pcolor(array_res)
    
    heatmap = ax.pcolor(array_res, cmap=cm.get_cmap("RdBu"), vmin=0.0, vmax=1.0)
      
    # put the major ticks at the middle of each cell
    ax.set_xticks(numpy.arange(array_res.shape[1])+0.5, minor=False)
    ax.set_yticks(numpy.arange(array_res.shape[0])+0.5, minor=False)
      
    ax.set_xticklabels(row_labels, minor=False,  rotation="vertical")
    ax.set_yticklabels(column_labels, minor=False)
    
    heatmap.axes.set_ylim(0, len(array_res))
    heatmap.axes.set_xlim(0, len(array_res))
    
    fig.colorbar(heatmap)
    
    plot.title("Grid search average scores")
    plot.xlabel("gamma")
    plot.ylabel("C")
    plot.gca().axis('tight')
    plot.savefig(plot_name)
    plot.show()
    print "saved plot", time.clock() - start_time

# from ml import vectorize
# from sklearn import preprocessing
# import math
start_time = time.clock()


print
print "Classifying new candidates as miRNA or not"
print " - loading data ...",

data_new = pickle.load( open("save_scaled_data_new.p", "rb"))
annotations_new = pickle.load( open("save_an_new.p", "rb"))
hp_candidates = pickle.load( open("save_hp_candidates_new.p", "rb"))


# for testing purposes:
low_confidence_data = pickle.load( open("save_low_confidence_data.p", "rb"))
low_confidence_reads = pickle.load( open("low_confidence_reads.p", "rb"))

print "..loaded"

# annotations = pickle.load( open("save_an.p", "rb"))
# data = pickle.load( open("save_scaled_data.p", "rb"))

train = data_new[1:]
train_annotations = annotations_new[1:]

test_all = data_new[0]
test_all_annotations = annotations_new[0]

all_examples = list(itertools.chain.from_iterable(data_new))
all_annotations = list(itertools.chain.from_iterable(annotations_new))



some_lc = [lc for lc, count in zip(low_confidence_data[0], low_confidence_reads) if count > 0]

print " - fixed data", len(train), len(train_annotations), time.clock() - start_time


folds = range(len(train))

# deep copy data for multiple threads, probably not needed
copy_train = [numpy.copy(train) for _ in folds]
copy_train_anno = [copy.deepcopy(train_annotations) for _ in folds]


# filters = [filter_miRNAs for _ in folds]
# filters = [filter_candidates for _ in folds]
# filters = [filter_dead for _ in folds]
filters = [None for _ in folds]

threads = len(train)
pool = Pool(threads)

zipzap = zip(copy_train, copy_train_anno, folds, filters)

print " - starting SVM grid search to find best C and gamma. Time used:", time.clock() - start_time

res = pool.map(one_grid, zipzap)
# res = map(one_grid, zipzap)
print " - finished SVM grids", time.clock() - start_time



best_scores, best_cs, best_gammas, results = zip(*res)

def log_avg(l): 
    '''   average of log values, then scaled back '''
    log_list = (math.log(x, 2.0) for x in l)
    log_val = sum(log_list) / len(l)
    avg_val = 2.0 ** log_val 
    return avg_val

def avg(l):
    return sum(l) / len(l)

best_c = log_avg(best_cs)
best_gamma = log_avg(best_gammas)

print 
print "Results", time.clock() - start_time
print "score?\t", avg(best_scores), log_avg(best_scores), best_scores
print "c-val\t", avg(best_cs), log_avg(best_cs), best_cs
print "gamma\t", avg(best_gammas), log_avg(best_gammas), best_gammas
print "\nlog avg c and gammas are used:"
print "\tC:", best_c
print "\tgamma:", best_gamma
print




# gamma_values = [10.0 ** i for i in xrange(-4,5)]
# c_values = [10.0 ** i for i in xrange(-4,5)]
# gamma_values = [10.0 ** i for i in xrange(-4,5)]
# c_values = [10.0 ** i for i in xrange(-4,5)]


# gamma_values = [4.0 ** i for i in xrange(-4,3)]
# c_values = [4.0 ** i for i in xrange(-3,4)]

gamma_values = [2.0 ** i for i in xrange(-18,-1)]
c_values = [2.0 ** i for i in xrange(-1,10)]

np_res = numpy.array(results)
array_res = numpy.mean(np_res, axis=0)

print unravel_index(array_res.argmax(), array_res.shape)

# print numpy.argmax(array_res, axis=0)
# _c, _gamma = numpy.argmax(array_res, axis=0)

_c, _gamma = unravel_index(array_res.argmax(), array_res.shape)
print "sum all best values:", _c, _gamma, c_values[_c], gamma_values[_gamma]





heatplot(array_res, c_values, gamma_values, filters)





train_folds = train
train_folds_annotations = train_annotations

train = list(itertools.chain.from_iterable(train))
train_annotations = list(itertools.chain.from_iterable(train_annotations))



#===============================================================================
#  first candidates, then HC mirnas, then dead miRNA: 
#  annotations are like [000000111111111000]
#===============================================================================
print "\n final testing"


test_miRNA = [c for c,a in zip(test_all, test_all_annotations) if a ]
test_miRNA_annotations = [1] * len(test_miRNA)


first_miRNA = test_all_annotations.index(1)
test_candidates = test_all[:first_miRNA]
test_candidates_annotations = test_all_annotations[:first_miRNA]

first_dead = test_all_annotations.index(0,first_miRNA)
test_dead = test_all[first_dead:]
test_dead_annotations = test_all_annotations[first_dead:]

test_not_c, test_not_c_annotations = remove_candidates(test_all, test_all_annotations, 5)

some_lc_an = [2]*len(some_lc)
# testwithlc = list(test_not_c) + some_lc
# testwithlc_an = list(test_not_c_annotations) + some_lc_an
testwithlc = list(test_all) + some_lc
testwithlc_an = list(test_all_annotations) + some_lc_an

print "adding low confidence", len(some_lc)

# t2, te_a = filter_miRNAs(test_all, test_all_annotations)
# 
# test_not_c += t2
# test_not_c_annotations += te_a

# print first_miRNA, first_dead, len(test_all_annotations)
# print len(test_candidates), len(test_miRNA), len(test_dead)
# print len(test_candidates) + len(test_miRNA) + len(test_dead)





#===============================================================================
# test candidates here
#===============================================================================




learner = svm.SVC(C=best_c,gamma=best_gamma, cache_size=500, probability=True)

learner.fit(train, train_annotations)

#===============================================================================
# #  NO NO NO training SVM using test data :(
# # # learner.fit(test_all, test_all_annotations) # NO NO NO NO NO NO
# # # learner.fit(test_not_c, test_not_c_annotations)
# # # learner.fit(testwithlc, testwithlc_an)
# # # learner.fit(all_examples, all_annotations)
#===============================================================================

print
print "------"
print sum(all_annotations), len(all_annotations)
print "classifying new candidates:"
print len(hp_candidates)

candidate_classes = learner.predict(hp_candidates)

print "nr of miRNAs:", sum(candidate_classes), candidate_classes

pickle.dump(candidate_classes, open("candidate_classified_miRNA.p", "wb"))



# def round_number(number, limit=0.99):
#     return 1 if number > limit else 0

candidate_scores = learner.predict_proba(hp_candidates)

print "sum scores:"
print sum(candidate_scores)
print candidate_scores
print

cc = list(candidate_classes)

csl = map(list, candidate_scores)
csl = list(csl)

print csl

cscores = zip(*csl)[1] 

print cscores

rounded_scores = [1 if s > 0.99 else 0 for s in cscores]

pickle.dump(rounded_scores, open("candidate_classified_99.p", "wb"))


# csc = list(candidate_scores[1])


assert len(cc) == len(cscores), (len(cc), len(cscores))

for s, c in sorted(zip(cscores, cc)):
    print c, s
print "----------"




# candidate_scores = [a for a,b in candidate_scores]
candidate_scores = list(candidate_scores)
candidate_scores = map(list, candidate_scores)

print len(candidate_scores), len(candidate_scores[0])

candidate_score_lists = zip(*candidate_scores)

print len(candidate_score_lists), len(candidate_score_lists[0])

# c_scored = [(t,c) for t,c in zip(candidate_scores, hp_candidates)]



sorted_cs = map(sorted, candidate_score_lists)

print len(candidate_scores), len(candidate_scores[0]), sorted_cs[0]

scaled123 = pickle.load( open("save_scaled_data.p", "rb"))


# plot.title("new miRNA candidates classification score")
# for s in sorted_cs:
#     plot.plot(s)
#     break
#     
# 
# plot.ylim(ymin=0,ymax=1)
# plot.show()


candidate_scores = learner.predict_proba(hp_candidates)
candidate_scores_list = map(list, candidate_scores)
candidate_score_trues = [x[1] for x in candidate_scores_list]

candidate_score_trues_sorted = sorted(candidate_score_trues)

# plot.title("new miRNA candidates classification score 2")
# plot.ylim(ymin=0,ymax=1)
# plot.plot(candidate_score_trues_sorted)
# plot.show()

#==============================================================================
# plot results using 9 folds + the best 
#==============================================================================

def score_once(param):
    
    train_, train_an_ = param
    
    learner = svm.SVC(C=best_c,gamma=best_gamma, cache_size=500, probability=True)
    learner.fit(train_, train_an_)
    _scores = learner.predict_proba(hp_candidates)
    _scores = map(list, _scores)
    _scores = [x[1] for x in _scores]
    
    return _scores


def sorted_by_first_list(several_lists):
    
    itemwise = zip(*several_lists)
    sorted_itemwise = sorted(itemwise)
    sorted_by_first = zip(*sorted_itemwise)
    return sorted_by_first

#===============================================================================
# loading data from mirdeep 
#===============================================================================

# train_folds
# train_folds_annotations



fold_scores = map(score_once, zip(train_folds, train_folds_annotations))


sd_scores = map(numpy.std, zip(*fold_scores))

fold_scores.insert(0, candidate_score_trues)

sorted_fold_scores = sorted_by_first_list(fold_scores)


# plot.title("new miRNA candidates classification score 3")
# plot.ylim(ymin=0,ymax=1)
# plot.plot(sorted_fold_scores[0])
# plot.show()


all_and_sd = sorted_by_first_list([candidate_score_trues, sd_scores])

plot.title("miRNA candidate classification score and SD")
plot.ylim(ymin=0,ymax=1)
# plot.plot(sorted_fold_scores[0])
plot.plot(all_and_sd[0], label="classification score")
plot.plot(all_and_sd[1], label="9-fold score SD")
plot.legend(loc="upper left")
plot.show()


# plot.title("new miRNA candidates classification score 5")
# plot.ylim(ymin=0,ymax=1)
# map(plot.plot, sorted_fold_scores)
# plot.show()
# 
# 
# each_sorted = map(sorted, fold_scores)
# plot.title("new miRNA candidates classification score 5")
# plot.ylim(ymin=0,ymax=1)
# map(plot.plot, each_sorted)
# plot.show()




classes_mirdeep = pickle.load( open("classes_mirdeep.p", "rb"))
scores_mirdeep = pickle.load( open("scores_mirdeep.p", "rb"))


fpr_mirdeep, tpr_mirdeep, _thresholds = metrics.roc_curve(classes_mirdeep, scores_mirdeep)
roc_auc = metrics.auc(fpr_mirdeep, tpr_mirdeep)

found_by_mirDeep = pickle.load( open("found_by_mirDeep.p", "rb"))




#===============================================================================
# test the precisi
#===============================================================================

print len( zip(test_all, found_by_mirDeep))
print len( test_all)
print len(found_by_mirDeep)


# roc_plot(learner, test_all, test_all_annotations)
# roc_plot(test_miRNA, test_miRNA_annotations)
# roc_plot(test_candidates, test_candidates_annotations)
# roc_plot(test_dead, test_dead_annotations)


score_all = learner.score(test_all, test_all_annotations)
score_miRNA = learner.score(test_miRNA, test_miRNA_annotations)
score_candidates = learner.score(test_candidates, test_candidates_annotations)
score_dead = learner.score(test_dead, test_dead_annotations)



print "final scores:"
print "\t total score:\t\t", score_all
print "\t HC miRNA score:\t", score_miRNA
print "\t candidate score:\t", score_candidates
print "\t dead score:\t\t", score_dead


print "\n"
print "\n-------------------"
# zipzap = zip(copy_train, copy_train_anno, folds, filters)

def all_folds_classification(param):
    
    (all_data, all_annotations, test_fold_nr, filter_test) = param

    test = all_data[test_fold_nr]
    test_annotations = all_annotations[test_fold_nr]

    
    test_hc, test_hc_an = filter_miRNAs(test, test_annotations)
    test_cand, test_cand_an = filter_candidates(test, test_annotations)
    test_dead, test_dead_an = filter_dead(test, test_annotations)

    
    p1 = list(all_data[:test_fold_nr])
    p2 = list(all_data[test_fold_nr+1:])
    
    train_lists = p1 + p2
    train_annotations_lists = all_annotations[:test_fold_nr] + all_annotations[test_fold_nr+1:]
    
    train = list(itertools.chain.from_iterable(train_lists))
    train_annotations = list(itertools.chain.from_iterable(train_annotations_lists))
    
    
    learner = None
    learner = svm.SVC(C=best_c, gamma=best_gamma, cache_size=500, probability=True)
    learner.fit(train, train_annotations)
    
    
    
    probs = learner.predict_proba(test)
    fpr, tpr, _thresholds = metrics.roc_curve(test_annotations, probs[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    
    print "area under curve:\t", roc_auc
    
    plot.plot(fpr, tpr)
    plot.title("ROC plot HC training folds")
    plot.xticks([0.1*x for x in range(0, 11)])
    plot.xlabel("False positive rate")
    plot.yticks([0.1*x for x in range(0, 11)])
    plot.ylabel("True positive rate")
    plot.ylim(ymin=0,ymax=1)
    # plot.savefig(filename)


    
    return (fpr, tpr)


res = map(all_folds_classification, zipzap)

plot.show()

# 
# map(plot.plot, res)
# # plot.plot(fpr, tpr)
# plot.title("ROC plot training folds")
# plot.xticks([0.1*x for x in range(0, 11)])
# plot.xlabel("False positive rate")
# plot.yticks([0.1*x for x in range(0, 11)])
# plot.ylabel("True positive rate")
# plot.ylim(ymin=0,ymax=1)
# # plot.savefig(filename)
# plot.show()


print "finished" #, time.clock() - start_time






















