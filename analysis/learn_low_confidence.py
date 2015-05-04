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
from sklearn import svm, metrics
import itertools

start_time = time.clock()

def roc_plot(test_data, test_data_annotations, name=""):
    probs = learner.predict_proba(test_data)
    fpr, tpr, _thresholds = metrics.roc_curve(test_data_annotations, probs[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    
    print 
    print "ROC PLOT:"
    print fpr, tpr
    print "area under curve:", roc_auc
    
    title = "ROC plot " + name
    
    f = name.replace(" ", "_")
    f = f.replace(",", "")
    filename = "results/ROC_" + f + ".png"

    plot.plot(fpr, tpr)
    plot.title(title)
    plot.xticks([0.1*x for x in range(0, 11)])
    plot.xlabel("False positive rate")
    plot.yticks([0.1*x for x in range(0, 11)])
    plot.ylabel("True positive rate")
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







print
print "Classifying low confidence miRNAs"
print " - loading data ...",

annotations = pickle.load( open("save_an.p", "rb"))
data = pickle.load( open("save_scaled_data.p", "rb"))
low_confidence_data = pickle.load( open("save_low_confidence_data.p", "rb"))
low_confidence_names = pickle.load( open("save_low_confidence_names.p", "rb"))

print "..loaded"

train_annotations = annotations[1:]
train = data[1:]
low_confidence_data = list(itertools.chain.from_iterable(low_confidence_data))
low_confidence_names = list(itertools.chain.from_iterable(low_confidence_names))

print " - fixed data", len(train), len(train_annotations), time.clock() - start_time


folds = range(len(train))

# deep copy data for multiple threads, probably not needed
copy_train = [numpy.copy(train) for _ in folds]
copy_train_anno = [copy.deepcopy(train_annotations) for _ in folds]


filters = [filter_miRNAs for _ in folds]
# filters = [filter_candidates for _ in folds]
# filters = [filter_dead for _ in folds]
# filters = [None for _ in folds]

threads = len(train)
pool = Pool(threads)

zipzap = zip(copy_train, copy_train_anno, folds, filters)

print " - starting SVM grid search to find best C and gamma", time.clock() - start_time

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
print "Results"
print "score?\t", avg(best_scores), log_avg(best_scores), best_scores
print "c-val\t", log_avg(best_cs),(avg(best_cs)), best_cs
print "gamma\t", log_avg(best_gammas), (avg(best_gammas)), best_gammas
print "\nlog avg c and gammas are used:"
print "\tC:", best_c
print "\tgamma:", best_gamma
print

np_res = numpy.array(results)
array_res = numpy.mean( np_res, axis=0 )


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
# fig, ax = plot.subplots()
# heatmap = ax.pcolor(array_res)
#  
# # put the major ticks at the middle of each cell
# ax.set_xticks(numpy.arange(array_res.shape[0])+0.5, minor=False)
# ax.set_yticks(numpy.arange(array_res.shape[1])+0.5, minor=False)
#  
# ax.set_xticklabels(row_labels, minor=False,  rotation="vertical")
# ax.set_yticklabels(column_labels, minor=False)
# 
# plot.xlabel("gamma")
# plot.ylabel("C")
# plot.show()
# plot.savefig(plot_name)
# print "saved plot", time.clock() - start_time


train = list(itertools.chain.from_iterable(train))
train_annotations = list(itertools.chain.from_iterable(train_annotations))



#===============================================================================
#  first candidates, then HC mirnas, then dead miRNA: 
#  annotations are like [000000111111111000]
#===============================================================================
print "\nFinal testing"
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


#===============================================================================
#  test low confidence here:
#===============================================================================

all_examples = list(itertools.chain.from_iterable(data))
all_annotations = list(itertools.chain.from_iterable(annotations))

test_all = data[0]
test_all_annotations = annotations[0]

# lc_learner = svm.SVC(C=best_c,gamma=best_gamma, cache_size=500, probability=True)
# lc_learner.fit(all_examples, all_annotations)
lc_learner = svm.SVC(C=best_c,gamma=best_gamma, cache_size=500, probability=True)
lc_learner.fit(test_all, test_all_annotations)

print
print "---------------"
print "low confidence data:"
print len(low_confidence_data), len(low_confidence_names)
print 


lc_class = lc_learner.predict(low_confidence_data)


print "low confidence miRNA classified as mirna (score over 0.469)"
print "\tnr", sum(lc_class),
print "\tfraction:",sum(lc_class) / (len(low_confidence_data)*1.0)
print

lc_0_1 = lc_learner.predict_proba(low_confidence_data)


print "sums of probabilities: [(not mirna), (mirna)]"
print sum(lc_0_1)


lc_0_1 = [(a,b) for a,b in lc_0_1]


overhalf = [ (mi, c) for (_not_mi, mi), c in zip(lc_0_1, lc_class) if c == 1]
underhalf = [ (mi, c) for (_not_mi, mi), c in zip(lc_0_1, lc_class) if c != 1]


results_lc = [ (mi, c, n) for (_not_mi, mi), c, n in zip(lc_0_1, lc_class, low_confidence_names)]

# print sum(overhalf)
print 
print "classified as True"
print len(overhalf),
print sorted(overhalf)
print
print "classified as False"
print len(underhalf),
print sorted(underhalf, reverse=True)
print


plot_res = sorted( [b for a,b in lc_0_1] )
plot.title("LC miRNA classification score")
plot.plot(plot_res)
plot.yticks([0.1*x for x in range(0, 11)])
plot.xlabel("Low confidence miRNAs")
plot.ylabel("Score")
plot.savefig("results/LC_classification.png")
plot.show()


filepathname = "results/LC_classification.txt"

sorted_res = sorted(results_lc, reverse=True)



print
print "writing to file..."
with open(filepathname, "w") as write_lc:
    
    for (mirna_score, is_miRNA, miRNA_name) in sorted_res:
        
        line1 = miRNA_name[1:] + "\n"
        line2 = "\tclass:\t" + str(is_miRNA) + "\n"
        line3 = "\tscore:\t" +str(mirna_score) + "\n"
        line4 = "\n"
        
        write_lc.write(line1)
        write_lc.write(line2)
        write_lc.write(line3)
        write_lc.write(line4)
        
print "... wrote results to:", filepathname
        
    
#===============================================================================
# finished testing low confidence
#===============================================================================





roc_plot(test_all, test_all_annotations, "HC, candidates and dead miRNA")
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


#     learner = svm.SVC(probability=True, cache_size=500)
#     learner.fit(train, train_annotations)
#     
#     
#     # roc plot 123
#
#     probs = learner.predict_proba(test_all)
#     fpr, tpr, _thresholds = metrics.roc_curve(test_all_annotations, probs[:,1])
#          
#     roc_auc = metrics.auc(fpr, tpr)
#     print "area under curve:", roc_auc
#      
#     plot.plot(fpr, tpr)
#     plot.show()


print "finished" #, time.clock() - start_time


