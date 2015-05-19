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

hp_candidates = pickle.load( open("save_hp_candidates_new.p", "rb"))

# data = pickle.load( open("save_scaled_data_new.p", "rb"))
# annotations = pickle.load( open("save_an_new.p", "rb"))



data = pickle.load( open("save_scaled_data.p", "rb"))
annotations = pickle.load( open("save_an.p", "rb"))
low_confidence_data = pickle.load( open("save_low_confidence_data.p", "rb"))
low_confidence_names = pickle.load( open("save_low_confidence_names.p", "rb"))

print "..loaded"

train_annotations = annotations[1:]
train = data[1:]
low_confidence_data = list(itertools.chain.from_iterable(low_confidence_data))
low_confidence_names = list(itertools.chain.from_iterable(low_confidence_names))



hp_candidates = list(itertools.chain.from_iterable(hp_candidates))


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

def learn_candidates(d):
    train_stuff, annotation_stuff = d

    lc_learner = svm.SVC(C=best_c,gamma=best_gamma, cache_size=500, probability=True)
    lc_learner.fit(train_stuff, annotation_stuff)
    lc_class = lc_learner.predict_proba(low_confidence_data)
    
    lc_class = map(list, lc_class)
    lc_class = [x[1] for x in lc_class] # using only score for is miRNA

    score_all = lc_learner.score(test_all, test_all_annotations)
    score_miRNA = lc_learner.score(test_miRNA, test_miRNA_annotations)
    score_candidates = lc_learner.score(test_candidates, test_candidates_annotations)
    score_dead = lc_learner.score(test_dead, test_dead_annotations)
    
    phc = lc_learner.predict(test_miRNA)
    plc = lc_learner.predict(low_confidence_data)
    
    print "----------"
    print lc_learner.predict_proba(test_miRNA)
    print sum(phc), len(phc), phc
    print
    print lc_learner.predict_proba(low_confidence_data)
    print sum(plc), len(plc), plc
    
    print
    print "final scores:"
    print "\t total score:\t\t", score_all
    print "\t HC miRNA score:\t", score_miRNA
    print "\t candidate score:\t", score_candidates
    print "\t dead score:\t\t", score_dead
    print
    return lc_class


def learn_new_stuff(d):
    train_stuff, annotation_stuff = d

    lc_learner = svm.SVC(C=best_c,gamma=best_gamma, cache_size=500, probability=True)
    lc_learner.fit(train_stuff, annotation_stuff)
    lc_class = lc_learner.predict_proba(hp_candidates)
    
    lc_class = map(list, lc_class)
    lc_class = [x[1] for x in lc_class] # using only score for is miRNA

    score_all = lc_learner.score(test_all, test_all_annotations)
    score_miRNA = lc_learner.score(test_miRNA, test_miRNA_annotations)
    score_candidates = lc_learner.score(test_candidates, test_candidates_annotations)
    score_dead = lc_learner.score(test_dead, test_dead_annotations)
    
    phc = lc_learner.predict(test_miRNA)
    plc = lc_learner.predict(hp_candidates)
    
    print "----------"
    print lc_learner.predict_proba(test_miRNA)
    print sum(phc), len(phc), phc
    print
    print lc_learner.predict_proba(hp_candidates)
    print sum(plc), len(plc), plc
    
    print
    print "final scores:"
    print "\t total score:\t\t", score_all
    print "\t HC miRNA score:\t", score_miRNA
    print "\t candidate score:\t", score_candidates
    print "\t dead score:\t\t", score_dead
    print
    return lc_class


def placement_scoring(result_lists):
    print
    
    
    placement_lists = map(single_placement, result_lists)
    
    placements = zip(*placement_lists)
    
    
    mean_placements = map(numpy.mean, placements)
    var_placements = map(numpy.std, placements)
#     var_placements = map(numpy.var, placements)
    
    
    for y in range(10):
        print placements[y]
        print mean_placements[y]
        print var_placements[y]
        print
    
    
    all_res = result_lists[0]
    
    
    if True:
#     if False:
        z = zip(mean_placements, var_placements, placements, all_res)
        sz = sorted(z)
        sorted_mean, sorted_stdev, sorted_placements, sorted_res_all = zip(*sz)
        
        
        print "mean vs allscored"
        plot.plot(sorted_mean, sorted_res_all)
        plot.show()

        print "mean", sorted_mean[:20]
        print sorted_placements[0]
        print sorted_res_all[:20]
        plot.plot(sorted_mean[:20])
        plot.show()
        
        print "stdev"
        plot.plot(sorted_stdev)
        plot.show()
        
        plot.plot(sorted_res_all)
        plot.show()
        
        fig, ax1 = plot.subplots()
           
        ax1.plot(sorted_mean)
    #     ax1.set_yscale("symlog")
        ax1.set_ylabel("average classification position 10-fold", color="b")
           
        ax2 = ax1.twinx()
        ax2.plot(sorted_stdev, "r")
        ax2.set_ylabel("stdev classification position 10-fold", color="r")
            
        for tl in ax2.get_yticklabels():
            tl.set_color("r")
    
        plot.show()
    
    
    
    
    
    
    
    
    z = zip(all_res, mean_placements, var_placements, placements)
    sz = sorted(z)
    sorted_res_all, sorted_mean, sorted_stdev, sorted_placements = zip(*sz)
    
    
    
    
#     z = zip(mean_placements, var_placements, placements, all_res)
#     sz = sorted(z)
#     sorted_mean, sorted_stdev, sorted_placements, sorted_res_all = zip(*sz)
    

    for y in range(10):
        print sorted_res_all[y]
        print sorted_mean[y]
        print sorted_stdev[y]
        print sorted_placements[y]
        print
    
    
    fig, ax1 = plot.subplots()
    fig.subplots_adjust(right=0.75)
      
    ax1.plot(sorted_res_all)
#     ax1.set_yscale("symlog")
    ax1.set_ylabel("classification score all", color="b")
      
    ax2 = ax1.twinx()
    ax2.plot(sorted_mean, "r")
    ax2.set_ylabel("average classification position 10-fold", color="r")
    
    
    for tl in ax2.get_yticklabels():
        tl.set_color("r")
        
    ax3 = ax1.twinx()
    ax3.plot(sorted_stdev, "g")
#     ax3.set_yscale("symlog")
    ax3.spines['right'].set_position(('axes', 1.2))
    ax3.set_frame_on(True)
    ax3.patch.set_visible(False)
    
    ax3.set_ylabel("variance classification position 10-fold", color="g")
    for tl in ax3.get_yticklabels():
        tl.set_color("g")
        
    
      
#     plot.title("Low confidence: reads vs classification score")
    plot.show()

#     
#     fig, ax1 = plot.subplots()
#       
#     ax1.plot(sorted_res_all)
# #     ax1.set_yscale("symlog")
#     ax1.set_ylabel("classification score all", color="b")
#       
#     ax2 = ax1.twinx()
#     ax2.plot(sorted_mean, "r")
#     ax2.set_ylabel("average classification position 10-fold", color="r")
#        
#     for tl in ax2.get_yticklabels():
#         tl.set_color("r")
#         
#     ax3 = ax1.twinx()
#     ax3.plot(sorted_stdev, "g")
# #     ax3.set_yscale("symlog")
#     ax3.set_ylabel("variance classification position 10-fold", color="g")
#       
#     for tl in ax3.get_yticklabels():
#         tl.set_color("g")
#       
# #     plot.title("Low confidence: reads vs classification score")
#     plot.show()
#      
#     
#     plot.plot(sorted_res_all)
# #     plot.plot(sorted_mean)
#     plot.show()
#     
#     
#     plot.plot(sorted_res_all, sorted_mean, 'rx')
#     plot.show()
#     
#     
#     plot.plot(sorted_res_all, sorted_stdev)
#     plot.show()
#     plot.plot(sorted_stdev, sorted_mean)
#     plot.show()
    
    assert 0

def single_placement(result_list):
    ''' returns the placement (by score), relative to the original list '''

#     add numbers find position for all results
    positions = range(len(result_list))
    
    pairs = zip(result_list, positions)
    sorted_pairs = sorted(pairs)
    
    _sorted_result, original_positions = zip(*sorted_pairs) # only need the position
    
    
    reverse_pos = reversed(positions)
    position_pairs = zip(original_positions, reverse_pos)
    score_positions = sorted(position_pairs)
    
    _original_pos, placements = zip(*score_positions)
    
    return placements
    
    
    
    
    


#     first sort by first entry

#     then sort by 


def predict_folds(data_list, annotation_list):
    
    print "10fold"
    
    res = map(learn_candidates, zip(data_list,annotation_list))
    
    all_d = list(itertools.chain.from_iterable(data_list))
    all_a = list(itertools.chain.from_iterable(annotation_list))
    
    best_res = learn_candidates((all_d, all_a))
    
    
    
    res.insert(0, best_res) # insert best results first 123
    
    
    print "making placement scores 123"
    placement_scoring(res)
    
    pickle.dump(res, open("lc_scores_10.p", "wb"))
    print len(res)
    
    
    
    sorted_by_best = sorted(zip(*res)) # all results sorted by result using all folds
    
    

    
    
    lc_reads = pickle.load( open("low_confidence_reads.p", "rb"))
    
    zip_best_reads = zip(best_res, lc_reads)
    sorted_reads = sorted(zip_best_reads)
    
    print len(sorted_reads)
    reads_sorted = zip(*sorted_reads)[1]
#     reads_sorted_logscaled = map(math.log, reads_sorted)
    
    mean_results = map(numpy.mean, sorted_by_best)
    stdev_results = map(numpy.var, sorted_by_best)
    
    
    
    print len(sorted_by_best), 1550
    
    res_sorted = zip(*sorted_by_best)
    
    best_res = res_sorted[0]
    print len(mean_results), len(stdev_results), len(best_res)
    
    plot.plot(mean_results, label="mean value")
    plot.plot(stdev_results, label="variance")
    plot.plot(best_res, label="score using all data")
#     plot.plot(reads_sorted_logscaled, label="reads, log scaled")
    plot.legend(loc='upper left')
    plot.show()
    
    print len(res_sorted), 11
    print len(res_sorted[0]), 1550
    
    
    
    
    
    
    
    
    
    
    
    
#     res = list(res)
#     for i in range(len(res)):
#         
#         res = list(res)
#         x = res.pop(0)
#         res.append(x)
#         sorted_by_best = sorted(zip(*res)) # all results sorted by result using all folds
#         
#         res = zip(*sorted_by_best)
#         
#         for r in reversed(res):
#             plot.plot(r)   
#         plot.savefig("results/lc_sorted_" + str(i) + ".png") 
#         plot.show()
#         
#     
#     
#     all_sorted = map(sorted, res_sorted)
#     
#     
#     for res in reversed(all_sorted):
#         plot.plot(res)    
#     plot.savefig("results/lc_all_sorted.png")
#     plot.show()
        

   

# predict_folds(data, annotations, learn_new_stuff)
predict_folds(data, annotations)

# assert 0

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

score_list = [b for a,b in lc_0_1]

print max(score_list)
print score_list


plot_res = sorted(score_list)
plot.title("LC miRNA classification score")
plot.plot(plot_res)
plot.yticks([0.1*x for x in range(0, 11)])
plot.xlabel("Low confidence miRNAs")
plot.ylabel("Score")
plot.savefig("results/LC_classification.png")
plot.show()



read_param = [x[0] for x in low_confidence_data]


val_tups = zip(read_param, score_list)

val_tups = sorted(val_tups)

read_param, score_list = zip(*val_tups)

plot.plot(read_param)
plot.plot(score_list)
plot.title("LC miRNA classification score")
# plot.yticks([0.1*x for x in range(0, 11)])
plot.xlabel("Low confidence miRNAs")
plot.ylabel("Score")
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
print "c-val\t", log_avg(best_cs),(avg(best_cs)), best_cs
print "gamma\t", log_avg(best_gammas), (avg(best_gammas)), best_gammas

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


