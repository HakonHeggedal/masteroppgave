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
from matplotlib import cm
from misc.correlation_plots import plot_correlation




FEATURE_NAMES = []
FEATURE_NAMES.append("mapped sequences")
# FEATURE_NAMES.append("mapped sequences rpm count")
FEATURE_NAMES.append("hairpin energy")
# FEATURE_NAMES.append("hairpin + 10nt energy")
FEATURE_NAMES.append("entropy nucleotides")
FEATURE_NAMES.append("entropy structure")
FEATURE_NAMES.append("heterogenity 5p start")
FEATURE_NAMES.append("heterogenity 5p end")
FEATURE_NAMES.append("heterogenity 3p start")
FEATURE_NAMES.append("heterogenity 3p end")
FEATURE_NAMES.append("candidate quality")
# FEATURE_NAMES.append("nr of bindings hairpin+10nt")
FEATURE_NAMES.append("overhang level outer")
FEATURE_NAMES.append("overhang outer")
FEATURE_NAMES.append("overhang level inner")
FEATURE_NAMES.append("overhang inner")
# FEATURE_NAMES.append("bulge factor")
FEATURE_NAMES.append("short seq align score")
FEATURE_NAMES.append("short vs long seqs")
# FEATURE_NAMES.append("small sequences aligned 5p")
# FEATURE_NAMES.append("small sequences aligned 3p")
FEATURE_NAMES.append("loop size")
FEATURE_NAMES.append("folds 5p")
FEATURE_NAMES.append("folds 3p")
FEATURE_NAMES.append("folds before")
FEATURE_NAMES.append("folds after")
FEATURE_NAMES.append("has hairpin structure")

# lc_10  = pickle.load( open("lc_scores_10.p", "rb"))
# lc_10_all = pickle.load( open("lc_scores_10_w_cand.p", "rb"))
# lc_names = pickle.load( open("save_low_confidence_names.p", "rb"))
# 
# low_confidence_reads = pickle.load( open("low_confidence_reads.p", "rb"))
# 
# 
# lc_names = list(itertools.chain.from_iterable(lc_names))
# 
# lc_score_all = lc_10_all[0]
# lc_score_nonhp = lc_10[0]
# 
# 
# 
# 
# print "before", len(lc_names)
# lc_names = [lc for lc, count in zip(lc_names, low_confidence_reads) if count > 0]
# print "after", len(lc_names)
# 
# 
# sorted_score_pairs_all = sorted(zip(lc_score_all, lc_names))
# scores_all, names_all = zip(*sorted_score_pairs_all)
# 
# 
# sorted_score_pairs_nonhp = sorted(zip(lc_score_nonhp, lc_names))
# scores_nonhp, names_nonhp = zip(*sorted_score_pairs_nonhp)
# 
# 
# LC_100_best_all = names_all[-100:]
# LC_100_best_nonhp = names_nonhp[-100:]
#  
# LC_100_worst_all = names_all[:100]
# LC_100_worst_nonhp = names_nonhp[:100]
# 
# print LC_100_best_all
# print LC_100_best_nonhp
# print LC_100_worst_all
# print LC_100_worst_nonhp
# 
# 
# pickle.dump(LC_100_best_all, open("LC_100_best_all.p", "wb"))
# pickle.dump(LC_100_best_nonhp, open("LC_100_best_nonhp.p", "wb"))
# pickle.dump(LC_100_worst_all, open("LC_100_worst_all.p", "wb"))
# pickle.dump(LC_100_worst_nonhp, open("save_scaled_data_new.p", "wb"))


# # plot.plot(scores_all[-100:])
# 
# plot.plot(scores_nonhp[-202:])
# # plot.plot(scores_all[:100])
# 
# plot.show()



# assert 0


# lc_10  = pickle.load( open("lc_scores_10.p", "rb"))
# 
# print len(lc_10)
# print len(lc_10[-1])
# assert 0

# lc_reads = pickle.load( open("low_confidence_reads.p", "rb"))
# 
# 
# # lc_reads = lc_reads[:100]
# 
# lc_reads = sorted(lc_reads, reverse=True)
# 
# lc_pos = range(1, len(lc_reads)+1)
# 
# plot.plot(lc_pos, lc_reads, label="reads")
# 
# 
# xticvals = [x for x in range(0, len(lc_reads)-100, 200)]
# xticvals.remove(800)
# maxel = len(lc_reads)+1
# zeroel = lc_reads.index(0.0)+1
# xticvals.append(maxel)
# xticvals.append(zeroel)
# xticvals = sorted(xticvals)
# print xticvals, maxel, zeroel
# 
# 
# 
# 
# plot.xlabel("LC miRNA")
# plot.ylabel("Reads")
# 
# plot.xlim(0, len(lc_reads)+1)
# plot.yscale("symlog")
# plot.xticks(xticvals)
# # plot.legend(loc='upper right')
# plot.title("LC miRNA reads")
# plot.show()
# 
# assert 0


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
    plot.ylim(ymin=0, ymax=1)
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
    
    plot.title("Grid search C and gamma")
    plot.xlabel("gamma")
    plot.ylabel("C")
    plot.gca().axis('tight')
    
    plot.savefig(plot_name)
    plot.show()
    print "saved plot", time.clock() - start_time



def compute_kernel_param(train, train_annotations, filters):

    folds = range(len(train))
    
    # deep copy data for multiple threads, probably not needed
    copy_train = [numpy.copy(train) for _ in folds]
    copy_train_anno = [copy.deepcopy(train_annotations) for _ in folds]
    
    

    
    threads = len(train)
    pool = Pool(threads)
    
    zipzap = zip(copy_train, copy_train_anno, folds, filters)
    
    print " - starting SVM grid search to find best C and gamma"
    
    res = pool.map(one_grid, zipzap)
    # res = map(one_grid, zipzap)
    print " - finished SVM grids"
    
    
    best_scores, best_cs, best_gammas, results = zip(*res)
    
    return best_scores, best_cs, best_gammas, results


print
print "Classifying low confidence miRNAs"
print " - loading data ...",

hp_candidates = pickle.load( open("save_hp_candidates_new.p", "rb"))

# data = pickle.load( open("save_scaled_data_new.p", "rb"))
# annotations = pickle.load( open("save_an_new.p", "rb"))

low_confidence_reads = pickle.load( open("low_confidence_reads.p", "rb"))

data = pickle.load( open("save_scaled_data.p", "rb"))
annotations = pickle.load( open("save_an.p", "rb"))
low_confidence_data = pickle.load( open("save_low_confidence_data.p", "rb"))
low_confidence_names = pickle.load( open("save_low_confidence_names.p", "rb"))

data_new = pickle.load( open("save_scaled_data_new.p", "rb"))
annotations_new = pickle.load( open("save_an_new.p", "rb"))

print "..loaded"

print len(data)
print 0
d = list(itertools.chain.from_iterable(data))
d2 = list(itertools.chain.from_iterable(data_new))
print len(d)
print len(d[0])
print "---"
# assert 0

plot_correlation(d, FEATURE_NAMES, False)
plot_correlation(d, FEATURE_NAMES, False)


#===============================================================================
# print "removing candidates from data"
# 
# data = data_new
# annotations = annotations_new
#===============================================================================


# # print map(sum, annotations)
# # # print filter_dead(data[0], )
# # 
# # # print len([filter_dead(d,a) for d,a in zip(data,annotations) ])
# print map(len, ([filter_dead(d,a)[1] for d,a in zip(data,annotations) ]))
# print map(len, ([filter_candidates(d,a)[1] for d,a in zip(data,annotations) ]))
# print map(len, ([filter_miRNAs(d,a)[1] for d,a in zip(data,annotations) ]))
# # print len([filter_miRNAs(d,a) for d,a in zip(data,annotations) ])
#  
#  
# assert 0



train_annotations = annotations[1:]
train = data[1:]
low_confidence_data = list(itertools.chain.from_iterable(low_confidence_data))
low_confidence_names = list(itertools.chain.from_iterable(low_confidence_names))



hp_candidates = list(itertools.chain.from_iterable(hp_candidates))













print " - fixed data", len(train), len(train_annotations), time.clock() - start_time


folds = range(len(train))

# filters = [filter_miRNAs for _ in folds]
# filters = [filter_candidates for _ in folds]
# filters = [filter_dead for _ in folds]
filters = [None for _ in folds]

best_scores, best_cs, best_gammas, grid_results = compute_kernel_param(train, train_annotations, filters)

print map(len, [best_scores, best_cs, best_gammas, grid_results])
pickle.dump(grid_results, open("svm_grid_results.p", "wb"))
pickle.dump(best_scores, open("svm_grid_best_scores.p", "wb"))
pickle.dump(best_cs, open("svm_grid_best_cs.p", "wb"))
pickle.dump(best_gammas, open("svm_grid_best_gammas.p", "wb"))

grid_results = pickle.load( open("svm_grid_results.p", "rb"))
best_scores = pickle.load( open("svm_grid_best_scores.p", "rb"))
best_cs = pickle.load( open("svm_grid_best_cs.p", "rb"))
best_gammas = pickle.load( open("svm_grid_best_gammas.p", "rb"))

print map(len, [best_scores, best_cs, best_gammas, grid_results])
print "------------++++++++"



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

np_res = numpy.array(grid_results)
array_res = numpy.mean( np_res, axis=0 )


# gamma_values = [10.0 ** i for i in xrange(-4,5)]
# c_values = [10.0 ** i for i in xrange(-4,5)]
# gamma_values = [10.0 ** i for i in xrange(-4,5)]
# c_values = [10.0 ** i for i in xrange(-4,5)]


# gamma_values = [4.0 ** i for i in xrange(-4,3)]
# c_values = [4.0 ** i for i in xrange(-3,4)]



gamma_values = [2.0 ** i for i in xrange(-18,-1)]
c_values = [2.0 ** i for i in xrange(-1,10)]


gamma_name = "_gamma_" + str(gamma_values[0]) + "-" + str(gamma_values[-1])
c_name = "C_" + str(c_values[0]) + "-" + str(c_values[-1])
f_name = "" if not filters[0] else filters[0].__name__
plot_name = "plot/Grid_avg_" + f_name + c_name + gamma_name + ".png"


row_labels = gamma_values
column_labels = c_values


gamma_text = ["2^"+str(i) for i in xrange(-18,-1)]

heatplot(array_res, c_values, gamma_text, filters)


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



print "C", best_c
print "gamma", best_gamma
assert 0


#===============================================================================
#  test low confidence here:
#===============================================================================

print len(low_confidence_data)
low_confidence_data = [lc for lc, count in zip(low_confidence_data, low_confidence_reads) if count > 0]
print len(low_confidence_data)

def learn_candidates(d):
    
    train_stuff, annotation_stuff = d


    lc_learner = svm.SVC(C=best_c,gamma=best_gamma, cache_size=500, probability=True)
    lc_learner.fit(train_stuff, annotation_stuff)
    lc_class = lc_learner.predict_proba(low_confidence_data)
    
    lc_class = map(list, lc_class)
    lc_class = [x[1] for x in lc_class] # using only score for is miRNA

#     score_all = lc_learner.score(test_all, test_all_annotations)
#     score_miRNA = lc_learner.score(test_miRNA, test_miRNA_annotations)
#     score_candidates = lc_learner.score(test_candidates, test_candidates_annotations)
#     score_dead = lc_learner.score(test_dead, test_dead_annotations)
#     
#     phc = lc_learner.predict(test_miRNA)
#     plc = lc_learner.predict(low_confidence_data)
#     
#     print "----------"
#     print lc_learner.predict_proba(test_miRNA)
#     print sum(phc), len(phc), phc
#     print
#     print lc_learner.predict_proba(low_confidence_data)
#     print sum(plc), len(plc), plc
#     
#     print
#     print "final scores:"
#     print "\t total score:\t\t", score_all
#     print "\t HC miRNA score:\t", score_miRNA
#     print "\t candidate score:\t", score_candidates
#     print "\t dead score:\t\t", score_dead
#     print
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

    
    
    placement_lists = map(single_placement, result_lists)
    
    placements = zip(*placement_lists)
    
    
    mean_placements = map(numpy.mean, placements)
    var_placements = map(numpy.std, placements)
#     var_placements = map(numpy.var, placements)
    
    print
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
        
        
#         print "mean vs allscored"
#         plot.plot(sorted_mean, sorted_res_all)
#         plot.title("10 fold mean vs all data score")
#         plot.show()
# 
#         print "mean", sorted_mean[:20]
#         print sorted_placements[0]
#         print sorted_res_all[:20]
#         plot.plot(sorted_mean[:20])
#         plot.title("20 best avg results")
#         plot.show()
#         
#         print "stdev"
#         plot.plot(sorted_stdev)
#         plot.title("stdev")
#         plot.show()
#         
#         plot.plot(sorted_res_all)
#         plot.title("all folds sorted by mean placement result")
#         plot.show()
        
        fig, ax1 = plot.subplots()
        
        
        ax1.plot(sorted_stdev, "r")
        ax1.set_ylabel("st.dev. classification position 9-folds + all folds", color="r")
        ax1.set_xlabel("candidates")
        
#         ax1.yaxis.tick_right()
        for tl in ax1.get_yticklabels():
            tl.set_color("r")
           
        ax2 = ax1.twinx()

        ax2.plot(sorted_mean, "b")
        for tl in ax2.get_yticklabels():
            tl.set_color("b")
    #     ax2.set_yscale("symlog")
        ax2.set_ylabel("average classification position 9-folds + all folds", color="b")
        plot.title("Average classification position ")
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
    ax1.set_ylabel("classification score all data", color="b")
    ax1.set_xlabel("LC miRNA")
      
    ax2 = ax1.twinx()
    ax2.set_axis_bgcolor('red')
    ax2.plot(sorted_mean, "r")
    
    ax2.set_ylabel("average classification position 9-fold", color="r")
    
    
    for tl in ax2.get_yticklabels():
        tl.set_color("r")
        
    ax3 = ax1.twinx()
    ax3.plot(sorted_stdev, "g", alpha=0.85)
#     ax3.set_yscale("symlog")
    ax3.spines['right'].set_position(('axes', 1.2))
    ax3.set_frame_on(True)
    ax3.patch.set_visible(False)
    
    ax3.set_ylabel("variance classification position 9-fold", color="g")
    for tl in ax3.get_yticklabels():
        tl.set_color("g")
        
    
    ax1.set_zorder(ax3.get_zorder()+2) # put ax1 in front of ax3
#     ax2.set_zorder(ax3.get_zorder()+1) # put ax1 in front of ax3
    ax1.patch.set_visible(False) # hide the 'canvas'
#     ax2.patch.set_visible(True) # hide the 'canvas'
    plot.title("LC scores with 9-fold position stability")
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
    
#     assert 0

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
    pickle.dump(res, open("lc_scores_10_w_cand.p", "wb"))
    print len(res)
    
    
    
    sorted_by_best = sorted(zip(*res)) # all results sorted by result using all folds
    
    

    
    
    lc_reads = pickle.load( open("low_confidence_reads.p", "rb"))
    
    zip_best_reads = zip(best_res, lc_reads)
    sorted_reads = sorted(zip_best_reads)
    
    print len(sorted_reads)
    reads_sorted = zip(*sorted_reads)[1]
#     reads_sorted_logscaled = map(math.log, reads_sorted)
    
    mean_results = map(numpy.mean, sorted_by_best[1:])
    stdev_results = map(numpy.var, sorted_by_best[1:])
    
    
    
    print len(sorted_by_best), 1550
    
    res_sorted = zip(*sorted_by_best)
    
    best_res = res_sorted[0]
    print len(mean_results), len(stdev_results), len(best_res)
    
    plot.plot(mean_results, label="average score")
    plot.plot(stdev_results, label="variance")
    plot.plot(best_res, label="score using all data")
#     plot.plot(reads_sorted_logscaled, label="reads, log scaled")
    plot.legend(loc='upper left')
#     plot.legend(loc="center right")
    plot.title("LC scores with 9-fold stability")
    plot.xlabel("LC miRNA")
    plot.ylabel("score")
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
plot.title("Low confidence scores, sorted")
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



roc_plot(test_all, test_all_annotations, "HC vs candidates and dead miRNA")
# roc_plot(test_miRNA, test_miRNA_annotations)
# roc_plot(test_candidates, test_candidates_annotations)
# roc_plot(test_dead, test_dead_annotations)


score_all = learner.score(test_all, test_all_annotations)
score_miRNA = learner.score(test_miRNA, test_miRNA_annotations)
score_candidates = learner.score(test_candidates, test_candidates_annotations)
score_dead = learner.score(test_dead, test_dead_annotations)


print "final scores:"
print "\t total score:\t\t", score_all, "\t", len(test_all_annotations)
print "\t HC miRNA score:\t", score_miRNA, "\t", len(test_miRNA_annotations)
print "\t candidate score:\t", score_candidates, "\t", len(test_candidates_annotations)
print "\t dead score:\t\t", score_dead, "\t", len(test_dead_annotations)
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


