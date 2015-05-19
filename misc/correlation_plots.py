'''
Created on 21. apr. 2015

@author: hakon
'''

from scipy.stats import spearmanr, pearsonr
from matplotlib import pyplot as plot
from matplotlib import cm
import numpy
import itertools
import pickle
import math
import matplotlib


def plot_pearson_correlation(candidates, feat_names):
    
    
    features = _get_property_list(candidates, feat_names)
    
    results = [[0]*len(feat_names) for _ in xrange(len(feat_names))]
    
    for i, feature1 in enumerate(features):
        for j, feature2 in enumerate(features):
            
            score, _pval = pearsonr(feature1, feature2)
            
            results[i][j] = score
    
    feat_names_space = [n.replace("_", " ") for n in feat_names]         
    
    plot_color_grid(results, feat_names_space, feat_names_space, "pearson")
    
    for r in results: print r



def plot_spearman_correlation(candidates, feat_names):
    
    features = _get_property_list(candidates, feat_names)

    results = [[0]*len(feat_names) for _ in xrange(len(feat_names))]
    
    for i, feature1 in enumerate(features):
        for j, feature2 in enumerate(features):
            
            score, _pval = spearmanr(feature1, feature2)
            
            results[i][j] = score
    
    feat_names_space = [n.replace("_", " ") for n in feat_names]         
    
    plot_color_grid(results, feat_names_space, feat_names_space, "spearman")
            
    for r in results: print r
    
    
def plot_correlation(vectors, feat_names=None, is_pearson=True):
    
    correlation_function = pearsonr if is_pearson else spearmanr
    correlation_function_name = "pearson" if is_pearson else "spearman"
    
    features = zip(*vectors)
    
    results = [[0]*len(features) for _ in xrange(len(features))]
    
    for i, feature1 in enumerate(features):
        for j, feature2 in enumerate(features):
            
            score, _pval = correlation_function(feature1, feature2)
            
            if math.isnan(score):
                print 
                print  "test"
                print i, j
                print feature1[:10]
                print feature2[:10]
                print score
#                 assert 0
            
            results[i][j] = score if not math.isnan(score) else 0
        print
    
    
    print len(features)
    print len(results)
    plot_color_grid(results, feat_names, feat_names, correlation_function_name)
    
    
    

def plot_color_grid(arraylike, row_labels, column_labels, name, vmin=-1, vmax=1):

    name_spaces = name.replace("_", " ")
    
    arraylike = numpy.array(arraylike)

    fig, ax = plot.subplots()
    
    plot.title("feature correlation " + name_spaces)
    heatmap = ax.pcolor(arraylike, cmap=cm.get_cmap("RdBu"), vmin=vmin, vmax=vmax)

    # put the major ticks at the middle of each cell
    ax.set_xticks(numpy.arange(arraylike.shape[0])+0.5, minor=False)
    ax.set_yticks(numpy.arange(arraylike.shape[1])+0.5, minor=False)
#     ax.set_xticks(numpy.arange(arraylike.shape[0]), minor=False)
#     ax.set_yticks(numpy.arange(arraylike.shape[1]), minor=False)
 
    if row_labels:
        ax.set_xticklabels(row_labels, minor=False,  rotation="vertical")
    if column_labels:
        ax.set_yticklabels(column_labels, minor=False)
     
    # figure does not know its own limits, so setting them here:
    heatmap.axes.set_ylim(0, len(arraylike))
    heatmap.axes.set_xlim(0, len(arraylike))
    
    fig.colorbar(heatmap)
    
    
    plot.show()
    plot.savefig("figures/"+name+".png")
    
    


def _get_property_list(candidates, feat_names):
    
    all_features = []
    
    for name in feat_names:
        feature = [getattr(c, name) for c in candidates]
        all_features.append(feature)
        
    return all_features

print matplotlib.__version__

