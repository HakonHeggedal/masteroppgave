'''
Created on 21. apr. 2015

@author: hakon
'''

from scipy.stats import spearmanr, pearsonr
from matplotlib import pyplot as plot
import numpy


def _get_properties(candidates, feat_names):
    
    all_features = []
    
    for name in feat_names:
        feature = [getattr(c, name) for c in candidates]
        all_features.append(feature)
        
    return all_features



def plot_pearson_correlation(candidates, feat_names):
    
    
    features = _get_properties(candidates, feat_names)
    
    results = [[0]*len(feat_names) for _ in xrange(len(feat_names))]
    
    for i, feature1 in enumerate(features):
        for j, feature2 in enumerate(features):
            
            score, _pval = pearsonr(feature1, feature2)
            
            results[i][j] = abs(score)
    
    feat_names_space = [n.replace("_", " ") for n in feat_names]         
    
    plot_color_grid(results, feat_names_space, feat_names_space, "pearson")
    
    for r in results: print r



def plot_spearman_correlation(candidates, feat_names):
    
    features = _get_properties(candidates, feat_names)

    results = [[0]*len(feat_names) for _ in xrange(len(feat_names))]
    
    for i, feature1 in enumerate(features):
        for j, feature2 in enumerate(features):
            
            score, _pval = spearmanr(feature1, feature2)
            
            results[i][j] = abs(score)
    
    feat_names_space = [n.replace("_", " ") for n in feat_names]         
    
    plot_color_grid(results, feat_names_space, feat_names_space, "spearman")
            
    for r in results: print r
    
    

    
    

def plot_color_grid(arraylike, row_labels, column_labels, name):
    
    
    arraylike = numpy.array(arraylike)
        #  create 
    _fig, ax = plot.subplots()
    _heatmap = ax.pcolor(arraylike)
     
    # put the major ticks at the middle of each cell
    ax.set_xticks(numpy.arange(arraylike.shape[0])+0.5, minor=False)
    ax.set_yticks(numpy.arange(arraylike.shape[1])+0.5, minor=False)
     
    ax.set_xticklabels(row_labels, minor=False,  rotation="vertical")
    ax.set_yticklabels(column_labels, minor=False)
    
#     plot.show()
    plot.savefig("figures/"+name+".png")


