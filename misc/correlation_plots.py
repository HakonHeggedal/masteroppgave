'''
Created on 21. apr. 2015

@author: hakon
'''

from scipy.stats import spearmanr, pearsonr
from matplotlib import pyplot as plot
from matplotlib import cm
import numpy


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
    
    
def plot_correlation(vectors, is_pearson=True):
    
    results = [[0]*len(vectors) for _ in xrange(len(vectors))]
    
    correlation_function = pearsonr if is_pearson else spearmanr
    correlation_function_name = "pearson" if is_pearson else "spearman"
    
    features = zip(*vectors)

    for i, feature1 in enumerate(features):
        for j, feature2 in enumerate(features):
            
            score, _pval = correlation_function(feature1, feature2)
            
            results[i][j] = score
    
    
    
    plot_color_grid(results, None, None, correlation_function_name)
    
    
    

def plot_color_grid(arraylike, row_labels, column_labels, name):
    
    
    arraylike = numpy.array(arraylike)
        #  create 
    fig, ax = plot.subplots()
    heatmap = ax.pcolor(arraylike, cmap=cm.get_cmap("RdBu"), vmin=-1, vmax=1)
     
    # put the major ticks at the middle of each cell
    ax.set_xticks(numpy.arange(arraylike.shape[0])+0.5, minor=False)
    ax.set_yticks(numpy.arange(arraylike.shape[1])+0.5, minor=False)
    
    if row_labels:
        ax.set_xticklabels(row_labels, minor=False,  rotation="vertical")
    if column_labels:
        ax.set_yticklabels(column_labels, minor=False)
    

    
    fig.colorbar(heatmap)
    
    plot.show()
#     plot.savefig("figures/"+name+".png")
    
    


def _get_property_list(candidates, feat_names):
    
    all_features = []
    
    for name in feat_names:
        feature = [getattr(c, name) for c in candidates]
        all_features.append(feature)
        
    return all_features


# # 
# a = [[-0.5, 0.0], [0.5, 1.0]]
# l = [1,2]
#   
# plot_color_grid(a, l, l, "only_test_delete_pls")



