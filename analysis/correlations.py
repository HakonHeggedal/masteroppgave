'''
Created on 5. mai 2015

@author: hakon
'''
import pickle
import itertools
from misc.correlation_plots import plot_correlation
from ml.vectorize import FEATURE_NAMES







scaled_data = pickle.load( open("save_scaled_data.p", "rb"))
print "plotting correlation 123"
one_list_data = list(itertools.chain.from_iterable(scaled_data))

print len(scaled_data)
print len(one_list_data)

print
print len(scaled_data[0])
print
print len(one_list_data[0]), one_list_data[0]

print len(zip(*one_list_data))
print len(zip(*one_list_data)[0])


# feat_list = zip(*one_list_data)

feature_names = FEATURE_NAMES

plot_correlation(one_list_data, feature_names, True)
plot_correlation(one_list_data, feature_names, False)


# plot_correlation(one_list_data)
