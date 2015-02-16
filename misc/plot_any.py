'''
Created on 16. feb. 2015

@author: hakon
'''


from matplotlib import pyplot
import numpy
import math
from scipy import stats



def plot(candidates, candidate_to_miRNAid, mirna_high_conf, name, isLog=False):
    """ plot any candidate feature given its name and candidate/miRNA classes """
    print
    print "plot" + name
    
    plot_name = name.replace("_", " ")
    outfile = name + ".png"
    
    candidate_only = []
    mirna_high = []
    mirna_low = []
    
    minval = 1000000
    maxval = -1000000
    for c in candidates:
        hashval = c.chromosome + c.chromosome_direction + str(c.pos_5p_begin)
        param = getattr(c, name) #.overhang_level_outer_10
        if isLog and param != 0:
            param = math.log(param)
        
        if param > maxval:
            maxval = param
        elif param < minval:
            minval = param
        
        if hashval in candidate_to_miRNAid:
            mi = candidate_to_miRNAid[hashval]
            
            if mi in mirna_high_conf:
                mirna_high.append(param)
            else:
                mirna_low.append(param)
        else:
            candidate_only.append(param)
    
    dens_cand = stats.kde.gaussian_kde(candidate_only)
    dens_high = stats.kde.gaussian_kde(mirna_high)
    dens_low = stats.kde.gaussian_kde(mirna_low)
  
    x = numpy.arange( minval, maxval, .1)
    y = numpy.arange( minval, maxval, .1)
    z = numpy.arange( minval, maxval, .1)
    
    ks_val, p_2s =  stats.ks_2samp(mirna_high, mirna_low)
    print "Kolmogorov-Smirnov test: miRNA high/low conf:"
    print "\tKS:\t\t\t", ks_val
    print "\ttwo sided p val:\t", p_2s

    pyplot.plot(x, dens_cand(x), "b")
    pyplot.plot(y, dens_high(y), "r")
    pyplot.plot(z, dens_low(z), "g")
    pyplot.savefig(outfile)
    pyplot.xlabel(plot_name)
    pyplot.ylabel("frequency")
    pyplot.show()
    


#     print candidate_only
#     print mirna_high
#     print mirna_low
    assert False
    
