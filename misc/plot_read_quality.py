'''
Created on 10. feb. 2015

@author: hakon
'''


from matplotlib import pyplot
import numpy
from scipy import stats

def plot(candidates, candidate_to_miRNAid, mirna_high_conf ):
    
    print
    print "plot read quality:"
    
    
    candidate_only = []
    mirna_high = []
    mirna_low = []
    minval = 1000000
    maxval = -1000000
    
    for c in candidates:
        hashval = c.chromosome + c.chromosome_direction + str(c.pos_5p_begin)
        param = c.quality
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
        
        
    print candidate_only
    print mirna_high
    print mirna_low
    
#     print minval, maxval
#     pyplot.hist(candidate_only)
    
    dens_cand = stats.kde.gaussian_kde(candidate_only)
    dens_high = stats.kde.gaussian_kde(mirna_high)
    dens_low = stats.kde.gaussian_kde(mirna_low)
#     
    x = numpy.arange( minval, maxval, .1)
    y = numpy.arange( minval, maxval, .1)
    z = numpy.arange( minval, maxval, .1)
    

    pyplot.plot(x, dens_cand(x))
    pyplot.plot(y, dens_high(y), "r")
    pyplot.plot(z, dens_low(z), "g")
    pyplot.savefig('quality.png')
    pyplot.xlabel("read quality")
    pyplot.ylabel("frequency")
    pyplot.show()