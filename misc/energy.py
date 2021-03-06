'''
Created on 6. feb. 2015

@author: hakon
'''

from matplotlib import pyplot
import numpy
from scipy import stats





def plot(candidates, candidate_to_miRNAid, mirna_high_conf ):
    
    print
    print "plot energy levels:"
    
    
    candidate_only = []
    mirna_high = []
    mirna_low = []
#     minval = 100000
#     maxval = -1000000
    for c in candidates:
        hashval = c.chromosome + c.chromosome_direction + str(c.pos_5p_begin)
        energy = c.hairpin_energy_10
#         if energy > maxval:
#             maxval = energy
#         elif energy < minval:
#             minval = energy
        
        if hashval in candidate_to_miRNAid:
            mi = candidate_to_miRNAid[hashval]
            
            if mi in mirna_high_conf:
                mirna_high.append(energy)
            else:
                mirna_low.append(energy)
        else:
            candidate_only.append(energy)
        
        
    print candidate_only
    print mirna_high
    print mirna_low
    
#     print minval, maxval
    
    
#     pyplot.hist(candidate_only)
    
    dens_cand = stats.kde.gaussian_kde(candidate_only)
    dens_high = stats.kde.gaussian_kde(mirna_high)
    dens_low = stats.kde.gaussian_kde(mirna_low)
#     
    x = numpy.arange( -100.0, 10.0, .1)
    y = numpy.arange( -100.0, 10.0, .1)
    z = numpy.arange( -100.0, 10.0, .1)

#  
    pyplot.plot(x, dens_cand(x))
    pyplot.plot(y, dens_high(y), "r")
    pyplot.plot(z, dens_low(z), "g")
    pyplot.savefig('energy.png')
    pyplot.xlabel("energy level")
    pyplot.ylabel("frequency")
    pyplot.show()
    
    pass

    assert False
    
