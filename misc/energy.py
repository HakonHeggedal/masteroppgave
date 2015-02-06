'''
Created on 6. feb. 2015

@author: hakon
'''

from matplotlib import pyplot
import numpy
from scipy import stats

def plot(candidates, candidate_to_miRNAid, mirna_high_conf):
    
    print
    print "plot energy levels:"
    
    
    candidate_only = []
    mirna_high = []
    mirna_low = []
    
    for c in candidates:
        hashval = c.chromosome + c.chromosome_direction + str(c.pos_5p_begin)
        energy = c.hairpin_energy_10
        
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
    
    
#     pyplot.hist(candidate_only)
    
    dens_cand = stats.kde.gaussian_kde(candidate_only)
    dens_high = stats.kde.gaussian_kde(mirna_high)
    dens_low = stats.kde.gaussian_kde(mirna_low)
#     
    x = numpy.arange(-100., 0, .1)
    y = numpy.arange(-100., 0, .1)
    z = numpy.arange(-100., 0, .1)

#  
    pyplot.plot(x, dens_cand(x))
    pyplot.plot(y, dens_high(y), "r")
    pyplot.plot(z, dens_low(z), "g")
    pyplot.savefig('energy.png')
    pyplot.show()


#     
#     density = stats.kde.gaussian_kde(data1)
#     
#     x = numpy.arange(0., 200, .1)
#     
#     # pyplot.show()
#     # pyplot.hold()
#     
#     data2 = [94, 2, 54, 4, 3, 13, 21, 4, 2, 2, 153, 4, 0, 2, 4, 0, 107, 0, 47, 7, 12, 0, 123, 2, 12, 1, 1, 16, 20, 2, 75, 10, 71, 53, 17, 12, 45, 0, 1, 3, 54, 18, 26, 78, 1, 3, 2, 54, 12, 21, 24, 1, 11, 22, 10, 6, 24, 13, 1, 3, 20, 111, 6, 2, 2, 42, 115, 131, 20, 1, 79, 119, 124, 1, 177, 130, 124, 9, 208, 14, 9, 8, 2, 2, 157, 25, 4, 48, 6, 1, 5, 46, 17, 6, 153, 116]
#     # lc_species = [0, 6, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 4, 0, 0, 0, 0, 3, 8, 0, 0, 5, 0, 0, 0, 7, 2, 0, 0, 0, 0, 0, 0, 0, 0, 8, 3, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 2, 13, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 3, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0]
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     density2 = stats.kde.gaussian_kde(data2)
#     y = numpy.arange(0., 200, .1)
#     
#     pyplot.plot(x, density(x))
#     pyplot.plot(y, density2(y), "r")
#     
#     pyplot.show()
    
    

    assert False
    
