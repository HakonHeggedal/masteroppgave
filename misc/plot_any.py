'''
Created on 16. feb. 2015

@author: hakon
'''


from matplotlib import pyplot
import numpy
import math
from scipy import stats



def plot(candidates, candidate_to_miRNAid, candidate_to_dead, mirna_high_conf, name, isLog=False):
    """ plot any candidate feature given its name and candidate/miRNA classes """
    print
    print "plot " + name
    
    name = name.replace(" ", "_")
    
    plot_name = name.replace("_", " ")
    outfile = name + ".png"
    
    candidate_only = []
    mirna_high = []
    mirna_low = []
    dead = []
    
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
            
#         if param > -10: continue
            
        if hashval in candidate_to_dead:
            dead.append(param)
        
        elif hashval in candidate_to_miRNAid:
            mi = candidate_to_miRNAid[hashval]
            
            if mi in mirna_high_conf:
                mirna_high.append(param)
            else:
                mirna_low.append(param)
        else:
            candidate_only.append(param)
    
#     print len(dead), dead
    
    dens_cand = stats.kde.gaussian_kde(candidate_only)
    dens_high = stats.kde.gaussian_kde(mirna_high)
    dens_low = stats.kde.gaussian_kde(mirna_low)
    dens_dead = stats.kde.gaussian_kde(dead)
  
    x = numpy.arange( minval, maxval, .1)
    y = numpy.arange( minval, maxval, .1)
    z = numpy.arange( minval, maxval, .1)
    ae = numpy.arange( minval, maxval, .1)
    
    ks_val, p_2s =  stats.ks_2samp(mirna_high, mirna_low)
    print "\tKolmogorov-Smirnov test: miRNA high/low conf:"
    print "\t\tKS: ", ks_val, "two sided p val: ", p_2s
    
    tval, prob = stats.ttest_ind(mirna_high, mirna_low, equal_var=True)
    print "\tStudents T-test (same variance): miRNA high/low conf:"
    print "\t\ttval: ",tval, "probability: ", prob 
    
    tval, prob = stats.ttest_ind(mirna_high, mirna_low, equal_var=False)
    print "\tWelchs T-test: miRNA high/low conf:"
    print "\t\ttval: ",tval, "probability: ", prob 

    pyplot.plot(x, dens_cand(x), "b")
    pyplot.plot(y, dens_high(y), "r")
    pyplot.plot(z, dens_low(z), "g")
    pyplot.plot(ae, dens_dead(ae), "m")
    pyplot.savefig(outfile)
    pyplot.xlabel(plot_name)
    pyplot.ylabel("frequency")
#     pyplot.show()
    pyplot.close()
    


#     print candidate_only
#     print mirna_high
#     print mirna_low
#     assert False
# FEATURES = ["hairpin_energy", "hairpin_energy_10", "hairpin_energy_40",
#                 "entropy_nucleotides", "entropy_structure", "heterogenity_5_begin",
#                 "heterogenity_5_end", "heterogenity_3_begin", "heterogenity_3_end",
#                 "quality", "bindings_max_10", "overhang_level_outer_10",
#                 "overhang_outer_10", "overhang_level_inner_10", "overhang_inner_10",
#                 "small_subs", "small_subs_5p", "small_subs_3p"]
# log_scaled = [False]*15 + [True]*3
# print len(log_scaled), log_scaled
# print len(FEATURES)
# 
# for feat_name, logs in zip(FEATURES, log_scaled):
#     print 123
