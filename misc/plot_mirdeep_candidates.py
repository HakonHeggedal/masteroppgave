'''
Created on 16. feb. 2015

@author: hakon
'''


from matplotlib import pyplot
import numpy
import math
from scipy import stats



def plot_candidate_results(candidates, candidate_to_miRNAid, candidate_to_dead, mirna_high_conf,
                           mirdeep_candidates, hp50, hp99, name, isLog=False):
    """ plot any candidate feature given its name and candidate/miRNA classes """
    
    print
    print "plot " + name
    
    name = name.replace(" ", "_")
    log_text = " log scaled" if isLog else ""
    
    plot_name = name.replace("_", " ")
    outfile = "figures/" + name + ".png"
    
    candidate_only = []
    mirna_high = []
    mirna_low = []
    dead = []
    mirdeep_c = []
    classed_50 = []
    classed_99 = []
    
    minval = 1000000
    maxval = -1000000
    for c in candidates:
        hashval = c.chromosome + c.chromosome_direction + str(c.hairpin_start)
        param = getattr(c, name) #.overhang_level_outer_10
        assert param is not None
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
            if hashval in mirdeep_candidates:
                mirdeep_c.append(param)
            if hashval in hp50: 
                classed_50.append(param)
            if hashval in hp99:
                classed_99.append(param)
            else:
                candidate_only.append(param)
            
    print "found by mirDeep:", len(mirdeep_c), mirdeep_c
    
#     print len(dead), dead
#     print filter(lambda x: x!= 1.0, candidate_only)
    dens_cand = stats.kde.gaussian_kde(candidate_only)
    dens_high = stats.kde.gaussian_kde(mirna_high)
    dens_low = stats.kde.gaussian_kde(mirna_low)
    dens_dead = stats.kde.gaussian_kde(dead)
    if len (mirdeep_c):
        dens_mirdeep = stats.kde.gaussian_kde(mirdeep_c)
    dens_50 = stats.kde.gaussian_kde(classed_50)
    dens_99 = stats.kde.gaussian_kde(classed_99)
    
  
    x = numpy.arange( minval, maxval, .1)
    y = numpy.arange( minval, maxval, .1)
    z = numpy.arange( minval, maxval, .1)
    ae = numpy.arange( minval, maxval, .1)
    oe = numpy.arange( minval, maxval, .1)
    c_50 = numpy.arange( minval, maxval, .1)
    c_99 = numpy.arange( minval, maxval, .1)
    
    ks_val, p_2s =  stats.ks_2samp(mirna_high, mirna_low)
    print "\tKolmogorov-Smirnov test: miRNA high/low conf:"
    print "\t\tKS: ", ks_val, "two sided p val: ", p_2s
    
    tval_same_var, prob_same = stats.ttest_ind(mirna_high, mirna_low, equal_var=True)
    print "\tStudents T-test (same variance): miRNA high/low conf:"
    print "\t\ttval: ",tval_same_var , "probability: ", prob_same 
    
    tval_diff_var, prob_diff = stats.ttest_ind(mirna_high, mirna_low, equal_var=False)
    print "\tWelchs T-test: miRNA high/low conf:"
    print "\t\ttval: ",tval_diff_var , "probability: ", prob_diff 
    
    
    fig, ax = pyplot.subplots()
    
    ax.plot(x, dens_cand(x), label="Candidates")
    ax.plot(y, dens_high(y), label="miRNA HC")
    ax.plot(z, dens_low(z), label="miRNA LC")
    ax.plot(ae, dens_dead(ae), label="dead miRNA")
    if len (mirdeep_c):
        ax.plot(oe, dens_mirdeep(oe), label="mirDeep2 new")
     
    ax.plot(c_50, dens_50(c_50), label="0.5 score")
    ax.plot(c_99, dens_99(c_99), label="0.99 score")
         
    print "doing magic"
    ax.set_xlabel(plot_name + log_text)
    ax.set_ylabel("frequency")
    print "doing more magic"
    box = ax.get_position()
    ax.set_position([box.x0-0.01, box.y0, box.width * 0.74, box.height])
    ax.legend(loc='upper left', bbox_to_anchor = (1.02, 1.0)) # make legend using this one
    
    print "should print now"
    pyplot.savefig(outfile)
    pyplot.show()
    pyplot.close()
#     pyplot.plot(x, dens_cand(x), "b", label="Candidates")
#     pyplot.plot(y, dens_high(y), "r", label="miRNA high confidence")
#     pyplot.plot(z, dens_low(z), "g", label="miRNA low confidence")
#     pyplot.plot(ae, dens_dead(ae), "m", label="dead miRNA")
#     if len (mirdeep_c):
#         pyplot.plot(oe, dens_mirdeep(oe), "k", label="dead miRNA")
#     
#     pyplot.plot(c_50, dens_dead(c_50), label="dead miRNA")
#     pyplot.plot(c_99, dens_dead(c_99), label="dead miRNA")
#         
#     pyplot.legend(loc='upper left') # make legend using this one
#     pyplot.xlabel(plot_name + log_text)
#     pyplot.ylabel("frequency")
#     pyplot.savefig(outfile)
#     pyplot.show()
#     pyplot.close()
    
    return (ks_val, p_2s, tval_same_var, prob_same, tval_diff_var, prob_diff)


fig, ax = pyplot.subplots()
 
ax.plot([4,2,3], label="Candidates")
ax.plot([1,0,3], label="miRNA HC")
ax.plot([1,1,3], label="miRNA LC")
ax.plot([1,6,3], label="dead miRNA")
ax.plot([1,2,2], label="mirDeep2 new")
ax.plot([1,2,5], label="0.99 scored")
ax.plot([3,3,3], label="classifier 0.99")
ax.set_ylabel("frequency")
ax.set_xlabel("frequency")
# pyplot.ax
# pyplot.subplots_adjust(right=0.8)
 
box = ax.get_position()
ax.set_position([box.x0-0.01, box.y0, box.width * 0.74, box.height])
 
ax.legend(loc='upper left', bbox_to_anchor = (1.02, 1.0))
pyplot.title("some random title heres he")
pyplot.show()
#     pyplot.close()









