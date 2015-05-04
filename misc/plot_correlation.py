'''
Created on 22. apr. 2015

@author: hakon
'''


from matplotlib import pyplot
import numpy


def plot_ttest(t_student, t_welch, lengths, is_maxlen=True):
    
    maxormin = "max" if is_maxlen else "min"
    
    smin = str(lengths[0])
    smax = str(lengths[-1])
    name = "figures/T-test" + smin + "_" + smax + maxormin + ".png"
    
    
    lengths = numpy.array(lengths)
    t_student = numpy.array(t_student)
    t_welch = numpy.array(t_welch)
    
    
    print len(lengths)
    print len(t_student)
    print len(t_welch)
    print lengths
    print
    print t_student
    print 
    print t_welch
    print
    
    pyplot.title("small sequences")
    pyplot.grid(True)
    pyplot.plot(lengths, t_student, label="Students T-test")
    pyplot.plot(lengths, t_welch, label="Welchs T-test")
    pyplot.xlabel(maxormin+" sequence length")
    pyplot.ylabel("T-Test score")
    pyplot.legend(loc='upper left')
    pyplot.ylim(ymin=0)
    
    # minorticks_on()
    pyplot.xticks(lengths)
    
    pyplot.show()
    pyplot.savefig(name)


def plot_kstest(ks_scores, lengths, is_maxlen=True):
    
    maxormin = "max" if is_maxlen else "min"
    
    smin = str(lengths[0])
    smax = str(lengths[-1])
    name = "figures/KS-test_" + smin + "_" + smax + maxormin + ".png"
    print len(lengths)
    print len(ks_scores)
    print lengths
    print
    print ks_scores
    
    lengths = numpy.array(lengths)
    ks_scores = numpy.array(ks_scores)
    
    pyplot.plot(lengths, ks_scores, label="KS-test" )
    pyplot.xlabel(maxormin+" sequence length")
    pyplot.ylabel("KS test score")
    pyplot.grid(True)
    pyplot.legend(loc='upper left')
    pyplot.xticks(lengths)
    pyplot.ylim(ymin=0)
    
    pyplot.show()
    pyplot.savefig(name)
    
    

# l = [8, 9, 10, 11, 12, 13, 14, 15, 16]
# 
# k = (0.32195773081201329, 0.39503893214682978, 0.48829810901001114, 0.60422691879866519, 0.62631813125695213, 0.66547274749721919, 0.66760845383759737, 0.6690989988876529, 0.65937708565072306)
# 
# plot_kstest(k, l, False)

# r = range(13, 18)
# 
# ks_val = (0.11948832035595103, 0.14335928809788656, 0.15668520578420464, 0.18002224694104552, 0.19942157953281425)
# t_var = (5.3311434702765537, 6.3255867814905367, 7.0479326684247665, 7.793978207440845, 8.6404347703693976)
# 
# t_sim = (4.2610009529504342, 5.0022385866068779, 5.5768320086580498, 6.0886191954143474, 6.6614277392142842)

# if __name__ == "__main__":
#     
#     plot_ttest(t_var, t_sim, r, True)
#     plot_kstest(ks_val, r, True)

# ks_p = (0.0016673991060482459, 7.3949393401517884e-05, 1.0154784126949554e-05, 2.0517035415567726e-07, 5.3051189184932447e-09)
# t_var_p = (1.0963486082783808e-07, 3.1602557009512655e-10, 2.558742607941345e-12, 1.0783462552345224e-14, 1.1950072524116202e-17)
# t_sim_p(2.6258094513164586e-05, 9.0472306098188822e-07, 4.9548024599765484e-08, 3.0557184883199246e-09, 1.0918613635302362e-10)
# 

