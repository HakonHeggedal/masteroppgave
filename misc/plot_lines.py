'''
Created on 22. apr. 2015

@author: hakon
'''


from matplotlib import pyplot


    

# # pyplot.plot( ks_p)
# # pyplot.plot( t_var)
# # pyplot.plot( t_var_p)
# 
# r = range(13, 18)
# 
# # pyplot.plot(r, ks_val, label="KS-test" )
# # pyplot.xlabel("max seq length")
# # pyplot.ylabel("KS test score")
# # pyplot.legend(loc='upper left')
# # pyplot.show()
# 
# 
# 
# pyplot.title("small sequences")
# pyplot.grid(True)
# pyplot.plot(r, t_var, label="Students T-test")
# pyplot.plot(r, t_sim, label="Welchs T-test")
# pyplot.xlabel("max sequence length")
# pyplot.ylabel("T-Test score")
# pyplot.legend(loc='upper left')
# 
# 
# # minorticks_on()
# pyplot.xticks(r)
# 
# pyplot.show()




def plot_ttest(t_student, t_welch, lengths, is_maxlen=True):
    
    maxormin = "max" if is_maxlen else "min"
    
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


def plot_kstest(ks_scores, lengths, is_maxlen=True):
    
    maxormin = "max" if is_maxlen else "min"
    
    pyplot.plot(lengths, ks_scores, label="KS-test" )
    pyplot.xlabel(maxormin+" sequence length")
    pyplot.ylabel("KS test score")
    pyplot.grid(True)
    pyplot.legend(loc='upper left')
    pyplot.xticks(lengths)
    pyplot.ylim(ymin=0)
    
    pyplot.show()
  
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

