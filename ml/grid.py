'''
Created on 19. feb. 2015

@author: hakon
'''

# used by pool.map
def one_grid(param):
    
    (all_data, all_annotations,test_fold_nr) = param
    
    c_values = [10.0 ** i for i in xrange(-3,4)]
    gamma_values = [10.0 ** i for i in xrange(-3,4)]
    
    
