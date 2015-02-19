'''
Created on 19. feb. 2015

@author: hakon
'''
import numpy


def compute_stem_start(candidates, mirna, hc_mirna):
    
    for c in candidates:
        mi = 0
        hashval = c.chromosome + c.chromosome_direction + str(c.pos_5p_begin)
        if hashval in mirna:
            mi = mirna[hashval]
            if mi not in hc_mirna: continue
        else:
            continue
        
        precursor = c.hairpin_fold_40
        entropies = c.bitpair_entropy_dict
        bitpairs = c.bitpair_probabilities
        
        j_five, j_three = _stem_junction_pos(precursor, entropies)
        
        print j_five, j_three
        print len(c.mapped_sequences)

        
#         print precursor
        
        print "found high confidence stuff:", mi
        for i in xrange(len(precursor)):
            
            print i,
            print precursor[i],
            print c.hairpin_padded_40[i],
            print entropies[i] if i in entropies else 0,
            print "  ", numpy.mean(bitpairs[i].values()) if i in bitpairs else 0.0,
            print "  ", bitpairs[i].items() if i in bitpairs else "",
            
            if i == j_five or i == j_three:
                print "<------ junction",
                
            if i == 40 or i == len(precursor)-40:
                print "<<<<<<< hairpin",
#             if i == 40 + c.pos_3p_end - c.pos_5p_begin:
#                 print 123,
                
            if i == 40 + c.pos_5p_end - c.pos_5p_begin:
                print "<-<-<-<-< bulge start",
                
            if i == 40 + c.pos_3p_begin - c.pos_5p_begin:
                print "<-<-<-<-< bulge end",
            
            print 
        
        print len(precursor), len(precursor)-80
        assert False

def _stem_junction_pos(precursor, entropy_dict):
        """get the coordinates of the stem-juction.
        That is the region of unpaired bases at the 5p and 3p end.
        This region tends to have a lot higher entropy than the rest,
        so we are gonna use that, since the forgi annonation is not reliable"""


        entropy = [0]*len(precursor)
        for i in range(len(precursor)):
            pos_entropy = entropy_dict[i] if i in entropy_dict else 0
            entropy[i] = pos_entropy

        #smooth
        smooth = 7
        entropy_smooth = [0]*len(entropy)
        
        entropy = entropy + [entropy[-1]] * smooth
        for i in range(len(entropy_smooth)):
            entropy_smooth[i] = numpy.mean(entropy[i:i+smooth+1])

        #(discrete) derivate
        derivative = []
        for i in range(len(entropy_smooth) - 1):
            derivative.append(entropy_smooth[i+1] - entropy_smooth[i])

        #2nd derivative
        derivative2 = []
        for i in range(len(derivative) - 1):
            derivative2.append(derivative[i+1] - derivative[i])

        j_five = 0
        j_three = 0
        skip_5 = False
        skip_3 = False
        for i in range(40, 0, -1):
            if derivative2[i]-derivative[i] > .09 and not skip_5:
                j_five = i
                skip_5 = True
            
            if derivative2[-i]-derivative[-i] > .09 and not skip_3:
                j_three = -i
                skip_3 = True
                
            if skip_5 and skip_3:
                break
            
        j_three = len(precursor) + j_three
        
        return j_five, j_three
    
    
