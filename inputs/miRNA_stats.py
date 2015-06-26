'''
Created on 9. des. 2014

@author: hakon
'''

from matplotlib import pyplot
import mirbase
import numpy
from inputs import miRNA, special_types

# import miRNAs and mature seqs.

    
hairpin_file = "hairpin.fa"
mature_seq_file = "mature.fa"
high_conf_file = "high_conf_hairpin.fa"
other_types = "mirTrons_other.txt"

hsa_to_hairpin, other_to_hairpin = mirbase.read_miRNA_fasta(hairpin_file)
hsa_to_mature, other_to_mature = mirbase.read_miRNA_fasta(mature_seq_file)
miRNA_high_conf = miRNA.read_high_confidence(high_conf_file) 


before = len(hsa_to_hairpin)
before2 = len(hsa_to_mature)

special_types.remove_mirTrons(hsa_to_hairpin, other_types)
special_types.remove_mirTrons(hsa_to_mature, other_types)



assert before != len(hsa_to_hairpin)
assert before2 != len(hsa_to_mature)



 
# miRNA_species = mirbase.similar_hairpins(hsa_to_hairpin, other_to_hairpin)
harpinID_to_mature, harpinID_to_matseqs = mirbase.combine_hairpin_mature(hsa_to_hairpin, hsa_to_mature)


print len(set(harpinID_to_matseqs) )
print len(set(harpinID_to_matseqs) - miRNA_high_conf)
print len(miRNA_high_conf)


matures = [0] * 30
hairpins = [0] * 200

for v in hsa_to_hairpin.itervalues():
    hairpins[len(v)] += 1

for v in hsa_to_mature.itervalues():
    matures[len(v)] += 1
    
matures = numpy.array(matures)
hairpins = numpy.array(hairpins)



mature_hc = [0] * 30
mature_lc = [0] * 30


#   mature seqs
for mirnaID, seq_set in harpinID_to_matseqs.iteritems():

    if mirnaID in miRNA_high_conf:
        for seq in seq_set:
            mature_hc[len(seq)] += 1

    else:
        for seq in seq_set:
            mature_lc[len(seq)] += 1
            
            

hp_hc = [0] * 161
hp_lc = [0] * 161
print

#     if mirnaID in miRNA_high_conf:            
for mirnaID, seq_set in harpinID_to_matseqs.iteritems():
    hairpin =  harpinID_to_mature[mirnaID]
    hairpin =  hsa_to_hairpin[mirnaID]
        
#     best_start = 1000
#     best_end = 0
    if len(seq_set) >= 2:
        
        best_start = min(map(lambda x: hairpin.find(x), seq_set))
        best_end = max(map(lambda x: hairpin.find(x)+len(x), seq_set))
    
    else:
        continue
        
        
#         for seq in seq_set:
#             assert seq in hairpin
#             
#             start = hairpin.find(seq)
#             end = start + len(seq)
#             
#             best_start = start if start < best_start else best_start
#             best_end = end if end > best_end else best_end

    hp_len = best_end - best_start 
#     print "got full hp", hp_len, best_start, best_end, len(hairpin), mirnaID
#     print "\t", seq_set, hairpin
    if hp_len > 33:
        if mirnaID in miRNA_high_conf:
            hp_hc[hp_len] += 1
        else:
            hp_lc[hp_len] += 1
            
#         assert 0
#         for seq in seq_set:
#             mature_hc[len(seq)] += 1

#     else:
#         for seq in seq_set:
#             mature_lc[len(seq)] += 1
    



# print [x for x in enumerate(matures)]
# print numpy.mean(matures)
# print numpy.median(matures)
# print "min ",next(i for i,x in enumerate(matures) if x != 0 )
# print "max ",next(len(matures)-i for i,x in enumerate(matures[::-1]) if x != 0 )
# 
# for i,x in enumerate(matures):
#     if x != 0:
#         print "min",i
#         break
#  
# for i in xrange(len(matures)-1, 0, -1):
#     if matures[i] != 0:
#         print "max", i
#         break
#  
# print matures[15], matures[16], matures[22], matures[28], matures[29]
# 
# print enumerate(hairpins)
# 
# print "min ",next(i for i,x in enumerate(hairpins) if x != 0 )
# print "max ",next(len(hairpins)-i-1 for i,x in enumerate(hairpins[::-1]) if x != 0 )
# 
# for i,x in enumerate(hairpins):
#     if x != 0:
#         print "min",i, hairpins[i]
#         break
#  
# for i in xrange(len(hairpins)-1, 0, -1):
#     if hairpins[i] != 0:
#         print "max", i, hairpins[i]
#         break
#     
# print hairpins[180], hairpins[181], hairpins[40], hairpins[41]


# derp = range(0, len(hairpins))
# derp = [float(x)-0.5 for x in derp]
# 
# pyplot.bar( derp, hairpins,width=1.0)
# pyplot.title("Length distribution for miRNA hairpins")
# pyplot.xlabel("length")
# pyplot.ylabel("frequency")
# pyplot.xticks(range(40,200,10))
# pyplot.grid(True)
# # pyplot.axis([15,30,0,1200])
# pyplot.savefig("mirna_hairpin_distr.pdf")
# pyplot.savefig("mirna_hairpin_distr.png")
# pyplot.show()



 
positions_all_mat = [float(x)-0.3 for x in range(len(matures))]
positions_hc = [float(x)-0.3 for x in range(len(mature_hc))]
positions_lc = [float(x)-0.3 for x in range(len(mature_lc))]
 
 
 
# pyplot.bar(derp,  matures, width=0.6, color="red")
pyplot.bar(positions_hc,  mature_hc, width=0.6, color="red", label="HC mature seqs")
# pyplot.bar(positions_lc,  mature_lc, bottom=mature_hc, width=0.6, color="green", label="LC mature seqs")
pyplot.title("Length distribution for miRNA mature sequences")
pyplot.xlabel("length")
pyplot.ylabel("frequency")
pyplot.grid(True)
pyplot.minorticks_on()
# pyplot.axis([16,28,0,1200])
# pyplot.xticks(range(16,28))
pyplot.locator_params(axis='x', nbins=20)
pyplot.yscale('symlog', nonposy='clip', basey=2)
pyplot.legend(loc='upper right')

pyplot.savefig("mirna_mature_distr.pdf")
pyplot.savefig("mirna_mature_distr.png")
pyplot.show()
pyplot.close()





print max(hp_hc)
print max(hp_lc)

# positions_all_mat = [float(x)-0.3 for x in range(len(matures))]
pos_hc = [float(x)-0.5 for x in range(len(hp_hc))]
pos_lc = [float(x)-0.5 for x in range(len(hp_lc))]



# pyplot.bar(derp,  matures, width=0.6, color="red")
pyplot.bar(pos_hc,  hp_hc, width=1.0, color="red", label="HC hairpins")
# pyplot.bar(pos_lc,  hp_lc, width=1.0, color="green", bottom=hp_hc, label="LC hairpins")
pyplot.title("Length distribution for miRNA hairpins")
pyplot.xlabel("length")
pyplot.ylabel("frequency")
pyplot.grid(True)
pyplot.minorticks_on()
# pyplot.axis([16,28,0,1200])
# pyplot.xticks(range(16,28))
# pyplot.locator_params(axis='x', nbins=20)
# pyplot.yscale('symlog', nonposy='clip', basey=2)
pyplot.legend(loc='upper right')

pyplot.savefig("mirna_hairpin_distr.pdf")
pyplot.savefig("mirna_hairpin_distr.png")
pyplot.show()




# pyplot.plot(hairpins)
# pyplot.show()