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

# print 12334567
# print len(hsa_to_hairpin), before
# print len(hsa_to_mature), before2

assert before != len(hsa_to_hairpin)
assert before2 != len(hsa_to_mature)

# print miRNA_high_conf
# print miRNA_high_conf.issubset(hsa_to_hairpin.keys())


 
# miRNA_species = mirbase.similar_hairpins(hsa_to_hairpin, other_to_hairpin)
harpinID_to_mature, harpinID_to_matseqs = mirbase.combine_hairpin_mature(hsa_to_hairpin, hsa_to_mature)


print len(set(harpinID_to_matseqs) )
print len(set(harpinID_to_matseqs) - miRNA_high_conf)
print len(miRNA_high_conf)
# print harpinID_to_matseqs
print
# print hairpinID_to_mature
# assert False

# calculate length distribution.


# for k,v in hsa_to_hairpin.iteritems():
#     print k,len(v)
# 
# for k,v in hsa_to_mature.iteritems():
#     print k,len(v)
#     

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

hcs = 0

for mirnaID, seq_set in harpinID_to_matseqs.iteritems():

    if mirnaID in miRNA_high_conf:
        for seq in seq_set:
            mature_hc[len(seq)] += 1
            hcs += 1
    else:
        for seq in seq_set:
            mature_lc[len(seq)] += 1
print "woot?"
print hcs 
print len(miRNA_high_conf)


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





derp = range(len(matures))
derp = [float(x)-0.3 for x in derp]

positions_hc = [float(x)-0.0 for x in range(len(mature_hc))]
positions_lc = [float(x)-0.3 for x in range(len(mature_lc))]


print matures
print derp

# pyplot.bar(derp,  matures, width=0.6, color="red")
pyplot.bar(positions_lc,  mature_lc, width=0.3, color="green")
pyplot.bar(positions_hc,  mature_hc, width=0.3, color="red")
pyplot.title("Length distribution for miRNA mature sequences")
pyplot.xlabel("length")
pyplot.ylabel("frequency")
pyplot.grid(True)
pyplot.minorticks_on()
# pyplot.axis([16,28,0,1200])
# pyplot.xticks(range(16,28))
pyplot.locator_params(axis='x', nbins=20)
pyplot.yscale('symlog', nonposy='clip', basey=2)
pyplot.savefig("mirna_mature_distr.pdf")
pyplot.savefig("mirna_mature_distr.png")
pyplot.show()




# pyplot.plot(hairpins)
# pyplot.show()