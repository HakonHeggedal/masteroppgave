'''
Created on 9. des. 2014

@author: hakon
'''

from matplotlib import pyplot
import mirbase
import numpy

# import miRNAs and mature seqs.

    
hairpin_file = "hairpin.fa"
mature_seq_file = "mature.fa"

hsa_to_hairpin, other_to_hairpin = mirbase.read_miRNA_fasta(hairpin_file)
hsa_to_mature, other_to_mature = mirbase.read_miRNA_fasta(mature_seq_file)


miRNA_species = mirbase.similar_hairpins(hsa_to_hairpin, other_to_hairpin)
hairpinID_to_mature = mirbase.combine_hairpin_mature(hsa_to_hairpin, hsa_to_mature)


# calculate length distribution.


for k,v in hsa_to_hairpin.iteritems():
    print k,len(v)
    

for k,v in hsa_to_mature.iteritems():
    print k,len(v)
    
    
matures = [0] * 30
hairpins = [0] * 200

for v in hsa_to_hairpin.itervalues():
    hairpins[len(v)] += 1

for v in hsa_to_mature.itervalues():
    matures[len(v)] += 1
    
matures = numpy.array(matures)
hairpins = numpy.array(hairpins)

print matures
print numpy.mean(matures)
print numpy.median(matures)
print "min ",next(i for i,x in enumerate(matures) if x != 0 )
print "max ",next(len(matures)-i for i,x in enumerate(matures[::-1]) if x != 0 )

for i,x in enumerate(matures):
    if x != 0:
        print "min",i
        break
 
for i in xrange(len(matures)-1, 0, -1):
    if matures[i] != 0:
        print "max", i
        break



print hairpins

print "min ",next(i for i,x in enumerate(hairpins) if x != 0 )
print "max ",next(len(hairpins)-i-1 for i,x in enumerate(hairpins[::-1]) if x != 0 )

for i,x in enumerate(hairpins):
    if x != 0:
        print "min",i, hairpins[i]
        break
 
for i in xrange(len(hairpins)-1, 0, -1):
    if hairpins[i] != 0:
        print "max", i, hairpins[i]
        break
    
print hairpins[180], hairpins[181], hairpins[40], hairpins[41]


derp = range(0, len(hairpins))
derp = [float(x)-0.5 for x in derp]

pyplot.bar( derp, hairpins,width=1.0)
pyplot.title("Length distribution for miRNA hairpins")
pyplot.xlabel("length")
pyplot.ylabel("frequency")
pyplot.grid(True)
# pyplot.axis([15,30,0,1200])
pyplot.savefig("mirna_hairpin_distr.pdf")
pyplot.savefig("mirna_hairpin_distr.png")
pyplot.show()





derp = range(len(matures))
derp = [float(x)-0.3 for x in derp]

pyplot.bar( derp, matures, width=0.6)
pyplot.title("Length distribution for mature miRNAs")
pyplot.xlabel("length")
pyplot.ylabel("frequency")
pyplot.grid(True)
# pyplot.axis([16,28,0,1200])
# pyplot.xticks(range(16,28))
pyplot.locator_params(axis='x', nbins=20)


pyplot.savefig("mirna_mature_distr.pdf")
pyplot.savefig("mirna_mature_distr.png")
pyplot.show()
# pyplot.





# pyplot.plot(hairpins)
# pyplot.show()