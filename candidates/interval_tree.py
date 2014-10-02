'''
Created on 1. okt. 2014

@author: hakon
'''
import intervaltree
from intervaltree import bio



GENOMENR_TO_FASTAFILE = {}
GENOMENR_TO_FASTAFILE["NC_000001.10"] = "chr1.fa"
GENOMENR_TO_FASTAFILE["NC_000002.11"] = "chr2.fa"
GENOMENR_TO_FASTAFILE["NC_000003.11"] = "chr3.fa"
GENOMENR_TO_FASTAFILE["NC_000004.11"] = "chr4.fa"
GENOMENR_TO_FASTAFILE["NC_000005.9"] = "chr5.fa"
GENOMENR_TO_FASTAFILE["NC_000006.11"] = "chr6.fa"
GENOMENR_TO_FASTAFILE["NC_000007.13"] = "chr7.fa"
GENOMENR_TO_FASTAFILE["NC_000008.10"] = "chr8.fa"
GENOMENR_TO_FASTAFILE["NC_000009.11"] = "chr9.fa"
GENOMENR_TO_FASTAFILE["NC_000010.10"] = "chr10.fa"
GENOMENR_TO_FASTAFILE["NC_000011.9"] = "chr11.fa"
GENOMENR_TO_FASTAFILE["NC_000012.11"] = "chr12.fa"
GENOMENR_TO_FASTAFILE["NC_000013.10"] = "chr13.fa"
GENOMENR_TO_FASTAFILE["NC_000014.8"] = "chr14.fa"
GENOMENR_TO_FASTAFILE["NC_000015.9"] = "chr15.fa"
GENOMENR_TO_FASTAFILE["NC_000016.9"] = "chr16.fa"
GENOMENR_TO_FASTAFILE["NC_000017.10"] = "chr17.fa"
GENOMENR_TO_FASTAFILE["NC_000018.9"] = "chr18.fa"
GENOMENR_TO_FASTAFILE["NC_000019.9"] = "chr19.fa"
GENOMENR_TO_FASTAFILE["NC_000020.10"] = "chr20.fa"
GENOMENR_TO_FASTAFILE["NC_000021.8"] = "chr21.fa"
GENOMENR_TO_FASTAFILE["NC_000022.10"] = "chr22.fa"
GENOMENR_TO_FASTAFILE["NC_000023.10"] = "chrX.fa"
GENOMENR_TO_FASTAFILE["NC_000024.9"] = "chrY.fa"
GENOMENR_TO_FASTAFILE["NC_001807.4"] = "chrM.fa"
MAX_HAIRPIN_LEN = 80
MIN_HAIRPIN_LEN = 46

#  0          1    2                                 3           4                           5                           6
# ['1-15830', '-', 'gi|224589818|ref|NC_000006.11|', '72113295', 'AGCTTCCAGTCGAGGATGTTTACA', 'IIIIIIIIIIIIIIIIIIIIIIII', '0']
def findall(unfixed_lines):
    fixed_lines = [line.strip().split("\t") for line in unfixed_lines]
    
    print len(fixed_lines)
    print fixed_lines[0]
    
    genome_tree = bio.GenomeIntervalTree()
    candidate_tree = bio.GenomeIntervalTre # only candidates here
    
    genomes = set()
    for prop in fixed_lines:
        
        is_forward_strand = prop[1] == "+"
        genomeref = prop[2]
        genomenr = genomeref.split("|")[3]
#         name_nr = genomenr.split("_")[1]
#         nr = name_nr.split(".")[0] # the genome number, use GENOMENR_TO_FASTAFILE to find genome file (hack)
        genome_offset = int(prop[3]) # offset into the genome, 0-indexed
        dna_sequence = prop[4]
        
        # test for candidate with this as 3' end
        outside = genome_offset + len(dna_sequence) - MAX_HAIRPIN_LEN - 1
        begin = genome_offset + len(dna_sequence) - MAX_HAIRPIN_LEN
        end = genome_offset + len(dna_sequence) - MIN_HAIRPIN_LEN
        
        three_hits = genome_tree[genomenr][begin, end]
        if three_hits:
            outside_hits = genome_tree[genomenr][outside]
            three_hits.difference_update(outside_hits)
            #TODO: add candidate(s) to candidate set
            
        
        # test for candidates with this as 5' end
        outside = genome_offset + MAX_HAIRPIN_LEN + 1
        begin = genome_offset + MIN_HAIRPIN_LEN
        end = genome_offset + MAX_HAIRPIN_LEN
        
        
        three_hits = genome_tree[genomenr][begin, end]
        if three_hits:
            outside_hits = genome_tree[genomenr][outside]
            three_hits.difference_update(outside_hits)
            #TODO: add candidate(s) to candidate set
        
        
        
        
        
        
        genomes.add(genomenr)
    
#     print "genomes: ", len(genomes)
#     for g in sorted(genomes):
#         print g
#     genomelist = [x for x in genomes.]
#     genomelist = genomelist.sort()
#     print "genomes??"
#     for x in genomelist:
#         print x
    
    
    