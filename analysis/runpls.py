# import os
# from bowtie import bowtie_get

import time
start_time = time.clock()
from candidates.microseqs import align_small_seqs
import numpy
from sklearn import svm


from inputs import merge
from inputs import mirbase
from inputs import miRNA
from genes import gene
from candidates import interval_tree_search
from candidates import heterogenity
from candidates import vienna
from candidates import tailing
from candidates import entropy
from candidates import quality
from candidates import structure
from candidates import overhang

from ml import vectorize


def _align_bowtie(bowtie_output_file, collapsed_seq_file):
    from subprocess import check_output
    import os # hack: setting path to bowtie index
    os.environ['BOWTIE_INDEXES'] = "/home/hakon/Skrivebord/h_sapiens_37_asm.ebwt"

    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", "-v 0", "-a", "-m 4",
                   human_index, collapsed_seq_file, bowtie_output_file]

    print check_output(bowtie_cmds)





def main():
#     start_time = time.clock()
    print "starting miRNA"
    
    print "merging collapsed files"

    fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed",
                    "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed",
                    "SRR207110.collapsed", "SRR207111.collapsed", "SRR207112.collapsed"]    
    
#     fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed",
#                     "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed"]
    fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed", "SRR207111.collapsed"]
#     fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed"]
#     fasta_files = ["SRR797062.collapsed"]
#     fasta_files =  ["SRR207110.collapsed", "SRR207111.collapsed", "SRR207112.collapsed"] 
#     fasta_file = "SRR797062.fa"
    
    dict_collapsed = merge.collapse_collapsed(fasta_files)
    
#     split small and larger sequences
#     write reads to file
    all_reads_file = "all.collapsed"
    reads, reads_count, small_reads, small_reads_count = merge.filter_seqeunces(dict_collapsed, 18)
    merge.write_collapsed(all_reads_file, reads, reads_count)

    print "reads:", len(reads), "small:", len(small_reads), len(small_reads_count)
    
    
#     aligning to genome using bowtie
    bowtie_output = "bowtie_out.map"
    _align_bowtie(bowtie_output, all_reads_file)
    print "finished bowtie in ", time.clock() - start_time, " seconds" 
    
#     read genome alignment from bowtie
    fixed_lines = [line.strip().split("\t") for line in open(bowtie_output)] 
    print "read positions in ", time.clock() - start_time, " seconds"
    
    
    
    hairpin_file = "hairpin.fa"
    mature_seq_file = "mature.fa"
    miRNA_file_name = "mirnas.fa"
    high_conf_file = "high_conf_hairpin.fa"
    miRNA_family_file = "miFam.dat"
    
    print "loading miRNA hairpins:"
    
    hsa_to_hairpin, other_to_hairpin = mirbase.read_miRNA_fasta(hairpin_file)
    hsa_to_mature, other_to_mature = mirbase.read_miRNA_fasta(mature_seq_file)
    
    miRNA_species = mirbase.similar_hairpins(hsa_to_hairpin, other_to_hairpin)
    
    hairpinID_to_mature = mirbase.combine_hairpin_mature(hsa_to_hairpin, hsa_to_mature)
    miRNA_high_conf = miRNA.read_high_confidence(high_conf_file)
    
    print len(miRNA_high_conf)
    print miRNA_high_conf.issubset(miRNA_species.keys())
    
    print "\nmifam\n"
    miRNA_fam = miRNA.read_family(miRNA_family_file)
    
#     assert False
    
    print len(miRNA_fam)
    for x in miRNA_species.keys():
        print x
        break
    for x in miRNA_fam.keys():
        print x
        break
    
    print len(set(miRNA_species.keys()) & set(miRNA_fam.keys()))
#     assert False
    
    
#     for k,v in miRNA_species.iteritems():
#         if v > 0:
#             print k,v
    
    
#     assert False
#     get mirbase entries
#     miRNAid_to_hairpin = mirbase.read_miRNA("mature.fa", "hairpin.fa")
#     print "got all miRNAs from mirbase files", len(miRNAid_to_hairpin)


    
#     miRNA_species = mirbase.mirna_copies(miRNAid_to_hairpin.values(), "hairpin.fa")
    
#     run write micro rnas to file

    mirbase.write_miRNA(hsa_to_hairpin, miRNA_file_name)
    print "wrote human miRNAs to file", time.clock() - start_time, " seconds"
    
#     run bowtie to find miRNA positions
    miRNA_bowtie_output = "miRNA.map"
    _align_bowtie(miRNA_bowtie_output, miRNA_file_name)
    
    print "aligned miRNAs in", time.clock() - start_time, " seconds"
    
    miRNA_bowtie_hits = [line.strip().split("\t") for line in open(miRNA_bowtie_output)] 
    
    unique_mirna_hits = set([x[0] for x in miRNA_bowtie_hits])
    
    print "miRNA bowtie hits:", len(miRNA_bowtie_hits)
    print "unique miRNA hits:", len(unique_mirna_hits)

    

#     assert False
#     using sequence tree to find possible candidates
#     candidate_tree, sequence_tree, candidates, seq_to_candidates = interval_tree_search.find_candidates(fixed_lines)
    candidate_tree, sequence_tree, candidates, seq_to_candidates = interval_tree_search.find_candidates_2(fixed_lines)
    
    print "found candidates in ", time.clock() - start_time, " seconds"
    print "bowtie hits", len(fixed_lines)
    print "candidate tree", len(candidate_tree)
    print "candidates", len(candidates)
    print "sequence tree", len(sequence_tree)
    print "mapped seqs", len(candidates[0].all_mapped_sequences)
    
#     assert False
#     candidate_list = gene.find_all(candidates)

# 0            1   2[0] [1]      [2] [3]
# ['1-15830', '-', 'gi|224589818|ref|NC_000006.11|',
#         NC_000006.11

        
    gene.include_padding(candidates)
    print "padded all candidates in ", time.clock() - start_time, " seconds"
    
#     assert False

#     assert False
    print "\nrunning viennafold"
    vienna.energy_fold(candidates) # slow
    
    print "candidates", len(candidates)
    print "align miRNAs to other sequences"
    candidate_to_miRNA = interval_tree_search.align_miRNAs(miRNA_bowtie_hits,
                                                           hairpinID_to_mature,
                                                           candidate_tree,
                                                           sequence_tree,
                                                           miRNA_species,
                                                           miRNA_high_conf)
    
    print "aligning small seqs"
    align_small_seqs(candidates, small_reads, small_reads_count)
    print "finished aligning small seqs"
    
#     sequence_freq = reads.readcollapsed(fasta_file)
#     print len(sequence_freq)    
    
    
    # run and set vienna RNAfold + energy on all candidates

    
    
    not_mapped_reads = [structure.Sequence(i,n,read) for i,(read,n) in 
                        enumerate(zip(reads, reads_count))
                        if read not in seq_to_candidates]
#     
#     print len(not_mapped_reads)
#     print len(seq_to_candidates)
#     print len(not_mapped_reads) + len(seq_to_candidates)
#     print len(reads)
#     
#     A/U ends for all remaining candidates
    tailing.tailing_au(candidates, not_mapped_reads)
     

    overhang.get_alignment(candidates)

#     degree of entropy in structure and nucleotides
    entropy.entropy(candidates)
#      
#     heterogenity (position counting)
    heterogenity.heterogenity(candidates)
#      
#     candidate quality: nr of sequence hits / all candidate hits for given sequences
    quality.candidate_quality(candidates, seq_to_candidates)
    
    
    
    print
    print "finished features in ", time.clock() - start_time, " seconds"
    
    for k,v in candidate_to_miRNA.iteritems():
        print k,v
    
    miRNAs = []
    miRNA_annotations = [] 
    only_candidates = []
    
    mi_keys = set(candidate_to_miRNA.keys())
    candidate_keys = set([c.chromosome + str(c.pos_5p_begin) for c in candidates])
    
    print "miRNAs:", len(mi_keys)
    print "all:", len(candidate_keys)
    print "pls be gone:", len(candidate_keys - mi_keys)
    
    
    for candidate in candidates:
#         print candidate, candidate in candidate_to_miRNA
        hashval = candidate.chromosome + str(candidate.pos_5p_begin)
        if hashval in candidate_to_miRNA:
            miRNAs.append(candidate)
            mi = candidate_to_miRNA[hashval]
            print mi, miRNA_species[mi]
            miRNA_annotations.append(miRNA_species[mi])
        else:
            only_candidates.append(candidate)
            
    print len(candidate_to_miRNA)
    print len(only_candidates)
    
    for mi in miRNAs:
        aaa = mi.small_subs
        bbb = sum([ float(x.data[1].split("-")[1]) for x in mi.mapped_sequences])
        ccc = aaa * 1.0 / bbb 
        print ccc, "\t\t\t", aaa, bbb, ccc
    
    assert  False

    learn_miRNAs = vectorize.candidates_to_array(miRNAs)
    class_miRNAs = numpy.array(miRNA_annotations)
    candidate_array = vectorize.candidates_to_array(only_candidates)
    
    for mi, nr in zip(learn_miRNAs, class_miRNAs):
        print mi, nr
    
    
    print "mirnas",len(learn_miRNAs)
    print "mirna species", len(class_miRNAs)
    print "candidates", len(candidate_array)
    
    learner = svm.SVC(probability=True)
    learner.fit(learn_miRNAs, class_miRNAs)
    res = learner.predict(candidate_array)
    print res


    with open("profit.txt", "w") as outfile:
        for r in res:
            outfile.write(str(r) + "\n")
        

# if __name__ == "__main__":
#     main()
main()
print "finished everything in ", time.clock() - start_time, " seconds"







