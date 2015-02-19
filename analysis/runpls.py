# import os
# from bowtie import bowtie_get

import time
from ml.miRNA_group import split_candidates
from inputs import special_types


from ml.miRNA_conf_group import create_folds
import itertools

from candidates.microseqs import align_small_seqs
import numpy
from sklearn import svm, preprocessing

import random
from inputs import merge
from inputs import mirbase
from inputs import miRNA

import math

from genes import gene

from candidates import interval_tree_search, stem
from candidates import heterogenity
from candidates import vienna
from candidates import tailing
from candidates import entropy
from candidates import quality
from candidates import structure
from candidates import overhang

# from misc import energy, plot_tailing
# from misc import plot_read_quality
# from misc import plot_entropy_dna
# from misc import plot_entropy_struct
# from misc import plot_inner_level
# from misc import plot_overhang_outer
# from misc import plot_overhang_inner
# from misc import plot_max_bindings
# from misc import h_5_b
# from misc import h_5_e
# from misc import h_3_b
# from misc import h_3_e
# from misc import plot_reads
from misc import plot_any

from ml import vectorize, param_estimate
from ml.vectorize import feature_names

from multiprocessing import Pool



def _align_bowtie(bowtie_output_file, collapsed_seq_file):
    from subprocess import check_output
    import os # hack: setting path to bowtie index
    os.environ['BOWTIE_INDEXES'] = "/home/hakon/Skrivebord/h_sapiens_37_asm.ebwt"

    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", "-v 0", "-a", "-m 4",
                   human_index, collapsed_seq_file, bowtie_output_file]

    print check_output(bowtie_cmds)





def main():
    start_time = time.clock()
    print "starting miRNA analysis"
    


    fasta_files = ["SRR797059.collapsed", "SRR797060.collapsed", "SRR797061.collapsed",
                    "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed",
                    "SRR207110.collapsed", "SRR207111.collapsed", "SRR207112.collapsed",
                    "SRR207113.collapsed", "SRR207114.collapsed", "SRR207115.collapsed",
                    "SRR207116.collapsed", "SRR207117.collapsed", "SRR207118.collapsed",
                    "SRR207119.collapsed",]
    
    fasta_files = ["SRR797059.collapsed", "SRR797060.collapsed", "SRR797061.collapsed",
                    "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed",
                    "SRR207110.collapsed", "SRR207111.collapsed", "SRR207112.collapsed"]    
#     
    fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed",
                    "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed"]
#     fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed", "SRR207111.collapsed"]
#     fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed"]
#     fasta_files = ["SRR797062.collapsed"]
#     fasta_files =  ["SRR207110.collapsed", "SRR207111.collapsed", "SRR207112.collapsed"] 
#     fasta_file = "SRR797062.fa"

    hairpin_file = "hairpin.fa"
    mature_seq_file = "mature.fa"
    miRNA_file_name = "mirnas.fa"
    high_conf_file = "high_conf_hairpin.fa"
    miRNA_family_file = "miFam.dat"
    other_types = "mirTrons_other.txt"
    
    all_reads_file = "all.collapsed"
    bowtie_output = "bowtie_out.map"
    
    miRNA_bowtie_output = "miRNA.map"
    
    ml_folds = 10

    print "merging",len(fasta_files), "collapsed files" if len(fasta_files)>1 else ""
    dict_collapsed = merge.collapse_collapsed(fasta_files, min_len=10, min_count=2)
    
#     split small and larger sequences
#     write reads to file

    reads, reads_count, small_reads, small_reads_count = merge.filter_seqeunces(dict_collapsed, 18)
    merge.write_collapsed(all_reads_file, reads, reads_count)

    print "long reads:", len(reads), "small:", len(small_reads), len(small_reads_count)
    
    
#     aligning to genome using bowtie

    _align_bowtie(bowtie_output, all_reads_file)
    print "finished bowtie in ", time.clock() - start_time, " seconds" 
    
#     read genome alignment from bowtie
    fixed_lines = [line.strip().split("\t") for line in open(bowtie_output)] 
    print "read positions in ", time.clock() - start_time, " seconds"

    
    print "loading miRNA hairpins:"
    hsa_to_hairpin, other_to_hairpin = mirbase.read_miRNA_fasta(hairpin_file)
    hsa_to_mature, other_to_mature = mirbase.read_miRNA_fasta(mature_seq_file)
    
    special_types.remove_mirTrons(hsa_to_hairpin, other_types)
    special_types.remove_mirTrons(hsa_to_mature, other_types)
    miRNA_species = mirbase.similar_hairpins(hsa_to_hairpin, other_to_hairpin)
    
    hairpinID_to_mature = mirbase.combine_hairpin_mature(hsa_to_hairpin, hsa_to_mature)
    miRNA_high_conf = miRNA.read_high_confidence(high_conf_file)
#     assert False
    
    print "\nhigh confidence set:",len(miRNA_high_conf),
    print miRNA_high_conf.issubset(miRNA_species.keys())
    
    print "\nreading miRNA family info (mifam)"
    miRNA_fam = miRNA.read_family(miRNA_family_file)

#     print len(set(miRNA_species.keys()) & set(miRNA_fam.keys()))
    
#     run write micro rnas to file

    mirbase.write_miRNA(hsa_to_hairpin, miRNA_file_name)
    print "\nwrote human miRNAs to file", time.clock() - start_time, " seconds"
    
#     run bowtie to find miRNA positions

    _align_bowtie(miRNA_bowtie_output, miRNA_file_name)
    
    print "aligned miRNAs in", time.clock() - start_time, " seconds"
    
    miRNA_bowtie_hits = [line.strip().split("\t") for line in open(miRNA_bowtie_output)] 
    
    unique_mirna_hits = set([x[0] for x in miRNA_bowtie_hits])
    
    print "miRNA bowtie hits:", len(miRNA_bowtie_hits)
    print "unique miRNA hits:", len(unique_mirna_hits)
    


#     using sequence tree to find possible candidates
#     candidate_tree, sequence_tree, candidates, seq_to_candidates = interval_tree_search.find_candidates(fixed_lines)
    candidate_tree, sequence_tree, candidates, seq_to_candidates = interval_tree_search.find_candidates_2(fixed_lines)
    

    print "\n\tfound candidates in ", time.clock() - start_time, " seconds"
    print "\tbowtie hits", len(fixed_lines)
    print "\tcandidate tree", len(candidate_tree)
    print "\tcandidates", len(candidates)
    print "\tsequence tree", len(sequence_tree)
    print "\tmapped seqs", len(candidates[0].all_mapped_sequences)


# 0            1   2[0] [1]      [2] [3]
# ['1-15830', '-', 'gi|224589818|ref|NC_000006.11|',
#         NC_000006.11

    print "\naligning miRNAs to sequences"
    candidate_to_miRNA = interval_tree_search.align_miRNAs(miRNA_bowtie_hits,
                                                           hairpinID_to_mature,
                                                           candidate_tree,
                                                           candidates,
                                                           sequence_tree,
                                                           seq_to_candidates,
                                                           miRNA_species,
                                                           miRNA_high_conf)
    
#     assert False
    print "\npadding all miRNA and Candidates"
    gene.include_padding(candidates)
    print "padded all candidates in ", time.clock() - start_time, " seconds"
    

    print "\nrunning viennafold"
    vienna.energy_fold2(candidates) # slow
    
    print len(candidates)


    stem.compute_stem_start(candidates, candidate_to_miRNA, miRNA_high_conf)
    #stats out here:
#     plot_any.plot(candidates, candidate_to_miRNA, miRNA_high_conf, "hairpin_energy_10")
#     assert False
#     energy.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
    
    # create mirna groups for classification
#     print "nr of candidates + miRNAS:", len(candidates)

#     training_lists, test_lists = split_candidates(candidates, candidate_to_miRNA, miRNA_fam, ml_folds)
    
#     print "finished making training/testsets ", time.clock()-start_time, "seconds"                             
#     assert False
    annotated_data, annotations, _unknown_data = create_folds(candidates, candidate_to_miRNA, miRNA_high_conf, miRNA_fam, ml_folds)
    

    align_small_seqs(candidates, small_reads, small_reads_count)
    
    
    not_mapped_reads = [structure.Sequence(i,n,read) for i,(read,n) in 
                        enumerate(zip(reads, reads_count))
                        if read not in seq_to_candidates]
    
#     
#     print len(not_mapped_reads)
#     print len(seq_to_candidates)
#     print len(not_mapped_reas) + len(seq_to_candidates)
#     print len(reads)
#     
#     A/U ends for all remaining candidates
    tailing.tailing_au(candidates, not_mapped_reads)
#     plot_tailing.plot(candidates, candidate_to_miRNA, miRNA_high_conf)

    overhang.get_alignment(candidates)
    
#     plot_max_bindings.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
#     plot_overhang_inner.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
#     plot_overhang_outer.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
#     plot_inner_level.plot(candidates, candidate_to_miRNA, miRNA_high_conf)

#     degree of entropy in structure and nucleotides
    entropy.entropy(candidates)
    
#     plot_entropy_dna.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
#     plot_entropy_struct.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
    
#      
#     heterogenity (position counting)
    heterogenity.heterogenity(candidates)
#     h_5_b.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
#     h_5_e.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
#     h_3_b.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
#     h_3_e.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
#      
#     candidate quality: nr of sequence hits / all candidate hits for given sequences
    quality.candidate_quality(candidates, seq_to_candidates)
#     plot_read_quality.plot(candidates, candidate_to_miRNA, miRNA_high_conf)
    

    
    print
    print "finished features in ", time.clock() - start_time, " seconds"

    


    
#     annotations = list(itertools.chain.from_iterable(annotations))
#     all_annotated = list(itertools.chain.from_iterable(annotated_data))
    
    vector_data = [vectorize.candidates_to_array(d) for d in annotated_data]
    single_vector_data = list(itertools.chain.from_iterable(vector_data) )
    
#     print len(vector_data)
    
    scaler = preprocessing.StandardScaler().fit(single_vector_data)
    
    scaled_vectors = [scaler.transform(d) for d in vector_data]

#     learner = svm.SVC(probability=True, cache_size=500)
    
    threads = ml_folds
    pool = Pool(threads)
    
    print "nr of threads:", threads
    
    param = zip([scaled_vectors]*ml_folds, [annotations]*ml_folds, range(ml_folds) )
    res_lists = pool.map(param_estimate.param_estimate_fold, param)
    
    
    feat_list = zip(*res_lists)
    for name, res in zip(feature_names, feat_list):
        m = numpy.mean(res)
        s = numpy.std(res)
        d = abs(m) - s
        print m, "\t", s, "\t", d, "\t", name
#         print sum(res)/len(res),"\t", name, res
    
    
    assert False
        

    # feature selection:
    # only one fold first:
    
    
    derp = annotated_data
    
    derp = map(vectorize.candidates_to_array, derp)
    
    
    
    test = annotated_data[0]
    test_annotations = annotations[0]

    train = list(itertools.chain.from_iterable(annotated_data[1:]))
    train_annotations = list(itertools.chain.from_iterable(annotations[1:]))

    train = vectorize.candidates_to_array(train)
    test = vectorize.candidates_to_array(test)

    test = preprocessing.scale(test)
    train = preprocessing.scale(train)

    # parameter selection:
    
    learner = svm.SVC(probability=True, cache_size=500)
    learner.fit(train, train_annotations)
    
    base_score = learner.score(test, test_annotations)
    
    print base_score

    
    features = len(train[0])
    
    scores = [0]*features
    
    for i in xrange(features):
        
        train_removed = numpy.delete(train,i, 1)
        test_removed = numpy.delete(test,i, 1)
        
        learner = svm.SVC(probability=True, cache_size=500) 
        learner.fit(train_removed, train_annotations)
        removed_score = learner.score(test_removed, test_annotations)
        scores[i] = base_score - removed_score
        
    print "each feature removed:"
    print scores
    
    for score, name in sorted(zip(scores, feature_names)):
        print  score,"\t", name


#     def score_features(train, train_annotations, test, test_annotations):
#         pass
#         
#         
# 
#     
#     
#     for i in xrange(len(train[0])):
        

        

#     print train[0]
#     print
#     train = numpy.delete(train,0, 1)
#     print train[0]

#     scaler = preprocessing.StandardScaler().fit(test)
    
#     test = preprocessing.normalize(test) #shit
#     train = preprocessing.normalize(train)
    
#     # flatten lists (for easy testing purposes only)
#     annotated_data = list(itertools.chain.from_iterable(annotated_data))
#     annotations = list(itertools.chain.from_iterable(annotations))
#     test_data_candidates = list(itertools.chain.from_iterable(unknown_data))
# 
#     # create vectors
#     annotated_data = vectorize.candidates_to_array(annotated_data)
#     annotations = numpy.array(annotations)
#     unknown_data = vectorize.candidates_to_array(test_data_candidates)
#     
#     # scale data -1, 1
#     annotated_data = preprocessing.scale(annotated_data)
#     unknown_data = preprocessing.scale(unknown_data)
#     
#     
#     learner = svm.SVR(probability=True, cache_size=1000)
# #     classer = svm.SVC(probability=True, cache_size=1000)
#     print "fit"
#     learner.fit(annotated_data, annotations)
# #     classer.fit(annotated_data, annotations)
# 
#     print "learn"
#     res = learner.predict(unknown_data)
# #     cls = learner.predict(unknown_data)
#     
#     print res
#     print max(res)
#     print sum(res)
#     print len(res)
#     
    
    
#     for val, c in sorted(zip(res, test_data_candidates), reverse=True):
#         hashval = c.chromosome + c.chromosome_direction + str(c.pos_5p_begin)
#         print val, "  \t", c.hairpin_energy_10, candidate_to_miRNA[hashval]


#     assert False
#     for mi, nr in zip(learn_miRNAs, class_miRNAs):
#         print mi, nr
    
    
#     print "mirnas",len(learn_miRNAs)
#     print "mirna species", len(class_miRNAs)
#     print "candidates", len(candidate_array)
#     
    
#     learner = svm.SVC(probability=True, cache_size=1000)
#     print "fit"
#     learner.fit(learn_miRNAs, class_miRNAs)
#     print "learn"
#     res = learner.predict(candidate_array)
#     
#     print max(res)
#     print res
# 
# 
#     with open("profit.txt", "w") as outfile:
#         for r in res:
#             outfile.write(str(r) + "\n")
        

# if __name__ == "__main__":
#     main()
main()
# print "finished everything in ", time.clock() - start_time, " seconds"




#     miRNAs = []
#     miRNA_annotations = []
#     only_candidates = []
#     
# #     mi_keys = set(candidate_to_miRNA.keys())
# #     candidate_keys = set([c.chromosome + str(c.pos_5p_begin) for c in candidates])
# #     
# #     print "miRNAs:", len(mi_keys)
# #     print "all:", len(candidate_keys)
# #     print "pls be gone:", len(candidate_keys - mi_keys)
#     
#     
#     for candidate in candidates:
# #         print candidate, candidate in candidate_to_miRNA
#         hashval = candidate.chromosome + candidate.chromosome_direction + str(candidate.pos_5p_begin)
#         if hashval in candidate_to_miRNA:
#             miRNAs.append(candidate)
#             mi = candidate_to_miRNA[hashval]
# #             print mi, miRNA_species[mi]
#             miRNA_annotations.append(miRNA_species[mi])
#         else:
#             only_candidates.append(candidate)
#             
#     print len(candidate_to_miRNA)
#     print len(only_candidates)



#     learn_one = training_lists[0]
#     train_one = training_lists[1] + training_lists[2] + training_lists[3] + training_lists[4]
#     print len(learn_one), len(train_one)
#     
#     learn_one_v = vectorize.candidates_to_array(learn_one)
#     train_one_v = vectorize.candidates_to_array(train_one)
#     
#     learn_one_s = preprocessing.scale(learn_one_v)
#     train_one_s = preprocessing.scale(train_one_v)
#     
#     one_learner = svm.OneClassSVM() #nu=0.1, kernel="rbf", gamma=0.1)
#     one_learner.fit(train_one_s)
#     res = one_learner.predict(learn_one_s)
#     
#     print max(res)
#     print len(res)
#     print sum(res)
#     print res
#     
#     # nr of high confidence among -1s
#     
#     hc_plus = 0
#     hc_minus = 0
#     for val, c in zip(res, learn_one):
#         hashval = c.chromosome + c.chromosome_direction + str(c.pos_5p_begin)
#         
#         if hashval not in candidate_to_miRNA: 
#             continue
#         
#         mi = candidate_to_miRNA[hashval]
#         if mi in miRNA_high_conf:
#             print val, mi
#         
#         if mi in miRNA_high_conf:
#             hc_plus += 1
#         else:
#             hc_minus += 1
#             
#     print "high confidence:"
#     print "+", hc_plus
#     print "-", hc_minus
#     print "total:", hc_plus * 1.0 / (hc_plus + hc_minus)


