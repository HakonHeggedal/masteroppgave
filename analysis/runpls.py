
import time
from inputs import special_types


from ml.miRNA_conf_group import create_folds
# import itertools

from candidates.microseqs import align_small_seqs, small_seq_stats,\
    length_distribution
# import numpy
from sklearn import svm, preprocessing, metrics


from inputs import merge
from inputs import mirbase
from inputs import miRNA


from genes import gene

from inputs import dead_mirna

from candidates import interval_tree_search, interval_tree_miRNA, interval_tree_dead
from candidates import hairpin
# from candidates import stem
from candidates import heterogenity
from candidates import vienna
from candidates import tailing
from candidates import entropy
from candidates import quality
from candidates import structure
from candidates import overhang


from misc import plot_any

from ml import vectorize



# from matplotlib import pyplot

import pickle
from candidates.hairpin import hairpin_stats
# from misc import
# from misc import 
from misc.plot_correlation import plot_kstest, plot_ttest
# from misc import 
from inputs.merge import one_large_fasta
from inputs.mirbase import human_only_file, not_human_file, write_human_hairpins
from misc.correlation_plots import plot_pearson_correlation, plot_spearman_correlation,\
    plot_correlation
import itertools
from other_stuff.mirna_split import fix_miRNA_training_test
from other_stuff.mirdeep_results import mirdeep_make_roc_data

def _align_bowtie(bowtie_output_file, collapsed_seq_file):
    from subprocess import check_output
    import os 
    # hack: setting path to bowtie index:
    os.environ['BOWTIE_INDEXES'] = "/home/hakon/Skrivebord/h_sapiens_37_asm.ebwt"

    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", "-v 0", "-a", "-m 3",
                   human_index, collapsed_seq_file, bowtie_output_file]

    print check_output(bowtie_cmds)



def main():
    start_time = time.clock()
    print "starting miRNA analysis"
    

    
#     fasta2 = ["Demux.SRhi10002.Adipocyte", "Demux.SRhi10002.Alveolar", "Demux.SRhi10002.Amniotic",
#              "Demux.SRhi10002.Dendritic1", "Demux.SRhi10002.Dendritic2", "Demux.SRhi10002.Endothelial",
#              "Demux.SRhi10002.Fibroblast1", "Demux.SRhi10002.Fibroblast2", "Demux.SRhi10002.Fibroblast3",
#              "Demux.SRhi10002.Intestinal", "Demux.SRhi10002.Meningeal", "Demux.SRhi10002.Mesenchymal",
#              "Demux.SRhi10002.Osteoblast", "Demux.SRhi10002.Pericytes", "Demux.SRhi10002.Renal",
#              "Demux.SRhi10002.Sebocyte1", "Demux.SRhi10002.Sebocyte2", "Demux.SRhi10002.SmoothBrachiocephalic",
#              "Demux.SRhi10002.SmoothProstate", "Demux.SRhi10002.SmoothSubclavian", "Demux.SRhi10002.SmoothUterine"]
#     
#     fasta3 = ["Demux.SRhi10003.Adipocyte", "Demux.SRhi10003.Amniotic%20Epithelial", "Demux.SRhi10003.amniotic%20membrane",
#               "Demux.SRhi10003.Endothelial0", "Demux.SRhi10003.Endothelial1", "Demux.SRhi10003.Endothelial2",
#               "Demux.SRhi10003.Fibroblast1", "Demux.SRhi10003.Fibroblast2", "Demux.SRhi10003.Fibroblast3",
#               "Demux.SRhi10003.Keratinocyte", "Demux.SRhi10003.Mesenchymaladipose", "Demux.SRhi10003.Mesenchymalbone",
#               "Demux.SRhi10003.Osteoblast", "Demux.SRhi10003.Pancreatic", "Demux.SRhi10003.Peripheral",
#               "Demux.SRhi10003.Prostate", "Demux.SRhi10003.Renal", "Demux.SRhi10003.Sertoli",
#               "Demux.SRhi10003.Skeletal", "Demux.SRhi10003.SmoothBrain", "Demux.SRhi10003.SmoothPulmonary",
#               "Demux.SRhi10003.SmoothUmbilical"]
#     
#     fasta2 = ["hg19/"+n for n in fasta2]
#     fasta3 = ["hg19/"+n for n in fasta3]
# 
#     fasta4 = ["hg19/Demux.SRhi10004."+str(i) for i in range(1,23)]
#     fasta5 = ["hg19/Demux.SRhi10005."+str(i) for i in range(1,24)]
#     fasta_files.extend(fasta2)
#     
#     fasta2.extend(fasta3)
#     fasta2.extend(fasta4)
#     fasta2.extend(fasta5)
  
  
#     
#     fasta_files = ["SRR797059.collapsed", "SRR797060.collapsed", "SRR797061.collapsed",
#                     "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed",
#                     "SRR207110.collapsed", "SRR207111.collapsed", "SRR207112.collapsed"]    
# 
#     fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed",
#                     "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed"]
#     fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed", "SRR207111.collapsed"]
#     fasta_files = ["SRR797060.collapsed", "SRR797061.collapsed"]
  
#     fasta_files =  ["SRR207110.collapsed", "SRR207111.collapsed", "SRR207112.collapsed"] 
#     fasta_file = "SRR797062.fa"



  
    fasta_files = ["SRR797062.collapsed"] #  small file for fast testing
    fasta_files_large_folder = ["fastas/Demux.SRhi." + str(i) + ".collapsed"  for i in range(296)]
      
  
    fasta_files_small = ["SRR797059.collapsed", "SRR797060.collapsed", "SRR797061.collapsed",
                    "SRR797062.collapsed", "SRR797063.collapsed", "SRR797064.collapsed",
                    "SRR207110.collapsed", "SRR207111.collapsed", "SRR207112.collapsed",
                    "SRR207113.collapsed", "SRR207114.collapsed", "SRR207115.collapsed",
                    "SRR207116.collapsed", "SRR207117.collapsed", "SRR207118.collapsed",
                    "SRR207119.collapsed"]
      
  
    fasta_files = fasta_files_large_folder
#     fasta_files = fasta_files_large_folder[:40]
  
    hairpin_file = "hairpin.fa"
    mature_seq_file = "mature.fa"
    miRNA_file_name = "mirnas.fa"
    high_conf_file = "high_conf_hairpin.fa"
    miRNA_family_file = "miFam.dat"
    other_types = "mirTrons_other.txt"
#     dead_mirnas = "miRNA.dead"
    dead_mirnas = "dead_list"
    dead_mirna_hairpins = "dead_hairpins.txt"
    dead_mirna_bowtie_file = "dead_hairpin_bowtie.fa"
    dead_mirna_bowtie_out = "dead_hairpin_locations.map"
      
      
    all_reads_file = "all.collapsed"
    bowtie_output = "bowtie_out.map"
      
    miRNA_bowtie_output = "miRNA.map"
    ml_folds = 10
    
    is_new_run = True
#     is_new_run = False
    
    #===========================================================================
    # #  making data for mirdeep2
    #===========================================================================
    # not_human_file(mature_seq_file, "other_matures.fa")
    # assert 0
    # human_only_file(mature_seq_file, "human_matures.fa")
    #  
    #  
    # hsa_to_hairpin, other_to_hairpin = mirbase.read_miRNA_fasta(hairpin_file)
    #  
    # write_human_hairpins(hsa_to_hairpin, "human_hairpins.fa")
    # assert 0
    #  
    #  
    # assert 0
    #  
    # one_large_fasta(fasta_files, "all_fasta.fa")
    # assert 0
    #   
    #  
    #===========================================================================
    
    
    if is_new_run:
        
        id_to_dead_hp, id_to_dead_mature = dead_mirna.get_hairpin(dead_mirnas)
    #     assert False
        
        print "merging",len(fasta_files), "collapsed files" if len(fasta_files)>1 else ""
            
        dict_collapsed = merge.collapse_collapsed(fasta_files, min_len=8, min_count=2)
            
    #     split small and larger sequences
    #     write reads to file
            
        reads, reads_count, small_reads, small_reads_count = merge.filter_seqeunces(dict_collapsed, 18)
        merge.write_collapsed(all_reads_file, reads, reads_count)
        
        print "long reads:", len(reads), "small:", len(small_reads), len(small_reads_count)
        print "fraction small / all:", len(small_reads)*1.0 / (len(small_reads) + len(reads))

            
            
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
            
        hairpinID_to_mature, harpinID_to_matseqs = mirbase.combine_hairpin_mature(hsa_to_hairpin, hsa_to_mature)
        miRNA_high_conf = miRNA.read_high_confidence(high_conf_file)
         
         
#         pickle.dump(harpinID_to_matseqs, open("harpinID_to_matseqs.p", "wb"))
#         pickle.dump(hsa_to_hairpin, open("hsa_to_hairpin.p", "wb"))
            
        print "\nhigh confidence set:", len(miRNA_high_conf),
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
            
            
        print "\nDead mirnas"
        mirbase.write_dead_mirna(id_to_dead_hp, dead_mirna_bowtie_file)
        _align_bowtie(dead_mirna_bowtie_out, dead_mirna_bowtie_file)
        dead_miRNA_hits = [line.strip().split("\t") for line in open(dead_mirna_bowtie_out)]
            
        print "TOTAL dead mirnas:", len(dead_miRNA_hits)
        
        
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
        candidate_to_miRNA = interval_tree_miRNA.align_miRNAs(miRNA_bowtie_hits,
                                                               hairpinID_to_mature,
                                                               harpinID_to_matseqs,
                                                               candidate_tree,
                                                               candidates,
                                                               sequence_tree,
                                                               seq_to_candidates,
                                                               miRNA_species,
                                                               miRNA_high_conf)
        
        mirdeep_make_roc_data(candidate_tree, candidate_to_miRNA, miRNA_high_conf)
        
        candidate_to_dead = interval_tree_dead.align_dead_miRNAs(dead_miRNA_hits,
                                                                 id_to_dead_hp,
                                                                 id_to_dead_mature,
                                                                 candidate_tree,
                                                                 candidates,
                                                                 sequence_tree,
                                                                 seq_to_candidates)
        

            
        print "\npadding all miRNA and Candidates"
        gene.include_padding(candidates)
        print "padded all candidates in ", time.clock() - start_time, " seconds"
            
        
        print "\nrunning vienna rnafold"
        vienna.energy_fold2(candidates)
          
        print "finished vienna folding"
          
    #     align_small_seqs(candidates, small_reads, small_reads_count)
    #     small_seq_stats(candidates)
    #     
    #     plot_any.plot(candidates, candidate_to_miRNA, candidate_to_dead, miRNA_high_conf, "ratio_short_long_5p", False )
    #     assert 0
        
        print "saving 123"
        
        pickle.dump(candidate_tree, open("candidate_tree.p", "wb"))
        pickle.dump(candidates, open("candidates_pre.p", "wb"))
        pickle.dump(candidate_to_miRNA, open("candidate_to_miRNA.p", "wb"))
        pickle.dump(miRNA_high_conf, open("miRNA_high_conf.p", "wb"))
           
        print "saving 234"
           
        pickle.dump(candidate_to_dead, open("candidate_to_dead.p", "wb"))
        pickle.dump(miRNA_fam, open("miRNA_fam.p", "wb"))
        pickle.dump(small_reads, open("small_reads.p", "wb"))
        pickle.dump(small_reads_count, open("small_reads_count.p", "wb"))
        pickle.dump(seq_to_candidates, open("seq_to_candidates.p", "wb"))
        
        pickle.dump(reads, open("reads.p", "wb"))
        pickle.dump(reads_count, open("reads_count.p", "wb"))
        
          
        print "saved 456"

    print "loading miRNAs"
    harpinID_to_matseqs = pickle.load( open("harpinID_to_matseqs.p", "rb"))
    hsa_to_hairpin = pickle.load( open("hsa_to_hairpin.p", "rb"))
    
    print "loading tree..."

    candidate_tree = pickle.load( open("candidate_tree.p", "rb"))
    
    print "loading picled stuff ...", time.clock() - start_time
    candidate_to_miRNA = pickle.load( open("candidate_to_miRNA.p", "rb"))
    miRNA_high_conf = pickle.load( open("miRNA_high_conf.p", "rb"))
    
    
    mirdeep_make_roc_data(candidate_tree, candidate_to_miRNA, miRNA_high_conf)
    
    candidates = pickle.load( open("candidates_pre.p", "rb"))
    
    candidate_to_dead = pickle.load( open("candidate_to_dead.p", "rb"))
    miRNA_fam = pickle.load( open("miRNA_fam.p", "rb"))
    
    
    small_reads = pickle.load( open("small_reads.p", "rb"))
    small_reads_count = pickle.load( open("small_reads_count.p", "rb"))
    seq_to_candidates = pickle.load( open("seq_to_candidates.p", "rb"))
    
    reads = pickle.load( open("reads.p", "rb"))
    reads_count = pickle.load( open("reads_count.p", "rb"))
    
    

    

    print "loaded back", time.clock() - start_time
    
    
    

    
    annotated_data, annotations, low_confidence_data = create_folds(candidates, candidate_to_miRNA, candidate_to_dead, miRNA_high_conf, miRNA_fam, ml_folds)
    
    
    fix_miRNA_training_test(annotated_data, annotations, low_confidence_data, hsa_to_hairpin, harpinID_to_matseqs, candidate_to_miRNA)
    
    assert 0
#     length_distribution(small_reads, small_reads_count)
#     
#     assert 0
    # overhang calculated using fold seq.
    overhang.get_alignment(candidates)
    

    # calculating hairpin stats (length + overhang using pairing prob.)
    hairpin_stats(candidates, candidate_to_miRNA, miRNA_high_conf)

    
    def _is_miRNA(c):
        hashval = c.chromosome+c.chromosome_direction+str(c.hairpin_start)
        return hashval in candidate_to_miRNA
        
    
    def _is_hc(c):
        hashval = c.chromosome+c.chromosome_direction+str(c.hairpin_start)
        if hashval in candidate_to_miRNA:
            if candidate_to_miRNA[hashval] in miRNA_high_conf:
                return True
        return False
    
    def _is_dead(c):
        hashval = c.chromosome+c.chromosome_direction+str(c.hairpin_start)
        return hashval in candidate_to_dead
    
    def get_miRNAid(c):
        hashval = c.chromosome+c.chromosome_direction+str(c.hairpin_start)
        return candidate_to_miRNA[hashval]
        
    
    
    print "all candidates+miRNA+other:", len(candidates)
    print "\twith hairpin struct:\t", len([c for c in candidates if c.has_hairpin_struct])
    
    _mirnas = [c for c in candidates if _is_miRNA(c)]
#     _mirnas2 = [c for c in candidates if c.miRNAid != None]
    _mirna_hc = [c for c in candidates if _is_hc(c)]
    
    _mirna_lc = [c for c in _mirnas if not _is_hc(c)]
    
    _mirna_dead = [c for c in candidates if _is_dead(c)]
    
    _not_mirnas = [c for c in candidates if not _is_miRNA(c) and not _is_dead(c)]
#     _not_mirnas2 = [c for c in candidates if c.miRNAid == None]
    
    print "mirnas:", len(_mirnas), len(candidate_to_miRNA)
    print "\twith hairpin struct:\t", len([m for m in _mirnas if m.has_hairpin_struct])
    
    print "candidates:", len(_not_mirnas)
    print "\twith hairpin struct:\t", len([m for m in _not_mirnas if m.has_hairpin_struct])
    
    print "HC mirnas:", len(_mirna_hc), len(miRNA_high_conf)
    print "\twith hairpin struct:\t", len([m for m in _mirna_hc if m.has_hairpin_struct])
    
    print "LC mirnas:", len(_mirna_lc)
    print "\twith hairpin struct:\t", len([m for m in _mirna_lc if m.has_hairpin_struct])
    


    
    
#     fail_hc = [c for c in _mirna_hc if not c.has_hairpin_struct]
#     fail_lc = [c for c in _mirnas if not c.has_hairpin_struct and not _is_hc(c)]
#     
#     hairpin.hairpin_stats(fail_hc, candidate_to_miRNA, miRNA_high_conf)
#     hairpin.hairpin_stats(fail_lc, candidate_to_miRNA, miRNA_high_conf)
#     hairpin.hairpin_stats(_mirnas, candidate_to_miRNA, miRNA_high_conf)

    
    


    
    not_mapped_reads = [structure.Sequence(i,n,read) for i,(read,n) in 
                        enumerate(zip(reads, reads_count))
                        if read not in seq_to_candidates]
    
    


#     aligning small sequences against hairpins
    align_small_seqs(candidates, small_reads, small_reads_count)
    
    
    lc_10 = pickle.load( open("lc_scores_10.p", "rb"))
    lc_names = pickle.load( open("save_low_confidence_names.p", "rb"))
    
    
    
    small_seq_stats(_mirna_lc, lc_10, lc_names, candidate_to_miRNA)
    print lc_names[0]
    assert 0
    
    
    
    small_seq_stats(candidates)    
#     small_seq_stats(_mirna_hc) # for testing only
    



#     A/U ends for all remaining candidates
    tailing.tailing_au_fast(candidates, not_mapped_reads)
#     tailing.tailing_au_simple(candidates, not_mapped_reads)

    
#     tailing.tailing_au(candidates, not_mapped_reads)
    

#     degree of entropy in structure and nucleotides
    entropy.entropy(candidates)
    

#      
#     heterogenity (position counting)
    heterogenity.heterogenity(candidates)

#      
#     candidate quality: nr of sequence hits / all candidate hits for given sequences
    quality.candidate_quality(candidates, seq_to_candidates)
    
    
#     plotting all features

    
       
#     FEATURES = ["hairpin_energy", "hairpin_energy_10", "hairpin_energy_40",
#                 "entropy_nucleotides", "entropy_structure", "heterogenity_5_begin",
#                 "heterogenity_5_end", "heterogenity_3_begin", "heterogenity_3_end",
#                 "quality", "bindings_max_10", "overhang_level_outer_10",
#                 "overhang_outer_10", "overhang_level_inner_10", "overhang_inner_10",
#                 "bulge_factor"]
#      
#      
#     log_scaled = [False]*16 + [True]*3
#     print log_scaled
#     for feat_name, logs in zip(FEATURES, log_scaled):
#       
#         plot_any.plot(candidates, candidate_to_miRNA, candidate_to_dead,
#                       miRNA_high_conf, feat_name, logs )
    
#     removed_nonvalues = [c for c in candidates if c.ratio_short_long_5p != 1.0]
#     
#     print "nonvalues left", len( [c for c in candidates if c.ratio_short_long_5p == 1.0])
#     
#     plot_any.plot(removed_nonvalues, candidate_to_miRNA, candidate_to_dead,
#                       miRNA_high_conf, "ratio_short_long_5p", isLog=True )
    
    short_correlate_13_17_max = []
    maxlen_range = range(13, 18)
    
    for i in maxlen_range:
        res = plot_any.plot(candidates, candidate_to_miRNA, candidate_to_dead,
                  miRNA_high_conf, "short_seq_align_10_"+str(i), isLog=False )
        
        short_correlate_13_17_max.append(res)
       
    for l in short_correlate_13_17_max:
        print l 
    
    print
    for l in zip(*short_correlate_13_17_max):
        print l
        
    (ks_val, p_2s_ks, t_student, p_student, t_welch, p_welch) = zip(*short_correlate_13_17_max)
    
    
    plot_kstest(ks_val, maxlen_range, True)
    plot_ttest(t_student, t_welch, maxlen_range, True)
    

    short_correlate_8_17_min = []
    minlen_range = range(8, 18)
    for i in minlen_range:
        res = plot_any.plot(candidates, candidate_to_miRNA, candidate_to_dead,
                  miRNA_high_conf, "short_seq_align_" + str(i) + "_17", isLog=False )
        
        short_correlate_8_17_min.append(res)
       
    for l in short_correlate_8_17_min:
        print l 
    
    print
    for l in zip(*short_correlate_8_17_min):
        print l
        
    (ks_val, p_2s_ks, t_student, p_student, t_welch, p_welch) = zip(*short_correlate_8_17_min)
    
    plot_kstest(ks_val, minlen_range, False)
    plot_ttest(t_student, t_welch, minlen_range, False)


    FEATURES = [
#                 "short_seq_align", # 0 until best val is found
#                 "ratio_short_long",
#                 "ratio_short_long_logval",
                "leading_au",
                "tailing_au",
                "overhang_inner",
                "overhang_outer",
                "loop_size",
                "folds_5p",
                "folds_3p",
                "folds_before",
                "folds_after", 
                ]
    
#     plot_pearson_correlation(candidates, FEATURES)
#     plot_spearman_correlation(candidates, FEATURES)



#     for feat_name in FEATURES:
#       
#         plot_any.plot(candidates, candidate_to_miRNA, candidate_to_dead,
#                       miRNA_high_conf, feat_name )
          
    
    print
    print " features finished:", time.clock() - start_time, " seconds"


#     
# 
# 
#     
# #     annotations = list(itertools.chain.from_iterable(annotations))
# #     all_annotated = list(itertools.chain.from_iterable(annotated_data))
#     
#     vector_data = [vectorize.candidates_to_array(d) for d in annotated_data]
#     single_vector_data = list(itertools.chain.from_iterable(vector_data) )
#     
# #     print len(vector_data)
#     
#     scaler = preprocessing.StandardScaler().fit(single_vector_data)
#     
#     scaled_vectors = [scaler.transform(d) for d in vector_data]
# 
# #     learner = svm.SVC(probability=True, cache_size=500)
#     
#     threads = ml_folds
#     pool = Pool(threads)
#     
#     print "nr of threads:", threads
#     
#     param = zip([scaled_vectors]*ml_folds, [annotations]*ml_folds, range(ml_folds) )
#     res_lists = pool.map(param_estimate.param_estimate_fold, param)
#     
#     
#     feat_list = zip(*res_lists)
#     for name, res in zip(feature_names, feat_list):
#         m = numpy.mean(res)
#         s = numpy.std(res)
#         d = abs(m) - s
#         print m, "\t", s, "\t", d, "\t", name
# #         print sum(res)/len(res),"\t", name, res
#     
#     
#     assert False
#         

    # feature selection:
    # only one fold first:
    
    
#     derp = annotated_data
#     derp = map(vectorize.candidates_to_array, derp)
    
    
    
#     test = annotated_data[0]
#     test_annotations = annotations[0]
# 
#     train = list(itertools.chain.from_iterable(annotated_data[1:]))
#     train_annotations = list(itertools.chain.from_iterable(annotations[1:]))
# 
#     train = vectorize.candidates_to_array(train)
#     test = vectorize.candidates_to_array(test)
#     
#     test = preprocessing.scale(test)
#     train = preprocessing.scale(train)
# 
#     # parameter selection:
#     
#     learner = svm.SVC(probability=True, cache_size=500)
#     learner.fit(train, train_annotations)
#     
#     
#     # roc plot 123
#     probs = learner.predict_proba(test)
#     print "probabilities", probs
#     fpr, tpr, _thresholds = metrics.roc_curve(test_annotations, probs[:,1])
#     
#     roc_auc = metrics.auc(fpr, tpr)
#     print "area under curve:", roc_auc
#     
#     pyplot.plot(fpr, tpr)
#     pyplot.show()
#     
#     base_score = learner.score(test, test_annotations)
#     print "base score: ", base_score
#     assert False
#     # roc plot 234
#      



#===============================================================================
# saving data for use in classification
#===============================================================================
    
    
    def is_good_or_miRNA(c):
        return _is_miRNA(c) or c.has_hairpin_struct or _is_dead(c)
    
    
    def is_bad_or_miRNA(c):
        return _is_miRNA(c) or _is_dead(c) or not c.has_hairpin_struct

    
    def is_good_candidate(c):
        return c.has_hairpin_struct and not _is_dead(c) and not _is_miRNA(c)
    
    
    
#     # removed bad candidates 
#     removed_bad_candidates = [c for c in candidates if is_good_or_miRNA(c)] # classifying LC, not using low quality miRNA
#     annotated_data, annotations, low_confidence_data = create_folds(removed_bad_candidates, candidate_to_miRNA, candidate_to_dead, miRNA_high_conf, miRNA_fam, ml_folds)
   
    # all candidates
    annotated_data, annotations, low_confidence_data = create_folds(candidates, candidate_to_miRNA, candidate_to_dead, miRNA_high_conf, miRNA_fam, ml_folds)
    
    # removed the good candidates (which are classified)
    print "candidates?"    
    not_good_candidates = [c for c in candidates if is_bad_or_miRNA(c)] # classifying bad candidates
    annotated_data_new, annotations_new, _low_confidence_data = create_folds(not_good_candidates, candidate_to_miRNA, candidate_to_dead, miRNA_high_conf, miRNA_fam, ml_folds)
    
    
    

    
    low_confidence_names = [map(get_miRNAid, l) for l in low_confidence_data]
    
    
    
    def sum_reads(c):
        reads = 0.0
        for el in c.mapped_sequences:
            name = el.data[1]
            reads += float(name.split("-")[1])
        
        return reads
        
            
    
    low_confidence_reads = [sum_reads(c) for fold in low_confidence_data for c in fold ]
    print low_confidence_reads
    
    pickle.dump(low_confidence_reads, open("low_confidence_reads.p", "wb"))
    assert 0
    
    
    # data for LC classification: 
    vector_data = map(vectorize.candidates_to_array, annotated_data)
    scaled_data = map(preprocessing.scale, vector_data)
    
    
        # low confidence miRNA
    low_confidence_data = map(vectorize.candidates_to_array, low_confidence_data)
    low_confidence_data = map(preprocessing.scale, low_confidence_data)
    
    
    # data for new miRNA classification:
    vector_data_new = map(vectorize.candidates_to_array, annotated_data_new)
    scaled_data_new = map(preprocessing.scale, vector_data_new)
    
    
        # good candidates (not miRNA, has hairpin struct)
    hp_candidates = [c for c in candidates if is_good_candidate(c)]
    print "hp_candidates", len(hp_candidates)
    hp_candidates = vectorize.candidates_to_array(hp_candidates)
    hp_candidates = map(preprocessing.scale, hp_candidates)
    
    print
    print "saving data ...",
    #===========================================================================
    # data used for classifying LOW CONFIDENCE:
    #===========================================================================
    pickle.dump(scaled_data, open("save_scaled_data.p", "wb")) # candidates, HC, DEAD
    pickle.dump(annotations, open("save_an.p", "wb")) # annotations for data over [00001111000]
    pickle.dump(low_confidence_data, open("save_low_confidence_data.p", "wb")) # low confidence miRNA
    pickle.dump(low_confidence_names, open("save_low_confidence_names.p", "wb")) # mirbase names, like >hsa-mir-516a-1
#     pickle.dump(annotated_data, open("save_da.p", "wb")) # extra stuff not needed
    
    #===========================================================================
    # data used for classifying new miRNA:
    #===========================================================================
    pickle.dump(scaled_data_new, open("save_scaled_data_new.p", "wb")) # non-hp candidates, HC, DEAD
    pickle.dump(annotations_new, open("save_an_new.p", "wb")) # annotations for data over [001111000]
    pickle.dump(hp_candidates, open("save_hp_candidates_new.p", "wb")) # annotations for data over [001111000]
    #TODO: position in genome for all candidates probably...
    
    
    print "...saved."
    print "Now loading back for testing ...",
    
    annotations = pickle.load( open("save_an.p", "rb"))
    scaled123 = pickle.load( open("save_scaled_data.p", "rb"))
    low_confidence_data = pickle.load( open("save_low_confidence_data.p", "rb"))
    low_confidence_names = pickle.load( open("save_low_confidence_names.p", "rb"))
#     annotated_data = pickle.load( open("save_da.p", "rb"))

    scaled_data_new = pickle.load( open("save_scaled_data_new.p", "rb"))
    annotations_new = pickle.load( open("save_an_new.p", "rb"))
    hp_candidates = pickle.load( open("save_hp_candidates_new.p", "rb"))


     
    print "... done", (len(scaled123))
    print "finished:", time.clock() - start_time, "seconds"
    
    

#    
#     
# if __name__ == "__main__":
#     main()
main()
# print "finished everything in ", time.clock() - start_time, " seconds"






