# import os
# from bowtie import bowtie_get
from candidates import interval_tree_search
import time
from genes import gene
from subprocess import check_output
import reads
from candidates.vienna import energy_fold


def main():
    start_time =time.clock()
    print "starting"
    fasta_file = "SRR797062.fa"
    bowtie_output = "test.map"
    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", "-v 0", "-a", "-m 10",  human_index, fasta_file, bowtie_output]
    bowtie_str = "bowtie -f -v 0 -a -m 10 h_sapiens_37_asm SRR797062.fa test.map"
    bowtie_str = "bowtie -f -v 0 h_sapiens_37_asm SRR797062.fa test.map"
    
    
    
#     print "running bowtie"
#     outs = check_output(bowtie_str, shell=True).strip()
# #     outs = check_output(["bowtie", human_index, "-c", "ACGTACGT", "ts.txt"], shell=True).strip()
#     print "run commands2"        
# #     print answers
# #     print "errors", errors
#     parts = outs.split("\t")
#     print len(outs), outs
#     print len(parts), parts
#      
    
#     bowtie_get.runcommand(bowtie_cmds)
#     bowtie_get.runcommand(["bowtie","-f", "h_sapiens_37_asm", "SRR797062.fa", "test.map"])
    
    print "finished bowtie in ", time.clock() - start_time, " seconds" 
#     unfixed_lines = open(bowtie_output).readlines() # read
#     fixed_lines = [line.strip().split("\t") for line in unfixed_lines] 
    fixed_lines = [line.strip().split("\t") for line in open(bowtie_output)] 
    
    print "read positions in ", time.clock() - start_time, " seconds" 
    candidates, sequences = interval_tree_search.find_candidates(fixed_lines)
    
    print len(candidates)
    print len(sequences)
    
    print "found candidates in ", time.clock() - start_time, " seconds"
    
    candidate_list = gene.find_all(candidates)
    
    print "candidates:", len(candidate_list)
    print "found all loki in ", time.clock() - start_time, " seconds"
    
    
    sequence_freq = reads.readcollapsed(fasta_file)
    
    print len(sequence_freq)
    
    energy_fold(candidate_list)


# if __name__ == "__main__":
#     main()

main()