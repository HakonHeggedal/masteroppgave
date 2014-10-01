import os
from bowtie import bowtie_get
from candidates import interval_tree


def main():
    fasta_file = "SRR797062.fa"
    bowtie_output = "test.map"
    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie", "-f", human_index, fasta_file, bowtie_output]
    
    bowtie_get.runcommand(bowtie_cmds)
#     bowtie_get.runcommand(["bowtie","-f", "h_sapiens_37_asm", "SRR797062.fa", "test.map"])
    
    unfixed_lines = open(bowtie_output).readlines() # read 
    interval_tree.findall(unfixed_lines)

    





















# if __name__ == "__main__":
#     main()

main()