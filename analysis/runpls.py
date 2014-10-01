from bowtie import bowtie_get
from analysis import bowtie



def main():
    fasta_file = "SRR797062.fa"
    bowtie_output = "test.map"
    human_index = "h_sapiens_37_asm"
    bowtie_cmds = ["bowtie","-f", human_index ,fasta_file , bowtie_output ]
    
#     bowtie_get.runcommand(bowtie_cmds)
    
    






















if __name__ == "__main__":
    main()