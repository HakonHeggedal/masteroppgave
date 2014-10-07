'''
Created on 6. okt. 2014

@author: hakon
'''
import os

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


def _find_loki(gene_name, begin, end, extra=40):
    ''' finds the position in the gene
        returns the sequence start to end, and the surroundings
    '''
    length = end - begin

    begin_extra = begin - extra
    end_extra = end + extra
    length_extra = end_extra - begin_extra
    
    with open(gene_name, "r") as gene:
        gene.seek(begin_extra)
        surroundings = gene.read(length_extra)
    
    loki = surroundings[extra:-extra]
    
    return loki, surroundings


def find_all(interval_trees, extra=40 ):
    
    p =  os.getcwd()
    os.chdir(os.pardir)
    os.chdir("genes")
    results = {} 
    
    for tree in interval_trees:
#         file_name = GENOMENR_TO_FASTAFILE[tree]
        file_name = tree + ".fa"
        
        print "reading from", file_name
        
        with open(file_name) as chr_file:
            
            _ = chr_file.readline() # first line is info, not used
            chr_lines = chr_file.readlines() # current genome
            line_len = len(chr_lines[0].strip())
            
            print line_len
            print len(chr_lines) * line_len, len(chr_lines)
            
            for interval in sorted(interval_trees[tree]):
                read_start = interval.begin - extra
                start_line = read_start // line_len
                start_pos = read_start % line_len

                read_end = interval.end + extra
                end_line = read_end // line_len
                end_pos = read_end % line_len
                
                print read_start, start_line
                print read_end, end_line
                
                seq = chr_lines[start_line][start_pos:].strip()
                for x in xrange(start_line+1, end_line-1):
                    print x
                    seq += chr_lines[x].strip()
                seq += chr_lines[end_line][:end_pos].strip() if start_line != end_line else ""
                seq = seq.upper()
                seq = seq.strip()
                seq = seq.strip("\n")
                print "sequence:", seq, len(seq)
                print interval.data
                print "interval.data[2] in seq", interval.data[2] in seq
                
#                 last = ""
# #                 print "can I haz?"
#                 for nr, line in enumerate(chr_lines):
#                     last = last.strip() + line.strip()
#                     last = last.upper()
# #                     print last, len(last)
# #                     
# #                     if nr == 10:
# #                         break
#                     if interval.data[2] in last:
#                         print "YES PLEASE", 
#                         print nr, last,
#                         print interval.data[2]
#                     
#                     last = line.strip()
#                 print last
                
#                 break
        break
    
    os.chdir(p)
        
        



        
        