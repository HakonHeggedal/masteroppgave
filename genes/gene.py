'''
Created on 6. okt. 2014

@author: hakon
'''
import os
from candidates import structure

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
    candidate_list = []
    fails = 0
    
    for tree in interval_trees:
#         file_name = GENOMENR_TO_FASTAFILE[tree]
        file_name = tree + ".fa"
        
        print "reading from", file_name
        
        with open(file_name) as chr_file:
            
            _ = chr_file.readline() # first line is info, not used
            chr_lines = chr_file.readlines() # current genome
            line_len = len(chr_lines[0].strip())
            
#             print line_len
#             print len(chr_lines) * line_len, len(chr_lines)
            
            for interval in sorted(interval_trees[tree]):
                read_start = interval.begin - extra
                start_line = read_start // line_len
                start_pos = read_start % line_len

                read_end = interval.end + extra
                end_line = read_end // line_len
                end_pos = read_end % line_len
                
#                 print
#                 print "new interval"
#                 print read_start, start_line
#                 print read_end, end_line
#                 print "size:", interval.end-interval.begin,
#                 print "\tinterval seqs:", interval.data[2], interval.data[4],
#                 print "\tinterval sizes:", len(interval.data[2]), len(interval.data[4])
#                 print "read size:", read_end - read_start
#                 print "full size should be:", interval.end-interval.begin + extra + extra
                
                #assembling sequence
                padded = chr_lines[start_line][start_pos:].strip()
                for x in xrange(start_line+1, end_line):
#                     print x
                    padded += chr_lines[x].strip()
                padded += chr_lines[end_line][:end_pos].strip() if start_line != end_line else ""
                
                padded = padded.upper()
                padded = padded.strip()
                padded = padded.strip("\n")
                
                hairpin = padded[extra:-extra]
                
#                 [strand, 5'name, 5'sequence, 3'name, 3'sequence]
#                 five = interval.data[2]
#                 three = interval.data[4]
#                 strand = interval.data[0]
#                 five_id = interval.data[1]
#                 three_id = interval.data[3]
                
#                 print "full sequence:", padded, len(padded)
#                 print interval.data
#                 print "interval.data[2] in padded", interval.data[2] in padded
#                 print "interval.data[4] in padded", interval.data[4] in padded
                
                if interval.data[2] not in padded or interval.data[4] not in padded:
                    fails += 1
#                 print len(hairpin), hairpin
#                 print interval.data[2], interval.data[4]
                
                canidate = structure.Candidate(hairpin, padded, interval.data)
                
                candidate_list.append(canidate)

#         break
    
    os.chdir(p)
    
    print "candidates not mapping to genome:", fails
    
    
    return candidate_list
        
        



        
        