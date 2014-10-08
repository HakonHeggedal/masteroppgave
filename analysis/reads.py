import os
# all RNA lines in the files
# rnas = [l for l in open("mature.fa") if l[0] is not ">"]

def readcollapsed(filename):
#     print os.getcwd()
    ''' returns a dict with sequence -> frequency for each sequence in fasta file '''
    sequences = []
#     print "trying to read", filename
    
    with open(filename, "r") as inputfile:
        count = 0
#         print "reading file:", filename
        for line in inputfile:
            line = line.strip()

            if line[0] == ">":
                count = int(line.split("-")[1])
#                 print count
            elif count > 0:
                sequences.append((line, count))
#                 print "seq",sequences
#                 freqdick[line] = count
#                 print freqdick[line], line
            else:
                count = 0
    
    
    
#     print len(sequences)
    return sequences
            

readcollapsed("SRR797062.fa")
