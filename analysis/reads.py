
# all RNA lines in the files
# rnas = [l for l in open("mature.fa") if l[0] is not ">"]

def readcollapsed(filename):
#     print os.getcwd()
    ''' returns a dict with sequence -> frequency for each sequence in fasta file '''
    sequences = []
#     print "trying to read", filename
    
    with open(filename, "r") as inputfile:
        name = ""
        
        for line in inputfile:
            line = line.strip()
            
            if line[0] == ">":
                name = line[1:].strip()
#                 print count
            elif len(name) > 0:
                sequences.append((name, line))
                name = ""

    
    
#     print len(sequences)
    return sequences
            

# print readcollapsed("SRR797062.fa")
