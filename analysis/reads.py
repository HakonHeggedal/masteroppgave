
# all RNA lines in the files
rnas = [l for l in open("mature.fa") if l[0] is not ">"]

def readcollapsed(filename):
    ''' returns a dict with sequence -> frequency for each sequence in fasta file '''
    freqdick = {}
    
    with open(filename, "r") as inputfile:
        
        for line in inputfile:
            line = line.strip()
            count = 0
            if line[0] == ">":
                count = int(line.split(-)[1])
            elif count > 0:
                freqdick[line] = count
                
    return freqdick
            
            