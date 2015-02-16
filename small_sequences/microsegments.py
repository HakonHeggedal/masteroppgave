'''compare miRNAs to small segments from same sequencing file.
'''

def compare(rnas, frequencies):
    ''' finds small segments (substrings) from rna-file, and aligns them with the rna
    input: list of rnas to be analysed,
    '''
    # 1: hash subparts of all rnas 
    rnadict = {}
    
    for nr, rna in enumerate(rnas):
        for i in xrange(len(rna)-3):
            subs = rna[i:i+4]
            if subs in rnadict:
                rnadict[subs].add((nr,i))
            else:
                rnadict[subs] = set([(nr,i)])
    
    # 2: find all substring-rnas
    substringdict = {}
    superstrings = set()
    
    for nr, rna in enumerate(rnas):
        print "rna: ", rna
        subs = rna[:4]
        
        if subs not in rnadict:
            continue
        
        if nr in superstrings:
            continue
        
        # all rnas that match this rna
        matching = set((x,p) for (x,p) in [rnadict[subs]] if nr != x)
        if len(matching) == 0:
            continue
        
        for i in xrange(4, len(rna), 4):
            subs = rna[i:i+4]
            if len(subs) < 4:
                i -= 4 - len(subs)
                subs = rna[i:i+4]
            
            new = rnadict[subs]
            matching = set([(rnanr, j) for (rnanr, j) in matching if (rnanr, j+i) in new])
            if len(matching) == 0:
                break

        if len(matching) > 0:
            new = [rnanr for (rnanr,_) in matching] 
            superstrings.update(new)
                
            substringdict[nr] = matching
    
    # 3: compare all substrings
    results = {}
    for nr, rna in enumerate(rna):
        if nr in superstrings:
            results[rna] = [0 * len(rna)]
    
    for nr, rna in enumerate(rna):
        if nr not in superstrings:
            # substring
            for nranr, pos_to_feat in substringdict:
                for i in xrange(pos_to_feat-len(rna), pos_to_feat, 1):
                    results[rnas[rnanr]] += frequencies[rnanr]
    
    return results
                
                
        
    


def _filtersimple(filename, minlen=4):
    ''' filter out useless small lines'''
    
    rnas = []
    freqs = []
    with open(filename) as lines:       
        for line in lines:    
            if line[0] == ">":
                freq = line.split("-")[-1].strip()
            else:
                rna = line.strip()
                if len(rna) >= minlen:
                    rnas.append(rna)
                    freqs.append(freq)
    
    return rnas, freqs
#     
#     
#     print "length should be the same:"
#     print len(rnas), len(freqs)
#     print 
#     for rna, freq in zip(rnas, freqs):
#         print rna, freq
# 
#     print
#     print "length should be the same:"
#     print len(rnas), len(freqs)
#     print 
        
        
rnas, freqs = _filtersimple("SRR797062.collapsed", 4)

print "compare"
compare(rnas, freqs)



