
import matplotlib.pyplot as plot
import SuffixTree.SubstringDict


def readcollapsed(filename):
    ''' returns a dict with sequence -> frequency for each sequence in fasta file '''
    sequences = []
    count = []
    
    with open(filename, "r") as inputfile:
        name = ""
        for line in inputfile:
            line = line.strip()
            if line[0] == ">":
                name = line.split("-")[1].strip()
            elif len(name) > 0 and len(line) > 10:
                sequences.append(line)
                
                count.append(name)
                name = ""

    return sequences, count


small_min = 5
small_max = 26

# import sequences
seqs, counts = readcollapsed("SRR797062.fa")
small_seqs = [x for x in seqs if len(x) >= small_min and len(x) < small_max]

print small_seqs[0:10]

print len(seqs), len(counts)

# hash sequence --> count
seq_to_count = {}

i= 1
for seq, count in zip(seqs, counts):
    seq_to_count[seq] = int(count)
    i+= 1
    
# print "IIII", i


# print len(seq_to_count)
# for (k,v) in seq_to_count.iteritems():
#     print k,v

print seqs[0] in seq_to_count



# create suffix Tree
suffixes = SuffixTree.SubstringDict()
for i, seq in enumerate(seqs):
    suffixes[seq] = i


# STATS:

# % of small sequences being substrings



# superstring distribution for substrings (how many strings has a small seq as substring)
# length distribution for substrings of most frequent superstrings
# distribution of superstrings for small sequences.



subs = 0
nosubs = 0
subs_scaled = 0
nosubs_scaled = 0


total_subs = len(small_seqs)

len_sub = [0]*small_max
len_sub_scaled = [0]*small_max
len_nosub = [0]*small_max
len_nosub_scaled = [0]*small_max



for seq in small_seqs:
#     print seq
    nr = len(suffixes[seq])
    if nr > 2:
        subs += 1
        subs_scaled += seq_to_count[seq]
        
        len_sub[len(seq)] += 1
        len_sub_scaled[len(seq)] += seq_to_count[seq]
    else:
        nosubs += 1
        nosubs_scaled += seq_to_count[seq]
        len_nosub[len(seq)] += 1
        len_nosub_scaled[len(seq)] += seq_to_count[seq]

tot_scaled = subs_scaled+nosubs_scaled
print subs, nosubs, total_subs, subs*1.0 / total_subs
print subs_scaled, nosubs_scaled, tot_scaled, subs_scaled*1.0 / tot_scaled

print [(x,y) for x,y in enumerate(len_sub)]

print [(x*1.0/(x+y), e) for (e,(x,y)) in enumerate(zip(len_sub,len_nosub)) if x > 0 or y > 0]









# print len(seqs)
# print seqs[0]
# print len(counts)
# print int(counts[0])
# 
# chr
# 
# chars = set("AGCTN")
# 
# for i, seq in enumerate(seqs):
#     for char in seq:
#         
#         
#         
#         if char not in chars:
#             print i, seq, char
# 
# print seqs[35]
# print seqs[36]
# print seqs[37]