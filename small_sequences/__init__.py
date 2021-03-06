
import matplotlib.pyplot as plot
import SuffixTree
# from SuffixTree import SuffixTree, STREE_DNA
import SuffixTree.SubstringDict
from SuffixTree import SubstringDict


def readcollapsed(filename):
    ''' returns a dict with sequence -> frequency for each sequence in fasta file '''
    sequences = []
    count = []
    
    with open(filename, "r") as inputfile:
        name = ""
        for line in inputfile:
            line = line.strip()
            if len(line) == 0:
                name = ""
            elif line[0] == ">":
                name = line.split("-")[1].strip()
            elif len(name) > 0 and len(line) > 10:
                sequences.append(line)
                count.append(name)
                name = ""

    return sequences, count



small_min = 12
small_max = 19
filename = "SRR797062.fa"
filename = "SRR207116.collapsed"

seqs, counts = readcollapsed(filename)
small_seqs = [x for x in seqs if len(x) >= small_min and len(x) < small_max]



print small_seqs[0:10]
print "nr of sequences:", len(seqs), len(counts)
print "nr of small seqs:", len(small_seqs)



# hash sequence --> count
seq_to_count = {}
i= 1
for seq, count in zip(seqs, counts):
    seq_to_count[seq] = int(count)
    i+= 1
    
# print "IIII", i

print seqs[0] in seq_to_count


seqs = [x for x in seqs if len(x) > small_max]
# create suffix Tree
print "adding to suffix tree"


suffixes = SuffixTree.SuffixTree()
# suffixes = SubstringDict()

print "ok?"
for i, seq in enumerate(seqs):
    suffixes.add(seq, i)
#     suffixes[seq] = i
#     print "!"
    
    if i%100000 == 0:
        print i, "of", len(seqs), i*100.0/len(seqs)
print "added to suffix tree"

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


print
print "analyzing small seqs"

for i, seq in enumerate(small_seqs):
    
    if i%1000 == 0:
        print
        print i, "of", len(small_seqs), i*100.0/len(small_seqs)
        print
    
    match_length, suffix_node, endpos = suffixes.match(seq)
#     print match_length, endpos, len(seq)
    
    
#     print "ok?"
#     nr = len(suffixes[seq]) # insanely slow...
    m = len(seq) is match_length
#     if len(seq) is not match_length:

        
        
    if len(seq) is match_length:
        subs += 1
        subs_scaled += seq_to_count[seq]
        
        len_sub[len(seq)] += 1
        len_sub_scaled[len(seq)] += seq_to_count[seq]
#         print "yes",

    else:
#     if nr <= 1:
        # not substring...
        nosubs += 1
        nosubs_scaled += seq_to_count[seq]
        len_nosub[len(seq)] += 1
        len_nosub_scaled[len(seq)] += seq_to_count[seq]
#         print "no",
    
    sh = ""
#     for s in seqs:
#         m_s = False
#         if seq in s and seq is not s:
# #             print "actually yes",
#             m_s = True
#             sh = s
#             break

#     if m_s and not m:
#         pass
# #         print "missing match:"
#         print match_length, endpos, len(seq)
#         print sh, seq
#         print
        
#         print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n\n\n!!!!!!!!!!!!"
#         
#     if m and not m_s:
#         print "---------------------\n\n\n"
        # nesten ikke lenger :)
        

tot_scaled = subs_scaled+nosubs_scaled


print "sub seq results"
print subs, nosubs, total_subs
print "percentage substrings",  subs*100.0 / total_subs

print "scaled values:"
print subs_scaled, nosubs_scaled, tot_scaled
print "% substring, scaled:", subs_scaled*100.0 / tot_scaled
 
print [(x,y) for x,y in enumerate(len_sub)]
 
relevant =  [(x*1.0/(x+y), e) for (e,(x,y)) in enumerate(zip(len_sub_scaled,len_nosub_scaled)) if x > 0 or y > 0]

print relevant
print
print "finished analysing small seqeunces"

print zip(*relevant)
[x,y] = [list(t) for t in zip(*relevant)]
plot.plot(y,x)

# plot.plot([1,2,3],[4,5,6])
plot.show()

# hits = 0
# alls = 0
# for i, sub in enumerate(small_seqs):
#     if i%1000 == 0:
#         print
#         print i, "of", len(small_seqs), i*100.0/len(small_seqs)
#      
#     for seq in seqs:
#         if sub in seq and sub != seq:
#             hits += 1
#             break
#  
# print "test", hits, hits*100.0 / len(small_seqs)


# 
# derpify pls
#         print "\t", len(seq), match_length, suffix_node.edgestr(), seq
#         print "pls?"
#         print suffix_node.num_children()
#         print suffix_node.num_leaves()
#         print suffix_node.children()
#         print "pls?"
# #         print suffix_node.find_child(seq[match_length+1])
#         break
# #         continue

