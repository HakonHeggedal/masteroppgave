import time
# import draw_graph


def position_frequency(string_list, symbol_set):
    
    longest = max(len(x) for x in string_list)
    results = [[0 for _ in xrange(len(symbol_set))] for _ in xrange(longest)]
    symbol_pos = {}
    
    for node1, s in enumerate(symbol_set):
        symbol_pos[s] = node1
        
    for string in string_list:
        for pos, letter in enumerate(string):
            results[pos][symbol_pos[letter]] += 1
            
    return results

def symbol_frequency(string_list):
    result = {}
    
    for s in string_list:
        for letter in s:
            result[letter] = result[letter] + 1 if letter in result else 1
            
    return result


def hashcluster(string_list, sublen = 3, offset = 1, maxlength = 10, filename = None, min_overlap = 15):
    """counts the overlap frequency between substrings of all strings in the string_list
    returns a list of dicts for each string.
    optional: write a "node - node - overlaps" file
    
    sublen : lenght of substring to be matched against each other
    offset: allowed offset between substrings (relative to original string)
    maxlength: max length of string
    filename: write to file if parameter is given
    min_overlap: threshold minimum nr of similar overlaps for writing to file
    """
    longest = min( max(len(x) for x in string_list), maxlength)
    
    maxrange = longest - sublen + 1
    posdicts = [{} for _ in xrange(maxrange)] # storing substrings
    print posdicts, len(posdicts)
    
    # add all substrings to dictionary in range (position +- offset)
    for nr, string in enumerate(string_list):
        print nr
        
        for node1 in xrange(len(string) - sublen + 1):
            subs = string[node1:node1+3]
            
            for j in xrange(-offset, offset+1):
                if 0 <= node1+j < maxrange:
                    if subs in posdicts[node1+j]:
#                         print "subs ", subs, " at position ", node1+j
                        posdicts[node1+j][subs].add(nr)
                    else:
#                         print "new ", subs, " at position ", node1+j
                        posdicts[node1+j][subs] = set([nr])
    
    print "creating result structure:", (time.clock() - starttime)

#     result = [[0 for _ in xrange(len(string_list))] for _ in xrange(len(string_list))]
    
    print "quadratic space, may possibly be out of memory..."
    results = [{} for _ in xrange(len(string_list))]
     
    print "filling out results", (time.clock() - starttime)
    for d in posdicts:
        for s in d.itervalues():
            for p1 in s:
                for p2 in s:
                    if p1 is not p2:
                        if p2 in results[p1]:
                            results[p1][p2] += 1
                        else:
                            results[p1][p2] = 1
    if filename == None:
        return results
    
    # printing resutls to file: node1, node2, overlap
    with open(filename, "w") as my_out_file:
        for node1,node1_dict in enumerate(a):
            print node1
            print [(node2,v) for (node2,v) in node1_dict.iteritems() if v >= min_overlap]
            if node1 < node2:
                for node2,v in node1_dict.iteritems():
                    if v >= min_overlap:
                        s = str(node2) +" "+ str(node1) +" "+ str(v) + "\n"
                        print s
                        my_out_file.write( s )
                        links.append((node1,node2))
    

def high_confindence_file(all_file, high_file, filename="high_confidence.txt"):
    all_rna = all_file
    high_rna = set()
    
#     with open(all_file) as lines:
#         ishuman = False        
#         for line in lines:    
#             if line[0] == ">":
#                 ishuman = line.find(">hsa") >= 0
#             elif ishuman:
#                 all_rna.append(line.strip())

    
    with open(high_file) as lines:
        ishuman = False        
        for line in lines:    
            if line[0] == ">":
                ishuman = line.find(">hsa") >= 0
            elif ishuman:
                high_rna.add(line.strip())
    
    print "finding confident: ", len(all_rna), len(high_rna)
    with open(filename, "w") as hc_file:
        for node1, rna in enumerate(all_rna):
            s = str(node1) + " 1\n" if rna in high_rna else str(node1) + " 0\n"
            hc_file.write(s)
    
        
#
# testing:
#


starttime = time.clock()
print "starting"
file_name = "mature.fa"
print "\nanalysing ", file_name 


rnas = [node1_dict.strip() for node1_dict in open(file_name) if node1_dict[0] is not ">"]
inv_rnas = [node1_dict.strip()[::-1] for node1_dict in open(file_name) if node1_dict[0] is not ">"]

rna_human = []
is_human = False
for line in open(file_name):
    if line[0] == ">":
        if line.find("Homo sapiens") > 0:
            is_human = True
    elif is_human:
        rna_human.append(line.strip())
        is_human = False    
    else:
        pass
        
rnas = rna_human[:]
# rnas = [node1_dict[::-1] for node1_dict in rnas]

freq_dict = symbol_frequency(rnas)
symbols = freq_dict.keys()
positions = position_frequency(rnas, symbols)

print "number of rnas:", len(rnas)
print "symbol frequency: ", freq_dict
print "most frequent positions:"
print symbols

for pos, freq in enumerate(positions):
    print freq, pos

a = hashcluster(rnas)
 
links = []

with open("human.txt", "w") as my_out_file:
    for node1,node1_dict in enumerate(a):
        print node1
        print [(node2,val) for (node2,val) in node1_dict.iteritems() if val > 13]

        for node2,val in node1_dict.iteritems():
            if node2 > node1 and val > 13:
                s = str(node1) +" "+ str(node2) +" "+ str(val) + "\n"
                print s
                my_out_file.write( s )
                links.append((node1,node2))

print len(a)

high_confindence_file(rna_human, "high_conf_mature.fa")

print "\nfinished in ", (time.clock() - starttime) ,"seconds"




