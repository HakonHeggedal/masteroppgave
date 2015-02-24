'''
Created on 23. feb. 2015

@author: hakon
'''

def _upper_DNA(seq):
    seq = seq.upper()
    return seq.replace("U", "T")


def get_hairpin(miRNA_file):

    miids = set()
    mi_id = None
    
    # find all IDs of all dead miRNAs
    with open(miRNA_file) as infile:
        for line in infile:
            if line.startswith("ID"):
                mi_id = line.split()[-1]
                if "hsa" in mi_id:
                    miids.add(mi_id)
    
    # brute force seach for miRNA hairpins
    hairpins = _brute_search_hairpin(miids)

    mature_seqs = _brute_search_mature( miids, "mature/mature", xrange(2,21,1))
#     m2 = _brute_search_mature( miids, "mature/maturestar", xrange(10,15,1)) # found none
    
    

    for key,item in mature_seqs.iteritems():
        item = list(item)[0]
        mature_seqs[key] = item
    
    def assert_contains((miid, mature)):
        assert mature in hairpins[miid]
#     
    map(assert_contains, mature_seqs.iteritems())
    

    
    def _print(x): print x
    print "dead mature seqs", len(mature_seqs)
    map(_print, mature_seqs.iteritems())
    print "dead hairpin seqs", len(hairpins)
    map(_print, hairpins.iteritems())
    
    print hairpins.viewkeys() - mature_seqs.viewkeys()
    
    return hairpins, mature_seqs
    


def _brute_search_hairpin(miRNA_ids):
    print "brute force seach for dead mirna"
    
    dead_to_hairpin = {}
    mi_id = None
    add_mi = False
    
    for i in xrange(2,21,1):
        name = "hairpins/hairpin" + str(i) + ".fa"
        
        with open(name) as infile:
#             print 123, name
            for line in infile:
                if line.startswith(">"):
                    mi_id = line[1:].split()[0]
                    add_mi = False
    #                 print mi_id
                    if mi_id in miRNA_ids:
#                         print "lol we got a hit 123", mi_id
                        if mi_id not in dead_to_hairpin:
                            add_mi = True
                elif add_mi:
                    part = line.strip()
                    print "\n123"
                    print part
                    part = _upper_DNA(part)
                    print part
                    if mi_id in dead_to_hairpin:
                        dead_to_hairpin[mi_id] += part
                    else:
                        dead_to_hairpin[mi_id] = part
#                     print "\t-", mi_id, dead_to_hairpin[mi_id]
                    
#     print len(miRNA_ids), len(dead_to_hairpin), miRNA_ids - dead_to_hairpin.viewkeys()
            
    return dead_to_hairpin


def _brute_search_mature( miRNA_ids2, file_prefix, filenrs): # =xrange(2,21,1)
    print "brute force seach for dead mirna mature seqs", file_prefix
    
    dead_to_sequence = {}
    
    add_mi = False
    
    for i in filenrs:
        name = file_prefix + str(i) + ".fa"
        
        with open(name) as infile:
            for line in infile:
                if line[0] == ">": add_mi = False
                
                if line.startswith(">hsa"):
                    line = line.strip()
                    
                    mi_id = line[1:].split()[0].lower()
                    if mi_id[-1] == "*":
                        mi_id = mi_id[:-1]
                    
                    if mi_id in miRNA_ids2:
                        add_mi = True
                    
                elif add_mi:
                    line = line.strip()
                    line = _upper_DNA(line)
                    if mi_id in dead_to_sequence:
                        dead_to_sequence[mi_id].add(line)
                    else:
                        dead_to_sequence[mi_id] = set([line])
#                     print "\t-", mi_id, dead_to_sequence[mi_id]
                    
    print len(miRNA_ids2), len(dead_to_sequence), miRNA_ids2 - dead_to_sequence.viewkeys()
            
    return dead_to_sequence

    
def _request_mirbase(name):
    import urllib2
    url = "http://www.mirbase.org/cgi-bin/get_seq.pl?acc=" + name
    try : response = urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        print e.code
    except urllib2.URLError as e:
        print e.reason
        print 123
    
    ans = response.readlines()
    if ans:
        return ans[2].strip()
    else:
        return ""


# print _request_mirbase("MI0005674")