'''
Created on 23. feb. 2015

@author: hakon
'''
import urllib2


def brute_search_hairpin(miRNA_ids):
    print "brute force seach for dead mirna"
    
    dead_to_hairpin = {}
    mi_id = None
    
    for i in xrange(2,21,1):
        name = "hairpins/hairpin" + str(i) + ".fa"
        print "\n", name
        
        with open(name) as infile:
#             print 123, name
            for line in infile:
                if line.startswith(">"):
                    mi_id = line.split()[1]
                    add_pls = False
    #                 print mi_id
                    if mi_id in miRNA_ids:
#                         print "lol we got a hit 123", mi_id
                        if mi_id not in dead_to_hairpin:
                            add_pls = True
                elif add_pls:
                    part = line.strip()
                    if mi_id in dead_to_hairpin:
                        dead_to_hairpin[mi_id] += part
                    else:
                        dead_to_hairpin[mi_id] = part
                    print "\t-", mi_id, dead_to_hairpin[mi_id]
                    
    print len(miRNA_ids), len(dead_to_hairpin), miRNA_ids - dead_to_hairpin.viewkeys()
            
    return dead_to_hairpin


def brute_search_mature(miRNA_ids, miRNA_id2, file_prefix="mature/mature", filenrs=xrange(2,21,1)):
    print "brute force seach for dead mirna mature seqs", file_prefix
    
    dead_to_sequence = {}
    mi_id = None
    
    for i in filenrs:
        name = file_prefix + str(i) + ".fa"
        print "\n", name
        
        with open(name) as infile:
#             print 123, name
            for line in infile:
                if line.startswith(">"):

                    if len(line.split()) < 2:
                        print line
                    mi_id = line.split()[1]
                    add_pls = False
    #                 print mi_id
                    if mi_id in miRNA_ids:
#                         print "lol we got a hit 123", mi_id
                        if mi_id not in dead_to_sequence:
                            add_pls = True
                elif add_pls:
                    line = line.strip()
                    if mi_id in dead_to_sequence:
                        dead_to_sequence[mi_id].append(line)
                    else:
                        dead_to_sequence[mi_id] = [line]
                    print "\t-", mi_id, dead_to_sequence[mi_id]
                    
    print len(miRNA_ids), len(dead_to_sequence), miRNA_ids - dead_to_sequence.viewkeys()
            
    return dead_to_sequence





# brute_search(None)

def get_hairpin(miRNA_file):
    from multiprocessing import Pool
    
    ishuman = False
    hairpin_list = []
    minames = set()
    total = 0
    name = ""
    
    with open(miRNA_file) as infile:
        
        for line in infile:
            if line.startswith("ID"):
                ishuman = "hsa" in line
            elif ishuman and line.startswith("FW"):
                name = line.split()[-1]
                minames.add(name)
                
    print minames
    
    hairpins = brute_search_hairpin(minames)
    m1 = brute_search_mature(minames, "mature/mature", xrange(2,21,1))
    m2 = brute_search_mature(minames, "mature/maturestar", xrange(10,15,1))
    
#     request_count = len(minames)
#     pool = Pool(request_count)
#     res = pool.map(_request_mirbase, list(minames))
#     res = [x for x in res if x != ""]
#     print res
    
#     hairpin = _request_mirbase(name)
#     print name, hairpin
#     total += 1
#     if hairpin:
#         hairpin_list.append(hairpin)
#         minames.append(name)

                    
#     print total 
#     print "\nprinting all found"
#     for name, hp in zip(minames, hairpin_list):
#         print name, hp
#         
#     print len(minames)
    
    
def _request_mirbase(name):
    
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