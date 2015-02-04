'''
Created on 4. feb. 2015

@author: hakon
'''

def remove_mirTrons(miRNA_set, filename):
    print "removing from miRNA set"
    with open(filename) as mirTrons:
        for special in mirTrons:
            special = ">" + special.strip()
            
            if special in miRNA_set:
                del(miRNA_set[special])
                print special, "is special"