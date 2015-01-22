'''
Created on 20. jan. 2015

@author: hakon
'''
def read_high_confidence(hc_file):
    
    high_confidence_miRNA = set()

    
    with open(hc_file) as infile:

        for line in infile:
            if line.startswith(">hsa"):
                mi_name = line.split(" ")[0]                
                high_confidence_miRNA.add(mi_name)
                
    return high_confidence_miRNA
                
                
                
def read_family(family_file):
    
    mirna_to_family = {}
    family_name = -1
    with open(family_file) as infile:
        
        for line in infile:
            line = line.strip()
            
            if line.startswith("ID"):
                family_name = line.split()[1]
                continue
            parts = line.split()
            if len(parts) == 3:
                mi_name = parts[2]
                if mi_name.startswith("hsa"):
                    mi_name = ">"+mi_name
                    mirna_to_family[mi_name] = family_name

                
    return mirna_to_family