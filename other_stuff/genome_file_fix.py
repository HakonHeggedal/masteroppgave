'''
Created on 30. apr. 2015

@author: hakon
'''

print 123
with open("hsa_37.fa") as read_file:
    with open("hsa_37_fix.fa", "w") as write_file:
        for line in read_file:
            line = line.strip()
            if line[0] == ">":
                line = line.split(" ")[0]
                
            write_file.write(line + "\n")


print 4321