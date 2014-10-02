x = "gi|224589820|ref|NC_000008.10|"



genomeref = "gi|224589820|ref|NC_000008.10|"
print genomeref
genomenr = genomeref.split("|")[3]
print genomenr
name_nr = genomenr.split("_")[1]
print name_nr
nr = name_nr.split(".")[0]
print nr
print int(nr)