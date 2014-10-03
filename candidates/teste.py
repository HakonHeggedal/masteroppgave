import intervaltree
from intervaltree import bio

i = intervaltree.IntervalTree()

i.addi(10, 20, "test1")
i.addi(15, 20, "test2")
i.addi(10, 16, "test3")

print 15 in i
print i[8:12]
print i[15] - set([])

print i.containsi(6, 11, "test1")


