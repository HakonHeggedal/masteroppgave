


# s = "012345"
# l = len(s)
# 
# offset = 0
# for x in xrange(-offset,offset+1):
#     print x

a = (1,2)
x = set([a])
b = (2,3)

x.add(b)

print x, len(x)

for el in x:
    print el

print 123