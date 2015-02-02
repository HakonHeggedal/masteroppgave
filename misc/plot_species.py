'''
Created on 21. jan. 2015

@author: hakon
'''
from matplotlib import pyplot
from matplotlib.ticker import ScalarFormatter



hc_distances = [13, 1, 12, 7, 0, 40, 6, 1, 2, 2, 10, 10, 0, 77, 10, 0, 4, 2, 18, 4, 7, 11, 46, 10, 4, 1, 3, 0, 6, 1, 0, 0, 5, 42, 2, 32, 2, 7, 1, 2, 7, 11, 6, 0, 1, 9, 0, 5, 10, 4, 10, 0, 3, 6, 5, 9, 22, 2, 3, 20, 6, 0, 23, 6, 0, 34, 24, 25, 31, 1, 0, 12, 4, 3, 27, 5, 34, 28, 27, 0, 3, 1, 2, 0, 2, 8, 1, 4, 37, 0, 3, 29, 1, 50, 10, 0, 10, 12, 10, 10, 2, 1, 1, 3, 2, 11, 3, 4, 1, 2, 11, 1, 2, 10, 9, 16, 5, 15, 3, 32, 3, 0, 0, 14, 1, 10, 2, 2, 0, 0, 6, 38, 5, 12, 6, 34, 10, 2, 0, 0, 7, 3, 17, 34, 0, 0]
hc_species = [4, 1, 3, 0, 2, 1, 4, 2, 2, 0, 5, 9, 0, 2, 1, 0, 0, 2, 3, 0, 2, 2, 7, 0, 8, 0, 1, 9, 6, 1, 0, 2, 0, 0, 1, 1, 1, 0, 5, 2, 7, 11, 9, 3, 2, 5, 6, 6, 8, 4, 4, 3, 9, 8, 11, 2, 0, 0, 4, 1, 0, 0, 1, 0, 0, 4, 0, 0, 0, 0, 0, 7, 2, 0, 0, 9, 1, 1, 0, 0, 2, 10, 2, 4, 0, 2, 0, 0, 0, 0, 3, 0, 0, 1, 3, 3, 2, 1, 3, 7, 1, 0, 0, 2, 5, 7, 0, 2, 3, 0, 0, 0, 2, 5, 4, 5, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 9, 7, 0, 0, 2, 0, 1, 7, 6, 0, 2, 1, 3, 3, 3, 0, 0, 0, 0, 0]


lc_distances = [94, 2, 54, 4, 3, 13, 21, 4, 2, 2, 153, 4, 0, 2, 4, 0, 107, 0, 47, 7, 12, 0, 123, 2, 12, 1, 1, 16, 20, 2, 75, 10, 71, 53, 17, 12, 45, 0, 1, 3, 54, 18, 26, 78, 1, 3, 2, 54, 12, 21, 24, 1, 11, 22, 10, 6, 24, 13, 1, 3, 20, 111, 6, 2, 2, 42, 115, 131, 20, 1, 79, 119, 124, 1, 177, 130, 124, 9, 208, 14, 9, 8, 2, 2, 157, 25, 4, 48, 6, 1, 5, 46, 17, 6, 153, 116]
lc_species = [0, 6, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 4, 0, 0, 0, 0, 3, 8, 0, 0, 5, 0, 0, 0, 7, 2, 0, 0, 0, 0, 0, 0, 0, 0, 8, 3, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 2, 13, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 3, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0]



print len(hc_distances),
print len(hc_species)
print len(lc_distances),
print len(lc_species)

max_dist1 = max(hc_distances)
max_dist2 = max(lc_distances)
max_dist = max(max_dist1, max_dist2)

max_spe1 = max(hc_species)
max_spe2 = max(lc_species)
max_spe = max(max_spe1, max_spe2)

doubles = []

i = 1

while i < 200:
    doubles.append(i)
    i <<= 1


pyplot.hist([hc_distances,lc_distances], 30, stacked=True)
pyplot.yscale("symlog",subsy=range(0,10))
pyplot.yticks(doubles)
pyplot.gca().set_xlim( -5, None )
pyplot.gca().yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
pyplot.grid(True)
pyplot.show()

# x_tics = range(-1, )
# 
# pyplot.plot(hc_distances, hc_species, "ro")
# pyplot.plot(lc_distances, lc_species, "go")
# pyplot.xscale("symlog")
# pyplot.yticks(range(-1, max_spe+1, 1))
# pyplot.xticks(range(-1, max_dist+1, 1))  #     (([-1, max_dist, -1, max_spe ]))
# pyplot.get
# 
# # pyplot.xscale("log")
# pyplot.grid(True)
# pyplot.show() 