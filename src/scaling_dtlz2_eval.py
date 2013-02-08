# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 10:05:03 2012

@author: marcusm
"""

import cPickle
import matplotlib.pyplot as plt
import multibench as mb
#from matplotlib.pylab import *

#for i in xrange(3,11):
for i in [3]:
    with open('scale_500') as pfile:
        spamA = cPickle.load(pfile)

    with open('scale_100') as pfile:
        spamB = cPickle.load(pfile)

    #boxplot for individual size 100
    statsA, statsB = spamA[0][1],spamB[0][1]
    res = 1
    maxstep = 30

    resultsA = zip(*[stat.ypsilon[:maxstep][::res] for stat in statsA])
    resultsB = zip(*[stat.ypsilon[:maxstep][::res] for stat in statsB])

    plt.figure()
    bpa = plt.boxplot(resultsA)
    bpb = plt.boxplot(resultsB)
    plt.setp(bpa['boxes'], color='black')
    plt.setp(bpa['fliers'], color='black')
    plt.setp(bpa['caps'], color='black')
    plt.setp(bpb['boxes'], color='green')
    plt.setp(bpb['fliers'], color='green')
    plt.setp(bpb['caps'], color='green')

    #figure()
    plt.title('fitness dimension %i' % i)
    plt.xticks(range(len(statsA[0].steps[:maxstep][::res])), statsA[0].steps[:maxstep][::res])
    plt.xlabel('evolution steps')
    plt.ylabel('convergence')
    plt.show()

