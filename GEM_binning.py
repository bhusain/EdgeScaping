# Preprocessing Script
# Benafsh Husain
# 2/16/18
# Description:  Program reads in GEM, deletes points <= 0, and populates
# an array of size 361 with a condensed version of the data

import pandas as pd
import matplotlib.pyplot as plt
import math
from random import randint


import numpy as np
import json

d = dict()

data = pd.read_csv('Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.GEM.txt',sep='\t')
data = data.transpose()

print(data.shape)

print(len(list(data.index)))

global_max = data.max().max()
global_min = data.min().min()

bins = round(global_max)
data = data.fillna(0)

print(bins)

tdata = []

total_data = np.zeros(20 * 20)


for n in range(0,10):
    i = randint(0, 70000)
    j = randint(0, 70000)
    if n%1000 == 0:
        print(n)
    if data.ix[:,i].max() > 0:
        if data.ix[:,j].max() > 0:
            k = np.zeros(20 * 20)
            for l,m in zip(data.ix[:,i],data.ix[:,j]) :
                if l > 3 and m > 3:
                    ind = (round(l) * 20) + round(m)
                    k[ind] = k[ind] + 1
            tdata.append(k)
            total_data = total_data + k

np.savetxt('test_prune.out', tdata, delimiter='\t')

np.savetxt('test_total_prune.out', total_data, delimiter='\t')
