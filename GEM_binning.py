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

data = data.fillna(0)

bins = int(round(global_max))
print(bins)

count = 0
for i in range(0, len(data)):
    gene1 = list(data.iloc[int(i)].values)
    only_pos_1 = [num for num in gene1 if num > 0]
    if len(only_pos_1) > 0:
        for j in range(i + 1, len(data)):
            gene2 = list(data.iloc[int(j)].values)
            only_pos_2 = [num for num in gene2 if num > 0]
            if len(only_pos_2) > 0:
                k = np.zeros(bins * bins)
                for l,m in zip(gene1,gene2):
                    if l > 0 and m > 0:
                        ind = int((round(l) * bins) + round(m))
                        #print(ind)
                        k[ind] = k[ind] + 1
                k = np.interp(k, (k.min(), k.max()), (0, +1))
                test = np.reshape(k, (bins,bins))
                plt.imshow(test)
                plt.show()
                print(i, j)
