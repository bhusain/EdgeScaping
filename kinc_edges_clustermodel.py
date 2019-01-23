import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
from random import randint
from matplotlib.lines import Line2D
from sklearn.cluster import KMeans
from scipy.spatial import distance
from sklearn.externals import joblib

data = pd.read_csv('Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.GEM.txt',sep='\t')
data = data.transpose()

data_edge = pd.read_csv('Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.sc.mcs30.md5.mmINF.th0.860100.coexpnet.edges.txt',sep='\t')
data_edge = data_edge.transpose()

clf = joblib.load('new_model_32')

centroid = clf.cluster_centers_
num_clusters = 32

print("Here")

tdata = []
gp_count = np.zeros((32,), dtype=int)
#14908
for i in range(0, 14908):
    if i%100 ==0:
        print(i)
    glist = []
    g1 = list(data_edge.iloc[0,i])
    g1 = ''.join(g1)
    for j in range(0, 2016):
        glist.append(data[g1][j])
    g2 = list(data_edge.iloc[1,i])
    g2 =''.join(g2)
    for k in range(0,2016):
        glist.append(data[g2][k])
    glist = np.asarray(glist)
    np.nan_to_num(glist, copy=False)
    glist[glist < 0] = 0
    glist = np.ndarray.tolist(glist)
    tdata.append(glist)


print("Now Here")
label = clf.predict(tdata)

print("Now Now Here")

for i in label:
    gp_count[i] = gp_count[i] + 1

print(gp_count)






