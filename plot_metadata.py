import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
from random import randint
from matplotlib.lines import Line2D
from sklearn.cluster import KMeans
from scipy.spatial import distance
from sklearn.externals import joblib
from scipy.spatial import ConvexHull
from pylab import *


data1 = pd.read_csv('Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.GEM.txt', sep='\t')

print(data1.shape)

xls = pd.ExcelFile('41598_2018_26310_MOESM5_ESM.xlsx')
df = pd.read_excel(xls, 'AllClinicalAnnotations')

#df = pd.read_csv('sorted_annot_modified_copy.txt', sep='\t')

print(df.shape)


global_max = data1.max().max()
global_min = data1.min().min()

data = data1.iloc[:, 1:]

labels = []

'''

p = "LGG"
n = "THCA"
r = "OV"
m = "BLCA"
o = "GBM"

for index, row in df.iterrows():
    #print(column[:25])
    #strng = df[df['file_id'].str.contains(column[:25])]['cancersubtype'].item()
    #strng = df[df['cancersubtype']].item()
    strng = row['cancersubtype']
    if strng == p:
        labels.append('y')
    if strng == n:
        labels.append('g')
    if strng == r:
        labels.append('b')
    if strng == m:
        labels.append('r')
    if strng == o:
        labels.append('c')

#strng = df[df['file_id'].str.contains(column[:25])].iloc[:, -1].item()
#print(str(df['file_id'].str.contains(column[:25])))
#print(strng)

legend_elements = [Line2D([0], [0], marker='o', color='y', label='LGG',
                          markerfacecolor='y', markersize=2),
                   Line2D([0], [0], marker='o', color='g', label='THCA',
                          markerfacecolor='g', markersize=2),
                   Line2D([0], [0], marker='o', color='b', label='OV',
                          markerfacecolor='b', markersize=2),
                   Line2D([0], [0], marker='o', color='r', label='BLCA',
                          markerfacecolor='r', markersize=2),
                   Line2D([0], [0], marker='o', color='c', label='GBM',
                          markerfacecolor='c', markersize=2)]


a = "notreported"
b = "stagei"
c = "stageii"
d = "stageiii"
e = "stageiv"
f = "stageiva"
g = "Solid Tissue Normal"

for index, row in df.iterrows():
    #print(column[:25])
    strng = row['tumor_stage']
    #strng2 = df[df['file_id'].str.contains(column[:25])]['classification_of_tumor'].item()
    strng2 = row['classification_of_tumor']
    if strng == a:
        labels.append('m')
    if strng == b:
        if strng2 == g:
            labels.append('g')
        else:
            labels.append('r')
    if strng == c:
        labels.append('b')
    if strng == d:
        labels.append('y')
    if strng == e:
        labels.append('c')
    if strng == f:
        labels.append('c')

legend_elements = [Line2D([0], [0], marker='o', color='m', label='NotReported',
                          markerfacecolor='m', markersize=2),
                   Line2D([0], [0], marker='o', color='c', label='Stage iv',
                          markerfacecolor='c', markersize=2),
                   Line2D([0], [0], marker='o', color='y', label='Stage iii',
                          markerfacecolor='y', markersize=2),
                   Line2D([0], [0], marker='o', color='b', label='Stage ii',
                          markerfacecolor='b', markersize=2),
                   Line2D([0], [0], marker='o', color='g', label='Normal',
                          markerfacecolor='g', markersize=2),
                   Line2D([0], [0], marker='o', color='r', label='Stage i',
                              markerfacecolor='r', markersize=2)]
'''
a = "Solid Tissue Normal"
b = "Primary Tumor"
c = "Recurrent Tumor"
d = "Metastatic"


for index, row in df.iterrows():
    #print(column[:25])
    strng = row['classification_of_tumor']
    if strng == a:
        labels.append('b')
    if strng == b:
        labels.append('m')
    if strng == c:
        labels.append('c')
    if strng == d:
        labels.append('y')

legend_elements = [Line2D([0], [0], marker='o', color='b', label='Normal',
                          markerfacecolor='b', markersize=2),
                   Line2D([0], [0], marker='o', color='m', label='Primary',
                          markerfacecolor='m', markersize=2),
                   Line2D([0], [0], marker='o', color='c', label='Recurrent',
                          markerfacecolor='c', markersize=2),
                   Line2D([0], [0], marker='o', color='y', label='Metastatic',
                          markerfacecolor='y', markersize=2)]

'''
a = "male"
b = "female"



for index, row in df.iterrows():
    #print(column[:25])
    strng = row['gender']
    if strng == a:
        labels.append('b')
    if strng == b:
        labels.append('r')


legend_elements = [Line2D([0], [0], marker='o', color='b', label='Male',
                          markerfacecolor='b', markersize=2),
                   Line2D([0], [0], marker='o', color='r', label='Female',
                          markerfacecolor='r', markersize=2)]
                          
'''
clf = joblib.load('./models/new_model_32')

centroid = clf.cluster_centers_
num_clusters = 32


fig = plt.figure()
for i in range(0,num_clusters):
    #plt.xlim(0, global_max)
    #plt.ylim(0, global_max)
    plt.xlabel('RNA seq expression')
    plt.ylabel('RNA seq expression')
    plt.scatter(centroid[i][:2016], centroid[i][2016:], s=2, c = labels, alpha=0.2)
    plt.legend(handles=legend_elements, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=4, mode="expand", borderaxespad=0.)
    fig.savefig( str(num_clusters) + '_cancer_tissue_type.png', dpi = 600)

plt.close()


'''
data = data1.transpose()

cluster = []
groups = []
labels = []

colors = ['r','b','y','g','c']

for i in list(data.index):
    #print(i[0])
    if i[0] == 'B':
        cluster.append(1)
        labels.append('r')
    if i[0] == 'O':
        cluster.append(2)
        labels.append('b')
    if i[0] == 'L':
        cluster.append(3)
        labels.append('y')
    if i[0] == 'T':
        cluster.append(4)
        labels.append('g')
    if i[0] == 'G':
        cluster.append(5)
        labels.append('c')



global_max = data.max().max()
global_min = data.min().min()

legend_elements = [Line2D([0], [0], marker='o', color='r', label='B',
                          markerfacecolor='r', markersize=2),
                   Line2D([0], [0], marker='o', color='b', label='O',
                          markerfacecolor='b', markersize=2),
                   Line2D([0], [0], marker='o', color='y', label='L',
                          markerfacecolor='y', markersize=2),
                   Line2D([0], [0], marker='o', color='g', label='T',
                          markerfacecolor='g', markersize=2),
                   Line2D([0], [0], marker='o', color='c', label='G',
                          markerfacecolor='c', markersize=2)]

clf = joblib.load('./models/new_model_32')

centroid = clf.cluster_centers_
num_clusters = 32

for i in range(0,num_clusters):
    #plt.xlabel(str(list(data)[i]))
    #plt.ylabel(str(list(data)[j]))
    #ax1.plot(x, y)
    #ax1.scatter(centroid[i][:2016], centroid[i][2016:], s=2, c = labels, alpha=0.4)
    #ax1.set_aspect('equal', 'box')
    plt.scatter(centroid[i][:2016], centroid[i][2016:], s=2, c = labels, alpha=0.2)
    plt.legend(handles=legend_elements, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=4, mode="expand", borderaxespad=0.)
    #plt.plot(gene_data[hull.vertices,0], gene_data[hull.vertices,1], alpha=0.0000002)
#plt.savefig(str(i) + '_both.eps', bbox_inches='tight', dpi=300)
    plt.show()

#plt.close()
'''
