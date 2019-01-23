import numpy as np
import pandas as pd
from random import randint
from sklearn.externals import joblib
import pickle
import time
from sklearn.metrics.cluster import adjusted_mutual_info_score


import numpy as np
import pandas as pd
from random import randint
from sklearn.externals import joblib
import pickle
import time
from sklearn.metrics.cluster import adjusted_mutual_info_score
from scipy.special import gamma,psi
from scipy import ndimage
from scipy.linalg import det
from numpy import pi

'''
data = pd.read_csv('Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.GEM.txt',sep='\t')
#data = data.transpose()

print(data.shape)

print(len(list(data.index)))

global_max = data.max().max()
global_min = data.min().min()

bins = round(global_max)
data = data.fillna(0)

print(bins)

from sklearn.externals import joblib


clf = joblib.load('./models/new_model_32')

#with open('max_genes.txt') as f:
#    lines = f.read().splitlines()

data1 = pd.DataFrame(data= None, columns=data.columns)
for i in range(0, len(lines)):
    data1.loc[lines[i]] = data.loc[lines[i]]

print(data1.shape)


data1.to_csv('pruned_GEM.txt', sep='\t')
'''

data = pd.read_csv('pruned_GEM.txt',sep='\t')
data = data.transpose()
data = data.fillna(0)
global_max = data.max().max()
global_min = data.min().min()


from sklearn.externals import joblib


clf = joblib.load('./models/new_model_32')

def mutual_information_2d(x, y, sigma=1, normalized=False):

    bins = (256, 256)

    jh = np.histogram2d(x, y, bins=bins)[0]

    # smooth the jh with a gaussian filter of given sigma
    ndimage.gaussian_filter(jh, sigma=sigma, mode='constant',
                                 output=jh)

    # compute marginal histograms
    jh = jh + EPS
    sh = np.sum(jh)
    jh = jh / sh
    s1 = np.sum(jh, axis=0).reshape((-1, jh.shape[0]))
    s2 = np.sum(jh, axis=1).reshape((jh.shape[1], -1))

    if normalized:
        mi = ((np.sum(s1 * np.log(s1)) + np.sum(s2 * np.log(s2)))
                / np.sum(jh * np.log(jh))) - 1
    else:
        mi = ( np.sum(jh * np.log(jh)) - np.sum(s1 * np.log(s1))
               - np.sum(s2 * np.log(s2)))

    return mi


g = open('distances_pruned.csv', 'a')
count = 0
count_name = 0
EPS = np.finfo(float).eps
for i in range(0,5870):
    name1 = data.ix[0, i]
    #print(name1)
    count_name = count_name + 1
    for j in range(i + 1, 5870):
        name2 = data.ix[0, j]
        t1data = []
        gene_data = data.ix[1:, i].tolist() + data.ix[1:, j].tolist()
        gene_data = np.asarray(gene_data)
        np.nan_to_num(gene_data, copy=False)
        #print(gene_data)
        #gene_data[gene_data < 0] = 0
        gene_data = np.ndarray.tolist(gene_data)
        t1data.append(gene_data)
        label1 = clf.predict(t1data)
        if label1 == 1 or label1 == 6 or label1 == 7 or label1 == 19 or label1 == 21 or label1 == 28:
            count = count + 1
            dist = np.linalg.norm(data.ix[1:, i] - data.ix[1:, j])
            #mi = adjusted_mutual_info_score(data.ix[1:, i], data.ix[1:, j])
            #mi = mutual_information_2d(data.ix[1:, i], data.ix[1:, j])
            #print(count, mi)
            g.write("%s,%s,%s\n" % (name1, name2, dist))
            if count%1000 ==0:
                print(count_name, count)

g.close()




'''

for i in range(0,70000):
    name = list(data)[i]
    if data.ix[:, i].max() > 0:
        gene_count = np.zeros(32)
        name = list(data)[i]
        all_genes_names.append(name)
        for p in range(0, 5000):
            j = randint(0, 70000)
            if data.ix[:,j].max() > 0 and i!=j:
                t1data = []
                gene_data = data.ix[:, i].tolist() + data.ix[:, j].tolist()
                gene_data = np.asarray(gene_data)
                # gene_pair = np.asarray(gene_data)
                np.nan_to_num(gene_data, copy=False)
                gene_data[gene_data < 0] = 0
                gene_data = np.ndarray.tolist(gene_data)
                t1data.append(gene_data)
                label1 = clf.predict(t1data)
                label = label1[0]
                gene_count[label] = gene_count[label] + 1

all_genes = []
all_genes_names = []
t0 = time.time()
f = open("gene_name.txt", 'w')
print("Here")
for i in range(0,70000):
    if i%1000 == 0:
        print(i)
    if data.ix[:, i].max() > 0:
        name = list(data)[i]
        f.write(str(name))
        f.write('\n')
f.close()
'''

