import scipy
import numpy as np
import pandas as pd

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.preprocessing import scale

from scipy.stats.stats import pearsonr 

from scipy.stats import invgamma 
from scipy.stats import beta
import matplotlib.pyplot as plt

import plotly
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.graph_objs import *


import os
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--Eigen",action= 'store_true',help= 'if given, multiplies PCs by respective eigenvalues across PCAs')

parser.add_argument("--Scale",action= 'store_true',help= 'if given, haplotype matrices are scaled prior to PCA.')

parser.add_argument("--range",default = '50,310',type= str,help = "min and max haplotype lengths to be tested, comma separated.")

parser.add_argument("--Npops",default = 10,type= int,help = "maximum number of pops to test at each step. must higher than 2.")

parser.add_argument("--ncomp",default = 5,type= int,help = "PCs to keep.")

parser.add_argument("--N",default = 100,type= int,help = "Number of haplotypes to be generated from each pop.")

parser.add_argument("--iter",default= 50,type= int,help= "Number of repeats per step")

args = parser.parse_args()

############ Functions

def pairwise_gen(x,y):
    miss= 0
    same= 0
    if len(x) != len(y):
        return 'vector lengths differ'
    else:
        for n in range(len(x)):
            if x[n] == y[n]:
                same += 1
        return 1 - same / (len(x) - miss)


def return_fsts2(freq_array):
    pops= range(freq_array.shape[0])
    H= {pop: [1-(freq_array[pop,x]**2 + (1 - freq_array[pop,x])**2) for x in range(freq_array.shape[1])] for pop in range(freq_array.shape[0])}
    Store= []

    for comb in it.combinations(H.keys(),2):
        P= [sum([freq_array[x,i] for x in comb]) / len(comb) for i in range(freq_array.shape[1])]
        HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
        per_locus_fst= [[(HT[x] - np.mean([H[p][x] for p in comb])) / HT[x],0][int(HT[x] == 0)] for x in range(len(P))]
        per_locus_fst= np.nan_to_num(per_locus_fst)
        Fst= np.mean(per_locus_fst)

        Store.append([comb,Fst])
    
    
    ### total fst:
    P= [sum([freq_array[x,i] for x in pops]) / len(pops) for i in range(freq_array.shape[1])]
    HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
    FST= np.mean([(HT[x] - np.mean([H[p][x] for p in pops])) / HT[x] for x in range(len(P))])
    
    return pd.DataFrame(Store,columns= ['pops','fst'])


def local_sampling_correct(data_now,n_comp):
    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(data_now)
    feats= pca.transform(data_now)
    
    N= 50
    bandwidth = estimate_bandwidth(feats, quantile=0.2)
    params = {'bandwidth': np.linspace(np.min(feats), np.max(feats),30)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    ## perform MeanShift clustering.
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False, min_bin_freq=5)
    ms.fit(feats)
    labels1 = ms.labels_
    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1))) if y != -1}

    ## Extract the KDE of each cluster identified by MS.
    Proxy_data= []

    for lab in label_select.keys():
        if len(label_select[lab]) < 3:
            continue
            
        Quanted_set= feats[label_select[lab],:]
        grid.fit(Quanted_set)

        kde = grid.best_estimator_
        Extract= kde.sample(N)
        Return= pca.inverse_transform(Extract)
        
        #Return= data_now[np.random.choice(label_select[lab],N),:]
        Proxy_data.extend(Return)
    
    Proxy_data= np.array(Proxy_data)
    pca2 = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(Proxy_data)
    var_comp= pca2.explained_variance_ratio_
    
    New_features= pca2.transform(data_now)# * var_comp
    return New_features, var_comp



def controled_fsts(vector_lib,Eigen,length_haps,Scale,N_pops,n_comp,Iter,N_sims):

    ## Population Sizes and labels
    bias_scheme= np.random.choice(range(25,200),N_pops,replace= False)
    unbiased_sheme= np.repeat(N_sims,N_pops)

    bias_labels= np.repeat(np.array([x for x in range(N_pops)]),bias_scheme)
    unbias_labels= np.repeat(np.array([x for x in range(N_pops)]),unbiased_sheme)

    ### store distances between centroids
    biased_pairwise= []
    unbiased_pairwise= []
    corrected_pairwise= []

    ### store PC projection:
    dist_PC_even= {x:[] for x in range(n_comp)}
    dist_PC_bias= {x:[] for x in range(n_comp)}
    dist_PC_corrected= {x:[] for x in range(n_comp)}
    
    ### store fsts
    fst_store= []

    ### store Pearson's r comparing gen_diffs and feature space diffs across scenarios
    biased_pears= []
    corrected_pears= []
    unbiased_pears= []

    ### triangular matrices extract.
    iu1= np.triu_indices(N_pops,1) # for centroid comparison

    iu_unbias= np.triu_indices(sum(unbiased_sheme),1)
    iu_bias= np.triu_indices(sum(bias_scheme),1)

    ### proceed.

    for rep in range(Iter):
        Pops= np.random.choice(vector_lib.shape[0],N_pops,replace= False)
        print('vectors selected: {}'.format(Pops))
        ########## FST

        freqs_selected= vector_lib[Pops,:length_haps]
        Pairwise= return_fsts2(freqs_selected)

        #fsts_compare = scale(Pairwise.fst)
        fsts_compare= Pairwise.fst
        fst_store.extend(fsts_compare)
        #########################################################
        ########### PCA ####################################
        #########################################################
        ############# unbiased sample

        #### generate data and perform PCA.
        data= []

        for k in range(N_pops):

            probs= vector_lib[Pops[k],:length_haps]

            m= unbiased_sheme[k]
            Haps= [[np.random.choice([1,0],p= [1-probs[x],probs[x]]) for x in range(length_haps)] for acc in range(m)]

            data.extend(Haps)

        data1= np.array(data)

        if Scale:
            data1= scale(data1)

        pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(data1)

        feat_unbias= pca.transform(data1)

        if Eigen:
            feat_unbias= feat_unbias * pca.explained_variance_ratio_

        ####### centroid comparison
        unbias_centroids= [np.mean(feat_unbias[[y for y in range(feat_unbias.shape[0]) if unbias_labels[y] == z],:],axis= 0) for z in range(N_pops)]
        unbias_centroids= np.array(unbias_centroids)

        unbias_pair_dist= pairwise_distances(unbias_centroids,metric= 'euclidean')
        unbias_pair_dist= unbias_pair_dist[iu1]

        #unbias_pair_dist= scale(unbias_pair_dist)
        unbiased_pairwise.extend(unbias_pair_dist)
        
        ## PC-wise centroid comparison
        for PC in range(unbias_centroids.shape[1]):
            unbias_PC_dist= pairwise_distances(unbias_centroids[:,PC].reshape(-1,1),metric= 'euclidean')
            unbias_PC_dist= unbias_PC_dist[iu1]
            #unbias_PC_dist= scale(unbias_PC_dist)
            dist_PC_even[PC].extend(unbias_PC_dist)

        #################################################
        ############## biased sample

        #### generate data and perform PCA
        data= []

        for k in range(N_pops):

            probs= vector_lib[Pops[k],:]

            m= bias_scheme[k]
            Haps= [[np.random.choice([1,0],p= [1-probs[x],probs[x]]) for x in range(length_haps)] for acc in range(m)]

            data.extend(Haps)

        data2= np.array(data)

        if Scale:
            data2= scale(data2)

        pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(data2)

        feat_bias= pca.transform(data2)

        if Eigen:
            feat_bias= feat_bias * pca.explained_variance_ratio_

        #### Centroid distances
        bias_centroids= [np.mean(feat_bias[[y for y in range(feat_bias.shape[0]) if bias_labels[y] == z],:],axis= 0) for z in range(N_pops)]
        bias_centroids= np.array(bias_centroids)

        bias_pair_dist= pairwise_distances(bias_centroids,metric= 'euclidean')
        bias_pair_dist= bias_pair_dist[iu1]
        #bias_pair_dist= scale(bias_pair_dist)
        biased_pairwise.extend(bias_pair_dist)
        
        ### PC-wise centroid comparison
        for PC in range(bias_centroids.shape[1]):
            bias_PC_dist= pairwise_distances(bias_centroids[:,PC].reshape(-1,1),metric= 'euclidean')
            bias_PC_dist= bias_PC_dist[iu1]
            #bias_PC_dist= scale(bias_PC_dist)
            dist_PC_bias[PC].extend(bias_PC_dist)

        ###############################################################"
        ################## bias correct
        ### perform MS correction on biased samples
        feat_correct,var_comp= local_sampling_correct(data2,n_comp)

        ### centroid Distances
        centroids= [np.mean(feat_correct[[y for y in range(feat_correct.shape[0]) if bias_labels[y] == z],:],axis= 0) for z in range(N_pops)]
        centroids= np.array(centroids)
        pair_dist= pairwise_distances(centroids,metric= 'euclidean')
        pair_dist= pair_dist[iu1]
        #pair_dist= scale(pair_dist)
        corrected_pairwise.extend(pair_dist)
        
        ### PC-wise centroid comparison
        for PC in range(centroids.shape[1]):
            corrected_PC_dist= pairwise_distances(centroids[:,PC].reshape(-1,1),metric= 'euclidean')
            corrected_PC_dist= corrected_PC_dist[iu1]
            #corrected_PC_dist= scale(corrected_PC_dist)
            dist_PC_corrected[PC].extend(corrected_PC_dist)
    
    
    t= np.array([
        unbiased_pairwise,
        biased_pairwise,
        corrected_pairwise
    ]).T
    
    ### Pearson's even vs. corrected
    corrected_PCs= [length_haps,N_pops]
    corrected_PCs.extend([pearsonr(dist_PC_even[x],dist_PC_corrected[x])[0] for x in dist_PC_even.keys()])
    
    ### Pearson's even vs. bias
    bias_PCs= [length_haps,N_pops]
    bias_PCs.extend([pearsonr(dist_PC_even[x],dist_PC_bias[x])[0] for x in dist_PC_even.keys()])
    
    ### Pearson's overall distances
    pearsons= [pearsonr(fst_store,t[:,x])[0] for x in range(t.shape[1])]
    Factors= [length_haps, N_pops]
    Factors.extend(pearsons)
    return Factors, bias_PCs, corrected_PCs




############ Frequency vectors

# Simulate frequency vectors. 
# We must first define the number of populations, the length of the haplotypes desired, and their respective population sizes

length_range= [int(x) for x in args.range.split(',')]

L= max(length_range)

import itertools as it
n= 25
steps_t= 15

# Vary a (beta distribution parameter).
a_range= np.linspace(1,2,steps_t)
a_set= [i for i in a_range for _ in range(n)]

# vary b.
b_range= np.linspace(0.1,.4,steps_t)
b_set= [i for i in b_range for _ in range(n)]

## length of haplotypes to extract.
L_set= [L] * n * 15


background= np.array([a_set,b_set,L_set]).T

vector_lib= []
for k in range(background.shape[0]):
    
    probs= beta.rvs(background[k,0], background[k,1], size=int(background[k,2]))
    probs[(probs > 1)]= 1
    
    
    vector_lib.append(probs)

vector_lib= np.array(vector_lib)


### PCA
## PCA on vectors simulated
n_comp = 3

pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(vector_lib)
features = pca.transform(vector_lib) * pca.explained_variance_ratio_


###################
###################



### Select frequency vectors and draw haplotypes.
Eigen = args.Eigen
Scale= args.Scale

print('Eigen= {}'.format(Eigen))
print('Scaled= {}'.format(Scale))

N_pops= args.Npops # Number of pops

n_comp= args.ncomp # components to keep following PCA

Iter= args.iter # repeats

N_sims= args.N # number of haplotypes to generate from each pop in the unbiased scenario.

#####
filename= 'Sampling_test_pearsons_' + str(args.ncomp) + 'PC_' + 'Eigen=' + str(Eigen) + '_Scaled=' + str(Scale) + '.txt'

print('to print to: ' + filename)

#####

length_step= 10

Pears= []

bias_PC_pears= []
corrected_PC_pears= []

for legos in range(length_range[0],length_range[1],length_step):
    
    for pop in range(3,N_pops):
        correlates, bias_pcs, corrected_pcs= controled_fsts(vector_lib,Eigen,legos,Scale,pop,n_comp,Iter,N_sims)
        
        print('Npops: {}; L: {}'.format(pop,legos))
        
        deviates= correlates
        
        bias_PC_pears.append(bias_pcs)
        corrected_PC_pears.append(corrected_pcs)
        Pears.append(deviates)


#os.makedirs(os.path.dirname(filename), exist_ok=True)

################################################################
############# PRINT ###########################################
################################################################

#### Genetic to feature space distances
Output= open(filename,'w')

Output.write('Size\tNpops\t' + '\t'.join(['unbiased','biased','corrected']))

Output.write('\n')

for liable in Pears:
    Output.write('\t'.join([str(x) for x in liable]))
    Output.write('\n')

Output.close()

#### Even to biased distances, PC-wise

filename= 'Sampling_BIAS_PC-wise_pears_' + str(args.ncomp) + 'PC_' + 'Eigen=' + str(Eigen) + '_Scaled=' + str(Scale) + '.txt'

Output= open(filename,'w')

Output.write('Size\tNpops\t' + '\t'.join(['PC' + str(x) for x in range(args.ncomp)]))

Output.write('\n')

for liable in bias_PC_pears:
    Output.write('\t'.join([str(x) for x in liable]))
    Output.write('\n')

Output.close()

#### Even to corrected distances, PC-wise

filename= 'Sampling_MScorr_PC-wise_pears_' + str(args.ncomp) + 'PC_' + 'Eigen=' + str(Eigen) + '_Scaled=' + str(Scale) + '.txt'
Output= open(filename,'w')


Output.write('Size\tNpops\t' + '\t'.join(['PC' + str(x) for x in range(args.ncomp)]))

Output.write('\n')

for liable in corrected_PC_pears:
    Output.write('\t'.join([str(x) for x in liable]))
    Output.write('\n')

Output.close()