import numpy as np
import spectralClustering as sc
import matplotlib.pyplot as plt
import networkx as nx
import pdb
from sklearn.preprocessing import normalize as nrmlz

# Debuggin and testing procedure for spectralClustering
#
#
#
#

data = None
clstr = None
truClass = None


def initialize(filename, nclusters=2, dims=2, eigengap=False,
               normalize=True, Fun=None):
    global data, clstr, truClass
    data = np.genfromtxt('../../DATASETS/'+filename, delimiter=' ', skiprows=1)
    dims = dims+1
    truClass = data[:, 0]
    data = data[:, 1:dims]

    if Fun is not None:
        data = Fun(data)
    if normalize:
        data = nrmlz(data, norm='l2', axis=0)

    clstr = sc.spCluster(data, nclusters, eigengap)


def refit(k, Draw=True):
    clstr.fitClusters(k)
    if Draw:
        plt.ion()
        plt.figure(1)
        plt.scatter(data[:, 0], data[:, 1], c=clstr.clusterIndices)
        plt.figure(2)
        drawAff()
        plt.draw()
    return clstr.clusterIndices


def drawAff():
    G = nx.from_numpy_matrix(clstr.affMatrix)
    nx.draw(G, dict(enumerate(clstr.data)), node_size=20)
