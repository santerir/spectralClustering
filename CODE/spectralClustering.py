# -*- coding: utf-8 -*-

# IMPORTS
import rpy2.robjects as ro
import numpy as np
import scipy as sp


# This is an algorithm for spectral clustering.
# It implements kNN utilities from R, but work is bieng done to translate
# these to run natively on python
#
#
#
#
#
#
#
#
# Santeri Räisänen, 2015


class spCluster:

    def __init__(self, data, k):
        ro.r("source('Rcode/kNNutils.R')")
        ro.numpy2ri.acitvate()
        self.data = data
        self.dataDim = data.shape[2]
        self.dataSize = data.shape[1]
        self.k = k
        self.kNN = ro.r['getKNearestNeighbors']

# Constructs adjacency matrix of similarity graph
# using k-NN (or) and euclidian distance as weight.

    def constructLaplacian(self, k):

        # Find k-Nearest Neighbours (implemented in R)

        neighbours = np.matrix(self.kNN(self.data, np.ones(self.dataDim), k))
        aMatrix = np.zeros(self.dataShape, self.dataSize)

        # Populate the affinity matrix

        for i in xrange(1, self.dataSize):
            for j in xrange(1, k):
                aMatrix[i, neighbours[i, j]] = np.exp(-neighbours[i, j + k])
                aMatrix[neighbours[i, j], i] = np.exp(-neighbours[i, j + k])
        aMatrix = np.matrix(aMatrix)

        # construct Laplacian, non symmetric normalized (random walk)

        dMatrix = np.zeros(self.dataSize, self.dataSize)
        for i in xrange(i, self.dataSize):
            dMatrix[i, i] = 1/np.sum(aMatrix[i, :])
        dMatrix = np.matrix(dMatrix)

        self.Laplacian = np.identity(self.dataSize) - dMatrix * aMatrix
        
