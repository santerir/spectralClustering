# -*- coding: utf-8 -*-

# IMPORTS
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as numpy2ri
import numpy as np
import scipy as sp
import sklearn.cluster as cl


# This is an algorithm for spectral clustering.
# It implements kNN utilities from R, but work is being done to translate
# these to run natively on python
#
#
#
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

    def __init__(self, data, nclusters=2, eigengap=False):
        ro.r("source('Rcode/kNNutils.R')")
        numpy2ri.activate()
        self.data = data
        self.eigengap = eigengap
        self.dataDim = data.shape[1]
        self.dataSize = data.shape[0]
        self.nclusters = nclusters
        self.kNN = ro.r['getKNearestNeighbors']

# Constructs adjacency matrix of similarity graph
# using k-NN (or) and euclidian distance as weight.

    def constructLaplacian(self, k):

        # Find k-Nearest Neighbours (implemented in R)

        neighbours = np.matrix(self.kNN(self.data, k))
        # Note that R uses indexing starting at 1
        neighbours[:, :k] -= 1
        tpl = (self.dataSize, self.dataSize)
        aMatrix = np.zeros(shape=tpl)

        # Populate the affinity matrix

        for i in xrange(0, self.dataSize):
            for j in xrange(0, k):
                aMatrix[i, neighbours[i, j]] =\
                    np.exp(-np.power(neighbours[i, j + k], 2))
                aMatrix[neighbours[i, j], i] =\
                    np.exp(-np.power(neighbours[i, j + k], 2))
        self.affMatrix = np.matrix(aMatrix)

        # Construct graph laplacian, unnormalized

        dMatrix = np.zeros(shape=tpl)
        for i in xrange(0, self.dataSize):
            dMatrix[i, i] = np.sum(aMatrix[i, :])
        dMatrix = np.matrix(dMatrix)

        self.Laplacian = dMatrix - aMatrix
        self.dMatrix = dMatrix

# Find k first  eigenvalues

    def findEigenvalues(self):
        if self.eigengap:
            w, v = sp.linalg.eigh(self.Laplacian, b=self.dMatrix,
                                  eigvals_only=False,
                                  turbo=True)
            gap = np.argmax(np.log(np.abs(w[:-1]-w[1:])))
            self.allEigValues = w
            v = v[:, :gap]
            w = w[:gap]

        else:
            w, v = sp.linalg.eigh(self.Laplacian, b=self.dMatrix,
                                  eigvals_only=False,
                                  eigvals=tuple([0, self.nclusters-1]))

        self.eigValues = w
        self.eigVectors = v

# Cluster eigenvectors by k-means

    def cluster(self):
        clstr = cl.KMeans(self.eigVectors.shape[1])
        self.clusterIndices = clstr.fit_predict(self.eigVectors)

    def fitClusters(self, k):
        self.constructLaplacian(k)
        print("Laplacian: Complete")
        self.findEigenvalues()
        print("Eigenvalues: Complete")
        self.cluster()

    def setData(self, data):
        self.data = data
