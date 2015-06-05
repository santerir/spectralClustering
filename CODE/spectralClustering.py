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
        self.k = k
        self.kNN = ro.r['getKNearestNeighbors']

    def constructAdjMatrix(self, k):
        neighbours = self.kNN(self.data, np.ones(self.dataDim), k)
