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

    def __init__(self, data):
        ro.r("source('Rcode/kNNutils.R')")
        ro.numpy2ri.acitvate()
        self.data = data
        self.kNN = ro.r['getKNearestNeighbors']

    def constructAdjMatrix(data, k): {}
