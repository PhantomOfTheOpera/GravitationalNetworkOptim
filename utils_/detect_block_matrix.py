import numpy as np
from numpy import linalg as LA
import copy
import graph_models
from itertools import product


def counts_disconnected_components(A: np.array, flw : list):
    
    flows = sorted(flw)

    for iteration, index in enumerate(flows):

        A = np.delete(A, index - iteration, 0)
        A = np.delete(A, index - iteration, 1)

    D = np.diag(np.sum(A, axis = 1))
    L = D - A
    w, v = LA.eigh(L)

    return w, v, len(np.where(np.isclose([0], w))[0])

def determine_sign_constant_vectors(matrix: np.array):

    indices = []
    for i in range(np.shape(matrix)[1]):

        if len(np.where(matrix[:, i] >= 0)[0]) == np.shape(matrix)[1] or len(np.where(matrix[:, i] <= 0)[0]) == np.shape(matrix)[1]:
            indices.append(i)

    return np.array(indices)

def get_indices_by_rows(A : np.array, rows : np.array):
    indices = []
    for row in rows:
        
        index = np.where(np.all(A == row, axis=1))
        indices.append(index[0][0])
    
    return indices

