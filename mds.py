import numpy as np
import networkx as nx
from numpy import linalg as LA
from numpy.linalg import matrix_power
import copy
import matplotlib.pyplot as plt
from utils_.dijkstra import find_all_paths

def brute_path_optimization(A : np.array):

    vertices_num = np.shape(A)[0]
    current_path_avg = np.inf
    current_center = 0

    for i in range(vertices_num):
        lengths, paths = find_all_paths(wmat = A, start=i)

        if np.mean(lengths) < current_path_avg:
            current_path_avg = np.mean(lengths)
            current_center = i

    return current_path_avg, current_center

def optimize_paths(A : np.array, start_vertex = 0):
    
    current_start = start_vertex

    while True:
        # print(current_start)
        lengths, paths = find_all_paths(wmat = A, start=current_start)

        total_start_length = np.sum(lengths)

        adjacent_vertices = np.where(~np.isclose([0], A[current_start]))[0]

        adj_lengths_minimum = np.inf
        adj_vertex_minimum = start_vertex

        for vertex in adjacent_vertices:
            # print(vertex)

            new_lengths = copy.deepcopy(lengths)
            for i in range(len(paths)):
                if vertex in paths[i]:
                    new_lengths[i] -=  A[current_start][vertex]
                else:
                    new_lengths[i] += A[current_start][vertex]

            if np.sum(new_lengths) < adj_lengths_minimum:
                adj_lengths_minimum = np.sum(new_lengths)
                adj_vertex_minimum = vertex
        
        if adj_lengths_minimum < total_start_length:
            current_start = adj_vertex_minimum
        else:
            break

    return current_start, np.mean(find_all_paths(wmat = A, start=current_start)[0])

def find_closest_vertex(G : np.array):
    centroid =  np.mean(G, axis = 0)
    distances = np.sum((G - centroid)**2, axis=1)
    return np.argmin(distances)


def _construct_inner_product_matrix_(Delta : np.array):
    """_summary_

    Args:
        Delta (np.array): dissimiliarity matrix

    Returns:
        np.array: Inner product matrix (n \cross n)

    """
    n = np.shape(Delta)[0]
    B = np.zeros((n ,n))
    third_term = np.sum(np.square(Delta)) / n**2

    for i in range(n):
        second_term = np.sum(np.square(Delta[i, :])) / n

        for j in range(n):
            first_term = np.sum(np.square(Delta[:, j])) / n

            B[i][j] = Delta[i][j] - first_term - second_term + third_term
    
    B *= -0.5

    return B

def calculate_embedding_matrix(B : np.array, m = 2):
    """_summary_

    Args:
        B (np.array): inner product matrix
        n (int, optional): Dimension of the embedding space. Defaults to 2.

    Returns:
        np.array:  Points coordinates matrix  (n \cross m)

    """
    w, v = LA.eigh(B)
    top_m_eigenval, V = w[-m:], v[:, -m:]
    Lambda = np.diag(top_m_eigenval)

    return V @ np.sqrt(Lambda)

def approximate_2_edge_paths_lengths(A : np.array):

    """_summary_

    Args:
        A (np.array): (n \cross n) Adjacency matrix
        n (int, optional): paths distance to approximate. Defaults to 2.

    """
    adj = copy.deepcopy(A)
    non_weighted_A = np.zeros(np.shape(adj))
    non_weighted_A[np.nonzero(adj)] = 1

    lengts_approx  = adj @ non_weighted_A + non_weighted_A @ adj
    lengts_approx[np.diag_indices_from(adj)] = 0


    paths_num = matrix_power(non_weighted_A, 2)
    paths_num[np.diag_indices_from(adj)] = 0
    
    for i, j in zip(*np.nonzero(lengts_approx)):
        lengts_approx[i][j] = lengts_approx[i][j] / paths_num[i][j]

    for i, j in zip(*np.where(adj == 0)):
            adj[i][j] = lengts_approx[i][j]

    return adj

def approximate_paths_lengths(A : np.array, m : int = 2):
    """_summary_

    Args:
        A (np.array): _description_
        m (int, optional): Number of 2-edge-length paths approximation, each multiplication approximates paths depth by power of 2. Defaults to 2. 
    """

    approx_matrix = A
    for i in range(m):
        approx_matrix = approximate_2_edge_paths_lengths(approx_matrix)
    
    return approx_matrix

def calculate_centroid_vertex(A: np.array, m = 2):

    A /= np.max(A)
    approx =  approximate_paths_lengths(A)

    for i, j in zip(*np.where(approx == 0)):
        if i != j:
            approx[i][j] = np.max([np.max(approx[i, :]), np.max(approx[:, j])])

    B = _construct_inner_product_matrix_(approx)

    X = calculate_embedding_matrix(B, m=2)

    ver = find_closest_vertex(X)

    return ver


if __name__ == "__main__":

    A = np.array([[0, 1, 0, 0, 0],
                  [1, 0, 3, 0, 0],
                  [0, 3, 0, 5, 0],
                  [0, 0, 5, 0, 7],
                  [0, 0, 0, 7, 0]], dtype=float)


    A = np.array([[0, 0, 1, 0, 0, 0, 0],
                [0, 0, 2, 0, 0, 0, 0],
                [1, 2, 0, 10, 0, 0, 0],
                [0, 0, 10, 0, 3, 4, 5],
                [0, 0, 0, 3, 0, 0, 0],
                [0, 0, 0, 4, 0, 0, 0],
                [0, 0, 0, 5, 0, 0, 0],
                ], dtype=float)

    A = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [1, 0, 2, 30, 0, 0, 0, 40, 0, 0],
                [0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 30, 0, 0, 3, 4, 5, 50, 0, 0],
                [0, 0, 0, 3, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 4, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 5, 0, 0, 0, 0, 0, 0],
                [0, 40, 0, 50, 0, 0, 0, 0, 3, 2],
                [0, 0, 0, 0, 0, 0, 0, 3, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 2, 0, 0]], dtype=float)
    
    A = np.array([[0, 1, 1, 0, 12,0],
                  [1, 0, 0, 0, 0, 5],
                  [1, 0, 0, 5,0, 0],
                  [0, 0, 5,0, 1, 0],
                  [12,0, 0, 1, 0, 1],
                  [0, 5,0, 0, 1, 0]], dtype=float)
    

    # A = np.array([[0, 0, 0, 0, 1, 0, 0, 10, 0, 0],
    #               [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    #               [0, 0, 0, 0, 0, 5, 2, 0, 0, 0],
    #               [0, 0, 0, 0, 0, 0, 2, 0, 0, 0],
    #               [1, 1, 0, 0, 0, 0, 0, 0, 15, 0],
    #               [0, 0, 5, 0, 0, 0, 3, 0, 0, 8],
    #               [0, 0, 2, 2, 0, 3, 0, 0, 0, 9],
    #               [10, 0, 0, 0, 0, 0, 0, 0, 2, 0],
    #               [0, 0, 0, 0, 15, 0, 0, 2, 0, 4],
    #               [0, 0, 0, 0, 0, 8, 9, 0, 4, 0]], dtype=float)
    
    print(calculate_centroid_vertex(A))
    A /= np.max(A)
    
    approx =  approximate_paths_lengths(A, m=2)
    
    
    for i, j in zip(*np.where(approx == 0)):
        if i != j:
            approx[i][j] = np.max([np.max(approx[i, :]), np.max(approx[:, j])])

    print(approx) 
    B = _construct_inner_product_matrix_(approx)

    X = calculate_embedding_matrix(B, m=2)

    print(X)
    ver = find_closest_vertex(X)
    print(ver)
    
    graph = nx.from_numpy_matrix(A, parallel_edges=False, create_using=None)
    nx.draw(graph, pos = X)
    nx.draw_networkx_labels(graph, X, font_size = 20, font_family = "sans-serif")
    plt.show()