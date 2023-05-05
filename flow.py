import numpy as np
from numpy import linalg as LA
import networkx as nx
from utils_.plotting import plot_path, _plot_clusters_
from copy import deepcopy
from graph_models import FlowGraph
import random
from itertools import combinations, product
from collections import defaultdict
from utils_.detect_block_matrix import counts_disconnected_components
from sklearn.cluster import KMeans

def construct_gravitaional_potential(A: np.array, flw : list):
    """_summary_

    Args:
        A (np.array): _description_
        flw (list): _description_

    Returns:
        _type_: _description_
    """
    
    w, v, disc_comp = counts_disconnected_components(A, flw)
    eigenvects = v[:, :disc_comp]
    d = np.linalg.norm(eigenvects, axis = 1)
    diag = np.diag(d)
    normed_eigenvects = np.linalg.inv(diag)@eigenvects
    
    kmeans = KMeans(n_clusters=disc_comp, random_state=0, n_init='auto').fit(normed_eigenvects)
    blocks_indexes = kmeans.labels_

    D = np.diag(np.sum(A, axis = 1))
   
    L = D - A
    flows = sorted(flw)
    for iteration, index in enumerate(flows):

        L = np.delete(L, index - iteration, 0)
        L = np.delete(L, index - iteration, 1)
        
    phi = np.zeros(len(L))

    for i in range(disc_comp):
        
        current_block_indeces = np.where(blocks_indexes == i)[0]
        # print(current_block_indeces)
        block_matrix = L[current_block_indeces, :][:, current_block_indeces]
        # print(block_matrix)
        w_block, v_block = LA.eigh(block_matrix)
        phi_ = np.abs(v_block[:, np.argmin(w_block[w_block > 0])])
        phi[current_block_indeces] = phi_

    for i in flows:
        phi = np.insert(phi, i, 0)
    # print(phi)
    blocks_length = list((len(blocks_indexes[np.where(blocks_indexes == i)]) for i in range(disc_comp)))
    return A, phi, blocks_length, blocks_indexes

def find_all_paths(A : np.array, current_vertex : int, phi : np.array, visited : list, path : list):
    """_summary_

    Args:
        A (np.array): _description_
        current_vertex (int): _description_
        phi (np.array): _description_
        visited (list): _description_
        path (list): _description_

    """
    
    # Mark the current node as visited and store in path
    visited[current_vertex]= True
    path.append(current_vertex)

    if phi[current_vertex] == 0:
        path.pop()
        yield path.copy()
        return 

        # If current vertex is not destination
        # Recur for all the vertices adjacent to this vertex
    adj_vertices = np.nonzero(A[current_vertex])
    # print(adj_vertices)
    try:
        phi_min_value = np.min(phi[adj_vertices])
        next_vertices = np.where(phi == phi_min_value)[0]
    except ValueError:
        next_vertices = np.empty([1])

    next_adjacent_vertices = np.intersect1d(adj_vertices, next_vertices)

    for i in next_adjacent_vertices:
        if not visited[i]:
            visited[i] = True
            if phi[i] == 0:
                path.append(i)

            yield from find_all_paths(A, i, phi, visited, path)
                
# Remove current vertex from path and mark it as unvisited
            visited[i]= False
            path.pop()
 
def optimization_step(A : np.array, phi : np.array, flow_vertices : list, search_size : int, non_flow_verices_number : int):
    """
    Args:
        A (np.array): _description_
        phi (np.array): _description_
    """
    ###
    # Для розрахування стокових вершин (станцій) для кожного клієнта


    vertex_stats = dict()
    non_flow_vert = np.where(phi != 0)[0]


    for vert in non_flow_vert:
        paths = list(find_all_paths(A, int(vert), phi, [False for i in range(search_size)], []))

        # for path in paths:
        #     from copy import deepcopy
        #     old_path = deepcopy(graph_painting[path[-1]])
        #     old_path.append(path)
        #     graph_painting[path[-1]] = deepcopy(old_path)

        # vertex_stats[vert] = {'final_vertices': tuple({path[-1] for path in paths}) , 'paths_len' : tuple(len(path) for path in paths)}
        if len(paths) != 0:
            vertex_stats[vert] = {'final_vertices': tuple({path[-1] for path in paths})}

    flow_stats, flow_distribution = _avg_busyness_(vertex_stats, phi)
    labels_list = np.ones(len(A))
    labels_list[list(flow_distribution.keys())] = 0

    for j in range(len(flow_distribution.keys())):
        dict_keys = list(flow_distribution.keys())
        dict_keys.pop(j)
        values = [flow_distribution[x] for x in dict_keys]
        nodes_union = set(set().union(*values))
        determined_nodes = list(flow_distribution[list(flow_distribution.keys())[j]] - nodes_union)
        labels_list[determined_nodes] = j + 2
        labels_list[list(flow_distribution.keys())[j]] = j + 2

    criteria = 0
    k = len(flow_vertices)
    n = non_flow_verices_number
    
    for b in flow_stats.values():
        criteria += (b/n - 1/k) ** 2
    business_metrics = np.sqrt(criteria) * np.sqrt(k) / (k - 1)
    business_metrics = round(business_metrics, 3)

    cluster_inter_metrics = _cluster_intersection_(flow_distribution, non_flow_verices_number)

    metrics_value = (business_metrics + cluster_inter_metrics) / 2

    return flow_stats, metrics_value, business_metrics, cluster_inter_metrics, labels_list

def one_step_permutation(graph : FlowGraph, rem_vert : list, j : int, visualization = False):

    A, phi, labels, blocks_indexes = construct_gravitaional_potential(graph.adjacency_matrix, rem_vert)
    ######
    value_ = labels.index(sorted(labels)[-1])
    value_indices = np.where(blocks_indexes == value_)
    # print(np.where(phi[value_indices] == np.max(phi[value_indices]))[0])
    ######
    # max_vertices = np.where(phi == np.max(phi))[0]
    max_vertices = np.where(phi[value_indices] == np.max(phi[value_indices]))[0]
    new_max = max_vertices[0]
    non_flow_vertex_number = graph.vertex_number - len(rem_vert)

    _, metrics, bus_metr, inter_metr, labels_list = optimization_step(A, phi, rem_vert, graph.vertex_number, non_flow_vertex_number)
    non_determined_nodes = np.where(labels_list == 1)
    labels_list[np.where(labels_list == 1)] = 0
    labels_list[non_determined_nodes] = max(labels_list) + 2


    if visualization:
        _plot_clusters_(j, labels_list, graph.graph, flow_nodes = rem_vert, positions = graph.pos,  title = f"Iteration {j}, metrics = {metrics}, business = {bus_metr}, intersection = {inter_metr}" , name = f"{j}")

    best_metr = np.inf
    busiest_flow = np.inf

    bus_metr = np.inf
    inter_metr = np.inf

    for flow, new_flow in product(rem_vert, max_vertices):
        to_remove = deepcopy(rem_vert)
        to_remove.remove(flow)
        to_remove += [new_flow]
        A_, phi_, _, __ = construct_gravitaional_potential(graph.adjacency_matrix, to_remove)
        _, metrics_, busi_, inters_, __  = optimization_step(A_, phi_, to_remove, graph.vertex_number, non_flow_vertex_number)
        if metrics_ <= best_metr:
            best_metr = metrics_
            busiest_flow = flow
            new_max = new_flow
            bus_metr = busi_
            inter_metr = inters_

    rem_vert.remove(busiest_flow)
    rem_vert += [new_max]
    
    return rem_vert, metrics

def _cluster_intersection_(flow_distribution : dict, non_flow_verices_number : int, method = 'avg'):
    """_summary_

    Args:
        flow_distribution (dict): _description_
    """

    pairs = list(combinations(flow_distribution.values(), 2))
    worst_case = 0
    intersections = []

    for x, y in pairs:
        joint_number = len(x.intersection(y))
        intersections.append(joint_number)
        if joint_number > worst_case:
            worst_case = joint_number
    if method == 'max':
        return round(worst_case / non_flow_verices_number, 3)
    # print((len(flow_distribution.values()) - 1) * (worst_case + 1))
    return round(np.sum((intersections - np.mean(intersections)) ** 2) ** 0.5 *(len(flow_distribution.values())) ** 0.5 / ((len(flow_distribution.values()) - 1) * (worst_case + 1))  , 3)

def _avg_busyness_(vertex_statistics : dict, phi : np.array):
    """

    Args:
        vertex_statistics (dict): _description_
        phi (np.array): _description_
    """

    flow_vert = np.where(phi == 0)[0]

    flow_stats = {vert: 0 for vert in flow_vert}
    flow_distribution =  defaultdict(set)

    for key in vertex_statistics:

        for flow in vertex_statistics[key]['final_vertices']:
            flow_stats[flow] += 1 / len(vertex_statistics[key]['final_vertices']) 
            flow_distribution[flow].add(key)
        

    for key in flow_stats:
        flow_stats[key] = np.abs(flow_stats[key])

    return flow_stats, flow_distribution

def check_loop(criterias_lst : list):

    if criterias_lst[-1] in criterias_lst[:-1]:
        return 1
    return 0
