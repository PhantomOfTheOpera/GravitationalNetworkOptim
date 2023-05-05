from flow import construct_gravitaional_potential, find_all_paths
import numpy as np
import random
from itertools import combinations
from collections import defaultdict
from mds import calculate_centroid_vertex, optimize_paths, brute_path_optimization
from flow import one_step_permutation, check_loop, optimization_step
from graph_models import FlowGraph
from copy import deepcopy

import sys
np.set_printoptions(threshold=sys.maxsize)

def final_clusterisation(A : np.array, flw : list):

    A, phi, _, __ = construct_gravitaional_potential(A, flw)

    vertex_stats = dict()

    non_flow_vert, flow_vert = np.where(phi != 0)[0], np.where(phi == 0)[0]

    for vert in non_flow_vert:
        paths = list(find_all_paths(A, int(vert), phi, [False for i in range(len(A))], []))
        vertex_stats[vert] = paths[0][-1]
    
    # get all servers' clients
    flow_distribution =  defaultdict(set)
    for key in vertex_stats:
        flow_distribution[vertex_stats[key]].add(key)
    
    cluster_indices = list()
    for server in flow_vert:
        cluster_indices.append([server] + list(flow_distribution[server]))

    return cluster_indices

def find_cluster_center(A : np.array, ind : list):
    """_summary_

    Args:
        A (np.array): _description_
        inds (list): _description_
    """
    if len(ind) == 1:
        return ind[0], set()
    
    A = A[ind, :][:, ind]
    approx_centr = calculate_centroid_vertex(A)

    real_center, avg_length = optimize_paths(A, start_vertex=approx_centr)

    return (ind[real_center], set(ind) - {ind[real_center]}, avg_length)

def optimize_server_location(A : np.array, groups_ind : list, compare = False):
    d_ = dict()

    for group in groups_ind:
        
        server, clients, avg_len = find_cluster_center(A, group)

        if compare:
            A_ = A[group, :][:, group]
            value = {"approximated" : round(avg_len, 3), "actual" : round(brute_path_optimization(A_)[0], 3)}
        else:
            value = {"approximated" : round(avg_len, 3)}
        d_[server] = [clients, value]

    return d_


def cluster_construction(graph: FlowGraph, A : np.array, k : int, num_iterations : int = 10):
        

    random.seed(2)
    rem_vert = random.sample(range(0,  A.shape[0]), k) 
    rem_vert = [28, 29, 30, 31, 32, 33]
    _,criteria_ = one_step_permutation(graph, deepcopy(rem_vert), 0, visualization=False)
    print(f"Initial Server locations: {rem_vert}, Initial criterion: {criteria_}")

    optimal_criteria_value = np.inf
    optimal_flows = rem_vert
    criterias_values = list()

    loops_number = 0
    
    for i in range(num_iterations):
        rem_vert_, criteria = one_step_permutation(graph, rem_vert, i, visualization=False)
        print(f"Server Nodes: {rem_vert_}, Optimization criterion: {criteria}")
        criterias_values.append(criteria) 

        if criteria < optimal_criteria_value:
            optimal_criteria_value = criteria
            optimal_flows = rem_vert_

        loops_number += check_loop(criterias_values)

        if loops_number > 2:
            break
    
    return optimal_flows


def brute_force_optimization(graph, A : np.array, k : int):

    possible_comb = list(combinations(list(range(0,  A.shape[0])), k))

    current_comb = possible_comb[0]
    current_value = np.inf

    for rem_vert in possible_comb:

        A, phi, labels, blocks_indexes = construct_gravitaional_potential(graph.adjacency_matrix, rem_vert)


        non_flow_vertex_number = graph.vertex_number - len(rem_vert)

        _, value, __, ___, ____ = optimization_step(A, phi, rem_vert, graph.vertex_number, non_flow_vertex_number)
        
        if value < 0.027:
            print(value)
            print(rem_vert)

        if value < current_value:
            current_value = value
            current_comb = rem_vert

    return current_value, current_comb


if __name__ == "__main__":
    from graph_models import unbalanced_3_2_1_tree, unbalanced_r_graph, random_powerlaw_tree, tutte_graph, chordal_cycle_graph, balanced_tree, weighted_balanced_tree, large_balanced_tree

    graph = unbalanced_3_2_1_tree
    # for unbalanced_3_2_1_tree --- good start [28, 29, 30, 31, 32, 33]

    adjacency_matrix = graph.adjacency_matrix
    # print(repr(adjacency_matrix))
    from utils_.detect_block_matrix import counts_disconnected_components
    print(counts_disconnected_components(adjacency_matrix, [])[-1])
    k = 3
    final_servers = cluster_construction(graph, adjacency_matrix, k)
    groups = final_clusterisation(adjacency_matrix, final_servers)
    print(groups)
    print(optimize_server_location(adjacency_matrix, groups, compare=True))

    # print(brute_force_optimization(graph, adjacency_matrix, k))
