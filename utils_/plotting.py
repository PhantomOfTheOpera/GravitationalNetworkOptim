import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from numpy import linalg as LA
import scipy as sp
from pylab import *


def plot_path(paths : list, G : nx.Graph, flow_nodes : list = [], positions : dict = {}, title : str ='Initial settings', name = '_'):

    edges = [(u, v) for (u, v, d) in G.edges(data=True)]
    pos = positions if positions else nx.spring_layout(G, k=7)
    # print([node for node in G])
    # print(f"flow_nodes : {flow_nodes}")
    color_map = ['blue' if node in flow_nodes else 'orange' for node in G]
    # positions for all nodes - seed for reproducibility
    plt.title(title, fontsize = 20)
    figure = plt.gcf() # get current figure
    figure.set_size_inches(16, 14)
    # nodes
    nx.draw_networkx_nodes(G, pos, node_size = 700, node_color = color_map)

    # edges
    nx.draw_networkx_edges(G, pos, edgelist = edges, width = 3, edge_color = "grey",  style = "dashed")

    for path, color in paths:
        #edgelist = path,
        nx.draw_networkx_edges(G, pos,  width = 3, alpha = 0.5, edge_color = color)

    labels = {node:node for node in G.nodes}
    # node labels
    nx.draw_networkx_labels(G, pos, font_size = 20, labels=labels, font_family = "sans-serif")
    # edge weight labels 

    # attrs = {}
    # for edge in edges:
    #     attrs[edge] = {"attr_name": f'{edge[0]} <--> {edge[1]}'}
    # nx.set_edge_attributes(G, attrs)
    # edge_labels = nx.get_edge_attributes(G, "attr_name")

    # nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=12, font_weight='bold')
    plt.savefig(f'exp_2/{name}.png', dpi=150)

    # plt.show()

def _plot_clusters_(k, labels : list, G : nx.Graph, flow_nodes : list = [], positions : dict = {}, title : str ='Initial settings', name = '_'):

    pos = positions
    fig = plt.figure()  
    plt.title(title, fontsize = 20)
    figure = plt.gcf() # get current figure
    figure.set_size_inches(16, 14)

    num_colors = len(set(labels))
    cmap = cm.get_cmap('Paired', num_colors)  
    diff_colors = [matplotlib.colors.rgb2hex(cmap(i)) for i in range(cmap.N)] 
    color_map = list(range(len(labels)))
    for i, label in enumerate(list(set(labels))):
        if label == 11:
            color = '#808080'
        else:
            color = diff_colors[i]
        for index in np.where(labels == label)[0]:
            color_map[index] = color

    # node_labels_ = {node : "F" if (5 * node[0]  + node[1]) in flow_nodes else f'' for node in G.nodes() }
    node_labels_ = {node : "F" if node in flow_nodes else f'' for node in G.nodes() }
    nx.draw(G, pos = pos, labels=node_labels_,  node_color=color_map, font_size = 35, font_family = "sans-serif")
    plt.savefig(f'exp/{name}.png', dpi=150)