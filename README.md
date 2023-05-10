# Gravitational Potemtial method for Network Optimization
This repository contatins code for the "Gravitational potential method and its application to network optimization" Bachelor Thesis.

It consists of five main parts:

* **main.py** and **optimization_for_diff_k.py** --- main modules, which perform the optimization for fixed number of servers and compare the results for different number of servers respectively.
* **flow.py** --- module, where the first part of the optimization is implemented.
* **mds.py** --- module, where the second part of the algorithm (finding the centroid node, which minimizes the paths' sum in the cluster) is implemented.
* **graph_models.py** --- module, where the requierd FlowGraph class is defined along with some models of graphs (balanced_tree, balanced_tree, random_powerlaw_tree. etc.)
* **utils_** --- directory with helper modeules 

## How to use

1. Create the desired graph model (instance of FlowGraph class, or select from already created in **graph_models.py**). 
The graph can be created from the edgelist as well as from the adjacency matrix (parameters ***edgelist***, ***adj_matrix***).

2. Pass the created graph to https://github.com/PhantomOfTheOpera/GravitationalNetworkOptim/blob/b02d6b1fb9040726bbac17bb2f079c9a87f98e9d/main.py#LL137C5-L137C32
