# Gravitational Potemtial method for Network Optimization
This repository contatins code for the "Gravitational potential method and its application to network optimization" Bachelor Thesis.

It consists of five main parts:

* **main.py** and **optimization_for_diff_k.py** --- main modules, which perform the optimization for fixed number of servers and compare the results for different number of servers respectively.
* **flow.py** --- module, where the first part of the optimization is implemented.
* **mds.py** --- module, where the second part of the algorithm (finding the centroid node, which minimizes the paths' sum in the cluster) is implemented.
* **graph_models.py** --- module, where the requierd FlowGraph class is defined along with some models of graphs (balanced_tree, balanced_tree, random_powerlaw_tree. etc.)
* **utils_** --- directory with helper modeules 

## How to use

0. Install the required packages with `python install -r requirements.txt`

1. Create the desired graph model (instance of FlowGraph class, or select from already created in **graph_models.py**). 
The graph can be created from the edgelist as well as from the adjacency matrix (parameters ***edgelist***, ***adj_matrix***).

2. Pass the created graph to https://github.com/PhantomOfTheOpera/GravitationalNetworkOptim/blob/991c1ad4fd0fa279eab8e47f2ccf04d71b145e1c/main.py#L137 or https://github.com/PhantomOfTheOpera/GravitationalNetworkOptim/blob/991c1ad4fd0fa279eab8e47f2ccf04d71b145e1c/optimization_for_diff_k.py#L10 respectively.

3. Run the desired file with `pyhton main.py` or  `pyhton optimization_for_diff_k.py`.

## Interpretation of outputs
### Module optimization_for_diff_k.py

The output is a table generated in LaTex (just compile the output in the terminal as a LaTex file). The table looks like this:
<p align="center">
  <img src="./img/exp_example.png" alt="Example of output" width="650">
</p>
### main.py

After defining the graph and setting the number of servers (`k`) the file generates output in format:

