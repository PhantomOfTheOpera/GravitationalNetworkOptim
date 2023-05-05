import random
import numpy as np
import graph_models
from flow import one_step_permutation, check_loop
from utils_.build_table import build_La_Tex_table



if __name__ == "__main__":

    graph = graph_models.large_balanced_tree
    k = 14
    table_lists = [['k', 'Random\_start', 'I_{avg}', 'Iteration\_number']]
    
    for j in range(14, 30):
        current_stats = [j]
        num_iterations = 10

        random.seed(41)
        rem_vert = random.sample(range(0, 360), k=j)
        print(rem_vert)
        current_stats.append(rem_vert)

        iterations_list = []
        optimal_criteria_value = np.inf
        optimal_flows = rem_vert
        criterias_values = list()

        loops_number = 0
        
        for i in range(num_iterations):
            
            rem_vert_, criteria = one_step_permutation(graph, rem_vert, i, visualization=False)
            criterias_values.append(criteria) 

            if criteria < optimal_criteria_value:
                optimal_criteria_value = criteria
                optimal_flows = rem_vert_

            loops_number += check_loop(criterias_values)

            if loops_number > 2:
                break

        current_stats.append(optimal_criteria_value)
        iterations_list.append(i)

        current_stats.append(iterations_list)
        table_lists.append(current_stats)

    print(criterias_values)    
    print(table_lists)
    build_La_Tex_table(table_lists)