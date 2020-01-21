import numpy as np
from collections import Counter

def calculate_topology_params(coalescent_tree):
    
    
    #colless
    colless = 0

    for nde in coalescent_tree.final_nodes:

        right = np.sum(i > nde.root_to_tip for i in coalescent_tree.sample_times)
        left = np.sum(i <= nde.root_to_tip for i in coalescent_tree.sample_times)


        difference = left - right

        colless += abs(difference)
        
    #At the moment, this includes the root in each step calculation  
    sackin = np.sum(coalescent_tree.total_steps)
    
    #WD_ratio
    
    depths = Counter(coalescent_tree.total_steps)
    
    max_depth = max(depths)
    max_width = depths.most_common()[0][1]
    
    WD_ratio = max_width/max_depth
    
    #delta_w
    
    tup_list = []

    for k,v in depth.items():
        tup = (k,v)
        tup_list.append(tup)

    sorted_depths = sorted(tup_list, key=lambda tup:tup[0])
    
    difference = 0

    for index, tup in enumerate(sorted_depths):
        if index > 0:
            new_difference =  abs(tup[1] - sorted_depths[index-1][1])
            if new_difference > difference:
                difference = new_difference
                
    delta_w = difference
    
    

