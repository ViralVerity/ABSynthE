import numpy as np
from collections import Counter

def calculate_topology_params(coalescent_tree):
    
    #colless and staircase measures
    colless = 0
    sum_ratios = []
    uneven = 0
    
    for nde in coalescent_tree.final_nodes:

        right = np.sum(i > nde.root_to_tip for i in coalescent_tree.sample_times)
        left = np.sum(i <= nde.root_to_tip for i in coalescent_tree.sample_times)

        difference = left - right

        colless += abs(difference)
        
        if right != left:
            uneven += 1

        if left<right:
            ratio = left/right
        elif right<left:
            ratio = right/left
        else:
            ratio = 1
        
    sum_ratios.append(ratio)
        
    staircase_1 = uneven/len(coalescent_tree.final_nodes)

    staircase_2 = np.mean(sum_ratios)

    #At the moment, this includes the root in each step calculation, which is correct as it should include the first branching
    sackin = np.sum(coalescent_tree.total_steps)
    
    #WD_ratio
    
    depths = Counter(coalescent_tree.total_steps)
    
    max_depth = max(depths)
    max_width = depths.most_common()[0][1]
    
    WD_ratio = max_width/max_depth
    
    #delta_w
    
    tup_list = []

    for k,v in depths.items():
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
    
    
    ##max_ladder and IL_nodes
    count_list = []

    for nde in coalescent_tree.final_nodes:

        ladder_count = 0

        go_up_ladder(nde, ladder_count, count_list)
        
    longest_ladder = max(count_list)
    
    ladder_likeness = longest_ladder/len(coalescent_tree.tips)
    
    return colless, sackin, WD_ratio, delta_w, longest_ladder, ladder_likeness, staircase_1, staircase_2

    
    
    
    
def go_up_ladder(nde, ladder_count, count_list):

    
    if (nde.new_children[0].type == "Ind" and nde.new_children[1].type == "Coal") or (nde.new_children[0].type == "Coal" and nde.new_children[1].type == "Ind"):
        
        ladder_count += 1
        
        for i in nde.new_children:
            if i.type == "Coal":
                go_up_ladder(i, ladder_count, count_list)
                
    else:
        
        count_list.append(ladder_count)
        
    return
        

