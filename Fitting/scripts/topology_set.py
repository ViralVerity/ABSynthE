import numpy as np
from collections import Counter
from collections import defaultdict
import datetime as dt

def calculate_topology_params(coalescent_tree, newick_string):
    
    #colless and staircase measures    
    node_to_all_children = defaultdict(list)
    
    for node in coalescent_tree.all_tips_nodes:
        
        if node.type == "coalescent":
            node_to_all_children[node].extend(node.node_children)
            for i in node.node_children:
                node_to_all_children[node].extend(node_to_all_children[i])
        else: #other one is an individual node and only included if the subtree contains a sample
            node_to_all_children[node] = []

    differences = []
    ratio_list = []
    uneven = 0
    for node in coalescent_tree.all_tips_nodes:
        if node.type == "coalescent":
            direct_children = node.node_children
            left = direct_children[0]
            right = direct_children[1]
            left_count = 0
            right_count = 0
            if node.type == "coalescent":
                for query in node_to_all_children[left]:
                    if query.type == "individual":
                        left_count += 1
            else:
                left_count += 1
                
            if node.type == "coalescent":
                for query in node_to_all_children[right]:
                    if query.type == "individual":
                        right_count += 1
            else:
                right_count += 1

            differences.append(abs(left_count - right_count))
            
            if left_count != right_count:
                uneven += 1
            
            if left_count == 0 and right_count == 0:
                print("zero on left and right count")
                print(node_to_all_children[left], node_to_all_children[right])
                print(left, right)
                print(len(coalescent_tree.all_tips_nodes))
            if left_count < right_count:
                ratio = left_count/right_count
            else:
                ratio = right_count/left_count

            ratio_list.append(ratio)
            
    colless = sum(differences)
    staircase_1 = uneven/len(coalescent_tree.final_nodes)
    staircase_2 = np.mean(ratio_list)

    if len(ratio_list) == 0:
        ind_count = 0
        for i in coalescent_tree.all_tips_nodes:
            if i.type == 'individual':
                ind_count += 1
        now = dt.datetime.now()
        
        with open(f"weird_run_{now}.txt", 'w') as fw:
            fw.write(f'{differences}\n')
            fw.write(f'len coalescent nodes: {len(coalescent_tree.final_nodes)}\n')
            fw.write(f'number of tips: {ind_count}\n')
            fw.write(f'total in all_tips_nodes: {len(coalescent_tree.all_tips_nodes)}\n')
            for i in coalescent_tree.all_tips_nodes:
                fw.write(f'{i.type}, children: {i.node_children}\n')
        
        with open(f"test_newick_{now}.txt", 'w') as fw:
            fw.write(newick_string)

    #Includes the root in each step calculation, which is correct as it should include the first branching
    sackin = np.sum(coalescent_tree.total_steps)
    
    #WD_ratio
    depths = Counter(coalescent_tree.total_steps)
    
    max_depth = max(depths)
    max_width = depths.most_common()[0][1]
    
    WD_ratio = max_width/max_depth
    
    #delta_w
    
    index = 0
    diffs = []
    for depth, count in depths.items():
        if index > 0:
            width_diff = abs(count - depths[index-1])
            diffs.append(width_diff)
        index += 1
    
    if len(diffs) > 0:
        delta_w = max(diffs)
    else:
        delta_w = np.nan
    
    ##max_ladder and IL_nodes
    count_list = []
    node_set = set()

    node_set = set()
    ladder_list = []
    for leaf in coalescent_tree.tips:
        go_up_ladder(coalescent_tree.root, leaf, node_set, [], ladder_list)

    max_ladder = max([len(i) for i in ladder_list])/len(coalescent_tree.tips)
    
    in_ladders = []
    for lst in ladder_list:
        for node in lst:
            if node.type == "coalescent":
                in_ladders.append(node)

    IL_nodes = len(in_ladders)/len(coalescent_tree.final_nodes)
    
    topology = [colless, sackin, WD_ratio, delta_w, max_ladder, IL_nodes, staircase_1, staircase_2]
    
    return topology


def go_up_ladder(root, node, node_set, ladder, ladder_list):
    
    if node == root:
        return
    
    sibling_nodes = [i for i in node.node_parent.node_children if i != node]
    if len(sibling_nodes) != 1:
        print(f'wrong len sibling nodes: {len(sibling_nodes)}')
        print(node.id, node.node_parent.id)
        for i in sibling_nodes:
            print(i, i.id)
            print(i.id)
    
    if node.type == "individual":
        if sibling_nodes[0].type != "individual":
            if node.node_parent not in node_set:
                ladder.append(node.node_parent)
                node_set.add(node.node_parent)
                go_up_ladder(root, node.node_parent, node_set, ladder, ladder_list)
            else:
                ladder_list.append(ladder)
                return
        else:
            ladder_list.append(ladder)
            return
    else:
        if sibling_nodes[0].type != "individual":
            if node.node_parent not in node_set:
                ladder.append(node.node_parent)
                go_up_ladder(root, node.node_parent,node_set, ladder, ladder_list)
                node_set.add(node.node_parent)
            else:
                ladder_list.append(ladder)
                return
        else:
            ladder_list.append(ladder)
            return
    

