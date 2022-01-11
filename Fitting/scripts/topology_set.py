import numpy as np
from collections import Counter
from collections import defaultdict
import datetime as dt


def find_node_to_all_children(node, node_to_all_children):
    
    if len(node.node_children) > 0:
        for child in node.node_children:
            node_to_all_children[node].append(child)
            if child.type == "coalescent":
                find_node_to_all_children(child, node_to_all_children)
    else:
        node_to_all_children[node] = []
        
    if node.node_parent:
        node_to_all_children[node.node_parent].extend(node_to_all_children[node])
            
    return node_to_all_children

def calculate_topology_params(coalescent_tree):
    
    #colless and staircase measures    
    node_to_all_children = defaultdict(list)
    
    node_to_all_children = find_node_to_all_children(coalescent_tree.root, node_to_all_children)

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
            if left.type == "coalescent":
                for query in node_to_all_children[left]:
                    if query.type == "individual":
                        left_count += 1
            else:
                left_count += 1
                
            if right.type == "coalescent":
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
                for i in node_to_all_children[left]:
                    print(i.type)
                for i in node_to_all_children[right]:
                    print(i.type)
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
        
        # with open(f"test_newick_{now}.txt", 'w') as fw:
        #     fw.write(newick_string)

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
    ladder_dict = defaultdict(list)
    go_down_ladder(coalescent_tree.root, None, ladder_dict)

    max_ladder = max([len(i) for i in ladder_dict.values()])/len(coalescent_tree.tips)
    
    in_ladders = []
    for lst in ladder_dict.values():
        for node in lst:
            in_ladders.append(node)

    IL_nodes = len(in_ladders)/len(coalescent_tree.final_nodes)
    
    topology = [colless, sackin, WD_ratio, delta_w, max_ladder, IL_nodes, staircase_1, staircase_2]
    
    return topology


def go_down_ladder(node, ladder_key, ladder_dict):
        
    if node.type == "coalescent":
        child1 = node.node_children[0]
        child2 = node.node_children[1]

        if child1.type != child2.type:
            if not ladder_key:
                ladder_key = node
            ladder_dict[ladder_key].append(node)
        else:
            ladder_key = None

        for child in node.node_children:
            go_down_ladder(child, ladder_key, ladder_dict)

    else:
        pass 
    

