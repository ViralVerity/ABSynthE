#take whole tree as input

import numpy as np
import datetime as dt

def calculate_branch_statistics(coalescent_tree):
    
    max_H = coalescent_tree.most_recent_date
    min_H = coalescent_tree.oldest_sample_date

    mean_lengths = np.mean(coalescent_tree.b_len_list)
    median_lengths = np.median(coalescent_tree.b_len_list)
    
    if len(coalescent_tree.b_len_list) >= 2:
        var_lengths = np.var(coalescent_tree.b_len_list)
    else:
        var_lengths = np.nan
        
    mean_internal = np.mean(coalescent_tree.internal_branches)
    median_internal = np.median(coalescent_tree.internal_branches)
    
    if len(coalescent_tree.internal_branches) >= 2:
        var_internal = np.var(coalescent_tree.internal_branches)
    else:
        var_internal = np.nan
        
    mean_external = np.mean(coalescent_tree.external_branches)
    median_external = np.median(coalescent_tree.external_branches)
    
    if len(coalescent_tree.external_branches) >= 2:
        var_external = np.var(coalescent_tree.external_branches)
    else:
        var_external = np.nan
        
    mean_ratio = mean_internal/mean_external
    median_ratio = median_internal/median_external
    
    if var_internal and var_external:
        var_ratio = var_internal/var_external
    else:
        var_ratio = np.nan
    
    branch_stats = [mean_lengths, median_lengths, var_lengths, mean_external, median_external, var_external, mean_internal, median_internal, var_internal, mean_ratio, median_ratio, var_ratio]
    

    if len(coalescent_tree.b_len_list) == 0 or len(coalescent_tree.external_branches) == 0 or len(coalescent_tree.internal_branches) == 0:
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
            fw.write(f'all branches: {coalescent_tree.b_len_list}\n')
            fw.write(f'external branches: {coalescent_tree.external_branches}\n')
            fw.write(f'internal branches: {coalescent_tree.internal_branches}\n')

        with open(f"test_newick_{now}.txt", 'w') as fw:
            fw.write(newick_string)

    return branch_stats