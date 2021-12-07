#take whole tree as input

import numpy as np

def calculate_branch_statistics(coalescent_tree):
    
    max_H = coalescent_tree.most_recent_date
    min_H = coalescent_tree.oldest_sample_date

    mean_lengths = np.mean(coalescent_tree.b_len_list)
    median_lengths = np.median(coalescent_tree.b_len_list)
    
    if len(coalescent_tree.b_len_list) >= 2:
        var = np.variance(coalescent_tree.b_len_list)
    else:
        var = None
        
    
    #These are done piecewise in the original over different time points
    mean_internal = np.mean(coalescent_tree.internal_branches)
    median_internal = np.median(coalescent_tree.internal_branches)
    
    if len(coalescent_tree.internal_branches) >= 2:
        var_internal = np.variance(coalescent_tree.internal_branches)
    else:
        var_internal = None
        
    mean_external = np.mean(coalescent_tree.external_branches)
    median_external = np.median(coalescent_tree.external_branches)
    
    if len(coalescent_tree.external_branches) >= 2:
        var_external = np.variance(coalescent_tree.external_branches)
    else:
        var_external = None
        
    
    mean_ratio = mean_internal/mean_external
    median_ratio = median_internal/median_external
    
    if var_internal and var_external:
        var_ratio = var_internal/var_external
    else:
        var_ratio = None
    
    tip_number = len(coalescent_tree.tips) 
    
    branch_stats = [tip_number, mean_lengths, median_lengths, var_lengths, mean_external, median_external, var_external, mean_internal, median_internal, var_internal, mean_ratio, median_ratio, var_ratio]
    
    return branch_stats