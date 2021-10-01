#take whole tree as input

import statistics

def calculate_branch_statistics(coalescent_tree):
    
    max_H = coalescent_tree.most_recent_date
    min_H = coalescent_tree.oldest_sample_date

    mean = statistics.mean(coalescent_tree.b_len_list)
    median = statistics.median(coalescent_tree.b_len_list)
    
    if len(coalescent_tree.b_len_list) >= 2:
        var = statistics.variance(coalescent_tree.b_len_list)
    else:
        var = None
        
    
    #These are done piecewise in the original over different time points
    int_mean = statistics.mean(coalescent_tree.internal_branches)
    int_median = statistics.median(coalescent_tree.internal_branches)
    
    if len(coalescent_tree.internal_branches) >= 2:
        int_var = statistics.variance(coalescent_tree.internal_branches)
    else:
        int_var = None
        
        
    
    ext_mean = statistics.mean(coalescent_tree.external_branches)
    ext_median = statistics.median(coalescent_tree.external_branches)
    
    if len(coalescent_tree.external_branches) >= 2:
        ext_var = statistics.variance(coalescent_tree.external_branches)
    else:
        ext_var = None
        
    
    mean_ratio = int_mean/ext_mean
    median_ratio = int_median/ext_median
    
    if int_var and ext_var:
        var_ratio = int_var/ext_var
    else:
        var_ratio = None
    
    tip_number = len(coalescent_tree.tips)
    
    return tip_number, max_H, min_H, mean, median, var, ext_mean, ext_median, ext_var, int_mean, int_median, int_var, mean_ratio, median_ratio, var_ratio