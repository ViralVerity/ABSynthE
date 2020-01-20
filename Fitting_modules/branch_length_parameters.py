#take whole tree as input
#get ones that haven't been removed

import statistics


def calculate_branch_statistics(coalescent_tree):
    
    max_H = coalescent_tree.most_recent_date
    min_H = coalescent_tree.oldest_sample_date

    mean = statistics.mean(coalescent_tree.b_len_list)
    median = statistics.median(coalescent_tree.b_len_list)
    var = statistics.variance(coalescent_tree.b_len_list)
    
    #These are done piecewise in the original over different time points
    int_mean = statistics.mean(coalescent_tree.internal_branches)
    int_median = statistics.median(coalescent_tree.internal_branches)
    int_var = statistics.variance(coalescent_tree.internal_branches)
    
    ext_mean = statistics.mean(coalescent_tree.external_branches)
    ext_median = statistics.median(coalescent_tree.external_branches)
    ext_var = statistics.variance(coalescent_tree.external_branches)
    
    mean_ratio = int_mean/ext_mean
    median_ratio = int_median/ext_median
    var_ratio = int_var/ext_var
    
    return max_H, min_H, mean, median, var, int_mean, int_median, int_var, ext_mean, ext_var, mean_ratio, median_ratio, var_ratio