from collections import defaultdict
from collections import OrderedDict
import numpy as np

def average_ltt_bins(ltt_dict, coalescent_times):
    
    full_len = max(coalescent_times)

    bin_size = full_len/20

    new_bin = 0
    new_bin_list = []
    new_tup_list = []
    
    for i in range(20):
        new_bin_list.append(new_bin)
        new_tup_list.append((new_bin, new_bin+bin_size))
        new_bin += bin_size
    
    #makes sure there's no rounding nonsense and the last bin ends with the end of the tree
    last_tup = new_tup_list[-1]
    start_final = last_tup[0]
    new_tup_list.pop()
    new_tup_list.append((start_final,full_len))

    new_lins = defaultdict(list)
    average_lins = {}

    for new_bin in new_tup_list:
        for old_bin in ltt_dict.keys():
            start1 = old_bin[0]
            start2 = new_bin[0]
            end1 = old_bin[1]
            end2 = new_bin[1]
            
            if start1 >= start2 and end1 <= end2:
                frac = (end1-start1)/bin_size
            elif start1 <= start2 and end1 < end2:
                frac = (end2-start2)/bin_size
            elif start1 >= start2 and end1 > end2 and start1 < end2:
                frac = (end2-start1)/bin_size
            elif start1 < start2 and end1 <= end2 and end1 > start2:
                frac = (end1-start2)/bin_size
            else:
                continue
            
            new_lins[new_bin].append(frac*ltt_dict[old_bin])
                    
    for new_bin, totals in new_lins.items():
        if len(totals) > 0:
            average_lins[new_bin] = np.mean(totals)
        else:
            average_lins[new_bin] = 0
        
    ltt_points = []
    for tup in average_lins.keys():
        ltt_points.append(tup[0])
    y_vals = list(average_lins.values())
    
    ltt_points.extend(y_vals)

    return ltt_points

    
def calculate_ltt_metrics(ltt_dict, coalescent_times, coalescent_tree):
    
    t_max_L = max(ltt_dict.items(), key = lambda k : k[1])[0][0]
    max_L = max(ltt_dict.items(), key = lambda k : k[1])[1]
        
    peak_y = max_L
    peak_x = t_max_L
    
    start_y = list(ltt_dict.values())[0]
    start_x = 0.0
    
    end_x = max(coalescent_times)
    end_y = list(ltt_dict.values())[-1]
    
    #peak and end could be the same
    if peak_x != end_x:
        slope_1 = (peak_y-start_y)/(peak_x-start_x)
        slope_2 = (peak_y-end_y)/(end_x-peak_x)
    else:
        slope_1 = (peak_y-start_y)/(peak_x-start_x)
        slope_2 = (peak_y-end_y)/1
    
    
    if slope_2 and slope_2 > 0:
        slope_ratio = slope_1/slope_2
    else:
        slope_ratio = np.nan

    sampling_times = []
    branching_times = []
    for node in coalescent_tree.final_nodes:
        if node.type == "individual":
            sampling_times.append(node.root_to_tip)
        elif node.type == "coalescent":
            branching_times.append(node.root_to_tip)
        else:
            print("Error: transmission node type still in the final nodes.\n")

    sampling_times = sorted(sampling_times)
    branching_times = sorted(branching_times)

    sampling_diffs = []
    for count,i in enumerate(sampling_times):
        try:
            next_one = sampling_times[count + 1]
            sampling_diffs.append(next_one - i)
        except IndexError:
            pass

    branching_diffs = []
    for count,i in enumerate(branching_times):
        try:
            next_one = branching_times[count + 1]
            branching_diffs.append(next_one - i)
        except IndexError:
            pass
        
    if len(sampling_diffs) > 0:
        mean_s_time = np.mean(sampling_diffs)
    else:
        mean_s_time = np.nan 
    if len(branching_diffs) > 0:
        mean_b_time = np.mean(branching_diffs)
    else:
        mean_b_time = np.nan
        
    ltt_metrics = [max_L, t_max_L, slope_1, slope_2, slope_ratio, mean_s_time, mean_b_time]
    
    return ltt_metrics
    
    
    
    
    

