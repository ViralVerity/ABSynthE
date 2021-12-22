from collections import defaultdict
from collections import OrderedDict

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
            elif start1 >= start2 and end1 > end2 and start1 < end2:
                frac = (end2-start1)/bin_size
            elif start1 < start2 and end1 <= end2 and end1 > start2:
                frac = (end1-start2)/bin_size
            else:
                continue
            
            new_lins[new_bin].append(frac*ltt_dict[old_bin])
                    
    for new_bin, totals in new_lins.items():
        average_lins[new_bin] = np.mean(totals)
        
    ltt_points = list(average_lins.keys())
    y_vals = list(average_lins.values())
    
    ltt_points.extend(y_vals)

    return ltt_points

    
def calculate_ltt_metrics(ltt_dict, coalescent_times):
    
    t_max_L = max(ltt_dict.items(), key = lambda k : k[1])[0][0]
    max_L = max(ltt_dict.items(), key = lambda k : k[1])[1]
        
    peak_y = max_L
    peak_x = t_max_L
    
    start_y = list(ltt_dict.values())[0]
    start_x = 0.0
    
    end_x = max(coalescent_times)
    end_y = list(ltt_dict.values())[-1]
    
    slope_1 = (peak_y-start_y)/(peak_x-start_x)
    slope_2 = (peak_y-end_y)/(end_x-peak_x)
    
    if slope_2 > 0:
        slope_ratio = slope_1/slope_2
    else:
        slope_ratio = None
    
    ltt_metrics = [max_L, t_max_L, slope_1, slope_2, slope_ratio]
    
    return ltt_metrics
    
    
    
    
    

