from collections import defaultdict
from collections import OrderedDict

def bin_sim(sim_ltt, obs_bins):
    
    intervals = []
    lineages = []
    bin_dict = defaultdict(list)
    sim_bins = {}



    for key, value in sim_ltt.items():
        start_interval = key[0]
        end_interval = key[1]
        lineage_number = value
        

        intervals.append((start_interval, end_interval))
        lineages.append(lineage_number)    

        
    

    for position, interval in enumerate(intervals):

        for bin_position,binn in enumerate(obs_bins):

            if interval[0] == 0 and binn[0] == 0.0:
                bin_dict[bin_position].append(position)

            if (interval[0] > binn[0] and interval[0] <= binn[1]) or (interval[1] > binn[0] and interval[1] <= binn[1]) or (interval[0] < binn[0] and interval[1] > binn[1]):

                bin_dict[bin_position].append(position)

        

    for binn_pos, interval_pos_list in bin_dict.items():

        if len(interval_pos_list) == 1:

            sim_bins[obs_bins[binn_pos]] = lineages[interval_pos_list[0]]


        if len(interval_pos_list) == 2:

            lineage_num_a = lineages[interval_pos_list[0]]
            lineage_num_b  = lineages[interval_pos_list[1]]

            interval_a = intervals[interval_pos_list[0]][1]

            bin_a = obs_bins[binn_pos][0]
            bin_b = obs_bins[binn_pos][1]


            weight_a = (interval_a - bin_a)
            weight_b = 1 - weight_a

            lineage_num = ((lineage_num_a*weight_a) + (lineage_num_b*weight_b))

            sim_bins[obs_bins[binn_pos]] = lineage_num


        elif len(interval_pos_list) > 2:

            last_position = len(interval_pos_list) - 1

            lineage_num_a = lineages[interval_pos_list[0]]
            lineage_num_b = lineages[interval_pos_list[last_position]]

            other_lins = [lineages[interval_pos_list[i]] for i in range(1, len(interval_pos_list)-1)]

            interval_a = intervals[interval_pos_list[0]][1]
            interval_b = intervals[interval_pos_list[last_position]][0]        

            bin_a = obs_bins[binn_pos][0]
            bin_b = obs_bins[binn_pos][1]

            weight_a = (interval_a - bin_a)
            weight_b = (bin_b - interval_b)

            other_weights = (1 - weight_a - weight_b)/(len(interval_pos_list)-2)

            final_other_lins = 0

            for i in other_lins:
                final_other_lins += (other_weights*i)

            lins = ((lineage_num_a*weight_a) + (lineage_num_b*weight_b) + final_other_lins)

            sim_bins[obs_bins[binn_pos]] = lins
            
            
    lineage_vector = list(sim_bins.values())
    
    return lineage_vector

    
def calculate_LTT_metrics(ltt_dict_input):
    
    ltt_dict = OrderedDict(ltt_dict_input)
    
    max_L = max(ltt_dict.values())
    
    t_max_L = max(ltt_dict, key = lambda k:ltt_dict[k])[0] #We'll get the start of that interval
    
    
    
    y1 = max_L
    y0 = next(iter(ltt_dict.values()))

    x1 = t_max_L
    x0 = 0.0
    
    y2 = next(reversed(ltt_dict.values()))

    x2 = next(reversed(ltt_dict.keys()))[1]


    slope_1 = abs((y2-y1)/(x2-x1))
    slope_2 = abs((y1-y0)/(x1-x0))
    
    slope_ratio = slope_1/slope_2
    
    return max_L, t_max_L, slope_1, slope_2, slope_ratio  
    
    
    
    
    

