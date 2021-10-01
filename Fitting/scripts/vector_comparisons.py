import branch_length_parameters as BL
import topology_set as TOP
import LTT_metrics as LTT


import json

import numpy as np


def normalise(vector):
    norm=np.linalg.norm(vector, ord=1)
    return vector/norm

def get_observed_SS():
    
    obs_tips = 371
    
    observed_dist = 93
    observed_ch = 30

    #Branch lengths#
    obs_max_H = 0.4272233755511325
    obs_min_H = 0.07927817007151033

    obs_mean_lengths = 0.024554999389390132
    obs_median_lengths = 0.01618327945212711
    obs_var_lengths = 0.0006457965315134023

    obs_mean_external = 0.030333644010762013
    obs_median_external = 0.021228056094926562
    obs_var_external = 0.0007945437958549986

    obs_mean_internal = 0.01874503420096475
    obs_median_internal = 0.011576011003413655
    obs_var_internal = 0.00043048336004913076

    obs_mean_ratio = 0.6179618312364396
    obs_median_ratio = 0.5453165825287359
    obs_var_ratio = 0.5417994102966885

    #Topology set#

    obs_colless = 56120
    obs_sackin = 8835
    obs_WD_ratio = 0.8235294117647058
    obs_delta_W = 9
    obs_max_ladder = 0.02425876010781671
    obs_IL_nodes = 0.46216216216216216
    obs_staircase_1 = 0.002702702702702703
    obs_staircase_2 = 0.11077844311377245

    #LTT set#

    max_L = 148.97678349186248
    t_max_L = 0.06408350633266988
    slope_1 = 404.73876858572703
    slope_2 = 1394.061280232578
    slope_ratio = 0.29033068655217403


    obs_LTT_points_pre = [59.64044861194819, 124.75396017627713, 133.9364509309589, 148.97678349186248, 148.33907671517255, 91.59093272776481, 70.07876927609597, 58.46111323028046, 50.24916890609602, 55.63545914760437, 72.6885652880981, 67.61114496403488, 41.31670189441371, 48.50878216596301, 40.355851127800555, 27.025321871739628, 19.004614160153153, 11.476947041818141, 2.0021308852811104, 2]
    
    LTT_bins = [(0, 0.021361168777556623), (0.021361168777556623, 0.042722337555113246), (0.042722337555113246, 0.06408350633266988), (0.06408350633266988, 0.08544467511022649), (0.08544467511022649, 0.10680584388778311), (0.10680584388778311, 0.12816701266533972), (0.12816701266533972, 0.14952818144289634), (0.14952818144289634, 0.17088935022045296), (0.17088935022045296, 0.19225051899800957), (0.19225051899800957, 0.2136116877755662), (0.2136116877755662, 0.2349728565531228), (0.2349728565531228, 0.25633402533067945), (0.25633402533067945, 0.27769519410823607), (0.27769519410823607, 0.2990563628857927), (0.2990563628857927, 0.3204175316633493), (0.3204175316633493, 0.3417787004409059), (0.3417787004409059, 0.36313986921846253), (0.36313986921846253, 0.38450103799601915), (0.38450103799601915, 0.40586220677357576), (0.40586220677357576, 0.4272233755511324)]
    
    obs_BL_pre = [obs_tips, obs_max_H, obs_min_H, obs_mean_lengths, obs_median_lengths, obs_var_lengths, obs_mean_external, obs_median_external, obs_var_external, obs_mean_internal, obs_median_internal, obs_var_internal, obs_mean_ratio, obs_median_ratio, obs_var_ratio]

    obs_top_pre = [obs_tips, obs_colless, obs_sackin, obs_WD_ratio, obs_delta_W, obs_max_ladder, obs_IL_nodes, obs_staircase_1, obs_staircase_2]

    obs_LTT_pre = [max_L, t_max_L, slope_1, slope_2, slope_ratio]
    
  
    obs_BL = normalise(obs_BL_pre)
    obs_top = normalise(obs_top_pre)
    obs_LTT = normalise(obs_LTT_pre)
    obs_LTT_points = normalise(obs_LTT_points_pre)
    
    return obs_BL, obs_top, obs_LTT, obs_LTT_points, LTT_bins, obs_tips, observed_dist, observed_ch


def compare_BL(obs_vectors, coalescent_tree):
    
    sim_BL_pre = BL.calculate_branch_statistics(coalescent_tree)
    sim_BL = normalise(sim_BL_pre)
        
    obs_BL = obs_vectors[0]
    
    branch_difference = np.linalg.norm(obs_BL - sim_BL)
    
    return branch_difference

def compare_topology(obs_vectors, coalescent_tree):
    
    sim_top_pre = TOP.calculate_topology_params(coalescent_tree)
    sim_top = normalise(sim_top_pre)
    
    obs_top = obs_vectors[1]
    
    topology_difference = np.linalg.norm(obs_top - sim_top)
    
    return topology_difference


def compare_LTT_stats(obs_vectors, coalescent_tree):
    
    sim_LTT_pre = LTT.calculate_LTT_metrics(coalescent_tree.lineages_through_time)
    sim_LTT = normalise(sim_LTT_pre)
    
    obs_LTT = obs_vectors[2]
    
    LTT_stat_difference = np.linalg.norm(obs_LTT - sim_LTT)
    
    return LTT_stat_difference


def compare_LTT_points(obs_vectors, coalescent_tree):
    
    LTT_bins = obs_vectors[4]
    
    sim_LTT_points_pre = LTT.bin_sim(coalescent_tree.lineages_through_time, LTT_bins)
        
    sim_LTT_points = normalise(sim_LTT_points_pre)
    
    obs_LTT_points = obs_vectors[3]
    
    if len(obs_LTT_points) == len(sim_LTT_points):
        LTT_point_difference = np.linalg.norm(obs_LTT_points - sim_LTT_points)
    
        return LTT_point_difference
    else:
        return
    
def get_tip_difference(obs_vectors, coalescent_tree):
    
    obs_tips = obs_vectors[5]
    sim_tips = len(coalescent_tree.tips)
    
    tip_difference = abs(obs_tips - sim_tips)
    
    return tip_difference
    









    


