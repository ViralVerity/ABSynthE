import branch_length_parameters as BL
import topology_set as TOP
import LTT_metrics as LTT


import json

import numpy as np

#Need to recalculate these when we get the new tree

def normalise(vector):
    norm=np.linalg.norm(vector, ord=1)
    return vector/norm

def get_observed_SS():
    
    obs_tips = 371

    #Branch lengths#
    obs_max_H = 0.4187139080715566 
    obs_min_H = 0.07898788067423779

    obs_mean_lengths = 0.024626764928399065
    obs_median_lengths = 0.016550288168320965
    obs_var_lengths = 0.0006479307331597084

    obs_mean_external = 0.030683040421764732
    obs_median_external = 0.022982686457671875
    obs_var_external = 0.0007669213909417984

    obs_mean_internal = 0.01853712152597904
    obs_median_internal = 0.012126822480269583
    obs_var_internal = 0.00045590389615045627

    obs_mean_ratio = 0.6041487828836514
    obs_median_ratio = 0.5276503468210313
    obs_var_ratio = 0.5944597471594775

    #Topology set#

    obs_colless = 55126
    obs_sackin = 6549
    obs_WD_ratio = 0.7666666666666667
    obs_delta_W = 9
    obs_max_ladder = 4
    obs_IL_nodes = 0.010958904109589041
    obs_staircase_1 = 0.0027472527472527475
    obs_staircase_2 = 0.6222222222222222

    #LTT set#

    max_L = 141.79283007931681
    t_max_L = 0.08374278161431133
    slope_1 = 417.3279994541845
    slope_2 = 1058.599752831301
    slope_ratio = 0.3942264282019818


    obs_LTT_points_pre = [53.14274216100124, 110.8591972315072, 134.32012378033238, 136.88603187626754, 141.79283007931681, 99.92958023626413, 75.66925004968692, 62.702353762861016, 51.64930119445268, 52.110751732458105, 72.74663279140115, 79.16184842239603, 51.048797530829425, 49.81492749361463, 44.18875202397806, 19.49287400488278, 12.4945282724165, 7.666291528846369, 2, 2]
    
    LTT_bins = [(0, 0.02093569540357783), (0.02093569540357783, 0.04187139080715566), (0.04187139080715566, 0.0628070862107335), (0.0628070862107335, 0.08374278161431133), (0.08374278161431133, 0.10467847701788915), (0.10467847701788915, 0.125614172421467), (0.125614172421467, 0.14654986782504484), (0.14654986782504484, 0.16748556322862268), (0.16748556322862268, 0.18842125863220052), (0.18842125863220052, 0.20935695403577836), (0.20935695403577836, 0.2302926494393562), (0.2302926494393562, 0.25122834484293405), (0.25122834484293405, 0.27216404024651186), (0.27216404024651186, 0.2930997356500897), (0.2930997356500897, 0.3140354310536675), (0.3140354310536675, 0.3349711264572453), (0.3349711264572453, 0.3559068218608231), (0.3559068218608231, 0.37684251726440093), (0.37684251726440093, 0.39777821266797875), (0.39777821266797875, 0.41871390807155656)]
    
    obs_BL_pre = [obs_max_H, obs_min_H, obs_mean_lengths, obs_median_lengths, obs_var_lengths, obs_mean_external, obs_median_external, obs_var_external, obs_mean_internal, obs_median_internal, obs_var_internal, obs_mean_ratio, obs_median_ratio, obs_var_ratio]

    obs_top_pre = [obs_colless, obs_sackin, obs_WD_ratio, obs_delta_W, obs_max_ladder, obs_IL_nodes, obs_staircase_1, obs_staircase_2]

    obs_LTT_pre = [max_L, t_max_L, slope_1, slope_2, slope_ratio]
    
  
    obs_BL = normalise(obs_BL_pre)
    obs_top = normalise(obs_top_pre)
    obs_LTT = normalise(obs_LTT_pre)
    obs_LTT_points = normalise(obs_LTT_points_pre)
    
    return obs_BL, obs_top, obs_LTT, obs_LTT_points, LTT_bins, obs_tips


def compare_BL(obs_vectors, coalescent_tree):
    
    sim_BL_pre = BL.calculate_branch_statistics(coalescent_tree)
    sim_BL = normalise(sim_BL_pre)
    
    print(sim_BL)
    
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
    









    


