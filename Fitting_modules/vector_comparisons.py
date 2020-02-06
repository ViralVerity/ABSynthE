import branch_length_parameters as BL
import topology_set as TOP

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


    ##
    
    #LTT points#
    
    ##
    
    obs_BL_pre = [obs_max_H, obs_min_H, obs_mean_lengths, obs_median_lengths, obs_var_lengths, obs_mean_external, obs_median_external, obs_var_external, obs_mean_internal, obs_median_internal, obs_var_internal, obs_mean_ratio, obs_median_ratio, obs_var_ratio]

    obs_top_pre = [obs_colless, obs_sackin, obs_WD_ratio, obs_delta_W, obs_max_ladder, obs_IL_nodes, obs_staircase_1, obs_staircase_2]

    obs_LTT_pre = []
    
    obs_LTT_points_pre = []
    
    obs_BL = normalise(obs_BL_pre)
    obs_top = normalise(obs_top_pre)
    obs_LTT = normalise(obs_LTT_pre)
    obs_LTT_points = normalise(obs_LTT_points_pre)
    
    return obs_BL, obs_top, obs_LTT, obs_LTT_points, obs_tips


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


##Need to do the code for calculating these things##

def compare_LTT_stats(obs_vectors, coalescent_tree):
    
    sim_LTT_pre = LTT.calculate_LTT_params(coalescent_tree)
    sim_LTT = normalise(sim_LTT_pre)
    
    obs_LTT = obs_vectors[2]
    
    LTT_stat_difference = np.linalg.norm(obs_LTT - sim_LTT)
    
    return LTT_stat_difference

def compare_LTT_points(obs_vectors, coalescent_tree):
    
    sim_LTT_points_pre = LTT.get_LTT_points(coalescent_tree)
    sim_LTT_points = normalise(sim_LTT_points_pre)
    
    obs_LTT_points = obs_vectors[3]
    
    LTT_point_difference = np.linalg.norm(obs_LTT_points - sim_LTT_points)
    
    return LTT_point_difference
    
def get_tip_difference(obs_vectors, coalescent_tree):
    
    obs_tips = obs_vectors[4]
    sim_tips = len(coalescent_tree.tips)
    
    tip_difference = abs(obs_tips - sim_tips)
    
    return tip_difference
    









    


