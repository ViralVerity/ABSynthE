import branch_length_parameters as bl
import topology_set as top
import LTT_metrics as ltt

import observed_summary_stats

import json

import numpy as np


##don't think any of these ever get called - where does the simulated get calculated?
def compare_BL(obs_vectors, coalescent_tree):
    
    sim_bl_pre = bl.calculate_branch_statistics(coalescent_tree)
    sim_bl = normalise(sim_bl_pre)
        
    obs_bl = obs_vectors[0]
    
    branch_difference = np.linalg.norm(obs_bl - sim_bl)
    
    return branch_difference

def compare_topology(obs_vectors, coalescent_tree):
    
    sim_top_pre = top.calculate_topology_params(coalescent_tree)
    sim_top = normalise(sim_top_pre)
    
    obs_top = obs_vectors[1]
    
    topology_difference = np.linalg.norm(obs_top - sim_top)
    
    return topology_difference


def compare_LTT_stats(obs_vectors, coalescent_tree):
    
    sim_ltt_pre = ltt.calculate_ltt_metrics(coalescent_tree.lineages_through_time)
    sim_ltt = normalise(sim_ltt_pre)
    
    obs_ltt = obs_vectors[2]
    
    ltt_stat_difference = np.linalg.norm(obs_ltt - sim_ltt)
    
    return LTT_stat_difference


def compare_LTT_points(obs_vectors, coalescent_tree):
    
    sim_ltt_points_pre = ltt.average_ltt_bins(coalescent_tree.lineages_through_time, coalescent_tree.coalescent_times)
    sim_ltt = normalise(sim_ltt_points_pre)
    
    obs_ltt_points = obs_vectors[3]
    
    ltt_point_difference = np.linalg.norm(obs_LTT_points - sim_LTT_points)
    
    return ltt_point_difference
  
    









    


