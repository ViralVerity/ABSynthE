import branch_length_parameters as bl
import topology_set as top
import LTT_metrics as ltt

import observed_summary_stats

import json

import numpy as np


def normalise(vector):
    norm=np.linalg.norm(vector, ord=1)
    return vector/norm

def get_observed_SS():
    
    branch_set, topology_set, ltt_set, ltt_points, dist, ch = import_observed_stats()
    
    obs_bl = normalise(branch_set)
    obs_top = normalise(topology_set)
    obs_ltt = normalise(ltt_set)
    obs_ltt_points = normalise(ltt_points)
    
    return obs_bl, obs_top, obs_ltt, obs_ltt_points, obs_tips, observed_dist, observed_ch


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
    
    sim_LTT_points_pre = ltt.bin_sim(coalescent_tree.lineages_through_time, LTT_bins)
        
    obs_LTT_points = obs_vectors[3]
    
    LTT_point_difference = np.linalg.norm(obs_LTT_points - sim_LTT_points)
    
    return LTT_point_difference
  
    
def get_tip_difference(obs_vectors, coalescent_tree): #not sure at the moment that this should be separate - I think t shoul dbe included in the normalised vector
    
    obs_tips = obs_vectors[5]
    sim_tips = len(coalescent_tree.tips)
    
    tip_difference = abs(obs_tips - sim_tips)
    
    return tip_difference
    









    


