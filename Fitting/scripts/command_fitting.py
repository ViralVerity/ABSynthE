#!/usr/bin/env python3
from absynthe import __version__

import argparse

from collections import defaultdict
from multiprocessing.pool import ThreadPool

from run_model_fitting import *

import absynthe.set_up.file_functions as file_funcs
from absynthe.set_up.make_contact_dicts import *
import absynthe.set_up.index_functions as index_funcs
import absynthe.set_up.distribution_functions as dist_funcs


@profile
def simulate_epidemic_all(parameters):
    
    parameters["output_directory"] = "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/all/"
    parameters["output_ltt"] = True
    
    a_sim, dist_mvmt, ch_mvmt = run_model(parameters, "all", 47) #calls the run_model_fitting.py script
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps

def simulate_epidemic_ltt(parameters):
    
    parameters["output_directory"] = "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/ltt/"
    parameters["output_ltt"] = True
    
    a_sim, dist_mvmt, ch_mvmt = run_model(parameters, "ltt", 7) #calls the run_model_fitting.py script
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps 
    
    
def simulate_epidemic_ltt_points(parameters):
    
    parameters["output_directory"] = "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/ltt_points/"
    parameters["output_ltt"] = True
    
    a_sim, dist_mvmt, ch_mvmt = run_model(parameters, "ltt_points", 22) #calls the run_model_fitting.py script 
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps
    
 
def simulate_epidemic_bl(parameters):
    
    parameters["output_directory"] = "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/branch/"
    parameters["output_ltt"] = False
    
    a_sim, dist_mvmt, ch_mvmt = run_model(parameters, "branch", 13) #calls the run_model_fitting.py script 
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps
    
def simulate_epidemic_top(parameters):
    
    parameters["output_directory"] =  "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/topology/"
    parameters["output_ltt"] = False
    
    a_sim, dist_mvmt, ch_mvmt = run_model(parameters, "topology", 9) #calls the run_model_fitting.py script 
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps


def get_jumps(mvmt_dict): 
    
    total_jumps = 0
    
    for v in mvmt_dict.values(): 
        total_jumps += len(v) 
   
    return total_jumps
