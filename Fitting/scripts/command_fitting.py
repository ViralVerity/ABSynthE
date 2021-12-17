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


def common_part(config, a,b,c):
    
    config["number_model_iterations"] = 1
    config["log_every"] =  1
    config["cfr"] = 0.7
    
    config["sampling_percentage"] = 0.16
    config["sampling_scheme"] = "uniform"
        
    config["input_directory"] = "../../SLE_EBOV_input_files/"

    config["case_limit"] = 3000
    config["day_limit"] = 124 
    config["output_tree"] = True
    config["output_skyline"] = False
    
    config["overwrite"] = False
    config["verbose"] = True

    config["calculate_R0"] = False
    
    cwd = os.getcwd()
    thisdir = os.path.abspath(os.path.dirname(__file__))
    
    # config = file_funcs.make_directories(config)
    # config = file_funcs.make_summary_files(config)

    config["distributions"] = dist_funcs.define_distributions() 
    config["population_structure"] = file_funcs.parse_population_information(os.path.join(config["input_directory"], "population_config.yaml"))
    config = make_contact_dicts(config["input_directory"], config)
    
    config["capped"] = True
    
    config['a'] = a 
    config['b'] = b
    config['c'] = c
    
    config["calculate_R0"] = False
    
    return config

def simulate_epidemic_all(a,b,c):
    
    config = {}

    config["output_directory"] = "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/all/"
    config["output_ltt"] = True

    config = common_part(config, a,b,c)
    
    a_sim, dist_mvmt, ch_mvmt = run_model(config, "all", 47) #calls the run_model_fitting.py script
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps

def simulate_epidemic_ltt(a,b,c):
    
    config = {}

    config["output_directory"] = "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/ltt/"
    config["output_ltt"] = True

    config = common_part(config,a,b,c)
    
    a_sim, dist_mvmt, ch_mvmt = run_model(config, "ltt", 7) #calls the run_model_fitting.py script
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps 
    
    
def simulate_epidemic_ltt_points(a,b,c):
    
    config = {}

    config["output_directory"] = "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/ltt_points/"
    config["output_ltt"] = True

    config = common_part(config, a,b,c)
    
    a_sim, dist_mvmt, ch_mvmt = run_model(config, "ltt_points", 22) #calls the run_model_fitting.py script 
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps
    
 
def simulate_epidemic_bl(a,b,c):
    
    config = {}

    config["output_directory"] = "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/branch/"
    config["output_ltt"] = False

    config = common_part(config,a,b,c)
    
    a_sim, dist_mvmt, ch_mvmt = run_model(config, "branch", 13) #calls the run_model_fitting.py script 
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps
    
def simulate_epidemic_top(a,b,c):
    
    config = {}

    config["output_directory"] =  "/localdisk/home/s1732989/ABSynthE/Fitting/test_results/topology/"
    config["output_ltt"] = False

    config = common_part(config, a,b,c)
    
    a_sim, dist_mvmt, ch_mvmt = run_model(config, "topology", 9) #calls the run_model_fitting.py script 
    
    dist_jumps = get_jumps(dist_mvmt)
    ch_jumps = get_jumps(ch_mvmt)
    
    return a_sim, ch_jumps, dist_jumps


def get_jumps(mvmt_dict): 
    
    total_jumps = 0
    
    for v in mvmt_dict.values(): 
        total_jumps += len(v) 
   
    return total_jumps
