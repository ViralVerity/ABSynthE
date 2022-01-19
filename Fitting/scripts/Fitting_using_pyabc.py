from collections import defaultdict
import random
from multiprocessing.pool import ThreadPool
import time
import os
import numpy as np

from pyabc.sampler import ConcurrentFutureSampler
from concurrent.futures import ThreadPoolExecutor

from command_fitting import *
import observed_summary_stats

import tempfile
import pyabc

import argparse

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False,
    description="fitting")
    
    parser.add_argument("--summary-stats-set", dest="summary_stats_set") #one of four
    parser.add_argument("--threads", dest="threads", default=28)
    
    args = parser.parse_args(sysargs)

    summary_stats_set = args.summary_stats_set
    threads = int(args.threads)

    run_number = 1

    observed_SS = observed_summary_stats.get_observed_SS(summary_stats_set)
    
    if summary_stats_set == "all":
        summary_stats = observed_SS[0]
        function = simulate_pyabc_all
    elif summary_stats_set == "branch":
        summary_stats = list(observed_SS[0][0])
        function = simulate_pyabc_bl
    elif summary_stats_set == "topology":
        summary_stats = list(observed_SS[0][1])
        function = simulate_pyabc_top
    elif summary_stats_set == "ltt":
        summary_stats = list(observed_SS[0][2])
        function = simulate_pyabc_ltt
    elif summary_stats_set == "ltt_points":
        summary_stats = list(observed_SS[0][3])
        function = simulate_pyabc_ltt_points

    observed = {"a":summary_stats, "b":observed_SS[1], "c":observed_SS[2]}
    # observed = {"a":summary_stats, "b": observed_SS[1], "c":observed_SS[2], "d":observed_SS[3]}

    parameters = dict(a=(0.5,1), b=(0,0.5), c=(0,0.5))
    prior = pyabc.Distribution(**{key: pyabc.RV("uniform", x, y - x) for key, (x,y) in parameters.items()})

    pool = ThreadPoolExecutor(max_workers=threads)
    sampler = ConcurrentFutureSampler(pool)

    db_path = (f"sqlite:///{summary_stats_set}.db")

    # abc = pyabc.ABCSMC(function, prior, distance, sampler=sampler) 
    # abc_id = abc.new(db_path, observed)
    # abc.run(max_nr_populations=10, minimum_epsilon=0.1)

    abc_continued = pyabc.ABCSMC(function, prior, distance, sampler=sampler)
    abc_continued.load(db_path, 4)
    abc_continued.run(max_nr_populations=7, minimum_epsilon=0.1)

    
# def simulate_pyabc_all(parameter):
#     result = simulate_epidemic_all(**parameter) 
#     return {"a":result[0], "b":result[1], "c":result[2], "d":result[3]} 
# def simulate_pyabc_ltt(parameter):
#     result = simulate_epidemic_ltt(**parameter) 
#     return {"a":result[0], "b":result[1], "c":result[2], "d":result[3]} 
# def simulate_pyabc_ltt_points(parameter):
#     result = simulate_epidemic_ltt_points(**parameter) 
#     return {"a":result[0], "b":result[1], "c":result[2], "d":result[3]} 
# def simulate_pyabc_bl(parameter):
#     result = simulate_epidemic_bl(**parameter) 
#     return {"a":result[0], "b":result[1], "c":result[2], "d":result[3]} 
# def simulate_pyabc_top(parameter):
#     result = simulate_epidemic_top(**parameter) 
#     return {"a":result[0], "b":result[1], "c":result[2], "d":result[3]} 

def simulate_pyabc_all(parameter):
    result = simulate_epidemic_all(**parameter) 
    return {"a":result[0], "b":result[1], "c":result[2]} 
def simulate_pyabc_ltt(parameter):
    result = simulate_epidemic_ltt(**parameter) 
    return {"a":result[0], "b":result[1], "c":result[2]} 
def simulate_pyabc_ltt_points(parameter):
    result = simulate_epidemic_ltt_points(**parameter) 
    return {"a":result[0], "b":result[1], "c":result[2]} 
def simulate_pyabc_bl(parameter):
    result = simulate_epidemic_bl(**parameter) 
    return {"a":result[0], "b":result[1], "c":result[2]} 
def simulate_pyabc_top(parameter):
    result = simulate_epidemic_top(**parameter) 
    return {"a":result[0], "b":result[1], "c":result[2]} 

def normalise(vector):
    norm=np.linalg.norm(vector, ord=1)
    return vector/norm

def distance(x,y): #inputs are the dictionaries
    
    sim_a = x['a']
    obs_a = y['a']
    
    sim_b = x['b']
    obs_b = y['b']
    
    sim_c = x['c']
    obs_c = y['c']

    # sim_d = x['d']
    # obs_d = y['d']
    
    check = True if True in np.isnan(np.array(sim_a)) else False
    if not check:
        if not all(v == 0 for v in sim_a):
            mid_a_x = normalise(sim_a)
        else:
            mid_a_x = sim_a
        new_a_x = np.array(mid_a_x)
        new_a_y = np.array(obs_a)
    else: #only compare values for which there are values 
        indices = np.argwhere(np.isnan(sim_a))
        processing = np.array(sim_a)
        processing = processing[~np.isnan((processing))]        
       
        new_a_x = normalise(processing)
        # new_a_x = np.array(mid_a_x)
        
        count = 0
        processing_y = list(obs_a)
        for i in indices:
            processing_y.pop(i[0]-count)
            count += 1
        
        new_a_y = np.array(processing_y)

    if new_a_x.size == new_a_y.size: 
        dist_a = np.linalg.norm(new_a_x - new_a_y)
    else:
        print(f"error - not the same len for a {(len(new_a_x))} against {len(new_a_y)}\n")
        print(sim_a)
        print(new_a_x)
        print(new_a_y)
        return
    
    if any([np.isnan(x) for x in new_a_x]): #keep this - still should be now nan in final part
        print(f"NaN in vector: {sim_a}, {new_a_x}")

    dist_b = np.linalg.norm(sim_b - obs_b)
    dist_c = np.linalg.norm(sim_c - obs_c)
    # dist_d = np.linalg.norm(sim_d - obs_d)
    
    final_b = dist_b/obs_b 
    final_c = dist_c/obs_c
    # final_d = dist_d/obs_d

    dist = dist_a + final_b + final_c
    # dist = dist_a + final_b + final_cs + final_d
    
    return dist
    

if __name__ == '__main__':
    main()
