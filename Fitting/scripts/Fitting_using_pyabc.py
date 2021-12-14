from collections import defaultdict
import random
from multiprocessing.pool import ThreadPool
import time
import os

from pyabc.sampler import ConcurrentFutureSampler
from concurrent.futures import ThreadPoolExecutor

import Tree_simulator_fitting as cts
from Simulate_epidemic_fitting import *
from movement_fitting import *

import observed_summary_stats

import tempfile
import pyabc

import argparse

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False,
    description=preamble(__version__))
    
    parser.add_argument("--results-path", dest="results_path")
    parser.add_argument("--summary-stats-set", dest="summary_stats_set") #one of four
    
    args.parser.parse_args(sysargs)

    results_path = args.results_path
    summary_stats_set = args.summary_stats_set

    run_number = 1
    try:
        os.mkdir(os.path.join(results_path, str(run_number)))
    except FileExistsError:
        pass

    observed_SS = get_observed_SS(summary_stats_set)
    
    if summary_stats_set == "all":
        summary_stats = observed_SS[0]
    elif summmary_stats_set == "branch":
        summary_stats = list(observed_ss[0][0])
    elif summmary_stats_set == "topology":
        summmary_stats = list(observed_ss[0][1])
    elif summmary_stats_set == "ltt":
        summmary_stats = list(observed_ss[0][2])
    elif summmary_stats_set == "ltt_points":
        summary_stats = list(observed_ss[0][3])

    observed = {"a":summary_stats, "b":observed_SS[6], "c":observed_SS[7]}
    
    parameters = dict(a=(0.5,1), b=(0,0.3), c=(0,0.5))
    prior = pyabc.Distribution(**{key: pyabc.RV("uniform", a, b - a) for key, (a,b) in parameters.items()})

    pool = ThreadPoolExecutor(max_workers=48)
    sampler = ConcurrentFutureSampler(pool)

    abc = pyabc.ABCSMC(simulate_pyabc, prior, distance, sampler=sampler) 

    db_path = (f"sqlite:///{summary_stats_set}.db")

    abc_id = abc.new(db_path, observed)

    history = abc.run(max_nr_populations=10, minimum_epsilon=0.3)
    
def simulate_pyabc(parameter):
    result = simulate_epidemic(**parameter) #this could be just call command? I dont' think so, 
    #but I think if I just have a script that then makes the config etc so we don't need to call it with args
    #I think it just runs the epidemic once and records the result. 
    return {"a":result[0], "b":result[1], "c":result[2]} #so this returns the sample as a dictionary



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
    
    if None not in sim_a_vector:
        mid_a_x = normalise(sim_a_vector)
        new_a_x = np.array(mid_a_x)
        new_a_y = np.array(obs_a_vector)
    else:
        indices = [i for i,x in enumerate(sim_a_vector) if x == None]
        
        processing = list(sim_a_vector)
        processing[:] = [x for x in processing if x != None]
       
        mid_a_x = normalise(processing)
        new_a_x = np.array(mid_a_x)
        
        count = 0
        processing_y = list(obs_a_vector)
        for i in indices:
            processing_y.pop(i-count)
            count += 1
        
        new_a_y = np.array(processing_y)

    if new_a_x.size == new_a_y.size: 
        dist_a = np.linalg.norm(new_a_x - new_a_y)
    else:
        print("error - not the same len for a" + str(len(new_a_x) + " " + str(len(new_a_y))))
        return
    
    dist_b = np.linalg.norm(sim_b - obs_b)
    dist_c = np.linalg.norm(sim_c - obs_c)
    
    final_b = dist_b/obs_b #not sure this is totally right - need to think about this more
    final_c = dist_c/obs_c

    dist = dist_a + final_b + final_c
    
    return dist
    










