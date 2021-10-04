from collections import defaultdict
import random
from multiprocessing.pool import ThreadPool
import time
import os

from pyabc.sampler import ConcurrentFutureSampler
from concurrent.futures import ThreadPoolExecutor

import Tree_simulator_fitting as cts
from Simulate_epidemic_fitting import *
from vector_comparisons import *
from movement_fitting import *


import tempfile
import pyabc

dropbox_path = "/disk2/home/s1732989/ABM/Fitting/"

#results_path = "LTT_ABCSMC/"
results_path = "top_ABCSMC/"
#results_path = "bl_ABCSMC/"

def normalise(vector):
    norm=np.linalg.norm(vector, ord=1)
    return vector/norm


#dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"
#results_path = "Looping models/Results/Fitting/LTT/"

run_number = 1

try:
    os.mkdir(os.path.join(dropbox_path, results_path, str(run_number)))
except FileExistsError:
    pass

observed_SS = get_observed_SS()

#new_LTT = list(observed_SS[3])
#top = list(observed_SS[1])
bl = list(observed_SS[0])

#observed = {"a":new_LTT, "b":observed_SS[7], "c":observed_SS[6]}
#observed = {"a":top, "b":observed_SS[7], "c":observed_SS[6]} 
observed = {"a":bl, "b":observed_SS[7], "c":observed_SS[6]}

def distance(x,y): #inputs are the dictionaries
    
    sim_a_vector = x['a']
    obs_a_vector = y['a']
    
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
        print("error - not the same len for a" + str(len(new_a_x) + " " + str(len(new_a_y))
        return
    
    dist_b = np.linalg.norm(sim_b - obs_b)
    dist_c = np.linalg.norm(sim_c - obs_c)
    
    final_b = dist_b/obs_b
    final_c = dist_c/obs_c

    dist = dist_a + final_b + final_c
    
    return dist
    


def simulate_pyabc(parameter):
    result = simulate_epidemic(**parameter)
    return {"a":result[0], "b":result[1], "c":result[2]} #so this returns the sample as a dictionary

parameters = dict(a=(0.5,1), b=(0,0.3), c=(0,0.5))

prior = pyabc.Distribution(**{key: pyabc.RV("uniform", a, b - a) for key, (a,b) in parameters.items()})


pool = ThreadPoolExecutor(max_workers=48)
sampler = ConcurrentFutureSampler(pool)

abc = pyabc.ABCSMC(simulate_pyabc, prior, distance, sampler=sampler)

#db_path = ("sqlite:///ltt.db")
db_path = ("sqlite:///topology.db")
#db_path = ("sqlite:///branch_lens.db")

abc_id = abc.new(db_path, observed)


history = abc.run(max_nr_populations=10, minimum_epsilon=0.3)







