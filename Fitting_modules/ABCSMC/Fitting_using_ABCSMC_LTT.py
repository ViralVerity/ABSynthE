#Do this separately for each parameter set

from collections import defaultdict
import random
from multiprocessing.pool import ThreadPool
import time
import os

import pyabc


import Tree_simulator_fitting as cts
from Simulate_epidemic_fitting import *
from vector_comparisons import *
from movement_fitting import *
from make_contact_dicts_chiefdom import *
import distribution_functions
 

print("Successfully importing modules")


##Setting things up for model running##
#dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"

#results_path_start = "Looping models/Results/Fitting/"

dropbox_path = "/localdisk/home/s1732989/ABM/Fitting/"

#Which parameter set are we using
results_path = "LTT/"
#results_path = "topology/"
#results_path = "branches/"

run_number = 1 

try:
    os.mkdir(os.path.join(dropbox_path, results_path, str(run_number)))
except FileExistsError:
    pass


##ABC setup###

observed_SS = get_observed_SS()

#For now, we have put this to 10% difference, can be expanded or contracted for accuracy/speed
#branch_threshold = 0.1
#top_threshold = 0.1
#LTT_stat_threshold = 0.1
#LTT_point_threshold = 0.1

#rejection_threshold_b = 9 
#rejection_threshold_c = 13 

accepted = []

def distance(observed_SS, sim):
    
    coalescent_tree = sim[0]
    dist = sim[1]
    ch = sim[2]
    
    obs_dist = observed_SS[7]
    obs_ch = observed_SS[8]
    
    a = compare_LTT_stats(observed_SS, coalescent_tree)
    b = get_jumps(obs_dist, dist)
    c = get_jumps(obs_ch, ch)
    
    distance = (a+b+c)/3
    

prior = pyabc.Distribution(a=RV("uniform", 0.5, 1), b=RV("uniform",0, 0.5), c=RV("uniform", 0, 0.15)
                     
                           
abc = pyabc.ABCSMC(model, prior, distance)
                           

                        
                           
                          
                           