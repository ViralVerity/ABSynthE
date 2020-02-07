
from collections import defaultdict
from Simulate_epidemic_fitting import *
from vector_comparisons import *
import random
from multiprocessing.pool import ThreadPool
import time

import sys
 
sys.path.insert(1, "../Simulation_scripts/")

#iteration_number_outside = 10000

print("Successfully importing modules")


import Tree_simulator as cts
import file_functions 
from make_contact_dicts_chiefdom import *

import distribution_functions


##Setting things up for model running##
dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"
results_path = "Looping models/Results/Fitting/test/"

size_file = open(dropbox_path + results_path + "epidemic_sizes.csv", 'w')

size_file.write("a_value,size,\n")

distributions = distribution_functions.define_distributions()

print("Defining contact structures")
contact_structure = make_contact_dicts(dropbox_path)


##ABC setup###
observed_SS = get_observed_SS()

iterations_per_value = 1 #So this might actually be only one, and we change a each time

rejection_threshold = 80
#rejection_threshold_other = ?? #play with these to get a good value. May be different for each set

rejection_threshold_b = 13 #For now - check with the new tree for these two. They are the HPDs
rejection_threshold_c = 7

accepted = []

run_number = 1 #iterate upwards as we go through a values

###Define rejection algorithm so that can be multithreaded###


def abc_algorithm(accepted):
    
    print("Starting ABC")
    
    count = 0
    N = 100
    
    while len(accepted) < N and count < 50000: 
      
        count += 1
        
        if count % 100 == 0:
            print("parameters tried = " + str(count))

        a = random.uniform(0.5,1) #Try this for now
      
        output = simulate_epidemic(a, iterations_per_value, distributions, contact_structure, size_file)
        
        if output: #So there'll only be output if the cases are between 1800 and 2800 already
            
            tree = output[0]
            
            difference = get_tip_difference(observed_SS, tree) #Just keeping this in to test, but can take it out in a bit
            
            if difference <= rejection_threshold:

                #Comment these out depending on what we are interested in for that run#
                #branch_diff = compare_BL(observed_SS, tree)
                #top_diff = compare_topology(observed_SS, tree)
                #LTT_stat_diff = compare_LTT_stats(observed_SS, tree)
                #LTT_point_diff = compare_LTT_points(observed_SS, tree)

                #if branch_diff <= branch_threshold:
                #if top_diff <= top_threshold:
                #if LTT_stat_diff <= LTT_stat_threshold:
                #if LTT_point_diff <= LTT_point_threshold:

                #print("accepted a value")
                accepted.append(a)
             

        
    return accepted


start = time.time()

pool = ThreadPool(8)

pool.map(abc_algorithm, (accepted,))  

print(accepted)


size_file.close()

end = time.time()

print("completed in " + str(end-start))


