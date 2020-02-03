
from collections import defaultdict
from Simulate_epidemic_fitting import *


import branch_length_parameters as BL
import topology_set as TOP


import sys
 
sys.path.insert(1, "../Simulation_scripts/")

#iteration_number_outside = 10000

print("Successfully importing modules")


import Tree_simulator as cts
import file_functions 
from make_contact_dicts import *

import distribution_functions



dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"
results_path = "Looping models/Results/Fitting/test/"

size_file = open(dropbox_path + results_path + "epidemic_sizes.csv", 'w')

size_file.write("a_value,size,\n")

distributions = distribution_functions.define_distributions() 

contact_structure = make_contact_dicts(dropbox_path)


iterations_per_value = 1 #So this might actually be only one, and we change a each time

a = ??

run_number = 1 #iterate upwards as we go through a values


tree = simulate_epidemic(a, iterations_per_value, distributions, contact_structure, size_file)


size_file.close()
