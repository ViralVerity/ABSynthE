
from collections import defaultdict
from Simulate_epidemic_fitting import *
from vector_comparisons import *
from movement_fitting import *
import random
from multiprocessing.pool import ThreadPool
import time

import Tree_sumulator_fitting as cts



import sys
 
sys.path.insert(1, "../Simulation_scripts/")

#iteration_number_outside = 10000

print("Successfully importing modules")

import file_functions 
from make_contact_dicts_chiefdom import *

import distribution_functions


##Setting things up for model running##
dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"
results_path = "Looping models/Results/Fitting/parameter_values/"

size_file = open(dropbox_path + results_path + "epidemic_sizes.csv", 'w')

size_file.write("a_value, b_value, c_value, size,\n")

distributions = distribution_functions.define_distributions()

print("Defining contact structures")
contact_structure = make_contact_dicts(dropbox_path)

district_list = ["bo", 'bombali', 'bonthe', 'kailahun', 'kambia', 'kenema', 'koinadugu', 'kono', 'moyamba', 'portloko', 'pujehun', 'tonkolili', 'westernarearural', 'westernareaurban']

ch_list = ['badjia', 'bagbo', 'bagbwe', 'baoma', 'bumpe-gao', 'gbo', 'jaiamabongor', 'kakua', 'komboya', 'lugbu', 'niawalenga', 'selenga', 'tikonko', 'valunia', 'wonde', 'bendu', 'bum', 'dema', 'imperri', 'jong', 'kpanda', 'kwamebaikrim', 'nongobabullom', 'sittia', 'sogbini', 'yawbeko', 'biriwa', 'bombalisebora', 'gbantikamaranka', 'gbendemungowahun', 'libeisaygahun', 'magbaiambandowahun', 'makarigbanti', 'safrokolimba', 'pakimasabong', 'sandaloko', 'sandatenraren', 'sellalimba', 'tambakha', 'brimaia', 'gbinledixing', 'magbema', 'mambolo', 'masungbala', 'samu', 'tonkolimba', 'dea', 'jaluahun', 'jawei', 'kissikama', 'kissiteng', 'kissitongi', 'luawa', 'malema', 'mandu', 'pejebongre', 'pejewest', 'penguia', 'upperbambara', 'yawei', 'dama', 'dodo', 'guara', 'goramamende', 'kanduleppiama', 'koya', 'languramaya', 'lowerbambara', 'malegohun', 'niawa', 'nomo', 'nongowa', 'simbaru', 'smallbo', 'tunkia', 'wando', 'dembeliasikunia', 'diang', 'folsaba', 'kasunko', 'mongo', 'neya', 'nieni', 'sengbe', 'sulima', 'wara-warabafodea', 'wara-warayagala', 'fiama', 'gbanekandor', 'gbane', 'gbense', 'goramakono', 'kamara', 'lei', 'mafindor', 'nimikoro', 'nimiyama', 'sandor', 'soa', 'tankoro', 'toli', 'bagruwa', 'banta', 'bumpeh', 'dasse', 'fakunya', 'kargboro', 'kaiyamba', 'kamajei', 'kongbora', 'kori', 'kowa', 'ribbi', 'timidale', 'upperbanta', 'burehkasseh', 'buyaromende', 'debia', 'kaffubullom', 'lokomasama', 'maforki', 'marampa', 'masimera', 'koya', 'sandamagbolontor', 'tmsafroko', 'barri', 'gallinesperri', 'kpaka', 'kpangakabonde', 'makpele', 'malen', 'manosakrim', 'kpangakrim', 'peje', 'sorogbema', 'sowa', 'yakemokepukumukrim', 'gbonkolenken', 'kafesimiria', 'kalansongoia', 'kholifamabang', 'kholifarowalla', 'kunike', 'kunikebarina', 'malalmara', 'sambaia', 'tane', 'yoni', 'westernurban', 'westernrural']


original_dist_mvmt = defaultdict(list)
original_ch_mvmt = defaultdict(list)

for item1 in district_list:
    for item2 in district_list:
        if item1 != item2:
            original_dist_mvmt[item1,item2] = []
            
for item3 in ch_list:
    for item4 in ch_list:
        if item3 != item4:
            original_ch_mvmt[item3, item4] = []
            
            

original_case_dict = {}
original_day_dict = defaultdict(list)
option_dict_districtlevel = defaultdict(list)
infected_individuals_set = set()
cdf_array = []
cdf_len_set = set()
#original_districts_present = []
#original_cluster_set = set()
original_trans_dict = defaultdict(list)
original_child_dict = defaultdict(list)
original_nodes = []
original_onset_times = []



epidemic_length = 148

for i in range(epidemic_length):
    original_day_dict[i] = []


##ABC setup###

observed_SS = get_observed_SS()
observed_dist = ??
observed_ch = ??

iterations_per_value = 1 #So this might actually be only one, and we change a each time

branch_threshold = ?
top_threshold = ?
LTT_stat_threshold = ?
LTT_point_threshold = ?




rejection_threshold_b = 13 #For now - check with the new tree for these two. They are the HPDs
rejection_threshold_c = 7

accepted = []

run_number = 1 




def abc_algorithm(accepted):
    
    print("Starting ABC")
    
    count = 0
    N = 1000
    
    ######FOR WHICH ONE WE'RE FITTING ON####
    LTT = True
    branches = False
    topology = False
    ######
    
    function = random.uniform
    
    while len(accepted) < N and count < 50000: 
      
        count += 1
        
        if count % 100 == 0:
            print("parameters tried = " + str(count))

        a = function(0.5,1) #Try this for now
        b = function(0,0.5) #Maybe make these much narrower, like 0 to 0.01
        c = function(0,0.5)
      
        output = simulate_epidemic(a, b, c, LTT, iterations_per_value, distributions, contact_structure, size_file, original_dist_mvmt, original_ch_mvmt, original_case_dict, original_day_dict, option_dict_district_level, infected_individuals_set, cdf_array, cdf_len_set, original_trans_dict, original_child_dict, original_nodes, original_onset_times)

        if output: #So there'll only be output if the cases are between 1800 and 2800 already
            
            if LTT:
                tree = output[0]
                LTT = output[1]
                district_mvmt = output[2]
                ch_mvmt = output[3]
            else:
                tree = output[0]
                district_mvmt = output[1]
                ch_mvmt = output[2]
            
            
            ch_difference = get_jumps(observed_ch, ch_mvmt)
            dist_difference = get_jumps(observed_dist, dist_mvmt)
                        
            if ch_difference <= rejection_b and dist_difference <= rejection_c:
                            
                if LTT:
                    LTT_stat_diff = compare_LTT_stats(observed_SS, tree)
                    LTT_point_diff = compare_LTT_points(observed_SS, tree)
                
                    if LTT_stat_diff <= LTT_stat_threshold and LTT_point_diff <= LTT_point_threshold:
                        
                        tup = (a,b,c)
                        
                        accepted.append(tup)
                        
                        accepted_file.write(f"{a}, {b}, {c}\n")
                        
                elif branches:
                    
                    branch_diff = compare_BL(observed_SS, tree)
                    
                    if branch_diff <= branch_threshold:
                        
                        tup = (a,b,c)
                        
                        accepted.append(tup)
                        
                        accepted_file.write(f"{a}, {b}, {c}\n")
                        
                elif topology:
                    
                    top_diff = compare_topology(observed_SS, tree)
                   
                    if top_diff <= top_threshold:
                        
                        tup = (a,b,c)
                        
                        accepted.append(tup)
                        
                        accepted_file.write(f"{a}, {b}, {c}\n")
                        
                        
            
    
    return accepted


start = time.time()

accepted_file = open(dropbox_path + results_path + str(run_number) + "/accepted_parameter_values.csv"

pool = ThreadPool(8)


pool.map(abc_algorithm, (accepted,))  

size_file.close()
                     
accepted_file.close()

end = time.time()

print("completed in " + str(end-start))


