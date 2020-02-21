from collections import defaultdict
import random
from multiprocessing.pool import ThreadPool
import time
import os


import Tree_simulator_fitting as cts
from Simulate_epidemic_fitting import *
from vector_comparisons import *
from movement_fitting import *
from make_contact_dicts_chiefdom import *
import distribution_functions
 

#iteration_number_outside = 10000

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


size_file = open(dropbox_path + results_path + str(run_number) + "/epidemic_sizes.csv", 'w')

size_file.write("a_value, b_value, c_value, size,\n")

distributions = distribution_functions.define_distributions()

print("Defining contact structures")
contact_structure = make_contact_dicts(dropbox_path)

district_list = ["bo", 'bombali', 'bonthe', 'kailahun', 'kambia', 'kenema', 'koinadugu', 'kono', 'moyamba', 'portloko', 'pujehun', 'tonkolili', 'westernarearural', 'westernareaurban']

ch_list = ['badjia', 'bagbo', 'bagbwe', 'baoma', 'bumpe-gao', 'gbo', 'jaiamabongor', 'kakua', 'komboya', 'lugbu', 'niawalenga', 'selenga', 'tikonko', 'valunia', 'wonde', 'bendu', 'bum', 'dema', 'imperri', 'jong', 'kpanda', 'kwamebaikrim', 'nongobabullom', 'sittia', 'sogbini', 'yawbeko', 'biriwa', 'bombalisebora', 'gbantikamaranka', 'gbendemungowahun', 'libeisaygahun', 'magbaiambandowahun', 'makarigbanti', 'safrokolimba', 'pakimasabong', 'sandaloko', 'sandatenraren', 'sellalimba', 'tambakha', 'brimaia', 'gbinledixing', 'magbema', 'mambolo', 'masungbala', 'samu', 'tonkolimba', 'dea', 'jaluahun', 'jawei', 'kissikama', 'kissiteng', 'kissitongi', 'luawa', 'malema', 'mandu', 'pejebongre', 'pejewest', 'penguia', 'upperbambara', 'yawei', 'dama', 'dodo', 'guara', 'goramamende', 'kanduleppiama', 'koya', 'languramaya', 'lowerbambara', 'malegohun', 'niawa', 'nomo', 'nongowa', 'simbaru', 'smallbo', 'tunkia', 'wando', 'dembeliasikunia', 'diang', 'folsaba', 'kasunko', 'mongo', 'neya', 'nieni', 'sengbe', 'sulima', 'wara-warabafodea', 'wara-warayagala', 'fiama', 'gbanekandor', 'gbane', 'gbense', 'goramakono', 'kamara', 'lei', 'mafindor', 'nimikoro', 'nimiyama', 'sandor', 'soa', 'tankoro', 'toli', 'bagruwa', 'banta', 'bumpeh', 'dasse', 'fakunya', 'kargboro', 'kaiyamba', 'kamajei', 'kongbora', 'kori', 'kowa', 'ribbi', 'timidale', 'upperbanta', 'burehkasseh', 'buyaromende', 'debia', 'kaffubullom', 'lokomasama', 'maforki', 'marampa', 'masimera', 'koya', 'sandamagbolontor', 'tmsafroko', 'barri', 'gallinesperri', 'kpaka', 'kpangakabonde', 'makpele', 'malen', 'manosakrim', 'kpangakrim', 'peje', 'sorogbema', 'sowa', 'yakemokepukumukrim', 'gbonkolenken', 'kafesimiria', 'kalansongoia', 'kholifamabang', 'kholifarowalla', 'kunike', 'kunikebarina', 'malalmara', 'sambaia', 'tane', 'yoni', 'westernurban', 'westernrural']

dist_keys = []
ch_keys = []

for item1 in district_list:
    for item2 in district_list:
        if item1 != item2:
            dist_keys.append((item1, item2))

for item3 in ch_list:
    for item4 in ch_list:
        if item3 != item4:
            ch_keys.append((item3, item4))


##ABC setup###

observed_SS = get_observed_SS()
observed_dist = 93
observed_ch = 30

iterations_per_value = 1 #So this might actually be only one, and we change a each time

#For now, we have put this to 10% difference, can be expanded or contracted for accuracy/speed
branch_threshold = 0.1
top_threshold = 0.1
LTT_stat_threshold = 0.1
LTT_point_threshold = 0.1

rejection_threshold_b = 9 
rejection_threshold_c = 13 

accepted = []


def abc_algorithm(accepted):
    
    accepted_file = open(dropbox_path + results_path + str(run_number) + "/accepted_parameter_values.csv", 'w')
    accepted_bc = open(dropbox_path + results_path + str(run_number) + "/accepted_b_and_c_values.csv", 'w')
   
    print("Starting ABC")
    
    count = 0
    count2 = 0
    N = 1000
    
    ######FOR WHICH ONE WE'RE FITTING ON - remember to check above and change results path####
    LTT = True
    branches = False
    topology = False
    ######
    
    function = random.uniform
    
    while len(accepted) < N and count < 50000000: 
        
        #These two bits are to log every 50 successful parameters
        if len(accepted)%50 == 1:
            count2 == 0
        if len(accepted) % 50 == 0 and N != 0 and count2 == 0:
            count2 += 1
            accepted_file.close()
            accepted_file = open(dropbox_path + results_path + str(run_number) + "/accepted_parameter_values.csv", 'a')
      
        
        count += 1
        
        if count % 100 == 0:
            print("parameters tried = " + str(count))

        a = function(0.5,1) #Try this for now
        b = function(0,0.2) #Maybe make these much narrower, like 0 to 0.1
        c = function(0,0.5) #I estimated it to be 0.07 so this gives a bit more
      
        output = simulate_epidemic(a, b, c, LTT, iterations_per_value, distributions, contact_structure, size_file, dist_keys, ch_keys)

        if output: #So there'll only be output if the cases are between 1800 and 2800 already
            
            tree = output[0]
            district_mvmt = output[1]
            ch_mvmt = output[2]
            
            
            ch_difference = get_jumps(observed_ch, ch_mvmt)
            dist_difference = get_jumps(observed_dist, district_mvmt)
                        
            if ch_difference <= rejection_threshold_b and dist_difference <= rejection_threshold_c:
                
                print("accepted b and c!")
                accepted_bc.write(f"{a}, {b}, {c}\n")
                
                #NB ALL LTT IS GOING BACKWARDS IN TIME
                if LTT:
                    LTT_stat_diff = compare_LTT_stats(observed_SS, tree)
                    LTT_point_diff = compare_LTT_points(observed_SS, tree)
                    
                    if LTT_point_diff: #This means that there has to be lineages in the all of the coalescent intervals
                
                        if LTT_stat_diff <= LTT_stat_threshold and LTT_point_diff <= LTT_point_threshold:

                            print("accepted a value!")
                            
                            tup = (a,b,c)
                            
                            print(str(tup))

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
                        
                        
      
                         
    accepted_file.close()
    accepted_bc.close()
    
    return accepted


start = time.time()


pool = ThreadPool(24)


pool.map(abc_algorithm, (accepted,))  

size_file.close()

end = time.time()

print("completed in " + str(end-start))


