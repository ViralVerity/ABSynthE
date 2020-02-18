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
 

iteration_number_outside = 10

print("Successfully importing modules")


##Setting things up for model running##
dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"

results_path_start = "Looping models/Results/Fitting/"

#dropbox_path = "/localdisk/home/s1732989/ABM/Fitting/"

#Which parameter set are we using
#results_path = "LTT/"
#results_path = "topology/"
#results_path = "branches/"
results_path = results_path_start + "b_and_c/"



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

#observed_SS = get_observed_SS()
observed_dist = 93
observed_ch = 30

#For now, we have put this to 10% difference, can be expanded or contracted for accuracy/speed

rejection_threshold_b = 9 
rejection_threshold_c = 13 



run_number = 1

try:
    os.mkdir(os.path.join(dropbox_path, results_path, str(run_number)))
except FileExistsError:
    pass

accepted = []

size_file = open(dropbox_path + results_path + str(run_number) + "/epidemic_sizes.csv", 'w')

size_file.write("a_value, b_value, c_value, size,\n")

accepted_c = open(dropbox_path + results_path + str(run_number) + "/accepted_c.csv", 'w')
accepted_b = open(dropbox_path + results_path + str(run_number) + "/accepted_b.csv", 'w')
accepted_both = open(dropbox_path + results_path + str(run_number) + "/accepted_both.csv", 'w')

accepted_b.write("a_value, accepted_b_value, rejected_c_value\n")
accepted_c.write("a_value, rejected_b_value, accepted_c_value\n")
accepted_both.write("a_value, b_value, c_value\n")

parameter_sets = []

with open("prelim_fit/epidemic_sizes.csv") as f:
    next(f)
    for l in f:
        toks = l.strip("\n").split(",")
        
        a = toks[0]
        b = toks[1]
        c = toks[2]
        
        parameter_sets.append((a,b,c))

def abc_algorithm(parameter_sets):
   
    print("Starting ABC")
    
    count = 0

    function = random.uniform
    
    for parameter_set in parameter_sets:
 
        count += 1
        
        if count/len(parameter_sets) % 10 == 0:
            print(str((count/len(parameter_sets))*100) + "% complete")

        a = float(parameter_set[0])
        b = float(parameter_set[1])
        c = float(parameter_set[2] )
        
        iterations_per_value = 10
        LTT = False
      
        output = simulate_epidemic(a, b, c, LTT, iterations_per_value, distributions, contact_structure, size_file, dist_keys, ch_keys)

        if output: #So there'll only be output if the cases are between 1800 and 2800 already
            
            district_mvmt = output[0]
            ch_mvmt = output[1]
            
            ch_difference = get_jumps(observed_ch, ch_mvmt)
            dist_difference = get_jumps(observed_dist, district_mvmt)
                        
            if ch_difference <= rejection_threshold_b and dist_difference <= rejection_threshold_c:
                
                accepted_both.write(f"{a}, {b}, {c}\n")
            
            elif ch_difference <= rejection_threshold_b and dist_difference > rejection_threshold_c:
                
                accepted_b.write(f"{a}, {b}, {c}\n")
                    
            elif ch_difference > rejection_threshold_b and dist_difference <= rejection_threshold_c:
                    
                accepted_c.write(f"{a}, {b}, {c}\n")
                
                
                
                
    accepted_c.close()
    accepted_b.close()          
    accepted_both.close()
    
    return accepted
                

    
start = time.time()

pool = ThreadPool(11)

pool.map(abc_algorithm, (parameter_sets,))  

size_file.close()

end = time.time()

print("completed in " + str(end-start))               
                
                
                
                
                
                
                
                