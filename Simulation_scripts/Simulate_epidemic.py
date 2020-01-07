iteration_number_outside = 50
iteration_count = -1

if iteration_count == -1:
    
    print("Successfully importing modules")
    
    #Check which of these are needed in this script
    import numpy as np
    from scipy import stats
    import scipy as sp
    import math
    import random
    from collections import defaultdict
    from collections import Counter
    import time
    from multiprocessing.pool import ThreadPool
    
    import Tree_simulator as cts
    import file_functions 
    from make_contact_dicts import *
    from individual_class import *
    from case_class import *
    import distribution_functions
    import index_functions
    from epidemic_function import *
    
    
    #import Fitting_functions as fit
    #import skygrid_prep
    
    
    #For the server
    #dropbox_path = "/localdisk/home/s1732989/ABM/"
    #results_path = "running_model/Results/no_caps/"
    
    
    dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"
    results_path = "Looping models/Results/testing_tidying/"
    
    run_number = 6
    
    capped = True
    case_limit = 250
    
    print("Defining parameters")
    
    try:
        file_functions.make_directories(dropbox_path, results_path, run_number)
    
    except FileExistsError:
        pass

    R0_output, size_output, most_recent_tip_file, length_output = file_functions.make_summary_files(dropbox_path, results_path, run_number)
    
    if capped:
        run_out_summary = file_functions.prep_runout_summary(dropbox_path, results_path, run_number)

    popn_size = 7092142
    epidemic_length = 1000
    cfr = 0.7
    sampling_percentage = 0.16
    
    district_list = ["bo", 'bombali', 'bonthe', 'kailahun', 'kambia', 'kenema', 'koinadugu', 'kono', 'moyamba', 'portloko', 'pujehun', 'tonkolili', 'westernarearural', 'westernareaurban']
    
    distributions = distribution_functions.define_distributions() 

    print("Importing dictionaries")
    
    contact_structure = make_contact_dicts(dropbox_path)

#Put these running functions in their own files
def run_model(iteration_number):
    iteration_count = -1
    
    for i in range(iteration_number):

        iteration_count += 1

        if iteration_count%5 == 0:
            write_file = True
        else:
            write_file = False

        if iteration_count%1 == 0:
            print(str(iteration_count) + " runs completed")

        original_dist_mvmt = defaultdict(list)

        for item1 in district_list:
            for item2 in district_list:
                if item1 != item2:
                    original_dist_mvmt[item1,item2] = []

        original_case_dict = {}
        original_day_dict = defaultdict(list)
        
        option_dict_countrylevel = defaultdict(list)
        option_dict_districtlevel = defaultdict(list)

        infected_individuals_set = set()

        cdf_array = []
        cdf_len_set = set()

        original_districts_present = []
        original_cluster_set = set()

        original_trans_dict = defaultdict(list)
        original_nodes = []
        
        original_onset_times = []
 


        for i in range(epidemic_length):
            original_day_dict[i] = []

        index_case_case, index_case_individual, original_case_dict, original_trans_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict = index_functions.make_index_case(contact_structure[0], cfr, distributions, original_case_dict, original_trans_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict)
        
        if write_file:

            info_file = file_functions.prep_info_file(dropbox_path, results_path, run_number, index_case_individual, iteration_count)
        
        susceptibles_left = True
        
        #Put run_epidemic in its own file if possible as it depends on other functions
        
        day_dict, case_dict, nodes, trans_dict, dist_mvmt, onset_times, districts_present, cluster_set, epidemic_capped = run_epidemic(0, original_day_dict, susceptibles_left , original_case_dict, original_trans_dict, infected_individuals_set, popn_size, option_dict_districtlevel, original_onset_times, original_nodes, original_cluster_set, cdf_len_set, cdf_array, original_districts_present, original_dist_mvmt, contact_structure, cfr, distributions, write_file, info_file, iteration_count, capped, epidemic_length, case_limit)

            
        remove_set = set()   
    
        #Removing cases that don't exist eg because the person was already infected, or because the parent had recovered/died
        for key, value in case_dict.items():
            if type(value) != Individual:
                remove_set.add(key)

        for item in remove_set:
            del case_dict[item]

        #Removing those people from the day dict
        for key, lst in day_dict.items():
            case_list = [item for item in lst if item not in remove_set] 
            day_dict[key] = case_list

        day_dict[0].append(index_case_case) #Put here so that it doesn't confuse the loop above because it has no parent AND otherwise it would get reassigned and stuff

        if epidemic_capped and not write_file: #ie if it's capped but not already being written
            
            runout_file = file_functions.prep_runout_file(dropbox_path, results_path, run_number, iteration_count)
         
            for indie in case_dict.values():

                day = trans_dict[indie.unique_id][1]
                symptoms = trans_dict[indie.unique_id][2]
                sampled = trans_dict[indie.unique_id][2]

                try:
                    runout_file.write(str(indie.unique_id) + "," + 
                    str(indie.parent.unique_id) +  "," +
                    str(indie.hh) +  "," +
                    str(indie.dist) + "," +
                    str(day) + "," +
                    str(symptoms) + "," +
                    str(sampled) +
                    "\n")
                except AttributeError:
                    runout_file.write(str(indie.unique_id) + "," + 
                    "NA" +  "," +
                    str(indie.hh) +  "," +
                    str(indie.dist) + "," +
                    str(day) + "," +
                    str(symptoms) + "," +
                    str(sampled) +
                    "\n")


        
        
        if write_file or epidemic_capped:

            tree_file, district_mvmt_file, skyline_file = file_functions.prep_other_files(dropbox_path, results_path, run_number, iteration_count)


        last_day = max(onset_times)

        if write_file or epidemic_capped:

            for key, value in dist_mvmt.items():
                if len(value) != 0:
                    district_mvmt_file.write(key[0] + "," + key[1] + "," + ",".join([str(i) for i in value]) + "\n")

            district_mvmt_file.close()

            result = cts.simulate_tree(trans_dict, nodes, sampling_percentage, last_day)
            
            if result:
                
                newick_string = result[0]
                skyline = result[1]
                tree = result[2]
                
                most_recent_tip_file.write(str(iteration_count) + "," + str(tree.most_recent_date) + "\n")

                tree_file.write(newick_string)

                tree_file.close()

                logpop_count = 0
                start_interval = 0.0
                
                for key, value in skyline.items():
                    logpop_count += 1
                    skyline_file.write(str(logpop_count) + ","
                                      + str(start_interval) + ","
                                      + str(key) + ","
                                      + str(value) + "\n")

                    start_interval = key

                skyline_file.close()
                
                if result[3]:
                    R0 = str(result[3])
                    R0_output.write(str(iteration_count) + "," + R0 + "\n")

        if write_file:
            info_file.close()
        if epidemic_capped and not write_file:
            runout_file.close()
        
        if epidemic_capped:
            run_out_summary.write(str(iteration_count) + "," + str(len(case_dict)) + "\n")
        
        length_output.write(str(iteration_count) + "," + str(last_day) + "\n")

        size_output.write(str(iteration_count) + "," + str(len(case_dict)) + "," + str(len(districts_present)) + "," + str(len(cluster_set)) + "\n")

    
            
            
pool = ThreadPool(8)

print("Running infection model")

pool.map(run_model,(iteration_number_outside,))

        
        
R0_output.close()
size_output.close()
length_output.close()
most_recent_tip_file.close()
                                  
if capped:
    run_out_summary.close()