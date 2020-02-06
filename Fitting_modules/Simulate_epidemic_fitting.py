
from multiprocessing.pool import ThreadPool
from collections import defaultdict

from epidemic_function_fitting import *
import index_functions_fitting

import sys
sys.path.insert(1, "../Simulation_scripts/")

import Tree_simulator as cts


def simulate_epidemic(a, iteration_number_outside, distributions, contact_structure, size_file):

    #pool = ThreadPool(4)

    #print("Running infection model")

    #pool.map(run_model,(a,iteration_number_outside, distributions, contact_structure, dropbox_path, results_path,))   

    capped = True

    result = run_model(a,iteration_number_outside, distributions, contact_structure, capped, size_file)

    return result

def run_model(a,iteration_number, distributions, contact_structure, capped, size_file):
    
    
    case_limit = 10000 #This may not be high enough, so we'll need to check if the epidemics are big enough
    popn_size = 7092142
    epidemic_length = 148
    cfr = 0.7
    sampling_percentage = 0.16
    

    district_list = ["bo", 'bombali', 'bonthe', 'kailahun', 'kambia', 'kenema', 'koinadugu', 'kono', 'moyamba', 'portloko', 'pujehun', 'tonkolili', 'westernarearural', 'westernareaurban']
        

    iteration_count = -1
    
    for i in range(iteration_number):
    
    ##Setting things up for running###
        iteration_count += 1

        original_dist_mvmt = defaultdict(list)

        for item1 in district_list:
            for item2 in district_list:
                if item1 != item2:
                    original_dist_mvmt[item1,item2] = []

        original_case_dict = {}
        original_day_dict = defaultdict(list)
        option_dict_districtlevel = defaultdict(list)
        infected_individuals_set = set()
        cdf_array = []
        cdf_len_set = set()
        original_districts_present = []
        original_cluster_set = set()
        original_trans_dict = defaultdict(list)
        original_child_dict = defaultdict(list)
        original_nodes = []
        original_onset_times = []
 

        for i in range(epidemic_length):
            original_day_dict[i] = []
        
        ###Making index case###
        index_case_case, index_case_individual, original_case_dict, original_trans_dict, original_child_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict = index_functions_fitting.make_index_case(contact_structure[0], cfr, distributions, original_case_dict, original_trans_dict, original_child_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict)
        
        
        susceptibles_left = True
        
        
        ###Run the epidemic###
        day_dict, case_dict, nodes, trans_dict, child_dict, dist_mvmt, onset_times, districts_present, cluster_set, epidemic_capped = run_epidemic(0, original_day_dict, susceptibles_left , original_case_dict, original_trans_dict, original_child_dict, infected_individuals_set, popn_size, option_dict_districtlevel, original_onset_times, original_nodes, original_cluster_set, cdf_len_set, cdf_array, original_districts_present, original_dist_mvmt, contact_structure, cfr, distributions, iteration_count, capped, epidemic_length, case_limit, a)

        
        remove_set = set()   
    
        ###Removing cases that don't exist eg because the person was already infected, or because the parent had recovered/died###
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

        ###Getting results and writing to file###

        last_day = max(onset_times)
        
        if len(case_dict) > 1800 and len(case_dict) < 2800: #Conditioning to speed up ABC
            
            size = len(case_dict)

            size_file.write(f"{a},{size}\n")

            result = cts.simulate_tree(trans_dict, child_dict, nodes, sampling_percentage, last_day)
            
            return result


        
        else:
            return

            


