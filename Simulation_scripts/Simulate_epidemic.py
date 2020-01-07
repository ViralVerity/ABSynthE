iteration_number_outside = 50
iteration_count = -1

if iteration_count == -1:
    
    print("Successfully importing modules")
    
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
    
    
    #import Fitting_functions as fit
    #import skygrid_prep
    
    
    #For the server
    #dropbox_path = "/localdisk/home/s1732989/ABM/"
    #results_path = "running_model/Results/no_caps/"
    
    
    dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"
    results_path = "Looping models/Results/testing_tidying/"
    
    run_number = 5
    
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
    
    inccdf, death_cdf, recovery_cdf = distribution_functions.define_distributions()
    #Make this a list called distributions or something

    print("Importing dictionaries")
    
    agent_location, dist_to_hh, hh_to_cluster, cluster_to_hh, hh_to_ppl, cluster_to_ppl, dist_to_ppl, district_distance, district_pops = make_contact_dicts(dropbox_path)
    #Make this contact_dicts or something

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

        index_case_case, index_case_individual, original_case_dict, original_trans_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict = index_functions.make_index_case(agent_location, cfr, inccdf, death_cdf, recovery_cdf, original_case_dict, original_trans_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict)
        
        if write_file:

            info_file = file_functions.prep_info_file(dropbox_path, results_path, run_number, index_case_individual, iteration_count)
        
        susceptibles_left = True
        
        #Put run_epidemic in its own file if possible as it depends on other functions
        
        day_dict, case_dict, nodes, trans_dict, dist_mvmt, onset_times, districts_present, cluster_set, epidemic_capped = run_epidemic(0, original_day_dict, susceptibles_left , original_case_dict, original_trans_dict, infected_individuals_set, popn_size, hh_to_ppl, cluster_to_hh, option_dict_districtlevel, district_distance, dist_to_ppl, original_onset_times, original_nodes, original_cluster_set, cdf_len_set, cdf_array, original_districts_present, original_dist_mvmt, write_file, info_file, iteration_count, capped, epidemic_length, case_limit)

            
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

def run_epidemic(start_day, day_dict, susceptibles_left , case_dict, trans_dict, infected_individuals_set, popn_size, hh_to_ppl, cluster_to_hh, option_dict_districtlevel, district_distance, dist_to_ppl, onset_times, nodes, cluster_set, cdf_len_set, cdf_array, districts_present, dist_mvmt, write_file, info_file, iteration_count, capped, epidemic_length, case_limit):
    
    epidemic_capped = False
    day_count = 0
    
    try: 
        for day, case_list in day_dict.items():    
            if not susceptibles_left:
                break
            day_count += 1
            if capped:
                if day_count > epidemic_length or len(case_dict) > case_limit:
                    epidemic_capped = True
                    break
            if len(case_list) != 0 and day >= start_day: #If there are new cases on this day
                for focal_case in case_list:
                    parent = case_dict[focal_case.parent]#Gets the individual object of parent (intialised last time) from the case dictionary using the case object

                    #May need to check that the new individual is coming out of this
                    assignment = focal_case.who_am_I(infected_individuals_set, popn_size, hh_to_cluster, dist_to_hh, cluster_to_ppl, hh_to_ppl, cluster_to_hh, option_dict_districtlevel, district_distance, dist_to_ppl, case_dict, parent, day, agent_location, cfr, inccdf, death_cdf, recovery_cdf) #Assign the current case id to an individual

                    if assignment == None and day != 0: #If individual is already infected
                        #remove_set.add(focal_case) #Case doesn't exist so must be removed from day/case dict
                        pass

                    elif assignment == False: #If there is no-one left
                        print("All members of the population infected")
                        #last_day = day
                        #print("last day = " + str(last_day))
                        susceptibles_left = False
                        return day_dict, case_dict, nodes, trans_dict, dist_mvmt, onset_times, dist_present, cluster_set
                    
                    #Test that this only comes up when type = Individual
                    else: #There is a successful assignation to a specific individual

                        onset_times.append(day)

                        focal_individual = case_dict[focal_case] #This is an Individual object got by who_am_I

                        focal_individual.parent = parent #Gives Individual the same parent as the Case object for recording

                        parent.children.append(focal_individual)

                        trans_dict[focal_individual.unique_id] = [focal_individual.parent.unique_id, day, (day+ focal_individual.incubation_day)]

                        nodes.append(focal_individual.unique_id)

                        if write_file == True:
                            info_file.write(str(focal_individual.unique_id) + "," + 
                                            str(focal_individual.parent.unique_id) +  "," +
                                            str(focal_individual.hh) +  "," +
                                            str(focal_individual.dist) + "," +
                                            str(day) + "," +
                                            str(day + focal_individual.incubation_day) + "," +
                                            str(day + focal_individual.incubation_day) +
                                            "\n")

                            if focal_individual.dist != focal_individual.parent.dist:
                                dist_mvmt[focal_individual.dist,focal_individual.parent.dist].append(day)


                        if focal_individual.dist not in districts_present:
                            districts_present.append(focal_individual.dist)

                        if focal_individual.comm not in cluster_set:
                            cluster_set.add(focal_individual.comm)

                        poss_case_dict = focal_individual.get_possible_cases() #Gives dict of contact_level: number of people
                        
                        for level, number in poss_case_dict.items():
                            for person in range(number):
                                day_inf_output = focal_individual.when_infected(day, person, cdf_len_set, cdf_array)[0]
                                
                                if day_inf_output == None:
                                    #print("Finished infection first")
                                    pass

                                else: #Need to test this as well - could return "Type"
                                    new_case = Case(len(case_dict), level, focal_case)
                                    case_dict[new_case] = None
                                    day_dict[day_inf_output].append(new_case)


    except RuntimeError:
        if not capped:
            print("Adding more days to" + str(iteration_count))
            original_length = len(day_dict)
            new_start = day #Should start again from when the error was thrown, so for now will recalculate all the infecteds for that day
            for i in range(1000):
                day_dict[original_length + i] = []

            run_epidemic(new_start, day_dict, susceptibles_left, case_dict, trans_dict, infected_individuals_set, popn_size, hh_to_ppl, cluster_to_hh, option_dict_districtlevel, district_distance, dist_to_ppl, onset_times, nodes, cluster_set, cdf_len_set, cdf_array, districts_present, dist_mvmt, write_file, info_file, iteration_count, capped, epidemic_length, case_limit)
        else:
            pass
        

    return day_dict, case_dict, nodes, trans_dict, dist_mvmt, onset_times, districts_present, cluster_set, epidemic_capped




             
            
            
pool = ThreadPool(8)

print("Running infection model")

pool.map(run_model,(iteration_number_outside,))

        
        
R0_output.close()
size_output.close()
length_output.close()
most_recent_tip_file.close()
                                  
if capped:
    run_out_summary.close()