def run_model(iteration_number):
    
    iteration_count = -1
    
    for i in range(iteration_number):
    
    ##Setting things up for running###
        iteration_count += 1

        if iteration_count%5 == 0:
            write_file = True
        else:
            write_file = False

        if iteration_count%10 == 0:
            print(str(iteration_count) + " runs completed")

        original_dist_mvmt = defaultdict(list)
        original_ch_mvmt = defaultdict(list)

        for item1 in district_list:
            for item2 in district_list:
                if item1 != item2:
                    original_dist_mvmt[item1,item2] = []
                    
        for item1 in ch_list:
            for item2 in ch_list:
                if item1 != item2:
                    original_ch_mvmt[item1, item2] = []

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
        index_case_case, index_case_individual, original_case_dict, original_trans_dict, original_child_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict = index_functions.make_index_case(contact_structure[0], cfr, distributions, original_case_dict, original_trans_dict, original_child_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict)
        
        if write_file:

            info_file = file_functions.prep_info_file(dropbox_path, results_path, run_number, index_case_individual, iteration_count)
        
        susceptibles_left = True
         
        ###Run the epidemic###
        day_dict, case_dict, nodes, trans_dict, child_dict, dist_mvmt, ch_mvmt, onset_times, districts_present, cluster_set, epidemic_capped = run_epidemic(0, original_day_dict, susceptibles_left , original_case_dict, original_trans_dict, original_child_dict, infected_individuals_set, popn_size, option_dict_districtlevel, original_onset_times, original_nodes, original_cluster_set, cdf_len_set, cdf_array, original_districts_present, original_dist_mvmt, original_ch_mvmt, contact_structure, cfr, distributions, write_file, info_file, iteration_count, capped, epidemic_length, case_limit, SDB_start, SDB_success)

        
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
        
        if epidemic_capped and not write_file: #ie if it's capped but not already being written
            
            runout_file = file_functions.prep_runout_file(dropbox_path, results_path, run_number, iteration_count)
         
            for indie in case_dict.values():

                day = trans_dict[indie.unique_id][1]
                symptoms = trans_dict[indie.unique_id][2]
                sampled = trans_dict[indie.unique_id][2]
                
                try:
                    runout_file.write(f"{indie.unique_id},{indie.parent.unique_id},{indie.hh},{indie.dist},{day},{symptoms},{sampled},\n")

                except AttributeError:
                    runout_file.write(f"{indie.unique_id},NA,{indie.hh},{indie.dist},{day},{symptoms},{sampled},\n")

        
        if write_file or epidemic_capped:
            tree_file, district_mvmt_file, ch_mvmt_file, skyline_file, ltt_file = file_functions.prep_other_files(dropbox_path, results_path, run_number, iteration_count)

        last_day = max(onset_times)

        if write_file or epidemic_capped:

            for key, value in dist_mvmt.items():
                if len(value) != 0:
                    district_mvmt_file.write(key[0] + "," + key[1] + "," + ",".join([str(i) for i in value]) + "\n")

            district_mvmt_file.close()
            
            for key, value in ch_mvmt.items():
                if len(value) != 0:
                    ch_mvmt_file.write(key[0] + "," + key[1] + "," + ",".join([str(i) for i in value]) + "\n")
                    
            ch_mvmt_file.close()
            
            result = cts.simulate_tree(trans_dict, child_dict, nodes, sampling_percentage, last_day)
            
            if result:
                
                newick_string = result[0]
                skyline = result[1]
                tree = result[2]
                lineages_through_time = result[5]
                
                most_recent_tip_file.write(str(iteration_count) + "," + str(tree.most_recent_date) + "\n")

                tree_file.write(newick_string)

                tree_file.close()

                logpop_count = 0
                #start_interval = 0.0
                
                for key, value in skyline.items():
                    logpop_count += 1
                    
                    skyline_file.write(f"{logpop_count},{key[0]},{key[1]},{value}\n")

                skyline_file.close()
                
                lineage_count = 0
                
                for k,v in lineages_through_time.items():
                    lineage_count += 1
                    ltt_file.write(f"{lineage_count},{k[0]},{k[1]},{v}\n")
                
                if result[3]:
                    R0 = str(result[3])
                    
                    R0_output.write(f"{iteration_count},{R0}\n")
                    

        if write_file:
            info_file.close()
        if epidemic_capped and not write_file:
            runout_file.close()
        
        if epidemic_capped:
            size = len(case_dict)
            
            run_out_summary.write(f"{iteration_count},{size}\n")
        
        size = len(case_dict)
        dists = len(districts_present)
        clusters = len(cluster_set)
        
        length_output.write(f"{iteration_count},{last_day}\n")

        size_output.write(f"{iteration_count},{size},{dists},{clusters}\n")

    
