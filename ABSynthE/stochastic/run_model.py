from collections import defaultdict
import absynthe.set_up.index_functions import *
import absynthe.stochastic.tree_simulator as tree_sim


def run_model(config):
    
    iteration_count = -1
    
    for i in range(config["number_model_iterations"]):
    
        ##writing to file and screen
        iteration_count += 1
        if iteration_count%config["log_every"] == 0:
            write_file = True
        else:
            write_file = False
        if iteration_count%10 == 0:
            sys.stdout.write(f'{iteration_count} runs completed')

        epidemic_config = index_functions.make_data_structures(config)        
        
        ###Making index case###
        epidemic_config = index_functions.make_index_case(config, epidemic_config)
        
        if write_file: #check that info file gets written to in the epidemic run like I think it does
            config["info_file"] = file_functions.prep_info_file(config["output_directory"]), index_case_individual, iteration_count)
        
        ###Run the epidemic###
        #this function needs work in terms of the config - also want to look at all these data_structures - part of config or something?
        config["susceptibles_left"] = True #needs to be external because run_epidemic is recursive
        day_dict, case_dict, nodes, trans_dict, child_dict, dist_mvmt, ch_mvmt, onset_times, districts_present, chiefdom_set, epidemic_stopped = run_epidemic(config, 0, day_dict, case_dict, trans_dict, child_dict, infected_individuals_set, option_dict_districtlevel, onset_times, nodes, chiefdom_set, cdf_len_set, cdf_array, districts_present, dist_mvmt, ch_mvmt, iteration_count)

    
        ###Removing cases that don't exist eg because the person was already infected, or because the parent had recovered/died###
        remove_set = set()   
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
        last_day = max(onset_times)

        ###Getting results and writing to file###
        # could separate these out and put them in the file functions file? Should probably just have a write to file option
        size = len(case_dict)
        dists = len(districts_present)
        chiefdoms = len(chiefdom_set)

        config["files"]["length_output"].write(f"{iteration_count},{last_day}\n")
        config["files"]["size_output"].write(f"{iteration_count},{size},{dists},{chiefdoms}\n")
        
        if epidemic_stopped:
            config["files"]["run_out_summary"].write(f"{iteration_count},{size}\n")
        
        if epidemic_stopped and not write_file: #ie if it's been stopped because it's reached the day or case cap but not already being written
            
            runout_file = file_functions.prep_runout_file(config["output_directory"], iteration_count)
            for indie in case_dict.values():

                day = trans_dict[indie.unique_id][1]
                symptoms = trans_dict[indie.unique_id][2]
                sampled = trans_dict[indie.unique_id][2]
                
                try:
                    runout_file.write(f"{indie.unique_id},{indie.parent.unique_id},{indie.hh},{indie.dist},{day},{symptoms},{sampled},\n")

                except AttributeError:
                    runout_file.write(f"{indie.unique_id},NA,{indie.hh},{indie.dist},{day},{symptoms},{sampled},\n")
            runout_file.close()

        
        if write_file or epidemic_stopped:
            tree_file, district_mvmt_file, ch_mvmt_file, skyline_file, ltt_file = file_functions.prep_other_files(config["output_directory"], iteration_count)
            
            for district_pair, count_list in dist_mvmt.items(): #what is the value here - is it counts or days that they're happening on?
                if len(count_list) != 0:
                    counts = ",".join([str(i) for i in value])
                    district_mvmt_file.write(f'{district_pair[0]},{district_pair[1]},{counts}"\n"')
            district_mvmt_file.close()
            
            for ch_pair, count_list in ch_mvmt.items():
                if len(count_list) != 0:
                    counts = ",".join([str(i) for i in value])
                    ch_mvmt_file.write(f'{ch_pair[0]},{ch_pair[1]},{counts}\n')      
            ch_mvmt_file.close()
            
            #all this tree sim still needs tidying
            #add in if statements for if skyline, if ltts
            result = tree_sim.simulate_tree(trans_dict, child_dict, nodes, sampling_percentage, last_day) 
            if result:
                newick_string = result[0]
                #skyline = result[1]
                tree = result[1]
                #lineages_through_time = result[5]
                
                connfig["files"]["most_recent_tip_file"].write(f'{iteration_count},{tree.most_recent_date}\n')

                #gets opened above - possibly not wanted? end pu with a lot of empty files?
                tree_file.write(newick_string)
                tree_file.close()

                #start_interval = 0.0
                #I think this is all commented out so that it doesn't slow down post-fitting things - should still be tidied up though.
                if config["make_skyline"]: 
                    logpop_count = 0
                    for key, value in skyline.items():
                        logpop_count += 1
                        skyline_file.write(f"{logpop_count},{key[0]},{key[1]},{value}\n")
                    skyline_file.close()
                
                if config["make_ltt"]: 
                    lineage_count = 0
                    for k,v in lineages_through_time.items():
                        lineage_count += 1
                        ltt_file.write(f"{lineage_count},{k[0]},{k[1]},{v}\n")
                    ltt_file.close()
                
                if result[2]:
                    R0 = str(result[2])
                    config["files"]["R0_output"].write(f"{iteration_count},{R0}\n")
                    

        if write_file: #is there a way to check if a file is open?
            config["info_file"].close()
            







            
            
