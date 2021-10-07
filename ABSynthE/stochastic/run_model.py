from collections import defaultdict
import sys
import absynthe.set_up.index_functions import *
import absynthe.stochastic.tree_simulator as tree_sim
from absynthe.stochastic.epidemic_function import *


def run_model(config):
    
    iteration_count = -1
    
    for i in range(config["number_model_iterations"]):
    
        ##writing to file and screen
        iteration_count += 1
        if iteration_count%config["log_every"] == 0:
            config["write_file"] = True
        else:
            config["write_file"] = False
        if iteration_count%10 == 0:
            sys.stdout.write(f'{iteration_count} runs completed')

        epidemic_config = index_functions.make_data_structures(config)        
        
        ###Making index case###
        epidemic_config = index_functions.make_index_case(config, epidemic_config)
        
        if write_file: #check that info file gets written to in the epidemic run like I think it does
            config["info_file"] = file_functions.prep_info_file(config["output_directory"]), index_case_individual, iteration_count)
        
        ###Run the epidemic###
        epidemic_config = run_epidemic(0, config, epidemic_config)
    
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
        
        if epidemic_stopped and not config["write_file"]: #ie if it's been stopped because it's reached the day or case cap but not already being written
            
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

        
        if config["write_file"] or epidemic_stopped:
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
            
def run_epidemic(start_day, config, epidemic_config):
    
    epidemic_config["epidemic_stopped"] = False
    day_count = 0
    
    try: #replace with an if statement
        for day, case_list in epidemic_config["day_dict"].items():    
            day_count += 1
            if config["capped"]:
                if day_count > config["day_limit"] or len(epidemic_config["case_dict"]) > config["case_limit"]: 
                    sys.stdout.write("Epidemic has reached day limit or case limit\n")
                    epidemic_config["epidemic_stopped"]  = True
                    #then need to assign final cases an individual if it's reached the case_limit, but not if it's reached the day limit. At the moment, they'll be removed
                    return epidemic_config
            if len(case_list) != 0 and day >= start_day: #If there are new cases on this day
                for focal_case in case_list:
                    parent = epidemic_config["case_dict"][focal_case.parent]#Gets the individual object of parent (intialised last time) from the case dictionary using the case object

                    #May need to check that the new individual is coming out of this
                    assignment = focal_case.who_am_I(parent, day, config, epidemic_config) #Assign the current case id to an individual

                    #using none and false different is not good - work on that
                    if assignment == None and day != 0: #If individual is already infected
                        pass
                    elif assignment == False: #If there is no-one left
                        sys.stdout.write("All members of the population infected\n")
                        epidemic_config["epidemic_stopped"] = True
                        return epidemic_config
                    
                    #Test that this only comes up when type = Individual
                    else: #There is a successful assignation to a specific individual

                        epidemic_config["onset_times"].append(day)

                        #adding to relevant structures to fill in information about the new case
                        focal_individual = epidemic_config["case_dict"][focal_case] #This is an Individual object got by who_am_I
                        focal_individual.parent = parent #Gives Individual the same parent as the Case object for recording
                        parent.children.append(focal_individual)

                        #when transmission actually happened - how is this used? might be better to have a dictionary.
                        epidemic_config["transmission_dict"][focal_individual.unique_id] = [focal_individual.parent.unique_id, day, (day+focal_individual.incubation_day)]
                        
                        #starting new instance for child dict - this could happen in the class definition?
                        epidemic_config["child_dict"][focal_individual.unique_id] = []
                        epidemic_config["child_dict"][focal_individual.parent.unique_id].append(focal_individual.unique_id)

                        epidemic_config["nodes"].append(focal_individual.unique_id)

                        if config["write_file"]:
                            config["info_file"].write(f'{focal_individual.unique_id},{focal_individual.parent.unique_id},{focal_individual.hh},{focal_individual.comm},{focal_individual.dist},{day},{day+focal_individual.incubation_day},{day + focal_individual.incubation_day}\n')

                            if focal_individual.dist != focal_individual.parent.dist:
                                epidemic_config["dist_mvmt"][focal_individual.dist,focal_individual.parent.dist].append(day)
                            if focal_individual.comm != focal_individual.parent.comm: #This is going to include between district movements too remember
                                epidemic_config["ch_mvmt"][focal_individual.comm,focal_individual.parent.comm].append(day)
                        
                        epidemic_config["districts_present"].add(focal_individual.dist)
                        epidemic_config["chiefdoms_present"].add(focal_individual.comm)

                        poss_case_dict = focal_individual.get_possible_cases() #Gives dict of contact_level: number of people
                        
                        for level, number in poss_case_dict.items():
                            for person in range(number):
                                day_inf_output = focal_individual.when_infected(day, person, cdf_len_set, cdf_array)[0]
                                
                                if not day_inf_output:
                                    #print("Finished infection first")
                                    pass

                                else: #Need to test this as well - could return "Type"
                                    new_case = Case(len(case_dict), level, focal_case)
                                    epidemic_config["case_dict"][new_case] = None
                                    epidemic_config["day_dict"][day_inf_output].append(new_case)


    except RuntimeError:
        if not config["capped"]:
            
            original_length = len(epidemic_config["day_dict"])
            new_start = day #Should start again from when the error was thrown, so for now will recalculate all the infecteds for that day
            for i in range(1000):
                epidemic_config["day_dict"][original_length + i] = []

            run_epidemic(new_start, epidemic_config)
        else:
            pass
        

    return epidemic_config






