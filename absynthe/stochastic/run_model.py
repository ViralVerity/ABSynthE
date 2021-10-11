from collections import defaultdict
import sys
import absynthe.set_up.index_functions as index_functions
import absynthe.set_up.file_functions as file_functions
import absynthe.stochastic.tree_simulator as tree_sim


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
            config["info_file"] = file_functions.prep_info_file(config["output_directory"], epidemic_config["index_case_individual"], iteration_count)
        
        ###Run the epidemic###
        epidemic_config = run_epidemic(0, config, epidemic_config)
    
        ###Removing cases that don't exist eg because the person was already infected, or because the parent had recovered/died###
        #NB if a case limit is set, it's possible this removal process will put the case count below the case limit.
        remove_set = set()   
        for case, assignment in epidemic_config["case_dict"].items():
            if not assignment: 
                remove_set.add(case)
        for case in remove_set:
            del epidemic_config["case_dict"][case]

        #Removing those people from the day dict
        for day, case_list in epidemic_config["day_dict"].items():
            new_case_list = [item for item in case_list if item not in remove_set] 
            epidemic_config["day_dict"][key] = new_case_list

        epidemic_config["day_dict"][0].append(epidemic_config["index_case_case"]) #Put here so that it doesn't confuse the loop above because it has no parent AND otherwise it would get reassigned and stuff
        last_day = max(epidemic_config["onset_times"])

        ###Getting results and writing to file###
        write_to_summary_files(config, epidemic_config, iteration_count, last_day)

        if epidemic_config["epidemic_stopped"] and not config["write_file"]:
            write_runout_file(config, epidemic_config)

        if epidemic_config["epidemic_stopped"] or config["write_file"]:
            record_individual_epidemic(config, epidemic_config)        

        if config["write_file"]: #is there a way to check if a file is open?
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
                    return epidemic_config
            
            if len(case_list) != 0 and day >= start_day: #If there are new cases on this day
                for focal_case in case_list:
                    parent = epidemic_config["case_dict"][focal_case.parent]#Gets the individual object of parent (intialised last time) from the case dictionary using the case object

                    #May need to check that the new individual is coming out of this
                    assignment, config, epidemic_config = focal_case.who_am_I(parent, day, config, epidemic_config) #Assign the current case id to an individual

                    if not assignment:
                        if assignment == 0 and day != 0: #If individual is already infected
                            pass
                        else: #If there is no-one left
                            sys.stdout.write("All members of the population infected\n")
                            epidemic_config["epidemic_stopped"] = True
                            return epidemic_config

                    #Test that this only comes up when type = Individual
                    else: #There is a successful assignation to a specific individual

                        epidemic_config["onset_times"].append(day)

                        #adding to relevant structures to fill in information about the new case
                        focal_individual = epidemic_config["case_dict"][focal_case] #Individual object got by who_am_I as value
                        focal_individual.parent = parent #Gives Individual the same parent as the Case object for recording
                        parent.children.append(focal_individual)

                        #when transmission actually happened - used to make the tree
                        epidemic_config["transmission_dict"][focal_individual.unique_id] = {"parent":focal_individual.parent.unique_id, "day_infected":day, "day_sampled":(day+focal_individual.incubation_day)} #sample day currently same as symptom onset
                        
                        #starting new instance for child dict - this could happen in the class definition?
                        epidemic_config["child_dict"][focal_individual.unique_id] = []
                        epidemic_config["child_dict"][focal_individual.parent.unique_id].append(focal_individual.unique_id)

                        epidemic_config["nodes"].append(focal_individual.unique_id)

                        if config["write_file"]:
                            config["info_file"].write(f'{focal_individual.unique_id},{focal_individual.parent.unique_id},{focal_individual.hh},{focal_individual.ch},{focal_individual.dist},{day},{day+focal_individual.incubation_day},{day + focal_individual.incubation_day}\n')

                            if focal_individual.dist != focal_individual.parent.dist:
                                epidemic_config["dist_mvmt"][focal_individual.dist,focal_individual.parent.dist].append(day)
                            if focal_individual.ch != focal_individual.parent.ch: #This is going to include between district movements too remember
                                epidemic_config["ch_mvmt"][focal_individual.ch,focal_individual.parent.ch].append(day)
                        
                        epidemic_config["districts_present"].add(focal_individual.dist)
                        epidemic_config["chiefdoms_present"].add(focal_individual.ch)

                        poss_case_dict = focal_individual.get_possible_cases() #Gives dict of contact_level: number of people
                        
                        for level, number in poss_case_dict.items():
                            for person in range(number):
                                day_inf_output = focal_individual.when_infected(day, person, epidemic_config["cdf_len_set"], epidemic_config["cdf_array"])[0]
                                
                                if not day_inf_output:
                                    #print("Finished infection first")
                                    pass

                                else: 
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



def write_to_summary_files(config, epidemic_config, iteration_count, last_day):

    size = len(epidemic_config["case_dict"])
    dists = len(epidemic_config["districts_present"])
    chiefdoms = len(epidemic_config["chiefdoms_present"])

    config["length_output"].write(f"{iteration_count},{last_day}\n")
    config["size_output"].write(f"{iteration_count},{size},{dists},{chiefdoms}\n")
    
    if epidemic_config["epidemic_stopped"]:
        config["run_out_summary"].write(f"{iteration_count},{size}\n")


def write_runout_file(config, epidemic_config):

    runout_file = file_functions.prep_info_file(config["output_directory"], epidemic_config["index_case_individual"],iteration_count)
        
    for individual in epidemic_config["case_dict"].values():
        #it's weird to use the transmission dict here
        day = epidemic_config["transmission_dict"][individual.unique_id]["day_sampled"]
        symptoms = epidemic_config["transmission_dict"][individual.unique_id]["day_sampled"] #for now they get sampled on the first day of symptoms
        sampled = epidemic_config["transmission_dict"][individual.unique_id]["day_sampled"] 
        
        if individual.parent:
            runout_file.write(f"{individual.unique_id},{individual.parent.unique_id},{individual.hh},{individual.dist},{day},{symptoms},{sampled},\n")
        
    runout_file.close()


def record_individual_epidemic(config, epidemic_config):

    district_mvmt_file, ch_mvmt_file = file_functions.prep_movement_files(config["output_directory"], iteration_count)
            
    for district_pair, count_list in epidemic_config["dist_mvmt"].items(): #what is the value here - is it counts or days that they're happening on?
        if len(count_list) != 0:
            counts = ",".join([str(i) for i in value])
            district_mvmt_file.write(f'{district_pair[0]},{district_pair[1]},{counts}"\n"')
    district_mvmt_file.close()
    
    for ch_pair, count_list in epidemic_config["ch_mvmt"].items():
        if len(count_list) != 0:
            counts = ",".join([str(i) for i in value])
            ch_mvmt_file.write(f'{ch_pair[0]},{ch_pair[1]},{counts}\n')      
    ch_mvmt_file.close()
    
    if config["output_tree"] or config["calculate_R0"] or config["output_ltt"] or config["output_skyline"]:
        result = tree_sim.simulate_tree(epidemic_config, config, last_day) 
        if result:
            coalescent_tree, newick_string, skyline, R0, those_sampled, ltt = result

            config["files"]["most_recent_tip_file"].write(f'{iteration_count},{tree.most_recent_date}\n')

            if config["output_tree"]:
                tree_file = open(os.path.join(config["output_directory"],"trees",f"tree_for_{iteration_count}.txt", 'w'))
                tree_file.write(newick_string)
                tree_file.close()

            if config["make_skyline"]: 
                skyline_file = open(os.path.join(output_directory,"skylines", f"skyline_for_{iteration_count}.csv", 'w'))
                skyline_file.write("number,start_interval,end_interval,logpopsize\n")

                logpop_count = 0
                for key, value in skyline.items():
                    logpop_count += 1
                    skyline_file.write(f"{logpop_count},{key[0]},{key[1]},{value}\n")
                skyline_file.close()
            
            if config["make_ltt"]: 
                ltt_file = open(os.path.join(output_directory,"lineages",f"ltt_for_{iteration_count}.csv", 'w')) 
                ltt_file.write("number,start,end,lineages\n")

                lineage_count = 0
                for k,v in lineages_through_time.items():
                    lineage_count += 1
                    ltt_file.write(f"{lineage_count},{k[0]},{k[1]},{v}\n")
                ltt_file.close()
            
            if config["calculate_R0"]:
                config["files"]["R0_output"].write(f"{iteration_count},{R0}\n")






