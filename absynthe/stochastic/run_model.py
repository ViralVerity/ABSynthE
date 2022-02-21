from collections import defaultdict
import sys
import os
import absynthe.set_up.index_functions as index_functions
import absynthe.set_up.file_functions as file_functions
import absynthe.stochastic.tree_simulator as tree_sim
import time

from absynthe.classes.case_class import Case
from absynthe.classes.individual_class import Individual

def run_model(config, iteration_count):
        
    ##writing to file and screen
    if iteration_count%config["log_every"] == 0:
        config["write_file"] = True
        print(f'{iteration_count} runs completed')
    else:
        config["write_file"] = False            

    epidemic_config = index_functions.make_data_structures(config)        
    
    ###Making index case###
    epidemic_config = index_functions.make_index_case(config, epidemic_config)
    
    if config["write_file"]:
        epidemic_config["info_file"] = file_functions.prep_info_file(config["output_directory"], epidemic_config["index_case_individual"], iteration_count)
    
    ###Run the epidemic###
    epidemic_config = run_epidemic(0, config, epidemic_config)

    ###Removing cases that don't exist eg because the person was already infected, or because the parent had recovered/died###
    #NB if a case limit is set, it's possible this removal process will put the case count below the case limit.
    new_case_dict = {k:v for k,v in epidemic_config["case_dict"].items() if v}
    epidemic_config['case_dict'] = new_case_dict

    #Removing those people from the day dict
    for day, case_list in epidemic_config["day_dict"].items():
        new_case_list = [item for item in case_list if item in epidemic_config['case_dict']] 
        epidemic_config["day_dict"][day] = new_case_list

    epidemic_config["day_dict"][0].append(epidemic_config["index_case_case"]) #Put here so that it doesn't confuse the loop above because it has no parent AND otherwise it would get reassigned and stuff
    last_day = max(epidemic_config["onset_times"])
    
    if config["day_limit"]:
        if config["day_limit"] < last_day:
            last_day = config["day_limit"]

    ###Getting results and writing to file###
    if epidemic_config["epidemic_stopped"] and not config["write_file"]: #special log file for when it runs out 
        file_functions.write_runout_file(config, epidemic_config, iteration_count)

    if epidemic_config["epidemic_stopped"] or config["write_file"]:
        R0,most_recent_date = record_individual_epidemic(iteration_count, config, epidemic_config, last_day)        

    if config["write_file"]: #is there a way to check if a file is open?
        epidemic_config["info_file"].close()

    result_dict = {}
    result_dict["iteration_count"] = iteration_count
    result_dict["length"] = last_day
    result_dict["cases"] = len(epidemic_config["case_dict"])
    result_dict["districts"] = len(epidemic_config["districts_present"])
    result_dict['chiefdoms'] = len(epidemic_config["chiefdoms_present"])
    result_dict['epidemic_stopped'] = epidemic_config['epidemic_stopped']
    result_dict["R0"] = R0
    result_dict['most_recent_date'] = most_recent_date
        
    return result_dict
            
def run_epidemic(start_day, config, epidemic_config):
    
    epidemic_config["epidemic_stopped"] = False
    start = time.time()
    # try: #replace with an if statement
    for day, case_list in epidemic_config["day_dict"].items():             
        if config["verbose"]:
            if day >= 100:
                if day%100 == 0:
                    diff = (time.time() - start)/60
                    print(f'{int(day/diff)} days per minute on day {day} and {int(len(epidemic_config["case_dict"])/diff)} cases per minute.')

        if config["capped"]:
            if config["day_limit"]:
                if day > config["day_limit"]: 
                    sys.stdout.write("Epidemic has reached day limit\n")
                    epidemic_config["epidemic_stopped"]  = True
                    return epidemic_config
            elif config["case_limit"]:
                if len(epidemic_config["case_dict"]) > config["case_limit"]:
                    sys.stdout.write("Epidemic has reached case limit\n")
                    epidemic_config["epidemic_stopped"]  = True
                    return epidemic_config
                
        if len(case_list) != 0 and day >= start_day: #If there are new cases on this day
            for focal_case in case_list:
                
                parent = epidemic_config["case_dict"][focal_case.parent]#Gets the individual object of parent (intialised last time) from the case dictionary using the case object
                result = focal_case.who_am_I(parent, day, config, epidemic_config) #Assign the current case id to an individual id

                if not result:
                    if result == 0: #If individual is already infected
                        pass
                    else: #If there is no-one left
                        print(result)
                        print(len(epidemic_config["infected_individuals_set"]))
                        print(config["population_structure"]["popn_size"])
                        sys.stdout.write("All members of the population infected\n")
                        epidemic_config["epidemic_stopped"] = True
                        return epidemic_config

                #Test that this only comes up when type = Individual
                else: #There is a successful assignation to a specific individual
                    assignment, config, epidemic_config = result
                    
                    focal_individual = Individual(assignment, parent, config["population_structure"]["agent_location"], config["cfr"], config["distributions"], day, epidemic_config)
                    epidemic_config["case_dict"][focal_case] = focal_individual
                    parent.children.append(focal_individual)

                    if config["write_file"]:
                        if focal_individual.death_state:
                            death = day+focal_individual.incubation_day+focal_individual.death_day
                            recovery = "NA"
                        else:
                            recovery = day+focal_individual.incubation_day+focal_individual.recovery_day
                            death = "NA"
                            
                        epidemic_config["info_file"].write(f'{focal_individual.unique_id},{focal_individual.parent.unique_id},{focal_individual.hh},{focal_individual.ch},{focal_individual.dist},{day},{day+focal_individual.incubation_day},{focal_individual.death_state},{death},{recovery},{day + focal_individual.incubation_day}\n')

                        if focal_individual.dist != focal_individual.parent.dist:
                            epidemic_config["dist_mvmt"][focal_individual.dist,focal_individual.parent.dist].append(day)
                        if focal_individual.ch != focal_individual.parent.ch: #This is going to include between district movements too remember
                            epidemic_config["ch_mvmt"][focal_individual.ch,focal_individual.parent.ch].append(day)
                    
                    ##get possible transmissions##
                    poss_case_dict = focal_individual.get_possible_cases(config) #Gives dict of contact_level: number of people
                    
                    for level, number in poss_case_dict.items():
                        for person in range(number):
                            day_inf_output = focal_individual.when_infected(day, person, epidemic_config["cdf_len_set"], epidemic_config["cdf_array"])[0]
                            
                            if not day_inf_output:
                                #print("Finished infection first")
                                continue

                            else:
                                if config["day_limit"]:
                                    if day_inf_output > config["day_limit"]:
                                        continue 

                                new_case = Case(len(epidemic_config["case_dict"]), level, focal_case)
                                epidemic_config["case_dict"][new_case] = None
                                epidemic_config["day_dict"][day_inf_output].append(new_case)


    # except RuntimeError: #couldn't work it out for internal, so for now just adding more days at the start
    #     if not config["capped"]:
    #         original_length = len(epidemic_config["day_dict"])
    #         new_start = day #Should start again from when the error was thrown, so for now will recalculate all the infecteds for that day
    #         for i in range(1000):
    #             epidemic_config["day_dict"][original_length + i] = []

    #         run_epidemic(new_start, epidemic_config)
    #     else:
    #         pass
        

    return epidemic_config



def record_individual_epidemic(iteration_count, config, epidemic_config, last_day):

    district_mvmt_file, ch_mvmt_file = file_functions.prep_movement_files(config["output_directory"], iteration_count)
            
    for district_pair, count_list in epidemic_config["dist_mvmt"].items(): #what is the value here - is it counts or days that they're happening on?
        if len(count_list) != 0:
            counts = ",".join([str(i) for i in count_list])
            district_mvmt_file.write(f'{district_pair[0]},{district_pair[1]},{counts}"\n"')
    district_mvmt_file.close()
    
    for ch_pair, count_list in epidemic_config["ch_mvmt"].items():
        if len(count_list) != 0:
            counts = ",".join([str(i) for i in count_list])
            ch_mvmt_file.write(f'{ch_pair[0]},{ch_pair[1]},{counts}\n')      
    ch_mvmt_file.close()
    
    R0 = None
    most_recent_date = None
    if config["output_tree"] or config["calculate_R0"] or config["output_ltt"] or config["output_skyline"]:
        result = tree_sim.simulate_tree(epidemic_config, config, last_day) 
        if result:
            coalescent_tree, newick_string, skyline, R0, those_sampled, lineages_through_time, coalescent_times = result
            most_recent_date = coalescent_tree.most_recent_date 
            
            if config["output_tree"]:
                tree_file = open(os.path.join(config["output_directory"],"trees",f"tree_for_{iteration_count}.txt"), 'w')
                tree_file.write(newick_string)
                tree_file.close()

            if config["output_skyline"]: 
                skyline_file = open(os.path.join(config["output_directory"],"skylines", f"skyline_for_{iteration_count}.csv"), 'w')
                skyline_file.write("number,start_interval,end_interval,logpopsize\n")

                logpop_count = 0
                for key, value in skyline.items():
                    logpop_count += 1
                    skyline_file.write(f"{logpop_count},{key[0]},{key[1]},{value}\n")
                skyline_file.close()
            
            if config["output_ltt"]: 
                ltt_file = open(os.path.join(config["output_directory"],"lineages",f"ltt_for_{iteration_count}.csv"), 'w') 
                ltt_file.write("number,start,end,lineages\n")

                lineage_count = 0
                for k,v in lineages_through_time.items():
                    lineage_count += 1
                    ltt_file.write(f"{lineage_count},{k[0]},{k[1]},{v}\n")
                ltt_file.close()
                
    return R0, most_recent_date
                






