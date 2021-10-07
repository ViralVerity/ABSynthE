from absynthe.classes.case_class import *
from absynthe.classes.individual_class import *
import sys

#needs tidying and integrating config for the args
def run_epidemic(start_day, config, epidemic_config):
    
    epidemic_config["epidemic_stopped"] = False
    day_count = 0
    
    try: #replace with an if statement
        for day, case_list in epidemic_config["day_dict"].items():    
            day_count += 1
            if config["capped"]:
                if day_count > config["day_limit"] or len(epidemic_config["case_dict"]) > config["case_limit"]: #but then will cases get assigned an individual? Might need to add in another function to do that
                    sys.stdout.write("Epidemic has reached day limit or case limit")
                    epidemic_config["epidemic_stopped"]  = True
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

            run_epidemic(new_start, day_dict, susceptibles_left, case_dict, trans_dict, child_dict, infected_individuals_set, popn_size, option_dict_districtlevel, onset_times, nodes, cluster_set, cdf_len_set, cdf_array, districts_present, dist_mvmt, ch_mvmt, contact_structure, cfr, distributions, write_file, info_file, iteration_count, capped, epidemic_length, case_limit)
        else:
            pass
        

    return day_dict, case_dict, nodes, trans_dict, child_dict, dist_mvmt, ch_mvmt, onset_times, districts_present, cluster_set, epidemic_stopped