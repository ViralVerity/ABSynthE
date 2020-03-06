from case_class import *
from individual_class import *

def run_epidemic(start_day, day_dict, susceptibles_left , case_dict, trans_dict, child_dict, infected_individuals_set, popn_size, option_dict_districtlevel, onset_times, nodes, cluster_set, cdf_len_set, cdf_array, districts_present, dist_mvmt, ch_mvmt, contact_structure, cfr, distributions, write_file, info_file, iteration_count, capped, epidemic_length, case_limit, SDB_start, SDB_success):
    
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
                    assignment = focal_case.who_am_I(infected_individuals_set, popn_size, option_dict_districtlevel, contact_structure, case_dict, parent, day, cfr, distributions, day, SDB_start, SDB_success) #Assign the current case id to an individual

                    if assignment == None and day != 0: #If individual is already infected
                        #remove_set.add(focal_case) #Case doesn't exist so must be removed from day/case dict
                        pass

                    elif assignment == False: #If there is no-one left
                        print("All members of the population infected")
                        susceptibles_left = False
                        return day_dict, case_dict, nodes, trans_dict, dist_mvmt, onset_times, dist_present, cluster_set
                    
                    #Test that this only comes up when type = Individual
                    else: #There is a successful assignation to a specific individual

                        onset_times.append(day)

                        focal_individual = case_dict[focal_case] #This is an Individual object got by who_am_I

                        focal_individual.parent = parent #Gives Individual the same parent as the Case object for recording

                        parent.children.append(focal_individual)

                        trans_dict[focal_individual.unique_id] = [focal_individual.parent.unique_id, day, (day+ focal_individual.incubation_day)]
                        
                        child_dict[focal_individual.unique_id] = []
                        
                        child_dict[focal_individual.parent.unique_id].append(focal_individual.unique_id)

                        nodes.append(focal_individual.unique_id)

                        if write_file == True:
                            info_file.write(f'{focal_individual.unique_id},{focal_individual.parent.unique_id}, {focal_individual.hh},{focal_individual.comm},{focal_individual.dist},{day},{day+focal_individual.incubation_day}, {day + focal_individual.incubation_day}\n')


                            if focal_individual.dist != focal_individual.parent.dist:
                                dist_mvmt[focal_individual.dist,focal_individual.parent.dist].append(day)
                            if focal_individual.comm != focal_individual.parent.comm: #This is going to include between district movements too remember
                                ch_mvmt[focal_individual.comm,focal_individual.parent.comm].append(day)


                        if focal_individual.dist not in districts_present:
                            districts_present.append(focal_individual.dist)

                        if focal_individual.comm not in cluster_set:
                            cluster_set.add(focal_individual.comm)

                        poss_case_dict = focal_individual.get_possible_cases() #Gives dict of contact_level: number of people
                        
                        for level, number in poss_case_dict.items():
                            for person in range(number):
                                #when the potential cases are infected
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
            
            original_length = len(day_dict)
            new_start = day #Should start again from when the error was thrown, so for now will recalculate all the infecteds for that day
            for i in range(1000):
                day_dict[original_length + i] = []

            run_epidemic(new_start, day_dict, susceptibles_left, case_dict, trans_dict, child_dict, infected_individuals_set, popn_size, option_dict_districtlevel, onset_times, nodes, cluster_set, cdf_len_set, cdf_array, districts_present, dist_mvmt, ch_mvmt, contact_structure, cfr, distributions, write_file, info_file, iteration_count, capped, epidemic_length, case_limit, SDB_start, SDB_success)
        else:
            pass
        

    return day_dict, case_dict, nodes, trans_dict, child_dict, dist_mvmt, ch_mvmt, onset_times, districts_present, cluster_set, epidemic_capped