from individual_class import *
from case_class import *
import random

def make_index_case(agent_location, cfr, distributions, original_case_dict, original_trans_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict):
    
    index_case_individual = Individual(random.choice(range(1382431,1908832)), agent_location, cfr, distributions) #These should be the IDs of the range in Kailahun

    index_case_individual.incubation_day = 0 #So that the first case is infectious on day one of the simulation

    index_case_case = Case(0, None, None)
    original_case_dict[index_case_case] = index_case_individual

    original_trans_dict[index_case_individual.unique_id] = ["NA", '0', index_case_individual.incubation_day]
    original_nodes.append(index_case_individual.unique_id)

    infected_individuals_set.add(index_case_individual.unique_id)

    original_districts_present.append(index_case_individual.dist)

    original_cluster_set.add(index_case_individual.comm)

    #poss_case_dict = get_possible_cases(index_case_individual)
    index_case_dict = {}
    index_case_dict["Hh"] = 5
    index_case_dict["Comm"] = 7
    index_case_dict["Dist"] = 2
    index_case_dict["Country"] = 0

    #print(poss_case_dict)

    for level, number in index_case_dict.items():
        for person in range(number):
            day_inf = 0
            new_case = Case(len(original_case_dict), level, index_case_case)
            original_case_dict[new_case] = None
            original_day_dict[day_inf].append(new_case)
            
           
            
    return index_case_case, index_case_individual, original_case_dict, original_trans_dict, original_nodes, infected_individuals_set, original_districts_present, original_cluster_set, original_day_dict