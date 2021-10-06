from absynthe.classes.individual_class import *
from absynthe.classes.case_class import *
import random


def make_data_structures(config):

    case_dict = {}
    trans_dict = defaultdict(list)
    child_dict = defaultdict(list)
    day_dict = defaultdict(list)

    option_dict_districtlevel = defaultdict(list) #needs a better name
    infected_individuals_set = set()
    cdf_array = []
    cdf_len_set = set()
    
    districts_present = []
    chiefdom_set = set()
    
    nodes = []
    onset_times = []

    if config["day_limit"]:
        for i in range(config["day_limit"]): 
            day_dict[i] = []
    else:
        for i in range(1000):
            day_dict[i] = []

    #setting up empty dictionaries
    dist_mvmt = defaultdict(list)
    ch_mvmt = defaultdict(list)

    for item1 in district_list:
        for item2 in district_list:
            if item1 != item2:
                dist_mvmt[item1,item2] = []
                
    for item1 in ch_list:
        for item2 in ch_list:
            if item1 != item2:
                ch_mvmt[item1, item2] = []

    return case_dict, day_dict, trans_dict, child_dict, option_dict_districtlevel, infected_individuals_set, cdf_array, cdf_len_set, districts_present, chiefdom_set, nodes, onset_times, dist_mvmt, ch_mvmt

def make_index_case(config, data_structures):

    case_dict, day_dict, trans_dict, child_dict, option_dict_districtlevel, infected_individuals_set, cdf_array, cdf_len_set, districts_present, chiefdom_set, nodes, onset_times, dist_mvmt, ch_mvmt = data_structures

    index_case_individual = Individual(random.choice(range(1169263,1214412)), config["agent_location"], config["cfr"], config["distributions"]) #These should be the IDs of the range in Kissi Teng, Kailahun

    index_case_individual.incubation_day = 0 #So that the first case is infectious on day one of the simulation

    index_case_case = Case(0, None, None)
    case_dict[index_case_case] = index_case_individual

    trans_dict[index_case_individual.unique_id] = ["NA", '0', index_case_individual.incubation_day]
    nodes.append(index_case_individual.unique_id)
    
    child_dict["NA"] = [index_case_individual.unique_id]

    infected_individuals_set.add(index_case_individual.unique_id)

    districts_present.append(index_case_individual.dist)

    cluster_set.add(index_case_individual.comm)

    #From Wauqier 2015 - specific to EBOV SLE
    index_case_dict = {}
    index_case_dict["Hh"] = 2
    index_case_dict["Comm"] = 7
    index_case_dict["Dist"] = 3
    index_case_dict["Country"] = 2

    for level, number in index_case_dict.items():
        for person in range(number):
            day_inf = 0
            new_case = Case(len(case_dict), level, index_case_case)
            case_dict[new_case] = None
            day_dict[day_inf].append(new_case)
               
    return index_case_case, index_case_individual, case_dict, trans_dict, child_dict, nodes, infected_individuals_set, districts_present, cluster_set, day_dict, dist_mvmt, ch_mvmt


