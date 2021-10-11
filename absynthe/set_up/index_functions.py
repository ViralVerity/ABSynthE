from absynthe.classes.individual_class import *
from absynthe.classes.case_class import *
import random
from collections import defaultdict
#stick this in make_contact_dicts
def make_data_structures(config):

    epidemic_config = {}

    epidemic_config["case_dict"] = {} #case object to individual object
    epidemic_config["transmission_dict"] = defaultdict(dict)
    epidemic_config["child_dict"] = {}
    epidemic_config["day_dict"] = {}

    epidemic_config["option_dict_district_level"] = defaultdict(list)
    epidemic_config["infected_individuals_set"] = set()
    epidemic_config["cdf_array"] = []
    epidemic_config["cdf_len_set"] = set()

    epidemic_config["districts_present"] = set() #was a list before, I don't think it needs to be
    epidemic_config["chiefdoms_present"] = set() 
    epidemic_config["nodes"] = []
    epidemic_config["onset_times"] = []

    epidemic_config["dist_mvmt"] = defaultdict(list)
    epidemic_config["ch_mvmt"] = defaultdict(list)


    if config["day_limit"]:
        for i in range(config["day_limit"]): 
            epidemic_config["day_dict"][i] = []
    else:
        for i in range(1000):
            epidemic_config["day_dict"][i] = []


    for item1 in config["population_structure"]["district_list"]:
        for item2 in config["population_structure"]["district_list"]:
            if item1 != item2:
                epidemic_config["dist_mvmt"][item1,item2] = []
                
    for item1 in config["population_structure"]["ch_list"]:
        for item2 in config["population_structure"]["ch_list"]:
            if item1 != item2:
                epidemic_config["ch_mvmt"][item1, item2] = []

    return epidemic_config

#put this in the individual class definition - add a index=False default
def make_index_case(config, epidemic_config):

    index_id = random.choice(range(1169263,1214412)) #These are the IDs of the range in Kissi Teng, Kailahun
    epidemic_config["index_id"] = index_id

    index_case_individual = Individual(index_id, None, config["population_structure"]["agent_location"], config["cfr"], config["distributions"], 0, epidemic_config) 
    index_case_individual.incubation_day = 0 #So that the first case is infectious on day one of the simulation
    index_case_case = Case(0, None, None)
    
    epidemic_config["case_dict"][index_case_case] = index_case_individual

    #From Wauqier 2015 - specific to EBOV SLE, it's the number of secondary cases from the index case at the funeral
    index_case_dict = {}
    index_case_dict["Hh"] = 2
    index_case_dict["Ch"] = 7
    index_case_dict["Dist"] = 3
    index_case_dict["Country"] = 2

    for level, number in index_case_dict.items():
        for person in range(number):
            day_inf = 0
            new_case = Case(len(epidemic_config["case_dict"]), level, index_case_case)
            epidemic_config["case_dict"][new_case] = None
            epidemic_config["day_dict"][day_inf].append(new_case)

    epidemic_config["index_case_individual"] = index_case_individual
    epidemic_config["index_case_case"] = index_case_case
               
    return epidemic_config


