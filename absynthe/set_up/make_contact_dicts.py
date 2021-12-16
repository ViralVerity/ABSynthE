import json
from collections import defaultdict
import os

##add the other function from index funcs in here

def make_contact_dicts(input_directory, config):
    
    agent_location = defaultdict(list)
    district_pops = {}

    #will need to be generalised if it's going to be useable for other levels of people - maybe adm1,adm2?
    with open(os.path.join(input_directory,"agent_location.txt")) as json_file:
        data = json.load(json_file)

    for key in data.keys():
        new_key = int(key)
        agent_location[new_key] = data[key]
    config["population_structure"]["agent_location"] = agent_location

    with open(os.path.join(input_directory,"household_to_people.txt")) as json_file:
        config["population_structure"]["hh_to_ppl"] = json.load(json_file)

    with open(os.path.join(input_directory,"chiefdom_to_people.txt")) as json_file:
        config["population_structure"]["ch_to_ppl"] = json.load(json_file)

    with open(os.path.join(input_directory,"district_to_people.txt")) as json_file:
        config["population_structure"]["dist_to_ppl"] = json.load(json_file)

    # with open(os.path.join(input_directory,"district_neighbouring.txt")) as json_file:
    #     config["population_structure"]["neighbouring_districts"] = json.load(json_file)
            
            
    return config