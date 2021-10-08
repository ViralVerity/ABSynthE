import json
from collections import defaultdict
import os

#Outputs are the same except that the cluster is now a chiefdom

def make_contact_dicts(input_directory, config):
    
    agent_location = defaultdict(list)
    district_pops = {}

    #will need to be generalised if it's going to be useable for other levels of people - maybe adm1,adm2?
    with open(os.path.join(input_directory,"agent_location.txt")) as json_file:
        data = json.load(json_file)

    for key in data.keys():
        new_key = int(key)
        agent_location[new_key] = data[key]
    config["agent_location"] = agent_location

    with open(os.path.join(input_directory,"District_to_hh.txt")) as json_file:
        config["dist_to_hh"] = json.load(json_file)

    with open((os.path.join(input_directory,"household_to_chiefdom.txt")) as json_file:
        config["hh_to_ch"] = json.load(json_file)

    with open((os.path.join(input_directory,"Chiefdom_to_hh.txt")) as json_file:
        config["ch_to_hh"] = json.load(json_file)

    with open((os.path.join(input_directory,"Hh_to_people.txt")) as json_file:
        config["hh_to_ppl"] = json.load(json_file)

    with open((os.path.join(input_directory,"chiefdom_to_people.txt")) as json_file:
        config["chiefdom_to_ppl"] = json.load(json_file)

    with open((os.path.join(input_directory,"District_to_ppl.txt")) as json_file:
        config["dist_to_ppl"] = json.load(json_file)

    with open((os.path.join(input_directory,"district_relative_distance.txt")) as json_file:
        config["district_distance"] = json.load(json_file)
        
    with open((os.path.join(input_directory,"dist_to_ch.txt")) as json_file:
        config["dist_to_ch"] = json.load(json_file)

  
    #Do we actually use this one? Drawn from census directly anyway
    with open((os.path.join(input_directory,"district_population.csv", 'r')) as f:
        for l in f:
            toks = l.strip("\n").split(",")
            district_pops[toks[0]] = int(toks[1])
    config["district_pops"] = district_pops
            
            
    return config