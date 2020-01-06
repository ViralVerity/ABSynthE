import json
from collections import defaultdict

def make_contact_dicts(dropbox_path):
    
    agent_location = defaultdict(list)
    district_pops = {}

    with open(dropbox_path + "Contact_structure/agent_location.txt") as json_file:
        data = json.load(json_file)

    for key in data.keys():
        new_key = int(key)
        agent_location[new_key] = data[key]


    with open(dropbox_path + "Contact_structure/District_to_hh.txt") as json_file:
        dist_to_hh = json.load(json_file)

    with open(dropbox_path + "Contact_structure/hh_to_cluster.txt") as json_file:
        hh_to_cluster = json.load(json_file)

    with open(dropbox_path + "Contact_structure/cluster_to_hh.txt") as json_file:
        cluster_to_hh = json.load(json_file)

    with open(dropbox_path + "Contact_structure/Hh_to_people.txt") as json_file:
        hh_to_ppl = json.load(json_file)

    with open(dropbox_path + "Contact_structure/cluster_to_ppl.txt") as json_file:
        cluster_to_ppl = json.load(json_file)

    with open(dropbox_path + "Contact_structure/District_to_ppl.txt") as json_file:
        dist_to_ppl = json.load(json_file)

    with open(dropbox_path + "Contact_structure/district_relative_distance.txt") as json_file:
        district_distance = json.load(json_file)

  

    with open(dropbox_path + "Contact_structure/district_population.csv", 'r') as f:
        for l in f:
            toks = l.strip("\n").split(",")
            district_pops[toks[0]] = int(toks[1])
            
            
    return agent_location, dist_to_hh, hh_to_cluster, cluster_to_hh, hh_to_ppl, cluster_to_ppl, dist_to_ppl, district_distance, district_pops