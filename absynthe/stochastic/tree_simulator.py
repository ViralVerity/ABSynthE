from collections import defaultdict
import random
import sys

from absynthe.classes.node_class import *
from absynthe.classes.tree_class import *


def sampling(transmission_dict, config, epidemic_len): 
    """Get who is sampled, inputs are list of individual ids and the sampled percentage"""
    not_enough_cases = False
    week_bins = []
    weeks_cases = defaultdict(list)
    those_sampled = set()
    
    total = len(transmission_dict) #Doesn't have NA in it because that can't be sampled
    number_of_people = round(total*config["sampling_percentage"])
    weeks = round(epidemic_len/7)
    
    if config["sampling_scheme"] == "uniform":
        if weeks > 0:
            samples_per_week = round(number_of_people/weeks)
        else:
            return those_sampled, not_enough_cases
        
        if samples_per_week < 1:
            not_enough_cases = True
            return those_sampled, not_enough_cases
        
        for i in range(weeks):
            weeks_cases[i] = []
    
        for person,person_dict in transmission_dict.items():
            week_number = int(person_dict["day_sampled"]/7)
            weeks_cases[week_number].append(person)
                    
            
        #Leaving this out because it's a specific conditioning on long epidemics - ie it has to have cases up until the very end. This is for fitting
        #week_cases[len(week_cases) - 1] = (140,148)
        #if len(weeks_cases[(140,148)]) == 0:
        #   not_enough_cases = True
        #  return
                    
        for week, lst in weeks_cases.items():
            if len(lst) >= samples_per_week:
                for i in random.sample(lst, k=samples_per_week):
                    those_sampled.add(i)
                
        more_needed = number_of_people - len(those_sampled)
        
        if more_needed > 0:
            try:
                for i in random.sample(transmission_dict.keys(), k=more_needed):
                    those_sampled.add(i)
            except ValueError:
                not_enough_cases = True 
                return those_sampled, not_enough_cases
    
    return those_sampled, not_enough_cases

def get_R0(node_dict): 
    """Count number of people in 4th and 3rd generations to calculate R0"""
    
    gen_3 = 0
    gen_4 = 0
    
    for nde in node_dict.values():
        if nde.generation == 3:
            gen_3 += 1
            gen_4 += len(nde.infections)
            
    if gen_3 != 0 and gen_4 != 0:
        R0 = gen_4/gen_3
    
        return R0

def plot_skyline(Ne_dict): #not really relevant in here

    for_plotting = {}
    for key, value in Ne_dict.items():
        
        new_key = key*-1
        
        for_plotting[new_key] = value

    names = list(for_plotting.keys())
    values = list(for_plotting.values())

    plt.step(names,values)


def simulate_tree(epidemic_config, config, epidemic_len): 
    """Function to simulate a coalescent tree"""
            
    R0 = None
    newick_tree = None
    ltt = None
    skyline = None

    transmission_dict = epidemic_config["transmission_dict"]
    if config["day_limit"]:
        if epidemic_len < config["day_limit"]:
            pass
        else:
            transmission_dict = defaultdict(dict)
            for person, person_dict in epidemic_config["transmission_dict"].items():
                if person_dict["day_sampled"] <= config["day_limit"]:
                    transmission_dict[person] = person_dict

    sampling_output = sampling(transmission_dict, config, epidemic_len)
    
    if not sampling_output:
        sys.stdout.write("Not enough cases to make a tree")
        return
    else:
        those_sampled, not_enough_cases = sampling_output
    
    if len(those_sampled) != 0:
   
        #Intialise as type individual nodes, make transmission tree and make subtrees
        
        index_case = epidemic_config["child_dict"]["NA"][0]
        node_dict = {}
        
        #this is recursive, so calling it once generates the whole node dictionary
        node(index_case, "individual", transmission_dict=epidemic_config["transmission_dict"], child_dict=epidemic_config["child_dict"], those_sampled=those_sampled, node_dict=node_dict)
        coalescent_tree = tree(tree_type="whole_tree",node_dict=node_dict, epidemic_len=epidemic_len)
        
        if config["calculate_R0"]:
            R0 = get_R0(node_dict)
        if config["output_tree"]:
            newick_string = coalescent_tree.to_newick(coalescent_tree.root, those_sampled)
        if config["output_ltt"] or config["output_skyline"]:
            skyline, lineages_through_time, coalescent_times = coalescent_tree.calculate_ne(those_sampled)
        
        return coalescent_tree, newick_string, skyline, R0, those_sampled, ltt
        
    else:
        # sys.stderr.write("No-one assigned for sampling in tree simulation\n")
        return

