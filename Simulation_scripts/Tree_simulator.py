from collections import defaultdict
from collections import OrderedDict
import random
import uuid
import numpy as np
from scipy import special

from node_class import *
from tree_class import *


###Probably can't go into a class yet - need to tidy up anyway though###
def sampling(trans_dict, sampled_percentage, epidemic_len): 
    """Get who is sampled, inputs are list of individual ids and the sampled percentage"""
    not_enough_cases = False
    week_bins = []
    weeks_cases = defaultdict(list)
    those_sampled = set()
    
    total = len(trans_dict) #Doesn't have NA in it because that can't be sampled
    
    number_of_people = round(total*sampled_percentage)

    weeks = epidemic_len/7
    
    samples_per_week = round(number_of_people/weeks)
    
    if samples_per_week < 1:
        not_enough_cases = True
        return not_enough_cases
    
    count = 0
    for i in range(round(weeks)):
        tup = (count, count + 6)
        
        week_bins.append(tup)
        
        count += 7
    
    week_bins[len(week_bins) - 1] = (140,148)
        
    for i in week_bins:
        weeks_cases[i] = []
    
    for key,value in trans_dict.items():
        date_sampled = value[2]
        for tup in week_bins:
            if date_sampled < tup[1] and date_sampled >= tup[0]:
                weeks_cases[tup].append(key)
                
    if len(weeks_cases[(140,148)]) == 0:
        not_enough_cases = True
        return not_enough_cases
                
    for week, lst in weeks_cases.items():
        if len(lst) >= samples_per_week:
            for i in random.sample(lst, k=samples_per_week):
                those_sampled.add(i)
            
    #print(those_sampled)
    
    
    more_needed = number_of_people - len(those_sampled)
    
    if more_needed > 0:
        try:
            for i in random.sample(trans_dict.keys(), k=more_needed):
                those_sampled.add(i)
        except ValueError:
            not_enough_cases = True
            return not_enough_cases
    
                         
    #those_sampled = set(random.sample(trans_dict.keys(), k = number_of_people))
    
    #print("number of sequences is " + str(len(those_sampled)))
    
    #TEST HERE that the number of sampled people is what it should be AND that they have all been marked as sampled
    
    return those_sampled, not_enough_cases



#With tree class
def update_coalescent_tree(subtree, whole_tree):
    """Update the coalescent tree with information from the subtree"""
    
    whole_tree.branch_lengths.update(subtree.branch_lengths)
    
    if subtree.sample_tip:
        whole_tree.nodes.add(subtree.sample_tip)

    for tip in subtree.transmission_tips:
        whole_tree.nodes.add(tip)

    for internal_node in subtree.coalescent_nodes:
        whole_tree.nodes.add(internal_node)
    
    
#With tree class
def connections(subtree, whole_tree, subtree_dict):
    """Stick the subtrees together by attaching transmission nodes to roots of next tree"""
    
    #Find the most recent tip of the whole coalescent tree
    if subtree.sample_tip: 
        if float(subtree.sample_tip.time_sampled) >= whole_tree.most_recent_date:
            whole_tree.most_recent_tip = subtree.sample_tip
            whole_tree.most_recent_date = float(subtree.sample_tip.time_sampled)
      
    
    #Deals with tips of the tree
    if len(subtree.transmission_tips) == 0:
        whole_tree.nodes.add(subtree.sample_tip) 
        whole_tree.branch_lengths.update(subtree.branch_lengths)
        pass
    
    else:

        for i in subtree.transmission_tips:
            connections(subtree_dict[i.id.id], whole_tree, subtree_dict)
    
    #Deals with the root
    if subtree.id.parent.id == "NA":
        
        whole_tree.root = subtree.root
        
        update_coalescent_tree(subtree, whole_tree)
                    
        for i in subtree.transmission_tips:
            connections(subtree_dict[i.id.id], whole_tree, subtree_dict)
        
        return subtree, whole_tree
        
    #Deals with all the nodes in between    
    else:

        subtree.root.tree_parent = subtree_dict[subtree.id.parent.id] #the subtree object

        subtree_parent = subtree.root.tree_parent

        right_node = [i for i in subtree_parent.transmission_tips if i.id == subtree.id][0] #the transmission node for this tree

        right_node.node_children.add(subtree.root)

        subtree.root.node_parent = right_node 
        
        update_coalescent_tree(subtree, whole_tree)
        
        return subtree, whole_tree


#Tree class
def remove_internals(whole_tree_node, whole_tree):
    """Remove internal nodes that joined subtrees together"""
    whole_tree_node.remove_func_called = True
    
    whole_tree_node.new_children = whole_tree_node.node_children.copy()
    
    if len(whole_tree_node.node_children) == 1:
        if whole_tree_node != whole_tree.root:
            
            #reassign next node's parent
            for child in whole_tree_node.new_children:
                child.node_parent = whole_tree_node.node_parent

                #Add next node to parent's child list
                whole_tree_node.node_parent.new_children.add(child)

                whole_tree.branch_lengths[child] = whole_tree.branch_lengths[child] + whole_tree.branch_lengths[whole_tree_node]

            whole_tree_node.removed = True

            #Remove self from parent's child list
            whole_tree_node.node_parent.new_children.remove(whole_tree_node)
        
        else:
            
            for child in whole_tree_node.new_children:
                child.node_parent = None
                whole_tree.branch_lengths[child] = 0.0
                whole_tree.root = child
                
            whole_tree_node.removed = True
            
    for i in whole_tree_node.node_children:
        whole_tree_node.for_loop_called = True
        remove_internals(i, whole_tree)
        
    return whole_tree

#Tree class
def get_tip_to_root(nde, whole_tree):
    """Traverses sampled coalescent tree to get node heights"""
  
    if nde.root_to_tip != 0.0:
        return nde.root_to_tip
    
    if nde == whole_tree.root:
        
        distance = 0.0
        
    else:
        
        distance = get_tip_to_root(nde.node_parent, whole_tree) + whole_tree.branch_lengths[nde]
        
    nde.root_to_tip = distance
    
    return nde.root_to_tip


#Tree class
def to_newick(nde, whole_tree, those_sampled):
    """Makes a newick string of the tree"""
    
    string = ",".join([to_newick(i, whole_tree, those_sampled) for i in nde.new_children if not i.removed])

    if len(nde.new_children) != 0:
        string = "(" + string

    if type(nde) == node:
        
        string += str(nde.id)
    else:

        string += ")"

    string += ":" + str(round(whole_tree.branch_lengths[nde], 2))
    
    if nde == whole_tree.root:
        string += ")"
        
    return string

#Tree class
def get_active_population(whole_tree):
    """Get active population at each coalescent interval"""
    
    coalescent_times = set()
    coal_list = []
    coalescent_intervals = defaultdict(tuple)
    
    active_population = defaultdict(list)

    sorted_dict = OrderedDict(sorted(whole_tree.node_heights.items(), key=lambda x:x[1]))

    for nde, height in sorted_dict.items():
       
        if type(nde) == coalescent_node:
            
            coalescent_times.add(height)
            coal_list.append(height)

    coalescent_times = sorted(coalescent_times)
    
    current_time = 0
    
    count = 0
    
    non_parent_set = set()
    
    for time in coalescent_times:
        
        count += 1
        
        coalescent_intervals[count] = (float(current_time),float(time))
        
        current_time = time
        
    
    for nde, height in sorted_dict.items():
        
        for number, times in coalescent_intervals.items():            
            if not nde.node_parent:
                non_parent_set.add(nde)
         
                if height < times[1] and whole_tree.node_heights[whole_tree.root] >= times[1]:
                    active_population[number].append(nde) 
            
            elif height < times[1] and nde.node_parent.height >= times[1]:

                active_population[number].append(nde)
                
   #  print("active pop")
#     for k,v in active_population.items():
#         print(k,v)
        
    if len(non_parent_set) > 1:
        print("NODES WITHOUT PARENTS" + str(len(non_parent_set)))
        
    #print(coalescent_intervals)
    
    return active_population, coalescent_intervals

#Tree class
def calculate_ne(whole_tree, those_sampled):
    """Get effective population sizes in each coalescent interval"""
    
    waiting_times = {}
    
    result = get_active_population(whole_tree)
    
    active_population = result[0]
    coalescent_intervals = result[1]
        
    Ne_dict = {}
      
    for key, value in coalescent_intervals.items():
        waiting_times[key] = value[1] - value[0]
     
        
    
    for key, value in waiting_times.items():
        
        tau = value/365
        
        #print(tau)
        
        if tau == 0:
            print("tau is zero here")
            print(key, value)
        
        lineages = len(active_population[key])

        count_weird_trees = 0

        if lineages == 1: #This should work because the tree can never start with one lineage
            try:
                Ne = Ne
            except UnboundLocalError:
                count_weird_trees += 1
                Ne = 0.000000001
                tree_file = open("error_trees" + str(count_weird_trees) + ".csv", 'w')
                to_newick(whole_tree.root, whole_tree, those_sampled)

        else:
            Ne = np.log(special.binom(lineages,2)) + np.log(tau) 
        
        new_key = (coalescent_intervals[key][0], coalescent_intervals[key][1])
        
        #print(new_key)
        
        #print(new_key, Ne) 
        
        Ne_dict[new_key] = Ne
    
    #print(Ne_dict)
    
    return Ne_dict, coalescent_intervals
    
#Tree class
def plot_skyline(Ne_dict):

    for_plotting = {}
    
    
    for key, value in Ne_dict.items():
        
        new_key = key*-1
        
        for_plotting[new_key] = value

    names = list(for_plotting.keys())
    values = list(for_plotting.values())

    plt.step(names,values)


def simulate_tree(trans_dict, nodes, sampling_proportion, epidemic_len): 
    """Function to simulate a coalescent tree"""
            
    sampling_output = sampling(trans_dict, sampling_proportion, epidemic_len)
    
    if type(sampling_output) == bool:
        not_enough_cases = sampling_output
    else:
        not_enough_cases = sampling_output[1]
        those_sampled = sampling_output[0]
    
    if not_enough_cases:
        return

    if len(those_sampled) != 0:
    
        node_dict = {}
        gen_3 = 0
        gen_4 = 0
        
        #Intialise as type individual nodes, and make transmission tree
        for person in nodes:
            node_class.node(person, "Ind", trans_dict, those_sampled, node_dict, gen_3, gen_4)
            
            if node_dict[person].transm_root:
                transm_root = node_dict[person]
                
        
        R0 = get_R0(gen_3, gen_4)

        subtree_dict_outside = {}

        big_tree = coalescent_tree()

        subtree_dict = get_subtrees(transm_root, subtree_dict_outside)[1]
        
        big_tree = connections(subtree_dict[transm_root.id], big_tree, subtree_dict)[1]

        #This is the maximum number of days the simulation can run for 
        smallest_root = epidemic_len 

        #Find out which subtree will contain the root for the coalescent tree - ie with the oldest sampled tip
        for item in subtree_dict.values():
            if item.sample_tip:
                if float(item.sample_tip.time_infected) < float(smallest_root):
                    smallest_root = item.sample_tip.time_infected
                    root_tree = item

        root_tree.coal_root = True

        remove_internals(big_tree.root, big_tree)
        
        for nde in big_tree.nodes:
            if type(nde) == node:
                if not nde.sampled:
                    print(nde)

        for nde in big_tree.nodes:
            if not nde.removed:
                big_tree.final_nodes.add(nde)
                get_tip_to_root(nde, big_tree)
                
                
            

        longest = big_tree.most_recent_tip.root_to_tip

        for nde in big_tree.final_nodes:
            nde.height = longest - nde.root_to_tip
            big_tree.node_heights[nde] = nde.height
            
            if not nde.remove_func_called:
                print(big_tree.root.node_children)
                #print("removal error")
                #if type(nde) == coalescent_node:
                 #   print(nde,nde.node_parent.node_children, nde.node_parent.for_loop_called)
                #if type(nde) == node:
                #    print(nde, nde.node_parent.node_children, nde.id, nde.node_parent.for_loop_called)
                #if type(nde) == transmission_node:
                 #   print(nde, nde.node_parent.node_children, nde.id.id, nde.node_parent.for_loop_called)
                if nde.node_parent.for_loop_called:
                    print(nde.node_parent, node)
             
                   
            


        newick_tree = to_newick(big_tree.root, big_tree, those_sampled)

        Ne_dict, coal_intervals = calculate_ne(big_tree, those_sampled)

        return newick_tree, Ne_dict, big_tree, R0, big_tree.most_recent_date, those_sampled, coal_intervals



    else:
        return

