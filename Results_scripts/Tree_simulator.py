from collections import defaultdict
from collections import OrderedDict
import random
import uuid
import numpy as np
from scipy import special

# ### Define node and tree classes

class node():
    
    def __init__(self, unique_id):
        
        self.id = unique_id
        self.sampled = False
        
        self.time_sampled = 0 #?

        self.children = set()
        
        self.sampled_children = set()
        
        self.to_root = []
        
        self.time_infected = 0 #?
        
        self.coalescent_time = [] #?
        
        self.height = 0.0 #?
        
        self.relative_height = 0 #?
        
        self.transm_root = False
        
        self.last = False
        
        self.removed = False
        
        self.node_children = set()
        
        self.root_to_tip = 0.0
        
        self.remove_func_called = False
        
        self.for_loop_called = False
        
        
class transmission_node():
    
    def __init__(self, unique_id):
        
        self.id = unique_id
        
        self.sampled = False
        
        self.time_sampled = 0 

        self.node_children = set()
        
        self.sampled_children = set()
        
        self.to_root = []
        
        self.time_infected = 0
        
        self.coalescent_time = []
        
        self.removed = False
        
        self.height = 0.0
        
        self.relative_height = 0
        
        self.root = False
        
        self.last = False
        
        self.root_to_tip = 0.0
        
        self.remove_func_called = False
        
        self.for_loop_called = False
        
class coalescent_node():
    
    def __init__(self, unique_id, height, lucky_pair, subtree):
        
        self.id = unique_id
        self.relative_height = height
        
        self.node_children = set()
        
        self.pair = lucky_pair
        
        self.transmission_node = False
        self.sampling_node = False #?
        self.internal_node = True #?
        
        self.last = False
        
        self.removed = False
        self.root_to_tip = 0.0
        
        self.height = 0.0
        self.remove_func_called = False
        
        self.subtree = subtree
        
        self.for_loop_called = False

class subtree():
    
    def __init__(self, person_tree):
        
        self.id = person_tree #will be the person it corresponds to
        
        self.person = person_tree
        
        self.transmission_tips = [] #Subtree objects
        
        self.sample_tip = None
        
        self.most_recent_tip = 0 #time wise when is the most recent, useful for the relative heights
        
        self.coalescent_nodes = [] #internal nodes as a list of coalescent_node objects
        
        self.root_time = 0
        
        self.relative_heights = {}
        
        self.branch_lengths = {}
        
        self.coal_root = False
        
    def get_branch_lengths(self):
        
        if not self.person.sampled and len(self.transmission_tips) == 1: #unsampled middle node 
            self.branch_lengths[self.transmission_tips[0]] = self.relative_heights[self.root] - self.relative_heights[self.transmission_tips[0]]
            self.branch_lengths[self.root] = 0.0
        
        elif self.person.sampled and len(self.transmission_tips) == 0: #sampled tip
            self.branch_lengths[self.sample_tip] = self.relative_heights[self.root] - self.relative_heights[self.sample_tip]
            self.branch_lengths[self.root] = 0.0
            
        else:
        
            for nde in self.transmission_tips:
                self.branch_lengths[nde] = self.relative_heights[nde.node_parent] - self.relative_heights[nde]

            for nde in self.coalescent_nodes:
                if not nde.last:
                    self.branch_lengths[nde] = self.relative_heights[nde.node_parent] - self.relative_heights[nde]
                else:
                    #should just be the root
                    self.branch_lengths[nde] = 0.0

            if self.person.sampled:
                self.branch_lengths[self.person] = self.relative_heights[self.sample_tip.node_parent] - self.relative_heights[self.sample_tip]

class coalescent_tree():
    
    def __init__(self):
        
        self.nodes = set()
        
        self.branch_lengths = {}
        
        self.node_heights = {}
        
        self.most_recent_date = 0.0
        
        self.final_nodes = set() #non-removed ones
        


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


def get_parent(trans_dict, person, node_dict):
    """Get parent as node object"""
    
    if trans_dict[person.id][0] in node_dict.keys():
        
        person.parent = node_dict[trans_dict[person.id][0]]
        person.parent.children.add(person)
    
    else: #if the parent is not yet in the node dict, make the parent. They shouldn't be made again because they'll now be in node_dict.keys()

        parent_input = trans_dict[person.id][0]

        person.parent = node(parent_input)

        if person.parent.id != "NA":
           
            input1 = trans_dict[person.parent.id][1] 
            input2 = trans_dict[person.parent.id][2] 

        
            uniform1 = np.random.uniform(0,1)
            
            person.parent.time_infected = float(trans_dict[person.parent.id][1]) + uniform1
            
            
            if input1 == input2:
                rnge = 1 - uniform1
                person.parent.time_sampled = float(person.parent.time_infected) + np.random.uniform(0,rnge)
                
            else:
                person.parent.time_sampled = float(trans_dict[person.parent.id][2])

        person.parent.children.add(person)

        node_dict[parent_input] = person.parent
        
    return person.parent
            

def initialise_person(input_person, trans_dict, those_sampled, node_dict):
    """Initialise the person as a node object"""
    
    #Making sure everyone is intitialised with useful information from the simulation output file
    if type(input_person) != node and input_person not in node_dict.keys(): #If person has not already been initialised

        person = node(input_person)
       
        input1 = trans_dict[person.id][1] 
        input2 = trans_dict[person.id][2] 

        uniform1 = np.random.uniform(0,1)

        person.time_infected = float(trans_dict[person.id][1]) + uniform1


        if input1 == input2:
            rnge = 1 - uniform1
            person.time_sampled = float(person.time_infected) + np.random.uniform(0,rnge)

        else:
            person.time_sampled = float(trans_dict[person.id][2]) 
            
        node_dict[input_person] = person
        
        person.parent = get_parent(trans_dict, person, node_dict)
        
            
    #Because we loop through every person, it's possible they've already been initialised
    elif input_person in node_dict.keys(): 
        person = node_dict[input_person]
    
    #Mostly for the recursion section, where the input is a node that we're interested in doing the next step with
    else:
        person = input_person
       
    
    return person

def get_root_list(input_person, trans_dict, those_sampled, node_dict): 
    """Get the full transmission tree for each individual"""
    
    person = initialise_person(input_person, trans_dict, those_sampled, node_dict)
        
    #Caching to save time
    if len(person.to_root) != 0:
        return(person.to_root)
        
    else:

        person.parent = get_parent(trans_dict, person, node_dict)
            
        #If the person is the root
        if person.parent.id == "NA": 
            path = []
            person.transm_root = True
        
        #Recur up the tree
        else:
            path = get_root_list(person.parent, trans_dict, those_sampled, node_dict) + [(person.parent)]

            
        person.to_root = path
        
        #Getting the sampling 
        if person.id in those_sampled:
            person.sampled = True
            for i in person.to_root:
                i.sampled_children.add(person) #So the sampled children is all sampled children downstream of the focal individual

        return path
    
def get_R0(node_dict):
    #Needs the gen node stuff and the get_root_list stuff. Let's generate along with the tree and just try and store it somewhere.
    gen_4 = []
    gen_3 = []
    
    for person in node_dict.values():
        
        if len(person.to_root) == 3:
        
            gen_4.append(person)
          
        elif len(person.to_root) == 2:
            
            gen_3.append(person)
            
        else:
            pass
        
    if len(gen_4) != 0 and len(gen_3) != 0:
        R0 = len(gen_4)/len(gen_3)
    
        return R0
    else:
        
        return
        


def get_subtrees(person, subtree_dict): #so this is the person's subtree we're working out
    
    """ Get each individual's virus tree, taking sampling into account"""   
    
    #If the person isn't sampled and that branch has no sampled children on it
    if not person.sampled and len(person.sampled_children) == 0:
        pass
   
    else:

        tree = subtree(person)

        subtree_dict[person.id] = tree

        tips = []
        lin_list = []
        current_height = 0 #Ie for each new subtree, we set this to zero

        #Assemble all tips involved to find most recent and get relative heights, and to get full lineage list for coalescence
        #Possibly get rid of relative height things, may not need them?

        #Transmission tips 
        if len(person.sampled_children) != 0:
            for i in person.children:
                if len(i.sampled_children) != 0 or i.sampled:

                    j = transmission_node(i) 

                    j.absolute_time = i.time_infected

                    tips.append(j.absolute_time)

                    lin_list.append(j)

                    tree.transmission_tips.append(j)

        #Sample tips
        if person.sampled: 
            tips.append(person.time_sampled)
            lin_list.append(person)
            tree.sample_tip = person

        tree.most_recent_tip = sorted([float(i) for i in tips])[::-1][0]

        tree.root_time = tree.most_recent_tip - float(person.time_infected)

        #Get relative heights for other tips in subtree
        for trans_node in tree.transmission_tips:
            trans_node.relative_height = tree.most_recent_tip - float(trans_node.absolute_time)
            tree.relative_heights[trans_node] = tree.most_recent_tip - float(trans_node.absolute_time)

        if person.sampled:
            person.relative_height = tree.most_recent_tip - float(person.time_sampled)
            tree.relative_heights[person] = tree.most_recent_tip - float(person.time_sampled)

        #Traverse the full transmission tree
        for i in person.children:
            if not i.sampled and len(i.sampled_children) == 0: #If there are no more samples on that branch
                pass
            else:
                get_subtrees(i, subtree_dict)


        #TEST that asks whether there is only one 0.0 in the heights

        if len(lin_list) == 1: #ie if they are a straight line (sampled tip or a non-sampled transmission)

            #no coalescent event needed and tree is just one line
            tree.root = coalescent_node(uuid.uuid1(), tree.root_time, person, tree) 

            tree.relative_heights[tree.root] = tree.root_time

            #NB it's [0] because it's a list of one element
            lin_list[0].node_parent = tree.root
            tree.root.node_children.add(lin_list[0])

            tree.coalescent_nodes.append(tree.root) 

            tree.get_branch_lengths()

            return tree, subtree_dict


        else: #if there are few things we need to coalesce eg transmission events, or sampling of an internal node  

            coalescent(lin_list, current_height, tree)

            #penultimate is assigned in coalescent
            tree.root = coalescent_node(uuid.uuid1(), tree.root_time, tree.penultimate, tree)

            tree.penultimate[0].node_parent = tree.root

            tree.root.node_children.add(tree.penultimate[0])

            tree.relative_heights[tree.root] = tree.root_time

            tree.root.last = True

            tree.coalescent_nodes.append(tree.root)

            tree.get_branch_lengths()

            return tree, subtree_dict

    
def sort_key(ele):
    """Small function used later to sort lists by relative height"""
    return ele.height

def coalescent(lineage_list, current_height, tree):
    """Get the coalescence times of any lineages that need to coalesce"""
    
    #Lineage list gets updated as lineages coalesce throughout the running of the function
    
    none_left = False
    
    if len(lineage_list) == 0:
        print("ERROR HERE, LIN LIST IS 0" + tree.id.id)
    
    if len(lineage_list) == 1: #one before the root of the subtree
        tree.penultimate = [lineage_list[0]]
        return 

    else:
        
        active_pop = [i for i in lineage_list if i.relative_height <= current_height]
        to_be_sampled = [i for i in lineage_list if i.relative_height > current_height]

        #More to be sampled? 
        if len(to_be_sampled) != 0:
            next_sample = next(iter(sorted(to_be_sampled,key=sort_key)))
        else:
            none_left = True

        if len(active_pop) == 1 and not none_left: 
        #If there's only one in the active population at this time point, but there's still more to be coalesced
            
            current_height = next_sample.relative_height
            
            coalescent(lineage_list, current_height, tree)

        #Got some stuff in active population to be coalescing
        else:
            
            j = len(active_pop) 

            populationDemographicModel = 1

            N0 = 1

            Ne = populationDemographicModel*N0 #At the moment, constant population size - this is inside a person

            #Wait time till next event is drawn therefore event happens at current height + tau

            tau = np.random.exponential(j*(j-1)/(2*Ne)) 

            #if it doesn't happen before the next lineage comes in
            if not none_left and current_height+tau > next_sample.relative_height: 
                
                current_height = next_sample.relative_height

                coalescent(lineage_list, current_height, tree)

            else:
                
                #Who is going to coalesce?
                lucky_pair = random.sample(active_pop, k=2)
                    
                #Have this in for the moment to force it happen before the root of the subtree
                if current_height + tau > tree.root_time:
                    tau = tree.root_time - current_height
                    #print("forced tau to zero in " + lucky_pair[0].subtree.id.id)
                    
            
                current_height += tau
                
                #ie the coalescent event of the pair selected above
                parent_node = coalescent_node(uuid.uuid1(), current_height, lucky_pair, tree)
                
                lucky_pair[0].node_parent = parent_node
                lucky_pair[1].node_parent = parent_node
                
                parent_node.node_children.add(lucky_pair[0])
                parent_node.node_children.add(lucky_pair[1])
                        
                tree.coalescent_nodes.append(parent_node)

                tree.relative_heights[parent_node] = current_height
                
                #Update lineage and see if it needs coalescing again
                updated_population = [i for i in lineage_list if i not in lucky_pair] + [parent_node]
                                
                coalescent(updated_population, current_height, tree)


def update_coalescent_tree(subtree, whole_tree):
    """Update the coalescent tree with information from the subtree"""
    
    whole_tree.branch_lengths.update(subtree.branch_lengths)
    
    if subtree.sample_tip:
        whole_tree.nodes.add(subtree.sample_tip)

    for tip in subtree.transmission_tips:
        whole_tree.nodes.add(tip)

    for internal_node in subtree.coalescent_nodes:
        whole_tree.nodes.add(internal_node)
    
    

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


def remove_extra_coalescent_events(whole_tree_node, whole_tree): #root of subtree should already have been removed
    
    if type(whole_tree_node) == coalescent_node:
        #print(whole_tree_node)
    
        if len(whole_tree_node.subtree.coalescent_nodes) > 2:
            print("testing for removal" + str(whole_tree_node))
            print("subtree is " + str(whole_tree_node.subtree.id.id))
            if whole_tree_node != whole_tree_node.subtree.penultimate[0]:
                print("keep " + str(whole_tree_node.subtree.penultimate[0]))
                print("removing " + str(whole_tree_node))
                print("subtree " + str(whole_tree_node.subtree.coalescent_nodes))
                for i in whole_tree_node.new_children:
                    i.node_parent = whole_tree_node.node_parent

            else:

                whole_tree.final_node_heights[whole_tree_node] = whole_tree_node.height
        else:
            
            whole_tree.final_node_heights[whole_tree_node] = whole_tree_node.height
            
    else:
        whole_tree.final_node_heights[whole_tree_node] = whole_tree_node.height
            
    
    
    return whole_tree.final_node_heights


def get_active_population(whole_tree):
    """Get active population at each coalescent interval"""
    
    coalescent_times = set()
    coal_list = []
    coalescent_intervals = defaultdict(tuple)
    
    active_population = defaultdict(list)
    
    #whole_tree.final_node_heights = {}
    
    #for nde in whole_tree.final_nodes:
        #whole_tree.final_node_heights = remove_extra_coalescent_events(nde, whole_tree)
    
    #sorted_dict = OrderedDict(sorted(whole_tree.final_node_heights.items(), key=lambda x:x[1]))  
    sorted_dict = OrderedDict(sorted(whole_tree.node_heights.items(), key=lambda x:x[1]))
#     print("heights = ")
    #print(sorted_dict)
    for nde, height in sorted_dict.items():
       
        if type(nde) == coalescent_node:
            
            coalescent_times.add(height)
            coal_list.append(height)
            
            #for nde2,height2 in sorted_dict.items():
             #   if height == height2 and nde != nde2:
              #      if type(nde) == coalescent_node and type(nde2) == coalescent_node:
               #         print(height, height2, nde, nde2, nde.subtree.id.id, nde2.subtree.id.id)
                #        print(nde.subtree.relative_heights, nde.subtree.root_time)
     
    
    # print("coal times = ")
#     print(sorted(coalescent_times))
#     print(coal_list)
    
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
    
    
def plot_skyline(Ne_dict):

    for_plotting = {}
    
    
    for key, value in Ne_dict.items():
        
        new_key = key*-1
        
        for_plotting[new_key] = value

    names = list(for_plotting.keys())
    values = list(for_plotting.values())

    plt.step(names,values)


def simulate_tree(trans_dict, nodes, sampling_proportion, epidemic_len): #change these arguments when putting in module for reading in
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

        #Traverse the tree and find out what's going on
        for person in nodes:
            lst = get_root_list(person, trans_dict, those_sampled, node_dict)
            if node_dict[person].transm_root:
                transm_root = node_dict[person]
                
        R0 = get_R0(node_dict)

        subtree_dict_outside = {}

        big_tree = coalescent_tree()

        subtree_dict = get_subtrees(transm_root, subtree_dict_outside)[1]
        
        big_tree = connections(subtree_dict[transm_root.id], big_tree, subtree_dict)[1]

        #This is the maximum number of days the simulation can run for 
        smallest_root = 148 

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

