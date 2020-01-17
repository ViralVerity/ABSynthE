import node_class as nc
import random
import uuid
import numpy as np
from collections import defaultdict
from collections import OrderedDict
from scipy import special

class tree():   
    
    def __init__(self, person_tree=None, subtree_dict=None, node_dict=None):
        
        #print("Making subtree for " + person_tree.id)
                
        self.nodes = set()
        self.tips = []
        self.branch_lengths = {}
        self.heights = {}
        
        if person_tree:
            self.subtree = True #a lil flag
            self.person = person_tree #will be the person as a node object it corresponds to if it's a subtree
            self.contains_sample = False      
            self.transmits = False
            
            self.absolute_tip_times = []            
                    
            self.sort_out_tips()

            if len(self.tips) != 0:
                
                self.coalescent(self.tips, 0.0)
            
                self.define_root()

                self.get_branch_lengths()
                

                        
        else:
            #Whole tree coalescent characteristics
    
            self.whole_tree = True

            self.most_recent_date = 0.0
            
            self.construct_tree(node_dict)
            
            self.final_nodes = self.nodes.copy()
            
            self.remove_internals(self.root)
            
            trans_count = 0
            coal_count = 0
            
            for i in self.nodes:
                if i.type == "Trans":
                    trans_count += 1
                if i.type == "Coal":
                    coal_count += 1
                    
                    
            if len(self.final_nodes) != (len(self.tips) - 1):
                print("error in final node number")
            if len(self.nodes) != (trans_count + coal_count):
                print("error in pre-removal number")
             
            for nde in self.final_nodes:
                
                if not nde.remove_func_called:
                    print("function not called on")
                    print("Node: " + str(nde))
                    print("Type: " + str(nde.type))
                    if nde.type == "Coal":
                        print("From subtree: " + str(nde.subtree.person.id))
                        print("Sampled? " + str(nde.subtree.contains_sample))
                    try:
                        print("Parent = " + str(nde.node_parent))
                    except AttributeError:
                        print("No parent")
                    print("Children = " + str(nde.node_children))
                
                elif nde.removed:
                    print("error in removing")
                    print("Node: " + str(nde))
                    print("Type: " + str(nde.type))
                    try:
                        print("Parent = " + str(nde.node_parent))
                    except AttributeError:
                        print("No parent")
                    print("Children = " + str(nde.node_children))
                
                elif nde.type == "Trans":
                    print("Trans node still in ")
                    print("Node: " + str(nde))
                    print("Type: " + str(nde.type))
                    try:
                        print("Parent = " + str(nde.node_parent))
                    except AttributeError:
                        print("No parent")
                    print("Children = " + str(nde.node_children))
                
                else:
                    self.get_tip_to_root(nde)
                    self.heights[nde] = self.most_recent_tip - nde.root_to_tip
                
            for tip in self.tips:
                self.heights[tip] = self.most_recent_tip - tip.time_sampled
            
            
    
    def sort_out_tips(self):
        
        focal_individual = self.person
        self.find_transmission_tips(focal_individual)
        
        self.find_sample_tips(focal_individual)
        
        if len(self.tips) != 0: 
            self.most_recent_tip = sorted([float(i) for i in self.absolute_tip_times])[::-1][0]
            self.root_time = self.most_recent_tip - float(focal_individual.time_infected)
        
        for tip in self.tips:
            tip.relative_height = self.most_recent_tip - tip.absolute_time
            self.heights[tip] = tip.relative_height 
            
    
    def find_transmission_tips(self, focal_individual):

        if len(focal_individual.sampled_infections) != 0:

            for case in focal_individual.infections:
                
                if len(case.sampled_infections) != 0 or case.sampled: 
                    
                    transmission_tip = nc.node(case, "Trans", infector=self.person, infectee=case) 
                    transmission_tip.absolute_time = case.time_infected
                    self.transmits = True
                    
                    self.absolute_tip_times.append(transmission_tip.absolute_time)
                    self.tips.append(transmission_tip)
    
    
    def find_sample_tips(self, focal_individual):
        
        if focal_individual.sampled: 
       
            self.absolute_tip_times.append(focal_individual.time_sampled)
            self.tips.append(focal_individual)
            
            self.contains_sample = True
          
            
    def define_root(self):
        
        self.root = nc.node(uuid.uuid1(), "Coal", height=self.root_time, children=[self.penultimate], subtree=self)
        
        self.heights[self.root] = self.root_time
                
        self.nodes.add(self.root)
        
        self.penultimate.node_parent = (self.root)
        
        self.branch_lengths[self.root] = 0.0


        
    def coalescent(self, lineage_list, current_height):    
        
        none_left = False
        
        if len(lineage_list) == 0:
            print("ERROR HERE, LIN LIST IS 0" + self.person.id)
        
        if len(lineage_list) == 1:
            self.penultimate = lineage_list[0]
            return
        
        else:
            
            active_pop = [i for i in lineage_list if i.relative_height <= current_height]
            to_be_sampled = [i for i in lineage_list if i.relative_height > current_height]
        
            #More to be sampled? 
            if len(to_be_sampled) != 0:
                next_sample = next(iter(sorted(to_be_sampled,key=lambda nde:nde.relative_height)))
            else:
                none_left = True
                
                
            if len(active_pop) == 1 and not none_left: 
            #If there's only one in the active population at this time point, but there's still more to be coalesced
            
                current_height = next_sample.relative_height
            
                self.coalescent(lineage_list, current_height)
            
            #Got some stuff in active population to be coalescing
            else:
                
                number_of_lineages = len(active_pop) 
                populationDemographicModel = 1
                N0 = 1
                Ne = populationDemographicModel*N0 #At the moment, constant population size - this is inside a person

                #Wait time till next event is drawn therefore event happens at current height + tau

                tau = np.random.exponential(number_of_lineages*(number_of_lineages-1)/(2*Ne)) 
               
                
                #if it doesn't happen before the next lineage comes in
                if not none_left and current_height+tau > next_sample.relative_height: 
                    
                    current_height = next_sample.relative_height
                    
                    self.coalescent(lineage_list, current_height)

            
                else:   
        
                    #Have this in for the moment to force it happen before the root of the subtree
                    if current_height + tau > self.root_time:
                        #print("root condition triggered")
                        tau = self.root_time - current_height
                   
                    #Who is going to coalesce?
                    lucky_pair = random.sample(active_pop, k=2)
                    #print("Coalescing " + str(lucky_pair))
                     
                    #print("current height = " + str(current_height))
                        
                    current_height += tau
                    
                    #print("tau = " + str(tau))
                        
                    #ie the coalescent event of the pair selected above
                    parent_node = nc.node(uuid.uuid1(), "Coal", height=current_height, children=lucky_pair, subtree=self)
                    #print("made a node" + str(parent_node))
                    #print("new node is " + str(parent_node.id) + " " + str(parent_node.relative_height))
                    
                    lucky_pair[0].node_parent = parent_node
                    lucky_pair[1].node_parent = parent_node
                    
                    self.nodes.add(parent_node)
                    self.heights[parent_node] = current_height
                    
                    updated_population = [i for i in lineage_list if i not in lucky_pair] + [parent_node]
                                
                    self.coalescent(updated_population, current_height)
    
    def get_branch_lengths(self): 
        
        for tip in self.tips:
            self.branch_lengths[tip] = self.heights[tip.node_parent] - self.heights[tip]

        for nde in self.nodes:
            #print("calculated branch lengths for " + str(nde))
            
            if nde != self.root:
                self.branch_lengths[nde] = self.heights[nde.node_parent] - self.heights[nde]
            else:
                self.branch_lengths[nde] = 0.0


    def update_coalescent_tree(self, subtree):
        """Update the coalescent tree with information from the subtree"""

        self.branch_lengths.update(subtree.branch_lengths)

        for tip in subtree.tips:
            if tip.type == "Ind" and subtree.contains_sample:
                self.tips.append(tip)
            else:
                self.nodes.add(tip)

        for nde in subtree.nodes:
            self.nodes.add(nde)

    
    def construct_tree(self, node_dict):
        
        for nde in node_dict.values():
            if nde.sampled or len(nde.sampled_infections) != 0:
              
                subtree = nde.subtree

                if subtree.person.index_case: 
                    self.root = subtree.root #If I want the subtree with the root in as well I can do it here

                if subtree.contains_sample:
                    if subtree.most_recent_tip >= self.most_recent_date:

                        #Will need to check that this is always a sample. In between it may be transmission but that's ok
                        self.most_recent_tip = subtree.most_recent_tip
                        self.most_recent_date = float(subtree.most_recent_tip)

                if subtree.transmits:
                    for tip in subtree.tips:
                        
                        if tip.type == "Trans":
                            
                            #print("trans tip = " + str(tip))
                            
                            donor_tree = tip.infector.subtree
                            recipient_tree = tip.infectee.subtree

                            recipient_tree.root.node_parent = tip
                            tip.node_children.append(recipient_tree.root) 


                self.update_coalescent_tree(subtree)
            
            
                
                
    def remove_internals(self, nde): 
        
        nde.remove_func_called = True
        
        nde.new_children = nde.node_children.copy()
        
        #print("node children at start of function")
        #print(nde.node_children)

        if len(nde.node_children) == 1: 
            if nde != self.root:
                
                parent = nde.node_parent 

                nde.removed = True
                parent.new_children.remove(nde)
                
                self.final_nodes.remove(nde)

                #Reassign parents and children to remove internal nodes
                for child in nde.node_children:
                    
                    #print(child, nde)
                    
                    child.node_parent = parent

                    parent.new_children.append(child)
                    
                    self.branch_lengths[child] = self.branch_lengths[child] + self.branch_lengths[nde]


            else:

                for child in nde.node_children:

                    child.node_parent = None
                    self.branch_lengths[child] = 0.0
                    self.root = child

                nde.removed = True
                
                self.final_nodes.remove(nde)

        #print("node children at end of function")
        #print(nde.node_children)  
        
        for i in nde.node_children:
            
            self.remove_internals(i)
                
        
        return self
    
    def get_tip_to_root(self, nde):

        try:
            if nde.root_to_tip != 0.0:
                return nde.root_to_tip
        except AttributeError:
            print("No root to tip for " + str(nde) + " " + str(nde.id) + " " + str(nde.type))

        if nde == self.root:
            distance = 0.0

        else:

            #print(nde, nde.type, nde.removed)
            #if nde == nde.subtree.root:
             #   print(nde)
            

            distance = self.get_tip_to_root(nde.node_parent) + self.branch_lengths[nde]

                     

        nde.root_to_tip = distance
        

        return nde.root_to_tip
    
    
    ###From here on we'll call in the main script not in the init function
    def to_newick(self, nde, those_sampled): #call on root of tree 
        
        string = ",".join([self.to_newick(i, those_sampled) for i in nde.new_children if not i.removed])

        if len(nde.new_children) != 0:
            string = "(" + string

        if nde.type == "Ind":

            string += str(nde.id)
        else:

            string += ")"

        string += ":" + str(round(self.branch_lengths[nde], 2))

        if nde == self.root:
            string += ")"

        return string
    
    
    def get_active_population(self):
        """Get active population at each coalescent interval"""
    
        coalescent_times = set()
        #coal_list = []
        coalescent_intervals = defaultdict(tuple)

        active_population = defaultdict(list)

        sorted_dict = OrderedDict(sorted(self.heights.items(), key=lambda x:x[1]))

        for nde, height in sorted_dict.items():

            if nde.type == "Coal":

                coalescent_times.add(height)
                #coal_list.append(height)

        coalescent_times = sorted(coalescent_times)

        current_time = 0

        count = 0

        non_parent_set = set() 

        for time in coalescent_times:

            count += 1

            coalescent_intervals[count] = (float(current_time),float(time))
            
            coalescent_points.add(time)

            current_time = time


        parent_index = 0

        for nde, height in sorted_dict.items():
            index = parent_index - 1 #the minus one so it can be in the same interval as its parent
            
            if not nde.node_parent:
                active_population[len(coalescent_points)].append(nde)
                non_parent_set.add(nde)

            for times in coalescent_points[parent_index:]:            

                index += 1
                
                if height < times and self.heights[nde.node_parent] >= times:

                    parent_index = index
                    
                    active_population[index + 1].append(nde)

                    break #because it can only appear once in the list
        if len(non_parent_set) > 1:
            print("NODES WITHOUT PARENTS" + str(len(non_parent_set)))

        return active_population, coalescent_intervals, sorted_dict
    
    
    def calculate_ne(self, those_sampled):
        """Get effective population sizes in each coalescent interval"""

        waiting_times = {}

        result = self.get_active_population()

        active_population = result[0]
        coalescent_intervals = result[1]
        sorted_dict = result[2]

        Ne_dict = {}

        for key, value in coalescent_intervals.items():
            waiting_times[key] = value[1] - value[0]

        for key, value in waiting_times.items():

            tau = value/365

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


            Ne_dict[new_key] = Ne


        return Ne_dict, coalescent_intervals, sorted_dict
    
    
    

    
        
        


        
        
        
                
                
                
                
                
                
                
                
                
                