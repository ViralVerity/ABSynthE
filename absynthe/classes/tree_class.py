import random
import uuid
import numpy as np
from collections import defaultdict
from collections import OrderedDict
from scipy import special

import absynthe.classes.node_class as node_class

class tree():   
    
    def __init__(self, tree_type=None, focal_node=None, subtree_dict=None, node_dict=None, epidemic_len=None):
        
        #print("Making subtree for " + person_tree.id)
        
        self.nodes = set()
        self.tips = []
        self.branch_lengths = {}
        self.heights = {}
        
        if tree_type == "subtree":
            self.subtree = True 
            self.focal_node = focal_node #The node corresponding to the person whose subtree it is 
            self.contains_sample = False      
            self.transmits = False
            
            self.absolute_tip_times = []            
                    
            self.sort_out_tips()

            if len(self.tips) != 0:

                self.coalescent(self.tips, 0.0)
                self.define_root()
                self.get_branch_lengths()
                        
        elif tree_type == "whole_tree":
            #Whole tree coalescent characteristics
    
            self.whole_tree = True

            self.most_recent_date = 0.0
            
            #For fitting summary statistics
            self.oldest_sample_date = epidemic_len
            self.b_len_list = []
            self.internal_branches = []
            self.external_branches = []
            self.sample_times = []
            self.total_steps = []
            
            self.construct_tree(node_dict)
            self.final_nodes = self.nodes.copy()
            self.remove_internals(self.root)
            
            
            transmission_node_count = 0
            coalescent_node_count = 0
            
            for node in self.nodes:
                if node.type == "transmission":
                    transmission_node_count += 1
                if node.type == "coalescent":
                    coalescent_node_count += 1
                    
            for nde in self.final_nodes:
                self.get_tip_to_root(nde)
                self.heights[nde] = self.most_recent_date - nde.root_to_tip
                
                #For fitting summary stats
                self.b_len_list.append(self.branch_lengths[nde])
                self.internal_branches.append(self.branch_lengths[nde])
            
            for tip in self.tips:
                self.heights[tip] = self.most_recent_date - tip.time_sampled
                
                #For fitting summary stats
                self.b_len_list.append(self.branch_lengths[tip])
                self.external_branches.append(self.branch_lengths[tip])
                self.sample_times.append(tip.time_sampled)
                tip.steps = tip.node_parent.steps + 1
                self.total_steps.append(tip.steps)
            
            self.all_tips_nodes = self.tips.copy()
            self.all_tips_nodes.extend(self.final_nodes)
                    
            ##tests
            # if len(self.final_nodes) != (len(self.tips) - 1):
            #     print("error in final node number")
            # if len(self.nodes) != (trans_count + coal_count):
            #     print("error in pre-removal number")
             
                    

    def sort_out_tips(self): #subtree
        
        self.find_transmission_tips(self.focal_node)
        self.find_sample_tips(self.focal_node)
        
        if len(self.tips) != 0: 
            self.most_recent_tip = sorted([float(i) for i in self.absolute_tip_times])[::-1][0]
            self.root_time = self.most_recent_tip - float(self.focal_node.time_infected)
        
        for tip in self.tips:
            tip.relative_height = self.most_recent_tip - tip.absolute_time
            self.heights[tip] = self.most_recent_tip - tip.absolute_time #relative height in subtree of tip
    
    def find_transmission_tips(self, focal_node): #subtree

        if len(focal_node.sampled_infections) != 0: #if any downstream children sampled
            for case_node in focal_node.infections: #case is another node object
                if len(case_node.sampled_infections) != 0 or case_node.sampled: 
                
                    transmission_tip = node_class.node(case_node, "transmission", infector=focal_node, infectee=case_node) 
                    transmission_tip.absolute_time = case_node.time_infected
                    
                    self.transmits = True
                    self.absolute_tip_times.append(transmission_tip.absolute_time)
                    self.tips.append(transmission_tip)
    
    
    def find_sample_tips(self, focal_node): #subtree
    #don't need to make a new node in here because it already exists when we make the full tree including unsampled tips
        if focal_node.sampled: 

            self.absolute_tip_times.append(focal_node.time_sampled)
            self.tips.append(focal_node)
            
            self.contains_sample = True
            self.sample_time = focal_node.time_sampled
          
            
    def define_root(self):
        
        self.root = node_class.node(uuid.uuid1(), "coalescent", height=self.root_time, children=[self.penultimate], subtree=self)
        self.heights[self.root] = self.root_time      
        self.nodes.add(self.root)
        self.penultimate.node_parent = (self.root)
        self.branch_lengths[self.root] = 0.0

    def coalescent(self, lineage_list, current_height):    
        
        none_left = False
        
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
                population_demographic_model = 1
                N0 = 1
                Ne = population_demographic_model*N0 #At the moment, constant population size - this is inside a person

                #Wait time till next event is drawn therefore event happens at current height + tau
                tau = np.random.exponential(number_of_lineages*(number_of_lineages-1)/(2*Ne)) 
               
                #if it doesn't happen before the next lineage comes in
                if not none_left and current_height+tau > next_sample.relative_height: 
                    
                    current_height = next_sample.relative_height
                    self.coalescent(lineage_list, current_height)

                else:   
                    #Have this in for the moment to force it happen before the root of the subtree
                    if current_height + tau > self.root_time:
                        tau = self.root_time - current_height
                   
                    #Who is going to coalesce?
                    lucky_pair = random.sample(active_pop, k=2)
                    current_height += tau
                                            
                    #ie the coalescent event of the pair selected above
                    parent_node = node_class.node(uuid.uuid1(), "coalescent", height=current_height, children=lucky_pair, subtree=self)
                    print(f"new parent node being assigned is: {parent_node}")
                    print(f"To: {lucky_pair[0]} and {lucky_pair[1]}")
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
            if tip.type == "individual" and subtree.contains_sample:
                self.tips.append(tip)
            else:
                self.nodes.add(tip)

        for nde in subtree.nodes:
            self.nodes.add(nde)

    
    def construct_tree(self, node_dict):
        
        for nde in node_dict.values():
            if nde.sampled or len(nde.sampled_infections) != 0:
                subtree = nde.subtree
                if subtree.focal_node.index_case: 
                    self.root = subtree.root #If I want the subtree with the root in as well I can do it here

                if subtree.contains_sample: #it will only ever contain one sample
                    if subtree.most_recent_tip >= self.most_recent_date:
                        self.most_recent_date = subtree.most_recent_tip
                        
                    if subtree.sample_time <= self.oldest_sample_date:
                        self.oldest_sample_date = subtree.sample_time

                if subtree.transmits:
                    for tip in subtree.tips:
                        if tip.type == "transmission":
                            
                            donor_tree = tip.infector.subtree
                            recipient_tree = tip.infectee.subtree

                            recipient_tree.root.node_parent = tip
                            tip.node_children.append(recipient_tree.root) 

                self.update_coalescent_tree(subtree)
                
    def remove_internals(self, nde): 
        
        nde.remove_func_called = True
        nde.old_children = nde.node_children.copy()

        if len(nde.old_children) == 1: 
            if nde != self.root:
                parent = nde.node_parent 
                nde.removed = True
                parent.node_children.remove(nde)
                
                self.final_nodes.remove(nde)

                #Reassign parents and children to remove internal nodes
                for child in nde.old_children:     
                    child.node_parent = parent
                    parent.node_children.append(child)
                    self.branch_lengths[child] = self.branch_lengths[child] + self.branch_lengths[nde]

            else:
                for child in nde.node_children:
                    child.node_parent = None
                    self.branch_lengths[child] = 0.0
                    self.root = child

                nde.removed = True
                self.final_nodes.remove(nde)

        
        for i in nde.old_children: #it has to go through all of them 
            self.remove_internals(i)
        
    
    def get_tip_to_root(self, nde):

        try:
            if nde.root_to_tip != 0.0:
                return nde.root_to_tip
        except AttributeError:
            print(f"No root to tip for {nde} {nde.id} {nde.type}")

        if nde == self.root:
            distance = 0.0
            steps = 0

        else:
            distance = self.get_tip_to_root(nde.node_parent) + self.branch_lengths[nde]
            steps = nde.node_parent.steps + 1
            
        nde.root_to_tip = distance
        nde.steps = steps
        
        return nde.root_to_tip
    
    
    ###From here on we'll call in the main script not in the init function
    def to_newick(self, nde, those_sampled): #call on root of tree 
        
        string = ",".join([self.to_newick(i, those_sampled) for i in nde.node_children if not i.removed])

        if len(nde.node_children) != 0:
            string = "(" + string
        if nde.type == "individual":
            string += str(nde.id)
        else:
            string += ")"

        string += ":" + str(self.branch_lengths[nde])

        if nde == self.root:
            string += ")"

        return string
    
    
    def get_active_population(self): 
        
        coalescent_times = []
        waiting_times = {}
        bins = []
        active_population = {}

        height_set = set()
        for i in self.heights.values():
            height_set.add(i)

        height_list = list(height_set)
        sorted_height_list = sorted(height_list)
        sorted_height_dict = {k: v for k, v in sorted(self.heights.items(), key=lambda item: item[1])}
        # prep = set(sorted_heights.values())
        # coalescent_times = set(prep)

        for height in sorted_height_list:
            if height != 0:
                coalescent_times.append(height)

        current_time = 0
        non_parent_set = set() 
        for time in coalescent_times:
            
            tup = (float(current_time),float(time))
            bins.append(tup)
            
            waiting_times[tup] = time - current_time
            active_population[tup] = 0
            
            current_time = time
                        
        active_population = {}

        for bin_pair in bins:
            active_population[bin_pair] = 0
            for node, height in sorted_height_dict.items():
                if not node.node_parent: #then it's the root
                    non_parent_set.add(node)
                    continue
                else:
                    start = height
                    end = self.heights[node.node_parent] 
                    if start <= bin_pair[0] and end > bin_pair[0]:
                        active_population[bin_pair] += 1

        if len(non_parent_set) > 1:
            print("NODES WITHOUT PARENTS" + str(len(non_parent_set)))

        
        return active_population, waiting_times, coalescent_times
    
    
    def calculate_ne(self, those_sampled):
        """Get effective population sizes in each coalescent interval for skyline"""

        active_population, waiting_times, coalescent_times = self.get_active_population() #active_population is also lineages through time
        Ne_dict = {}
        Ne = 0
        for times, tau in waiting_times.items():
            if tau <= 0:
                print("tau is not right here")
                print(times, tau)

            lineages = active_population[times]
            if lineages > 1:
                a = np.log(special.binom(lineages,2))
                b = np.log(tau)
                Ne = a + b
            else:
                Ne = Ne #ie from last time around - is that right?
            
            Ne_dict[times] = Ne


        #so this was if the tree starts with one lineage. Commenting out because that should never happen and it should break if it does
        # except UnboundLocalError:
        #     count_weird_trees += 1
        #     Ne = 0.000000001
        #     tree_file = open("error_trees" + str(count_weird_trees) + ".csv", 'w')
        #     to_newick(whole_tree.root, whole_tree, those_sampled)
        
        self.lineages_through_time = active_population
        self.coalescent_times = coalescent_times
                
        return Ne_dict, active_population, coalescent_times
    
    

    
        
        


        
        
        
                
                
                
                
                
                
                
                
                
                