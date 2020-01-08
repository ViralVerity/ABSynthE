from node_class import *
import random
import uuid
import numpy as np

class tree():
    #Need to work out when subtree init happens in the code - new loop? on individual node init?
    
    def __init__(self, person_tree=None, subtree_dict=None):
        
        self.nodes = []
        self.tips = []
        self.branch_lengths = {}
        self.heights = {}
        
        if person_tree:
            self.person = person_tree #will be the person as a node object it corresponds to if it's a subtree
        
            #Subtree characteristics
            self.transmission_tips = [] #Subtree objects (that it's transmitting to?) This is quite odd, may be able to tidy away
            #self.sample_tip = None    Might need this not sure - might be a point later if I ask if there's a sample
            #self.coalescent_nodes = [] #internal nodes as a list of coalescent_node objects 
            
            self.coal_root = False #is it the root of the whole tree (?)
            
            self.subtree = True #a lil flag
            
            subtree_dict[self.person] = self
            
            #From here, adding in things from original get_subtrees function so will need to coordinate with above variables
            
            self.absolute_tip_times = []
            #self.lin_list = []
            #self.current_height = 0
            
            self.sort_out_tips()
            self.coalescent(self.tips, 0.0)
            self.define_root()
            self.get_branch_lengths()
            
            return subtree_dict
            
        
        else:
            #Whole tree coalescent characteristics
            
            self.final_nodes = set() #non-removed ones
            
            self.whole_tree = True
    
    
    def sort_out_tips(self):
        
        focal_individual = self.person
        
        self.find_transmission_tips(focal_individual)
        self.find_sample_tips(focal_individual)
    
        self.most_recent_tip = sorted([float(i) for i in self.absolute_tip_times])[::-1][0]
        self.root_time = self.most_recent_tip - float(focal_individual.time_infected)
        
        for tip in self.tips:
            #Might not need both of these
            tip.relative_height = self.most_recent_tip - tip.absolute_time
            self.heights[tip] = tip.relative_height #Not sure about this being a dictionary
    
    def find_transmission_tips(self, focal_individual):
                
        if len(focal_individual.sampled_infections) != 0:
            for case in focal_individual.infections:
                if len(case.sampled_infections) != 0 or case.sampled:
                    
                    transmission_tip = node(case, "Trans") 
                    transmission_tip.absolute_time = case.time_infected
                    
                    self.absolute_tip_times.append(transmission_tip.absolute_time)
                    #self.lin_list.append(transmission_tip)    
                    self.tips.append(transmission_tip)
    
    
    def find_sample_tips(self, focal_individual):
        
        if focal_individual.sampled: 
       
            self.absolute_tip_times.append(focal_individual.time_sampled)
            #self.lin_list.append(focal_individual)
            self.tips.append(focal_individual)
            
            self.sample_tip = focal_individual #Do we need this?
          
            
    def define_root(self):
        
        self.root = node(uuid.uuid1(), "Coal", height=self.root_time, children=self.penultimate, subtree=self)
        
        self.heights[self.root] = self.root_time
        
        self.root.last = True #?
        
        self.nodes.append(self.root)
        
        self.penultimate.node_parent = (self.root)
        
            
       
        
    def coalescent(self, lineage_list, current_height):
        def sort_key(ele):
            """Small function used later to sort lists by relative height"""
            return ele.relative_height
        
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
                next_sample = next(iter(sorted(to_be_sampled,key=sort_key)))
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
                        tau = self.root_time - current_height
                   
                    #Who is going to coalesce?
                    lucky_pair = random.sample(active_pop, k=2)
                    print("Coalescing " + str(lucky_pair))
                     
                    #ie the coalescent event of the pair selected above
                    parent_node = node(uuid.uuid1(), "Coal", height=current_height, children=lucky_pair, subtree=self)
                    
                    print("new node is " + str(parent_node.id) + " " + str(parent_node.relative_height))
                    
                    lucky_pair[0].node_parent = parent_node
                    lucky_pair[1].node_parent = parent_node
                    
                    self.nodes.append(parent_node)
                    self.heights[parent_node] = current_height
                    
                    updated_population = [i for i in lineage_list if i not in lucky_pair] + [parent_node]
                                
                    self.coalescent(updated_population, current_height)
    
    
    #Need to tidy this
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