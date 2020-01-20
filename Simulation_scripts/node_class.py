import numpy as np
from tree_class import *

class node():
    
    def __init__(self, unique_id, node_type, trans_dict=None, child_dict=None, those_sampled=None, node_dict=None, height=None, children=None, subtree=None, infector=None, infectee=None):
        
        self.id = unique_id #Need to think about the IDs for the three types
        
        self.type = node_type #ie transmission, coalescent or individual - individual is a person.
        
        self.branch_len_calculated = False
        
        self.remove_func_called = False
        
        if self.type == "Ind": #ie if we're just looking at this person, not as a node in the tree
            self.infections = set() 
            self.sampled_infections = set()
                                    
            if self.id in those_sampled:
                self.sampled = True
            else:
                self.sampled = False
            
            parent_id = trans_dict[self.id][0]
            node_dict[self.id] = self
            
            if parent_id in node_dict.keys():
                self.parent = node_dict[parent_id]
            else:
                self.parent = None
                
            if parent_id == "NA":
                self.generation = 0
                self.index_case = True
            else:
                self.generation = self.parent.generation + 1
                self.index_case = False
                
                
            self.get_useful_info(trans_dict, those_sampled, node_dict) #Gets info like time course of infection
            self.find_children(trans_dict, child_dict, those_sampled, node_dict)
            
            #If the lineage is sampled
            if self.sampled or len(self.sampled_infections) != 0:
                
                self.root_to_tip = self.time_sampled
                
                self.subtree = tree(self, {})

        
        if self.type == "Trans":
            self.infector = infector
            self.infectee = infectee
        

        ##These three arguments are optional, because they're only relevant for coalescent nodes
        if height:
            self.relative_height = height
        if children:
            self.node_children = children
        elif not children:
            self.node_children = []
        if subtree:
            self.subtree = subtree
            

        self.removed = False
        
        #Test functions#
        self.remove_func_called = False
        self.for_loop_called = False
        
        if self.type == "Coal":
            self.root_to_tip = 0.0
        
        
    def get_useful_info(self, trans_dict, those_sampled, node_dict):

        input1 = trans_dict[self.id][1] 
        input2 = trans_dict[self.id][2] 

        uniform1 = np.random.uniform(0,1)

        if self.index_case:
            self.time_infected = 0.0
        else:
            self.time_infected = float(trans_dict[self.id][1]) + uniform1

        if input1 == input2:
            rnge = 1 - uniform1
            self.time_sampled = float(self.time_infected) + np.random.uniform(0,rnge)

        else:
            self.time_sampled = float(trans_dict[self.id][2]) + uniform1
            
        
        self.absolute_time = self.time_sampled #For use in the tree
    
    
    def find_children(self, trans_dict, child_dict, those_sampled, node_dict):
        
        secondary_infections = child_dict[self.id] #immediate transmission
        
        for child in secondary_infections:

            new_child = node(child, "Ind", trans_dict, child_dict, those_sampled, node_dict)

            self.infections.add(new_child)
          
            new_child.parent = self
            new_child.generation = self.generation + 1
                        
            if new_child.sampled:
                
                self.sampled_infections.add(new_child)

            self.sampled_infections = self.sampled_infections.union(new_child.sampled_infections)

            
    
    
    
        
        
        
        
        
        
        
        
        

