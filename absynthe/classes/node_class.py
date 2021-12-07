import numpy as np
from absynthe.classes.tree_class import *

class node():
    
    def __init__(self, unique_id, node_type, transmission_dict=None, child_dict=None, those_sampled=None, node_dict=None, height=None, children=None, subtree=None, infector=None, infectee=None):
        
        self.id = unique_id #Need to think about the IDs for the three types
        self.type = node_type #ie transmission, coalescent or individual - meaning for a transmission tree, a coalescent tree, or just as a person which is a tip if it's sampled
        
        self.removed = False
        
        if self.type == "individual": #if we're looking at the person - can be a tip or can be internal in the tree
            self.infections = set() 
            self.sampled_infections = set()
                                    
            if self.id in those_sampled:
                self.sampled = True
            else:
                self.sampled = False
            
            parent_id = transmission_dict[self.id]["parent"]
            node_dict[self.id] = self

            if parent_id:
                self.index_case = False
                if parent_id in node_dict.keys():
                    self.parent = node_dict[parent_id]
                    self.generation = self.parent.generation + 1
                else:
                    self.parent = None #if it's the first one of the siblings to go through      
            else:
                self.generation = 0
                self.index_case = True

            self.find_time_of_day_infected_sampled(transmission_dict, those_sampled, node_dict) 
            self.find_children(transmission_dict, child_dict, those_sampled, node_dict)
            
            #If the lineage is sampled or it's child is
            if self.sampled or len(self.sampled_infections) != 0:
                self.root_to_tip = self.time_sampled
                self.subtree = tree(tree_type = "subtree", focal_person=self)

        
        if self.type == "transmission":
            self.infector = infector
            self.infectee = infectee
        

        if self.type == "coalescent":
            self.relative_height = height
            self.subtree = subtree
            self.root_to_tip = 0.0 #because it's the root of the subtree I think

            if children:
                self.node_children = children
            elif not children:
                self.node_children = []
        
        #Test functions#
        self.remove_func_called = False
        self.for_loop_called = False
        self.branch_len_calculated = False
        self.remove_func_called = False
        
        return node_dict #?
        
        
    def find_time_of_day_infected_sampled(self, transmission_dict, those_sampled, node_dict):
        #this is to get the time of day someone was infected/sampled
        #it's in terms of years
        
        day_infected = transmission_dict[self.id]["day_infected"] 
        day_sampled = transmission_dict[self.id]["day_sampled"]

        uniform1 = np.random.uniform(0,1)

        if self.index_case:
            self.time_infected = 0.0
        else:
            self.time_infected = (day_infected + uniform1)/365

        if day_infected == day_sampled:
            rnge = time_infected - (day_sampled/365)
            self.time_sampled = (self.time_infected + np.random.uniform(0,rnge))
        else:
            self.time_sampled = (day_sampled + uniform1)/365
            
        self.absolute_time = self.time_sampled #For use in the tree
    
    
    def find_children(self, transmission_dict, child_dict, those_sampled, node_dict):
        
        self.secondary_infections = child_dict[self.id] #immediate transmission
        
        for child in secondary_infections:

            new_child = node(child, "individual", transmission_dict, child_dict, those_sampled, node_dict)

            self.infections.add(new_child)

            if not new_child.parent:
                new_child.parent = self
                new_child.generation = self.generation + 1
                        
            if new_child.sampled:
                self.sampled_infections.add(new_child)

            #to give all downstream sampled infections
            self.sampled_infections = self.sampled_infections.union(new_child.sampled_infections)

