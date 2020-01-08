import numpy as np

class node():
    
    def __init__(self, unique_id, node_type, trans_dict=None, those_sampled=None, node_dict=None, height=None, lucky_pair=None, subtree=None):
        
        self.id = unique_id
        self.sampled = False
        
        self.type = node_type #ie transmission, coalescent or individual - individual is a person.
        #NB before, node was each person I think, transmission_node and coalescent_node were separate classes
        #So really "node" before wasn't a node, it was a "node" on a transmission tree
        
        if self.type == "Ind": #ie if we're just looking at this person, not as a node in the tree
            self.get_useful_info(trans_dict, those_sampled, node_dict) #Gets info like time course of infection
        
        ##These three arguments are optional, because they're only relevant for coalescent nodes
        if height:
            self.height = height
        if lucky_pair:
            self.pair = lucky_pair
        if subtree:
            self.subtree = subtree
        
        self.children = set()
        self.sampled_children = set()
        self.node_children = set()
        
        self.to_root = []
        self.root_to_tip = 0.0

        self.transm_root = False
        self.last = False        
        self.removed = False
        
        #Test functions#
        self.remove_func_called = False
        self.for_loop_called = False
        
    def get_useful_info(self, trans_dict, those_sampled, node_dict):
        #Need to put the if statement that was in the original function to see if the person has already been initialised into the bigger code when we call __init__
        
        input1 = trans_dict[self.id][1] 
        input2 = trans_dict[self.id][2] 

        uniform1 = np.random.uniform(0,1)

        self.time_infected = float(trans_dict[self.id][1]) + uniform1

        if input1 == input2:
            rnge = 1 - uniform1
            self.time_sampled = float(self.time_infected) + np.random.uniform(0,rnge)

        else:
            self.time_sampled = float(trans_dict[self.id][2]) 
            
        node_dict[self.id] = self
        
        #self.get_parent(trans_dict, node_dict)
        #For now this is taken out, because it recur back to the root of the tree. May be helpful to combine functions, but could be less readable. But then it would get the path to root on init which would be nice
    
    
    def get_parent(self, trans_dict, those_sampled, node_dict):
        """Get parent as node object"""

        if trans_dict[self.id][0] in node_dict.keys():

            self.parent = node_dict[trans_dict[self.id][0]]
            self.parent.children.add(self)

        else: #if the parent is not yet in the node dict, make the parent. They shouldn't be made again because they'll now be in node_dict.keys()

            parent_id = trans_dict[self.id][0]
            
            if parent_id != "NA":
                print("making parent node")
                self.parent = node(parent_id, "Ind", trans_dict, those_sampled, node_dict)     
                node_dict[parent_id] = self.parent
                self.parent.children.add(self)


            
    def get_root_list(self, trans_dict, those_sampled, node_dict):
        
        #Caching to save time
        if len(self.to_root) != 0:
            return(self.to_root)

        else:
            self.get_parent(trans_dict, node_dict)

            #If the person is the root
            if self.parent.id == "NA": 
                path = []
                self.transm_root = True

            #Recur up the tree
            else:
                path = self.parent.get_root_list(trans_dict, those_sampled, node_dict) + [(self.parent)]
            #Not sure this will work here - may need to make this function a method of the tree or something

            self.to_root = path

            #Getting the sampling 
            if self.id in those_sampled:
                self.sampled = True
                for i in self.to_root:
                    i.sampled_children.add(self) #So the sampled children is all sampled children downstream of the focal individual

            return path
        
        
        
        
        
        
        
        
        
        
        
        

