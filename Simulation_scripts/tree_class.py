from node_class import *

class tree():
    
    def __init__(self, person_tree=None, subtree_dict=None):
        
        if person_tree:
            self.person = person_tree #will be the person as a node object it corresponds to if it's a subtree
        
            #Subtree characteristics
            self.transmission_tips = [] #Subtree objects (that it's transmitting to?) This is quite odd, may be able to tidy away
            #self.sample_tip = None    Might need this not sure - might be a point later if I ask if there's a sample
            self.coalescent_nodes = [] #internal nodes as a list of coalescent_node objects 
            self.relative_heights = {}
            self.branch_lengths = {}
            
            self.coal_root = False #is it the root of the whole tree (?)
            
            self.subtree = True #a lil flag
            
            subtree_dict[self.person] = self
            
            #From here, adding in things from original get_subtrees function so will need to coordinate with above variables
            
            self.absolute_tip_times = []
            self.tips = []
            self.lin_list = []
            self.current_height = 0
            
            self.sort_out_tips()
            
        
        else:
            #Whole tree coalescent characteristics
            self.nodes = set()
            self.branch_lengths = {}
            self.node_heights = {}
            self.most_recent_date = 0.0
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
            self.relative_heights[tip] = tip.relative_height
    
    def find_transmission_tips(self, focal_individual):
                
        if len(focal_individual.sampled_children) != 0:
            for child in focal_individual.children:
                if len(child.sampled_children) != 0 or child.sampled:
                    
                    transmission_tip = node(child, "Trans") 
                    transmission_tip.absolute_time = child.time_infected
                    
                    self.absolute_tip_times.append(transmission_tip.absolute_time)
                    self.lin_list.append(transmission_tip)    
                    self.tips.append(transmission_tip)
    
    
    def find_sample_tips(self, focal_individual):
        
        if focal_individual.sampled: 
       
            self.absolute_tip_times.append(focal_individual.time_sampled)
            self.lin_list.append(focal_individual)
            self.tips.append(focal_individual)
            
            self.sample_tip = focal_individual #Do we need this?
            
            
            
            
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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