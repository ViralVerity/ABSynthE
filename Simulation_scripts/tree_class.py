class tree():
    
    def __init__(self, person_tree=None):
        
        if person_tree:
            self.id = person_tree #will be the person it corresponds to if it's a subtree
            self.person = person_tree
        
            #Subtree characteristics
            self.transmission_tips = [] #Subtree objects
            self.sample_tip = None    
            self.most_recent_tip = 0 #time wise when is the most recent, useful for the relative heights 
            self.coalescent_nodes = [] #internal nodes as a list of coalescent_node objects 
            self.root_time = 0 
            self.relative_heights = {}
            self.branch_lengths = {}
            self.coal_root = False
            
            self.subtree = True #a lil flag
        
        else:
            #Whole tree coalescent characteristics
            self.nodes = set()
            self.branch_lengths = {}
            self.node_heights = {}
            self.most_recent_date = 0.0
            self.final_nodes = set() #non-removed ones
            
            self.whole_tree = True
    
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