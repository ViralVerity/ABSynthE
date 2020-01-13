import node_class as nc
import random
import uuid
import numpy as np

class tree():   
    
    def __init__(self, person_tree=None, subtree_dict=None, node_dict=None):
        
        #print("Making subtree for " + person_tree.id)
                
        self.nodes = []
        self.tips = []
        self.branch_lengths = {}
        self.heights = {}
        
        if person_tree:
            self.subtree = True #a lil flag
            self.person = person_tree #will be the person as a node object it corresponds to if it's a subtree
            self.contains_sample = False      
            self.transmits = False
            
            self.absolute_tip_times = []            
                  
            self.coal_root = False #is it the root of the whole tree (?) Maybe just have tree.root as a subtree object
            
            self.sort_out_tips()
            
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
        
        #if self.person.id == "1414058":
         #   print("root time = " + str(self.root_time) + " for " + str(self.person.id))
                
        self.nodes.append(self.root)
        
        self.penultimate.node_parent = (self.root)
        
        self.branch_lengths[self.root] = 0.0


        
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
                    
                    self.nodes.append(parent_node)
                    self.heights[parent_node] = current_height
                    
                    updated_population = [i for i in lineage_list if i not in lucky_pair] + [parent_node]
                                
                    self.coalescent(updated_population, current_height)
    
    def get_branch_lengths(self): 
        
        for tip in self.tips:
            self.branch_lengths[tip] = self.heights[tip.node_parent] - self.heights[tip]

        for nde in self.nodes:
            if nde != self.root:
                self.branch_lengths[nde] = self.heights[nde.node_parent] - self.heights[nde]


    def update_coalescent_tree(self, subtree):
        """Update the coalescent tree with information from the subtree"""

        self.branch_lengths.update(subtree.branch_lengths)

        for tip in subtree.tips:
            if tip.type == "Ind" and subtree.contains_sample:
                self.tips.append(tip)
            else:
                self.nodes.append(tip)

        for nde in subtree.nodes:
            self.nodes.append(nde)

    
    def construct_tree(self, node_dict):
        
        for nde in node_dict.values():
            
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
                        
                        donor_tree = tip.infector.subtree
                        recipient_tree = tip.infectee.subtree
                        
                        recipient_tree.root.node_parent = tip
                        tip.node_children.append(recipient_tree.root) 
            
                    
            self.update_coalescent_tree(subtree)
            
            
                
                
    def remove_internals(self, nde):
        

        nde.new_children = nde.node_children.copy()

        if len(nde.node_children) == 1: 
            if nde != self.root:
                
                parent = nde.node_parent 

                nde.removed = True
                parent.new_children.remove(nde)
                
                self.final_nodes.remove(nde)

                #Reassign parents and children to remove internal nodes
                for child in nde.node_children:

                    child.node_parent = parent

                    parent.new_children.append(child)

                    self.branch_lengths[child] = self.branch_lengths[child] + self.branch_lengths[nde]


            else:

                for child in nde.new_children:

                    child.node_parent = None
                    self.branch_lengths[child] = 0.0
                    self.root = child

                nde.removed = True
                
                self.final_nodes.remove(nde)


            
        for i in nde.node_children:
            self.remove_internals(i)
                
        
        return self
    
    

        
        
        
                
                
                
                
                
                
                
                
                
                