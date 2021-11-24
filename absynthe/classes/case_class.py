import sys
import random

class Case():
    
    def __init__(self, case_id, level, parent):

        self.children = []
        self.case_id = case_id
        self.level = level
        self.parent = parent
        
    def who_am_I(self, parent_individual, day, config, epidemic_config): #can add in day 0 specific stuff to each level if need be
        """Input is case object that has already been initialised.
        Finds out which of the parent's potential contacts are still susceptible"""

        if day == 0:
            if self.level == "Hh":
                poss_case = random.choice([person for person in config["population_structure"]["hh_to_ppl"][parent_individual.hh] if person != parent_individual.unique_id]) 
            elif self.level == "Ch":
                poss_case = random.choice([person for person in config["population_structure"]["ch_to_ppl"][parent_individual.ch] if person not in parent_individual.household_list]) 
            elif self.level == "Dist":
                poss_case = random.choice(person for person in config["population_structure"]["ch_to_ppl"]["Kissi_Tongi"]) #because we know which chiefdom those cases are in

        else:
            if len(epidemic_config["infected_individuals_set"]) == config["population_structure"]["popn_size"]: 
                return False

            if self.level == "Hh": #within household
                poss_case = random.choice([person for person in config["population_structure"]["hh_to_ppl"][parent_individual.hh] if person != parent_individual.unique_id]) 

            elif self.level == "Ch": #within chiefdom
                poss_case = random.choice([person for person in config["population_structure"]["ch_to_ppl"][parent_individual.ch] if person != parent_individual.unique_id]) 
            
            elif self.level == "Dist": #within district                
                poss_case = random.choice([person for person in config["population_structure"]["dist_to_ppl"][parent_individual.dist] if person != parent_individual.unique_id]) 

            elif self.level == "Country":        
                district = random.choice(config["population_structure"]["district_distance"][parent_individual.dist]) #not sure how to change this now that it's within
                poss_case = random.choice(config["population_structure"]["dist_to_ppl"][district])

            else:
                sys.stderr.write("ERROR no level assigned")

        #Is the person actually susceptible
        if poss_case not in epidemic_config["infected_individuals_set"]:
            return poss_case, config, epidemic_config
            
        elif day == 0: #So that there are actually 14 cases in the first transmission cluster
            self.who_am_I(parent_individual, day, config, epidemic_config)

        else:
            return 0
    
    
    
    
    
    
    
    