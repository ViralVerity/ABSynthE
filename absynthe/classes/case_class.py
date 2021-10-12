import sys
import random

class Case():
    
    def __init__(self, case_id, level, parent):

        self.children = []
        self.case_id = case_id
        self.level = level
        self.parent = parent
            
    def get_options_district(self, option_dict_district_level, config, parent_individual):

        if parent_individual.ch in option_dict_district_level.keys():
            dist_ppl_list = option_dict_district_level[parent_individual.ch]
        elif parent_individual.dist == "westernarearural" or parent_individual.dist == "westernareaurban":
            #not sure if this is the right behaviour
            dist_ppl_list = ([config["population_structure"]["hh_to_ppl"][hh] for hh in config["population_structure"]["ch_to_hh"][parent_individual.ch] 
                                                if hh != parent_individual.hh])
        else:
            dist_ppl_list = ([config["population_structure"]["ch_to_ppl"][chief] for chief in config["population_structure"]["dist_to_ch"][parent_individual.dist] if chief != parent_individual.ch])
            option_dict_district_level[parent_individual.ch] = dist_ppl_list

        return dist_ppl_list, option_dict_district_level


    def who_am_I(self, parent_individual, day, config, epidemic_config): 
        """Input is case object that has already been initialised.
        Finds out which of the parent's potential contacts are still susceptible"""
                
        if len(epidemic_config["infected_individuals_set"]) == config["population_structure"]["popn_size"]: 
            return False

        if self.level == "Hh": #within household
            poss_case = random.choice([person for person in config["population_structure"]["hh_to_ppl"][parent_individual.hh] if person != parent_individual.unique_id]) 

        elif self.level == "Ch": #within chiefdom, between household
            poss_case = random.choice(random.choice([config["population_structure"]["hh_to_ppl"][hh] for hh in config["population_structure"]["ch_to_hh"][parent_individual.ch] 
                                                if hh != parent_individual.hh]))
        
        elif self.level == "Dist": #within district, between chiefdom
                #NB if person in WAR and WAU, they're currently drawn from anywhere in that district for this level
                dist_ppl_list, option_dict_district_level = self.get_options_district(epidemic_config["option_dict_district_level"], config, parent_individual)
                poss_case = random.choice(random.choice(dist_ppl_list))
                epidemic_config['option_dict_district_level'] = option_dict_district_level

        elif self.level == "Country":        
            if day == 0:
                district = "kenema"
            else:
                district = random.choice(config["population_structure"]["district_distance"][parent_individual.dist])
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
    
    
    
    
    
    
    
    