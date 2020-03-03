from individual_class import *

class Case():
    
    def __init__(self, case_id, level, parent):

        self.children = []

        self.case_id = case_id

        self.level = level
        
        self.parent = parent

        #print("case ID " + str(self.case_id) + " is level " + str(self.level))

        if self.case_id != 0 and self.level == None:
            print("ERROR " + str(self.case_id))
            
            
    def get_options_district(self, option_dict_districtlevel, dist_to_ch, cluster_to_ppl, parent_individual):

        if parent_individual.comm not in option_dict_districtlevel.keys():

            dist_ppl_list = ([cluster_to_ppl[clust] for clust in dist_to_ch[parent_individual.dist] if clust != parent_individual.comm])

            option_dict_districtlevel[parent_individual.comm] = dist_ppl_list

        else:

            dist_ppl_list = option_dict_districtlevel[parent_individual.comm]
    
        
        return dist_ppl_list, option_dict_districtlevel


    
    def who_am_I(self, infected_individuals_set, popn_size, option_dict_districtlevel, contact_structure, case_dict, parent_individual, day, cfr, distributions, current_day, SDB_start, SDB_success): 
        """Input is case object that has already been initialised.
        Finds out which of the parent's potential contacts are still susceptible"""
        
        agent_location, dist_to_hh, hh_to_cluster, cluster_to_hh, hh_to_ppl, cluster_to_ppl, dist_to_ppl, dist_to_ch, district_distance,  district_pops = contact_structure
        
        if len(infected_individuals_set) == popn_size: 
            return False

        if self.level == "Hh":
            poss_case = random.choice([person for person in hh_to_ppl[parent_individual.hh] if person != parent_individual.unique_id]) 

        elif self.level == "Comm":
            
            poss_case = random.choice(random.choice([hh_to_ppl[hh] for hh in cluster_to_hh[parent_individual.comm] 
                                                if hh != parent_individual.hh]))
          
        elif self.level == "Dist":
            
            if parent_individual.dist == "westernarearural" or parent_individual.dist == "westernareaurban":
                return
            
            else:
            
                dist_ppl_list = self.get_options_district(option_dict_districtlevel, dist_to_ch, cluster_to_ppl, parent_individual)[0]

                poss_case = random.choice(random.choice(dist_ppl_list))

        elif self.level == "Country":        

            if day == 0:
                district = "kenema"
            else:
                district = random.choice(district_distance[parent_individual.dist])

            poss_case = random.choice(dist_to_ppl[district])

        else:
            print("ERROR no level assigned")


        #Is the person actually susceptible
        if poss_case not in infected_individuals_set:
            new_individual = Individual(poss_case, agent_location, cfr, distributions, day_arg=current_day, SDB_start_arg=SDB_start, SDB_success_arg=SDB_success)
            case_dict[self] = new_individual
            infected_individuals_set.add(poss_case)
            
        elif day == 0: #So that there are actually 14 cases in the first transmission cluster
            self.who_am_I(infected_individuals_set, popn_size, option_dict_districtlevel, contact_structure, case_dict, parent_individual, day, cfr, distributions)

        else:
            #print("Already infected")
            return

        #print("Time taken to define self = " + str(end-start))
        return poss_case, case_dict, infected_individuals_set
    
    
    
    
    
    
    
    