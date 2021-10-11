#going to need:

import numpy as np
import random
import absynthe.set_up.distribution_functions

class Individual(): 
    def __init__(self, unique_id, parent, agent_location, cfr, distributions, day, epidemic_config): 
        """Defines infection course parameters for individual"""
        
        inc_cdf = distributions["inc_cdf"]
        death_cdf = distributions["death_cdf"]
        recovery_cdf = distributions["recovery_cdf"]
        
        self.unique_id = unique_id
        self.parent = parent
        self.children = [] 

        self.hh = agent_location[self.unique_id][0]
        self.ch = agent_location[self.unique_id][1]
        self.dist = agent_location[self.unique_id][2]

        self.infection_day = day
        self.incubation_time(inc_cdf)
        self.death_prob(cfr)

        if self.death_state == True: 
            self.death_time(death_cdf)
            self.infectious_period = self.death_day + 7 #ebola specific
        else:
            self.recovery_time(recovery_cdf)
            self.infectious_period = self.recovery_day

        self.add_self_to_dicts(epidemic_config)


    def death_prob(self, cfr):
        death_poss = np.random.uniform(0, 1.0)
        if death_poss > cfr: #So they are still alive
            self.death_state = False
        else:
            self.death_state = True
        return self.death_state

    def incubation_time(self, inc_cdf):
        random_number = random.uniform(0,1)
        self.incubation_day = np.argmax(inc_cdf > random_number)
        return self.incubation_day

    def death_time(self, death_cdf):
        random_number = random.uniform(0,1)
        self.death_day = np.argmax(death_cdf > random_number)
        return self.death_day

    def recovery_time(self, recovery_cdf):
        #Can't recover before day 4 - taken from the NEJM paper figure
        random_number = random.uniform(recovery_cdf[3],1)
        self.recovery_day = np.argmax(recovery_cdf > random_number)
        return self.recovery_day


    def add_self_to_dicts(self, epidemic_config):

        epidemic_config["child_dict"][self.unique_id] = []
        if self.parent:
            epidemic_config["child_dict"][self.parent.unique_id].append(self.unique_id)
            epidemic_config["transmission_dict"][self.unique_id] = {"parent":self.parent.unique_id, "day_infected":self.infection_day, "day_sampled":(self.infection_day+self.incubation_day)} #sample day currently same as symptom onset
        else:
            if self.unique_id != epidemic_config["index_id"]:
                sys.stderr.write(f'{self.unique_id} has no parent and is not the index case.\n')
                sys.error(-1)
            epidemic_config["child_dict"]["NA"] = [self.unique_id]
            epidemic_config["transmission_dict"][self.unique_id] = {"parent":None, "day_infected":self.infection_day, "day_sampled":(self.infection_day+self.incubation_day)}

        epidemic_config["nodes"].append(self.unique_id)
        epidemic_config["districts_present"].add(self.dist)
        epidemic_config["chiefdoms_present"].add(self.ch)

        epidemic_config["onset_times"].append(self.infection_day)
        epidemic_config["infected_individuals_set"].add(self.unique_id)


    ###Running epidemic functions###
    def get_possible_cases(self):
        """Get the number of exposed secondary cases at each contact level"""
        #Not storing each individuals possible case dict due to memory concerns
        #If you want to get this, add self. in front of each poss_contact_dict mention

        poss_contact_dict = {}
        function = np.random.poisson
        lamb = np.random.gamma(0.37, 1.76) #lamb_m is 0.65 - that is from a paper somewhere, should be in the docs

        #######FROM ABCSMC FITTING PROCESS######
        a = 0.65
        b = 0.11
        c = 0.32
        ###############
        
        Hh_number = function(lamb)
        ch_number = function(a*lamb)
        dist_number = function(b*lamb)
        country_number = function(c*lamb)

        if Hh_number != None:
            poss_contact_dict["Hh"] = Hh_number
        else:
            poss_contact_dict["Hh"] = 0

        if ch_number != None:
            poss_contact_dict["Ch"] = ch_number
        else:
            poss_contact_dict["Ch"] = 0

        if dist_number != None:
            poss_contact_dict["Dist"] = dist_number
        else:
            poss_contact_dict["Dist"] = 0

        if country_number != None:
            poss_contact_dict["Country"] = country_number
        else:
            poss_contact_dict["Country"] = 0

        return poss_contact_dict  

    
    def when_infected(self, current_day, possible_case, cdf_len_set, cdf_array):
        """Gets when secondary cases are infected, if that happens before the end of self's infectious period"""
        if self.infectious_period not in cdf_len_set:
            cdf = distribution_functions.get_cdf(self.infectious_period)
            cdf_array.append(cdf)
            cdf_len_set.add(self.infectious_period)
        
        else:
            for item in cdf_array:
                if self.infectious_period == len(item):
                    cdf = item

        random_number = random.uniform(0,1)

        try:
            day = np.argmax(cdf > random_number)
            day_inf = day + current_day + self.incubation_day
            return day_inf, cdf_len_set
            
        except ValueError: #so this comes up if they finish their infection before the exposure happens - should have an if statement really
            return None, cdf_len_set, cdf_array

        



















