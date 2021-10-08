#going to need:

import numpy as np
import random
import absynthe.set_up.distribution_functions

class Individual(): 
    def __init__(self, unique_id, agent_location, cfr, distributions): 
        """Defines infection course parameters for individual"""
        
        inc_cdf = distributions["inc_cdf"]
        death_cdf = distributions["death_cdf"]
        recovery_cdf = distributions["recovery_cdf"]
        
        self.unique_id = unique_id
        self.children = [] 
        
        self.hh = agent_location[self.unique_id][0]
        self.ch = agent_location[self.unique_id][1]
        self.dist = agent_location[self.unique_id][2]

        self.incubation_time(inc_cdf)
        self.death_prob(cfr)

        if self.death_state == True: 
            self.death_time(death_cdf)
            self.infectious_period = self.death_day + 7
        else:
            self.recovery_time(recovery_cdf)
            self.infectious_period = self.recovery_day

    def death_prob(self, cfr):
        death_poss = np.random.uniform(0, 1.0)
        if death_poss > cfr: #So they are still alive
            self.death_state = False
        else:
            self.death_state = True
        return self.death_state

    def incubation_time(self, incc_df):
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

        if comm_number != None:
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

        



















