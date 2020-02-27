#going to need:

import numpy as np
import random
import distribution_functions

class Individual(): 
    def __init__(self, unique_id, agent_location, cfr, distributions): 
        """Defines infection course parameters for individual"""
        inccdf = distributions[0]
        death_cdf = distributions[1]
        recovery_cdf = distributions[2]
        
        self.unique_id = unique_id

        self.children = [] 

        self.hh = agent_location[self.unique_id][0]
        self.comm = agent_location[self.unique_id][1]
        #self.ch = agent_location[self.unique_id][2]
        self.dist = agent_location[self.unique_id][2]

        self.incubation_time(inccdf)
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

    def incubation_time(self, inccdf):

        random_number = random.uniform(0,1)
        self.incubation_day = np.argmax(inccdf > random_number)

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

        lamb = np.random.gamma(0.37, 1.76) #lamb_m is 0.65

        #######FROM ABCSMC FITTING PROCESS######
        a = 1
        b = 0.05 
        c = 0.06
        ###############
        
        Hh_number = function(lamb)
        comm_number = function(a*lamb)
        dist_number = function(b*lamb)
        country_number = function(c*lamb)

        if Hh_number != None:
            poss_contact_dict["Hh"] = Hh_number
        else:
            poss_contact_dict["Hh"] = 0

        if comm_number != None:
            poss_contact_dict["Comm"] = comm_number
        else:
            poss_contact_dict["Comm"] = 0

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

            try:
                day_inf = day + current_day + self.incubation_day
                return day_inf, cdf_len_set
            except TypeError:
                #Raise exception here instead of return?
                return "Type", cdf_len_set, cdf_array
            
        except ValueError:
            return None, cdf_len_set, cdf_array

        



















