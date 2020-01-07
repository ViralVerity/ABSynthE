#going to need:

import numpy as np
import random

class Individual(): 
        def __init__(self, unique_id, agent_location, cfr, inccdf, death_cdf, recovery_cdf): 

            self.unique_id = unique_id

            self.children = []
            self.exptimer = 0
            self.inftimer = 0

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
        
        
        def get_possible_cases(self):
            #Concern here is that now I'm storing the poss contact dict for every individual, which I don't need to do
            #Can I maybe remove it after I've used it? 
            #Tested - it def is storing it as an attribute of the individual
            
            self.poss_contact_dict = {}

            function = np.random.poisson

            lamb = np.random.gamma(0.37, 1.76) #lamb_m is 0.65

            a = 0.85
            b = a*0.5
            c = 0.07

            Hh_number = function(lamb)
            comm_number = function(a*lamb)
            dist_number = function(b*lamb)
            country_number = function(c*lamb)

            if Hh_number != None:
                self.poss_contact_dict["Hh"] = Hh_number
            else:
                self.poss_contact_dict["Hh"] = 0

            if comm_number != None:
                self.poss_contact_dict["Comm"] = comm_number
            else:
                self.poss_contact_dict["Comm"] = 0

            if dist_number != None:
                self.poss_contact_dict["Dist"] = dist_number
            else:
                self.poss_contact_dict["Dist"] = 0

            if country_number != None:
                self.poss_contact_dict["Country"] = country_number
            else:
                self.poss_contact_dict["Country"] = 0

            return self.poss_contact_dict  

