
from multiprocessing.pool import ThreadPool
from collections import defaultdict

from epidemic_function_fitting import *
import index_functions_fitting

import Tree_simulator_fitting as cts

from make_contact_dicts_chiefdom import *
import distribution_functions

import LTT_metrics as LTT_metrics
from branch_length_parameters import *
from topology_set import *

def normalise(vector):
    norm=np.linalg.norm(vector, ord=1)
    return vector/norm

distributions = distribution_functions.define_distributions()

#dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/"
dropbox_path = "/localdisk/home/s1732989/ABM/Fitting/"

size_file = open(dropbox_path + "ABCSMC/epidemic_size.csv", 'w')

print("Defining contact structures")
contact_structure = make_contact_dicts(dropbox_path)

district_list = ["bo", 'bombali', 'bonthe', 'kailahun', 'kambia', 'kenema', 'koinadugu', 'kono', 'moyamba', 'portloko', 'pujehun', 'tonkolili', 'westernarearural', 'westernareaurban']

ch_list = ['badjia', 'bagbo', 'bagbwe', 'baoma', 'bumpe-gao', 'gbo', 'jaiamabongor', 'kakua', 'komboya', 'lugbu', 'niawalenga', 'selenga', 'tikonko', 'valunia', 'wonde', 'bendu', 'bum', 'dema', 'imperri', 'jong', 'kpanda', 'kwamebaikrim', 'nongobabullom', 'sittia', 'sogbini', 'yawbeko', 'biriwa', 'bombalisebora', 'gbantikamaranka', 'gbendemungowahun', 'libeisaygahun', 'magbaiambandowahun', 'makarigbanti', 'safrokolimba', 'pakimasabong', 'sandaloko', 'sandatenraren', 'sellalimba', 'tambakha', 'brimaia', 'gbinledixing', 'magbema', 'mambolo', 'masungbala', 'samu', 'tonkolimba', 'dea', 'jaluahun', 'jawei', 'kissikama', 'kissiteng', 'kissitongi', 'luawa', 'malema', 'mandu', 'pejebongre', 'pejewest', 'penguia', 'upperbambara', 'yawei', 'dama', 'dodo', 'guara', 'goramamende', 'kanduleppiama', 'koya', 'languramaya', 'lowerbambara', 'malegohun', 'niawa', 'nomo', 'nongowa', 'simbaru', 'smallbo', 'tunkia', 'wando', 'dembeliasikunia', 'diang', 'folsaba', 'kasunko', 'mongo', 'neya', 'nieni', 'sengbe', 'sulima', 'wara-warabafodea', 'wara-warayagala', 'fiama', 'gbanekandor', 'gbane', 'gbense', 'goramakono', 'kamara', 'lei', 'mafindor', 'nimikoro', 'nimiyama', 'sandor', 'soa', 'tankoro', 'toli', 'bagruwa', 'banta', 'bumpeh', 'dasse', 'fakunya', 'kargboro', 'kaiyamba', 'kamajei', 'kongbora', 'kori', 'kowa', 'ribbi', 'timidale', 'upperbanta', 'burehkasseh', 'buyaromende', 'debia', 'kaffubullom', 'lokomasama', 'maforki', 'marampa', 'masimera', 'koya', 'sandamagbolontor', 'tmsafroko', 'barri', 'gallinesperri', 'kpaka', 'kpangakabonde', 'makpele', 'malen', 'manosakrim', 'kpangakrim', 'peje', 'sorogbema', 'sowa', 'yakemokepukumukrim', 'gbonkolenken', 'kafesimiria', 'kalansongoia', 'kholifamabang', 'kholifarowalla', 'kunike', 'kunikebarina', 'malalmara', 'sambaia', 'tane', 'yoni', 'westernurban', 'westernrural']

dist_keys = []
ch_keys = []

for item1 in district_list:
    for item2 in district_list:
        if item1 != item2:
            dist_keys.append((item1, item2))

for item3 in ch_list:
    for item4 in ch_list:
        if item3 != item4:
            ch_keys.append((item3, item4))
            
            
LTT_bins = [(0, 0.021361168777556623), (0.021361168777556623, 0.042722337555113246), (0.042722337555113246, 0.06408350633266988), (0.06408350633266988, 0.08544467511022649), (0.08544467511022649, 0.10680584388778311), (0.10680584388778311, 0.12816701266533972), (0.12816701266533972, 0.14952818144289634), (0.14952818144289634, 0.17088935022045296), (0.17088935022045296, 0.19225051899800957), (0.19225051899800957, 0.2136116877755662), (0.2136116877755662, 0.2349728565531228), (0.2349728565531228, 0.25633402533067945), (0.25633402533067945, 0.27769519410823607), (0.27769519410823607, 0.2990563628857927), (0.2990563628857927, 0.3204175316633493), (0.3204175316633493, 0.3417787004409059), (0.3417787004409059, 0.36313986921846253), (0.36313986921846253, 0.38450103799601915), (0.38450103799601915, 0.40586220677357576), (0.40586220677357576, 0.4272233755511324)]
            
            

            
def simulate_epidemic(a, b, c, distributions=distributions, contact_structure=contact_structure, district_keys=dist_keys, ch_keys=ch_keys, size_file=size_file):

    original_dist_mvmt = defaultdict(list)
    original_ch_mvmt = defaultdict(list)

    for pair in district_keys:
        original_dist_mvmt[pair] = []
    
    for pair2 in ch_keys:
        original_ch_mvmt[pair2] = []

    original_case_dict = {}
    original_day_dict = defaultdict(list)
    option_dict_district_level = defaultdict(list)
    infected_individuals_set = set()
    cdf_array = []
    cdf_len_set = set()
    original_trans_dict = defaultdict(list)
    original_child_dict = defaultdict(list)
    original_nodes = []
    original_onset_times = []

    epidemic_length = 148

    for i in range(epidemic_length):
        original_day_dict[i] = []

    capped = True
    
    #Change this here instead of above
    LTT = False
    
    #b = 0.02
    #c = 0.1
    
    iteration_number_outside = 1
    
    result = run_model(a, b, c, LTT, iteration_number_outside, distributions, contact_structure, capped, size_file, original_dist_mvmt, original_ch_mvmt, original_case_dict, original_day_dict, option_dict_district_level, infected_individuals_set, cdf_array, cdf_len_set, original_trans_dict, original_child_dict, original_nodes, original_onset_times)

    return result

def run_model(a, b, c, LTT, iteration_number, distributions, contact_structure, capped, size_file, original_dist_mvmt, original_ch_mvmt, original_case_dict, original_day_dict, option_dict_district_level, infected_individuals_set, cdf_array, cdf_len_set, original_trans_dict, original_child_dict, original_nodes, original_onset_times):
    
    
    case_limit = 10000 #This may not be high enough, so we'll need to check if the epidemics are big enough
    popn_size = 7092142
    epidemic_length = 148
    cfr = 0.7
    sampling_percentage = 0.16
    

    iteration_count = -1
    
    for i in range(iteration_number):
    
    ##Setting things up for running###
        iteration_count += 1
        
        ###Making index case###
        index_case_case, index_case_individual, original_case_dict, original_trans_dict, original_child_dict, original_nodes, infected_individuals_set, original_day_dict = index_functions_fitting.make_index_case(contact_structure[0], cfr, distributions, original_case_dict, original_trans_dict, original_child_dict, original_nodes, infected_individuals_set, original_day_dict)
        
        
        susceptibles_left = True
        
        
        ###Run the epidemic###
        day_dict, case_dict, nodes, trans_dict, child_dict, dist_mvmt, ch_mvmt, onset_times, epidemic_capped = run_epidemic(0, original_day_dict, susceptibles_left , original_case_dict, original_trans_dict, original_child_dict, infected_individuals_set, popn_size, option_dict_district_level, original_onset_times, original_nodes, cdf_len_set, cdf_array, original_dist_mvmt, original_ch_mvmt, contact_structure, cfr, distributions, iteration_count, capped, epidemic_length, case_limit, a, b, c)

        
        remove_set = set()   
    
        ###Removing cases that don't exist eg because the person was already infected, or because the parent had recovered/died###
        for key, value in case_dict.items():
            if type(value) != Individual:
                remove_set.add(key)

        for item in remove_set:
            del case_dict[item]

        #Removing those people from the day dict
        for key, lst in day_dict.items():
            case_list = [item for item in lst if item not in remove_set] 
            day_dict[key] = case_list

        day_dict[0].append(index_case_case) #Put here so that it doesn't confuse the loop above because it has no parent AND otherwise it would get reassigned and stuff

        ###Getting results and writing to file###

        last_day = max(onset_times)
        
        #if len(case_dict) > 1800 and len(case_dict) < 2800: #Conditioning to speed up ABC
            
        size = len(case_dict)

        size_file.write(f"{a}, {b}, {c}, {size}\n")


        tree = cts.simulate_tree(trans_dict, child_dict, nodes, sampling_percentage, last_day, LTT)

        #This is where we choose the SS group
        #sim_LTT_pre = LTT_metrics.calculate_LTT_metrics(tree.lineages_through_time)
        
        if tree:
            #sim_LTT_points_pre = LTT_metrics.bin_sim(tree.lineages_through_time, LTT_bins)
            #sim_LTT = normalise(sim_LTT_points_pre)
            
            bl = calculate_branch_statistics(tree)
            
            #top = calculate_topology_params(tree)


            ch_jumps = 0
            for v in ch_mvmt.values():
                ch_jumps += len(v)
            dist_jumps = 0
            for v in dist_mvmt.values():
                dist_jumps += len(v)
         


            #return sim_LTT, ch_jumps, dist_jumps
            #return top, ch_jumps, dist_jumps
            return bl, ch_jumps, dist_jumps
            
        
        
        else: #If no sampled cases
            return 0,0,0
        #return tree



            


