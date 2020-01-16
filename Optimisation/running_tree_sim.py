


import Tree_simulator_for_op as cts
from collections import defaultdict

##(trans_dict, nodes, sampling_proportion, epidemic_len)


dropbox_path = "/Users/s1743989/VirusEvolution Dropbox/Verity Hill/Agent_based_model/Looping models/"
results = "Results/no_caps/2/log_files/"




transm_dict = defaultdict(list)
nodes = []
times = []
child_dict = defaultdict(list)


with open(dropbox_path + results + "information_file_for_38.csv") as f:
    next(f)
    for l in f:
        toks = l.strip("\n").split(",")
        
        focal = toks[0]
    
        transm_dict[focal] = [toks[1], int(toks[4]), int(toks[5])]
        
        child_dict[toks[1]].append(focal)
        
        nodes.append(toks[0])
        
        times.append(int(toks[5]))
        
        nodes.append(focal)
            
        if focal not in child_dict.keys():
               child_dict[focal] = []
        
        
epidemic_len = max(times)


newick_tree, Ne_dict, big_tree, R0, big_tree.most_recent_date, those_sampled, coal_intervals = cts.simulate_tree(transm_dict, child_dict, nodes, 1, epidemic_len )


