    
def import_observed_stats():
    
    tips = 214
    
    dist = 100
    ch = 123
    
    #Branch lengths#
    max_H = 0.3859675963699254
    min_H = 0.11747444568536669

    mean_lengths = 0.035080050967117304
    median_lengths = 0.021939933283573243
    var_lengths = 0.0013770629364285776

    mean_external = 0.051805151396388485
    median_external = 0.04441308703686353
    var_external = 0.001384042360694057

    mean_internal = 0.018197166571532235
    median_internal = 0.0072969106972364806
    var_internal = 0.0008026179459713262

    mean_ratio = 0.35126171975246506
    median_ratio = 0.16429640865044878
    var_ratio = 0.5799085120262046
    
    branch_set = [tips, mean_lengths, median_lengths, var_lengths, mean_external, median_external, var_external, mean_internal, median_internal, var_internal, mean_ratio, median_ratio, var_ratio]

    #Topology set#

    colless = 2168
    sackin = 3204 #includes root
    WD_ratio = 0.75
    delta_W = 12
    max_ladder = 0.014018691588785047
    IL_nodes = 0.4507042253521127
    staircase_1 = 0.704225352112676
    staircase_2 = 0.5386550674701448
    
    topology_set = [tips, colless, sackin, WD_ratio, delta_W, max_ladder, IL_nodes, staircase_1, staircase_2]

    #LTT set#

    max_L = 95
    t_max_L = 0.16645019089776691
    slope_1 = 558.7257034575601
    slope_2 = 423.6566107364787
    slope_ratio = 1.318817384877529
    mean_s_time = 0.0012605312238711675
    mean_b_time = 0.0016789798767443802
    

    # LTT points # x first and then y
    ltt_points = [tips, 0, 0.019298379818496268, 0.038596759636992536, 0.057895139455488805, 0.07719351927398507, 0.09649189909248135, 0.11579027891097762, 0.1350886587294739, 0.15438703854797017, 0.17368541836646645, 0.19298379818496272, 0.212282178003459, 0.23158055782195527, 0.25087893764045155, 0.2701773174589478, 0.28947569727744404, 0.3087740770959403, 0.32807245691443654, 0.3473708367329328, 0.36666921655142903, 2, 2, 2, 3, 6, 14, 34, 74, 84, 90, 56, 55, 62, 54, 40, 41, 50, 54, 34, 25] 
    
    ltt_set = [tips, max_L, t_max_L, slope_1, slope_2, slope_ratio, mean_s_time, mean_b_time]
    
    return branch_set, topology_set, ltt_set, ltt_points, dist, ch