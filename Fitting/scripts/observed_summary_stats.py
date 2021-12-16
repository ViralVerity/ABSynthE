    
import numpy as np
    
def import_observed_stats(summary_stats_set):
    
    tips = 214
    
    dist = 100
    ch = 123
    # neighbouring = XX
    
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
    
    
    branch_set = [mean_lengths, median_lengths, var_lengths, mean_external, median_external, var_external, mean_internal, median_internal, var_internal, mean_ratio, median_ratio, var_ratio]
    
    if summary_stats_set != "all":
        branch_set.append(tips)
        
    #Topology set#

    colless = 2168
    sackin = 3204 #includes root
    WD_ratio = 0.75
    delta_W = 12
    max_ladder = 0.014018691588785047
    IL_nodes = 0.4507042253521127
    staircase_1 = 0.704225352112676
    staircase_2 = 0.5386550674701448
    
    topology_set = [colless, sackin, WD_ratio, delta_W, max_ladder, IL_nodes, staircase_1, staircase_2]
    
    if summary_stats_set != "all":
        topology_set.append(tips)

    #LTT set# backwards in time

    max_L = 95
    t_max_L = 0.21369863013683843
    slope_1 = 435.19230769260884
    slope_2 = 539.8534746773089
    slope_ratio = 0.8061304189117982
    mean_s_time = 0.0012605312238711675
    mean_b_time = 0.0016789798767443802

    ltt_set = [max_L, t_max_L, slope_1, slope_2, slope_ratio, mean_s_time, mean_b_time]
    if summary_stats_set != "all":
        ltt_set.append(tips)

    # LTT points # x first and then y
    ltt_points = [0, 0.019298379818496268, 0.038596759636992536, 0.057895139455488805, 0.07719351927398507, 0.09649189909248135, 0.11579027891097762, 0.1350886587294739, 0.15438703854797017, 0.17368541836646645, 0.19298379818496272, 0.212282178003459, 0.23158055782195527, 0.25087893764045155, 0.2701773174589478, 0.28947569727744404, 0.3087740770959403, 0.32807245691443654, 0.3473708367329328, 0.36666921655142903, 1.4391462194824414, 3.581938515240152, 4.420100920079058, 2.4914173412345284, 3.905970583926359, 2.7226772671032586, 3.3017143646432223, 3.3979011154953764, 5.0119232333320545, 4.8993542269646975, 3.784612056539124, 3.130287797569775, 2.4932891885501345, 0.9844615336004773, 0.951677090642562, 1.063827615278881, 1.2894395928177025, 1.0293869592741283] 
    if summary_stats_set != "all":
        ltt_points.append(tips)
    
    return branch_set, topology_set, ltt_set, ltt_points, dist, ch, tips

    # return branch_set, topology_set, ltt_set, ltt_points, dist, ch, neighbouring, tips

def get_observed_SS(summary_stats_set):
    
    branch_set, topology_set, ltt_set, ltt_points, observed_dist, observed_ch, tips = import_observed_stats(summary_stats_set)
    # branch_set, topology_set, ltt_set, ltt_points, observed_dist, observed_ch, neighbouring, tips = import_observed_stats(summary_stats_set)


    if summary_stats_set != "all":
        obs_bl = normalise(branch_set)
        obs_top = normalise(topology_set)
        obs_ltt = normalise(ltt_set)
        obs_ltt_points = normalise(ltt_points)
        output = obs_bl, obs_top, obs_ltt, obs_ltt_points
    else:
        sum_stats = []
        sum_stats.extend(branch_set)
        sum_stats.extend(topology_set)
        sum_stats.extend(ltt_set)
        sum_stats.extend(ltt_points)
        sum_stats.append(tips)

        print(sum_stats)

        sum_stats = normalise(sum_stats)
        output = sum_stats
    
    return output, observed_ch, observed_dist
    # return output, observed_ch, observed_dist, neighbouring

def normalise(vector):
    norm=np.linalg.norm(vector, ord=1)
    return vector/norm
