def get_LTT_SS(ltt_dict_input):
    
    ltt_dict = OrderedDict(ltt_dict_input)
    
    max_L = max(ltt_dict.values())
    
    t_max_L = max(ltt_dict, key = lambda k:ltt_dict[k])[0] #We'll get the start of that interval
    
    
    
    y1 = max_L
    y0 = next(iter(ltt_dict.values()))

    x1 = t_max_L
    x0 = 0.0
    
    y2 = next(reversed(ltt_dict.values()))

    x2 = next(reversed(ltt_dict.keys()))[1]


    slope_1 = abs((y2-y1)/(x2-x1))
    slope_2 = abs((y1-y0)/(x1-x0))
    
    slope_ratio = slope_1/slope_2
    
    return slope_1, slope_2, slope_ratio
    
    
    
    
    