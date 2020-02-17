def get_jumps(observed_mvmt, mvmt_dict):
    
    total_jumps = 0
    for v in mvmt_dict.values():

        total_jumps += len(v)
        
        
    difference = abs(observed_mvmt - total_jumps)
    
    percentage_difference = difference/observed_mvmt
        
    return percentage_difference


