import random
from multiprocessing.pool import ThreadPool
import time
from itertools import product


distribution = [0,1,2,3,4,5,6,7,8,9,10]


def test_multithread(input_lst):
    
    N = 10000
    
    while len(input_lst) < N:
        
        a = random.choice(distribution)
        
        if a % 2 == 0:
            
            input_lst.append(a)
            
            
            
    return input_lst
    
    
lst = []

start = time.time()

pool = ThreadPool(4)


end_lst = pool.map(test_multithread, (lst))   
 
end = time.time()
    
print(end-start)





lst2 = []

start = time.time()

test_multithread(lst2) 
 
end = time.time()
    
print(end-start)