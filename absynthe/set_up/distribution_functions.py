import scipy as sp
import numpy as np
from scipy import stats

#these are all ebov specific too - maybe need a directory with all ebov specific files and flags to point to them

def define_distributions():
    
    clinical_x = np.linspace(0, 40, 40)

    inc_shape = (8.5/7.6)**2
    inc_scale = (7.6**2)/8.5
    inc_cdf = sp.stats.gamma.cdf(clinical_x, inc_shape, loc = 0, scale = inc_scale)

    death_shape = (8.6/6.9)**2
    death_scale = (6.9**2)/8.6
    death_cdf = sp.stats.gamma.cdf(clinical_x, death_shape, loc = 0, scale = death_scale)

    recovery_shape = (15.2/6.2)**2
    recovery_scale = (6.2**2)/15.2
    recovery_cdf = sp.stats.gamma.cdf(clinical_x, recovery_shape, loc = 0, scale = recovery_scale)
    
    distribution_dict = {}
    distribution_dict["inc_cdf"] = inc_cdf
    distribution_dict["death_cdf"] = death_cdf
    distribution_dict["recovery_cdf"] = recovery_cdf

    return distribution_dict

def get_cdf(dim):

    x = np.linspace(0,dim, dim) #This is where difference between living/dead comes in
    mu = 3.1 
    sigma = 2.5

    a = shape = (mu/sigma)**2
    scale = (sigma**2)/mu

    cdf = sp.stats.gamma.cdf(x,a,loc=0,scale=scale)

    return cdf

