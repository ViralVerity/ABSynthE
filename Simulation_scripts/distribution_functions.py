import scipy as sp
import numpy as np

def define_distributions():
    
    clinical_x = np.linspace(0, 40, 40)

    incshape = (8.5/7.6)**2
    incscale = (7.6**2)/8.5
    inccdf = sp.stats.gamma.cdf(clinical_x, incshape, loc = 0, scale = incscale)

    death_shape = (8.6/6.9)**2
    death_scale = (6.9**2)/8.6
    death_cdf = sp.stats.gamma.cdf(clinical_x, death_shape, loc = 0, scale = death_scale)

    recovery_shape = (15.2/6.2)**2
    recovery_scale = (6.2**2)/15.2
    recovery_cdf = sp.stats.gamma.cdf(clinical_x, recovery_shape, loc = 0, scale = recovery_scale)
    
    
    return inccdf, death_cdf, recovery_cdf

def get_cdf(dim):

    x = np.linspace(0,dim, dim) #This is where difference between living/dead comes in
    mu = 3.1
    sigma = 2.5

    a = shape = (mu/sigma)**2
    scale = (sigma**2)/mu

    cdf = sp.stats.gamma.cdf(x,a,loc=0,scale=scale)

    return cdf

