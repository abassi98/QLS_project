import numpy as np


def monod(c,K):
    """
    
    Define the monod function
    
    Parameters
    ----------
    c : vector of dimension N_R. The resource concetrations of the system
    K : vector of dimension N_R. The half saturation constant
    
    Returns
    -------
    A vector of dimension N_R, which encapsualtes the dependence on the resource concetrations
    
    """
    return c/(c+K)
