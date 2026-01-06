import numpy as np

def get_partition(xt, xd, y, ori_sites):
    """
        Finds the partition function for the occupation of the origin region in the open conformation. 
    """
    D = (xd + xt*y + 1)**2.-4.*xt*(y-1.)
    if D<0:
        print("sqrt of negative value")
    sqrt_D = np.sqrt(D)
    lambda_plus = (xd + xt*y + 1)/2 + sqrt_D/2
    return lambda_plus**ori_sites

def open_probability(q_open, epsilon_cost):
    """
        based on the partition function and on the energetic cost of melting the double strands, 
        it calculates the probability of the 'open' (initiation prone) conformation. 
    """
    log_open=np.log(q_open)-epsilon_cost
    open_term = np.exp(log_open)
    return open_term/(open_term + 1)

def get_p_open(a_atp, a_adp, c_atp, c_adp, y, kori, ori_sites, epsilon_cost):
    """
        calculates the probability of open conformation starting from the concentrations of 
        proteins and occupied binding sites.  
    """
    x_t_open = (a_atp-c_atp)/kori
    x_d_open = (a_adp-c_adp)/kori
    q_open = get_partition(x_t_open, x_d_open, y, ori_sites)
    p_open = open_probability(q_open, epsilon_cost)
    return p_open

def get_firing_rate(p_open, k_max):
    """
        calculates the firing rate from the probability of open conformation. 
    """
    return k_max*p_open

def fr(a_atp, a_adp, c_atp, c_adp, y, kori, ori_sites, epsilon_cost, k_max):
    """
        calculates the firing rate starting from the concentrations of 
        proteins and occupied binding sites.  
    """
    p_open=get_p_open(a_atp, a_adp, c_atp, c_adp, y, kori, ori_sites, epsilon_cost)
    f_rate=get_firing_rate(p_open, k_max)
    return f_rate