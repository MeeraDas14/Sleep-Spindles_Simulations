

import matplotlib.pyplot as plt    # import matplotlib
import numpy as np                 # import numpy
import scipy as sp                 # import scipy
import math                        # import basic math functions
import random                      # import basic random number generator functions


# TRN cells Note: TRN and TC cells have the same Na and K channel kinetics: refer Austin Soplata 2017, Destexhe et al, 1993; Traub and Miles 1991

#channel opening rate is given by alpha and closing rate by beta values for all ion channels and is dependent on the membrane potential value 
# at each time point (V[i])

def alpha_n_TRN(V):
    '''
    input variables 
    V: membrane potential at every time point t
    
    output variables:
    
    alpha_n: the opening rate of the potassium gate as a function of membrane potential (V) at each time point t
    '''
    result = 0.032*(15-V)/(np.exp((15-V)/5)-1)
    return result

def beta_n_TRN(V):
    '''
    input variables 
    
    V: membrane potential at every time point t
    
    output variables:
    beta_n: the closing rate of the potassium gate as a function of membrane potential (V) at each time point t
    '''
    result= 0.5*np.exp((10-V)/40)
    return result

# Na+ rate functions for activation gate
def alpha_m_TRN(V):
    '''
    input variables 
    V: membrane potential at every time point t
    
    output variables:
    alpha_m: the opening rate of the sodium activation gate as a function of membrane potential (V) at each time point t
    '''
    result = 0.32*(13-V)/(np.exp((13-V)/4)-1)
    return result

def beta_m_TRN(V):
    '''
    input variables 
    V: membrane potential at every time point t
    
    output variables:
    beta_m: the closing rate of the sodium activation gate as a function of membrane potential (V) at each time point t
    '''
    result = 0.28*(V-40)/(np.exp((V-40)/5)-1)
    return  result

# Na+ rate functions for inactivation gate
def alpha_h_TRN(V):
    '''
    input variables 
    V: membrane potential at every time point t
    
    output variables:
    alpha_h: the opening rate of the sodium inactivation gate as a function of membrane potential (V) at each time point t
    '''
    result= 0.128*np.exp((17-V)/18)
    return result

def beta_h_TRN(V):
    '''
    input variables 
    V: membrane potential at every time point t
    
    output variables:
    beta_h: the closing rate of the sodium inactivation gate as a function of membrane potential (V) at each time point t
    '''
    result= 4/(1+np.exp((40-V)/5))
    return result
    

# infinity forms for K and Na channels in both TRN and TC

def n_inf(V):
    result= alpha_n_TRN(V)/(alpha_n_TRN(V)+beta_n_TRN(V))
    return result

def tau_n(V):
    result = 1/(alpha_n_TRN(V)+beta_n_TRN(V))
    return result

def m_inf(V):
    result=alpha_m_TRN(V)/(alpha_m_TRN(V)+beta_m_TRN(V))
    return result

def tau_m(V):
    result = 1/(alpha_m_TRN(V)+beta_m_TRN(V))
    return result

def h_inf(V):
    result=alpha_h_TRN(V)/(alpha_h_TRN(V)+beta_h_TRN(V))
    return result

def tau_h(V):
    result = 1/(alpha_h_TRN(V)+beta_h_TRN(V))
    return result


# T type calcium channels in TRN

#channel ACTIVATION is given by m with the decay constant tau and INACTIVATION by h 
# the activation and inactivation depend on membrane potential  ---> all adapted from Destexhe1994 and Soplata17

# ALTER THE half life value of V for different simulations : V+52 or V+80
def m_Ca_inf(V):
    result = 1/(1+np.exp(-(V+52)/7.4))
    return float(result)

def tau_m_Ca(V):
    result= 0.44 + 0.15/(np.exp((V+27)/10) +np.exp(-(V+102)/15))
    return result
    
def h_Ca_inf(V):
    result= 1/(1+np.exp((V+80)/5))
    return result

def tau_h_Ca(V):
    result= 22.7 + 0.27/(np.exp((V+48)/4) + np.exp(-(V+407)/50))
    return float( result)




# for the calcium dependent channels, the ACTIVATION and INACTIVATION is dependent on the 
# calcium concentration which is in turn dependent on I_T current

def m_CAN_inf(Ca_conc,n):
    '''input- Ca_conc- intracellular calcium concentration
                n - the value of power raised to the calcium concentration 
                as done in DEstexhe 1994
        out put- m_can_inf'''
    result=20*(np.power(Ca_conc,n))/(20*(np.power(Ca_conc,n))+0.002)
    return result

def tau_m_CAN(Ca_conc,n):
    result= 1/(20*(np.power(Ca_conc,n))+0.002)
    return result

def m_K_Ca_inf(Ca_conc,n):
    result=48*(np.power(Ca_conc,n))/(48*(np.power(Ca_conc,n))+0.03)
    return result

def tau_m_K_Ca(Ca_conc,n):
    result=1/(48*(np.power(Ca_conc,n))+0.03)
    return result


# T type calcium channels in TC

#channel ACTIVATION is given by m with the decay constant tau and INACTIVATION by h 
# the activation and inactivation depend on either membrane potential 

def m_Ca_TC_inf(V):
    result = 1/(1+np.exp(-(V+57)/6.2))
    return float(result)
    
def h_Ca_TC_inf(V):
    result= 1/(1+np.exp((V+81)/4))
    return result

def tau_Ca_TC_h(V):
    result= (30.8 +  (211.4+np.exp((V+113.2)/5))/(1+np.exp((V+84)/3.2)))/3.73
    return float( result)

# Hyperpolarization activated current --> Destexhe et al, 1993

def H_TC_inf(V):
    result = 1/(1+np.exp((V+69.8)/6.5))
    return result

def tau_H_S(V):
    result = np.exp((V+183.6)/15.24)
    return result

def tau_H_F(V):
    result = np.exp((V+158.6)/11.2)/(1+np.exp((V+75)/5.5))
    return result
    
