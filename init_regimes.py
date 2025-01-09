
import matplotlib.pyplot as plt    # import matplotlib
import numpy as np                 # import numpy
import scipy as sp                 # import scipy
import math                        # import basic math functions
import random
import sys
import seaborn as sns  
# sys.path.append('C:\Users\mahal\Thesis_code\oscillations-and-memory')
import gating_variables as gate



def init_isolated_TRN(T, regime_name,k):

    '''input arguments: T: duration of simulation'''
        
    g_Na_TRN=200   # Na+ channels (200)200
    g_K_TRN=20  # K+ channels ( 15,20,25) 20
    g_T_TRN=3# t channel

    if regime_name=="regime 1":
       
        g_L_TRN= 0.09  
        I_bg=0 

    elif regime_name=="regime 1_v2":
        
        g_L_TRN= 0.1 
        I_bg=0 

    elif regime_name=="regime 2":
        g_L_TRN=0.01 # 0.01 or 0.02
        I_bg = 0

    elif regime_name=="regime 3":
        g_L_TRN=0.05
        I_bg=0

    elif regime_name=="test leaky":
        g_L_TRN_list= [0.01, 0.03, 0.05, 0.06, 0.09, 0.1]
        g_L_TRN = g_L_TRN_list[k]
        I_bg=0
        # g_L_TRN= g_L_TRN_list[k]

    elif regime_name == "test gT":
        g_T_TRN_list= [1, 1.5, 1.75, 2, 2.5, 2.75, 3]
        g_T_TRN = g_T_TRN_list[k]
        g_L_TRN= 0.09
        I_bg=0

    elif regime_name == "gT=0":
        g_L_TRN= 0.09
        g_T_TRN = 0
        I_bg= 0



    elif regime_name =="gNa = 0":
       
        g_L_TRN= 0.09
        g_T_TRN = 3
        I_bg= 0
        g_Na_TRN = 0

    elif regime_name =="gK = 0":
        g_L_TRN= 0.09
        g_T_TRN = 3
        I_bg= 0
        g_K_TRN = 0

    elif regime_name == "gK=gNa=0":
        g_L_TRN= 0.09
        g_T_TRN = 3
        I_bg= 0
        g_K_TRN = 0
        g_Na_TRN = 0

        

    Ca_conc=np.zeros(len(T))
    Ca_conc[0]=0.002 #mM of calcium   ---> CHECK THIS INITIAL CONC.

    
    # Initiatlise ionic reversal potential (in mV)
    V_Na_TRN=50 # Na+ channels
    V_K_TRN=-100 # K+ channels
    V_L_TRN=-77 # leaky channels CHANGE TO -77 9/3
    V_T_TRN=120
    
    V_TRN=np.zeros(len(T))
    V_TRN[0]= -60 ## switch V_init to 0 from -85 (1 March 2022), switched on 2 March
    

    # current contributed by ion channels in each time point 
    I_K_TRN=np.zeros(len(T)) # Potassium
    I_Na_TRN=np.zeros(len(T)) # Sodium
    I_L_TRN=np.zeros(len(T))  # Leaky channels
    I_T_TRN=np.zeros(len(T))

    
    
    
    n_TRN=np.zeros(len(T))
    m_TRN=np.zeros(len(T))
    h_TRN=np.zeros(len(T))

    #initial rate values for K and Na rate equations
    n_TRN[0]=gate.alpha_n_TRN(V_TRN[0])/(gate.alpha_n_TRN(V_TRN[0])+gate.beta_n_TRN(V_TRN[0]))
    m_TRN[0]=gate.alpha_m_TRN(V_TRN[0])/(gate.alpha_m_TRN(V_TRN[0])+ gate.beta_m_TRN(V_TRN[0]))
    h_TRN[0]=gate.alpha_h_TRN(V_TRN[0])/(gate.alpha_h_TRN(V_TRN[0])+gate.beta_h_TRN(V_TRN[0]))

    #initial rate values for 
    m_Ca=np.zeros(len(T))
    h_Ca=np.zeros(len(T))
    m_Ca[0]=gate.m_Ca_inf(0)
    h_Ca[0]=gate.h_Ca_inf(0)

    return g_Na_TRN, g_K_TRN, g_T_TRN, g_L_TRN, I_bg, V_Na_TRN, V_K_TRN, V_T_TRN, V_L_TRN, V_TRN, I_K_TRN, I_Na_TRN, I_L_TRN, I_T_TRN, n_TRN, m_TRN, h_TRN, m_Ca, h_Ca



    
def init_isolated_TC(T, regime_name,k):
    
    '''input arguments: T: duration of simulation'''

    #if required to abolish either Na, K , Ca or H current then set the gmax value as o  
    g_Na_TC=90 
    g_K_TC=10  
    
    
    

    if regime_name=="regime 2":
        g_L_TC=0.01 
        g_KL_TC= 0
        # g_H=0.015 #0.015
        I_bg=0
        g_H=0.015
        g_T_TC=2

    elif regime_name=="spontaneous":
        g_L_TC=0.01 
        g_KL_TC= 0.0172
        # g_H=0.015
        I_bg=0
        g_H=0.015
        g_T_TC=2


    elif regime_name=="test leaky":
        g_L_TC_list= list([0.01, 0.03, 0.05, 0.07, 0.09, 0.1])
        g_L_TC = g_L_TC_list[k]
        # g_L_TC_list=list([0.3, 0.4, 0.5])
        # g_L_TC_list=list([0.03, 0.06, 0.09]) # used this for comparing gLmax effect on burst latency
        g_KL_TC= 0
        # g_H=0.015
        I_bg= 0
        g_H=0.015
        g_T_TC=2

    elif regime_name=="test gT1":
        g_T_TC_list= list([1.5,2.0])

        g_T_TC = 1.5
        g_L_TC=0.01 
        g_KL_TC= 0
        I_bg= 0
        g_H=0.015
    
    elif regime_name=="test gT2":
        g_T_TC_list= list([1.5,2.0])

        g_T_TC = 2.0
        g_L_TC=0.01 
        g_KL_TC= 0
        I_bg= 0
        g_H=0.015

    elif regime_name == "test_gT":
        g_T_TC_list= [0.5, 0.75, 1.5, 1.65, 1.8, 1.95, 2.0, 2.5]

        g_T_TC = g_T_TC_list[k]
        g_L_TC=0.01 
        g_KL_TC= 0
        I_bg= 0
        g_H=0.015
  
    elif regime_name=="test gH1":
        # g_H_TC_list= list([ 0, 0.01, 0.015, 0.02])
        g_H = 0
        g_L_TC=0.01 
        g_T_TC = 2
        g_KL_TC= 0
        I_bg= 0

    
    elif regime_name=="test gH2":
        # g_H_TC_list= list([ 0, 0.01, 0.015, 0.02])
        g_H = 0.01
        g_L_TC=0.01 
        g_T_TC = 2
        g_KL_TC= 0
        I_bg= 0

    elif regime_name=="test gH3":
        # g_H_TC_list= list([ 0, 0.01, 0.015, 0.02])
        g_H = 0.015
        g_L_TC=0.01 
        g_T_TC = 2
        g_KL_TC= 0
        I_bg= 0
    
    elif regime_name=="test gH4":
        # g_H_TC_list= list([ 0, 0.01, 0.015, 0.02])
        g_H = 0.02
        g_L_TC=0.01 
        g_T_TC = 2
        g_KL_TC= 0
        I_bg= 0

    elif regime_name== "test_gH":
        g_H_TC_list= [0, 0.005,0.01, 0.015, 0.017, 0.02, 0.025, 0.05, 0.1,0.2,0.5]
        g_H = g_H_TC_list[k]
        g_L_TC=0.01 
        g_T_TC = 2
        g_KL_TC= 0
        I_bg= 0

    elif regime_name == "gT=0":
        g_H = 0.015
        g_L_TC=0.01 
        g_T_TC = 0
        g_KL_TC= 0
        I_bg= 0

    
    elif regime_name == "gH=0":
        g_H = 0
        g_L_TC=0.01 
        g_T_TC = 2
        g_KL_TC= 0
        I_bg= 0

    elif regime_name == "gT and gH=0":
        g_H = 0
        g_L_TC=0.01 
        g_T_TC = 0
        g_KL_TC= 0
        I_bg= 0

    elif regime_name =="gNa = 0":
        g_H = 0.015
        g_L_TC=0.01 
        g_T_TC = 0
        g_KL_TC= 0
        I_bg= 0
        g_Na_TC = 0

    elif regime_name =="gK = 0":
        g_H = 0.015
        g_L_TC=0.01 
        g_T_TC = 0
        g_KL_TC= 0
        I_bg= 0
        g_K_TC = 0

    elif regime_name =="gH=gK=gNa=0":
        g_H = 0
        g_L_TC=0.01 
        g_T_TC = 0
        g_KL_TC= 0
        I_bg= 0
        g_K_TC = 0
        g_Na_TC = 0


    Ca_conc=np.zeros(len(T))
    Ca_conc[0]=0.002 #mM of calcium   ---> CHECK THIS INITIAL CONC.

    
    V_Na_TC=50
    V_K_TC=-100
    V_H=-43 #mV
    V_L_TC=-70 #(changed from -70 mV)
    V_T_TC = 120 # used this instead of using the dynamic equation in Austin Sopalata 2017
    V_KL_TC= -100

    # membrane potential at each time point
    V_TC=np.zeros(len(T))
    V_TC[0]=-60 #change between 0 to -60 or -65 to Vleak value i.e, -70 

    #current contributed by ion channels in each time point 
    I_K_TC=np.zeros(len(T)) # Potassium
    I_Na_TC=np.zeros(len(T)) # Sodium
    I_L_TC=np.zeros(len(T)) # k Leaky channels
    I_KL_TC= np.zeros(len(T)) # k Leaky channels
    I_H_TC=np.zeros(len(T))
    I_T_TC=np.zeros(len(T))
  

    # open probabilities of each channel gates at each time point
    n_TC=np.zeros(len(T))
    m_TC=np.zeros(len(T))
    h_TC=np.zeros(len(T))
    S_H_TC=np.zeros(len(T))
    F_H_TC=np.zeros(len(T))

     #initial rate values for K and Na rate equations
    n_TC[0]=gate.alpha_n_TRN(V_TC[0])/(gate.alpha_n_TRN(V_TC[0])+gate.beta_n_TRN(V_TC[0]))
    m_TC[0]=gate.alpha_m_TRN(V_TC[0])/(gate.alpha_m_TRN(V_TC[0])+gate.beta_m_TRN(V_TC[0]))
    h_TC[0]=gate.alpha_h_TRN(V_TC[0])/(gate.alpha_h_TRN(V_TC[0])+gate.beta_h_TRN(V_TC[0]))
    alpha_S_TC_0= gate.H_TC_inf(V_TC[0])/gate.tau_H_S(V_TC[0])
    beta_S_TC_0= (1-gate.H_TC_inf(V_TC[0]))/gate.tau_H_S(V_TC[0])
    S_H_TC[0]= alpha_S_TC_0/(alpha_S_TC_0 + beta_S_TC_0)

    alpha_F_TC_0= gate.H_TC_inf(V_TC[0])/gate.tau_H_F(V_TC[0])
    beta_F_TC_0= (1-gate.H_TC_inf(V_TC[0]))/gate.tau_H_F(V_TC[0])
    F_H_TC[0]= alpha_F_TC_0/(alpha_F_TC_0 + beta_F_TC_0)
    

    #initial rate values for 
    m_Ca_TC=np.zeros(len(T))
    h_Ca_TC=np.zeros(len(T))

    m_Ca_TC[0]=gate.m_Ca_TC_inf(0)
    h_Ca_TC[0]=gate.h_Ca_TC_inf(0)

    return g_Na_TC, g_K_TC, g_T_TC, g_L_TC, I_bg, g_H, g_KL_TC, V_Na_TC, V_K_TC, V_T_TC, V_L_TC, V_H, V_KL_TC, V_TC, I_K_TC, I_Na_TC, I_L_TC, I_T_TC,I_H_TC, I_KL_TC, n_TC, m_TC, h_TC, m_Ca_TC, h_Ca_TC, S_H_TC, F_H_TC



    

def init_synaptic(T, manipulation):

    # conductance values for GABAa, AMPA and GABAb. Note: here, we only vary GABAb as due to a slower rate, this affects the ;latency of bursting effectively compared to GABAa. 
    # Besides, both GABAa and AMPA are essential for burst generation and changing their values effects occurence of bursts
    
    
    if manipulation == "standard":
        # state variables for synaptic current from TRN to TC
        g_GABAb = 0.01          #  0.001 mS/cm^2 ---.  ''' increase GABAb to see if that allows bursting in TC without an external pulse ''' (Destexhe93 - 4 nS, or 0.013 mS/cm2)     
        g_GABAa = 0.069 #mS/cm^2 (from 0.069)
        g_AMPA = 0.4 #mS/cm^2  ((destexhe93- 0.02, or 1nS) 


    elif manipulation == "test GABAb_1":
        g_GABAb = 0.01
        g_GABAa=0.069 #mS/cm^2 (from 0.069)
        g_AMPA=0.4 #mS/cm^2  ((destexhe93- 0.02, or 1nS) 

    elif manipulation == "test GABAb_2":
        g_GABAb = 0.03
        g_GABAa = 0.069 #mS/cm^2 (from 0.069)
        g_AMPA = 0.4 #mS/cm^2  ((destexhe93- 0.02, or 1nS) 

    elif manipulation == "test GABAb_3":
        g_GABAb = 0.05
        g_GABAa = 0.069 #mS/cm^2 (from 0.069)
        g_AMPA = 0.4 #mS/cm^2  ((destexhe93- 0.02, or 1nS) 

    elif manipulation == "g_GABAb=0":
        g_GABAb=0
        g_GABAa=0.069 #mS/cm^2 (from 0.069)
        g_AMPA=0.4 #mS/cm^2  ((destexhe93- 0.02, or 1nS) 

    elif manipulation == "g_GABAa=0":
        g_GABAb=0.01
        g_GABAa=0 #mS/cm^2 (from 0.069)
        g_AMPA=0.4 #mS/cm^2  ((destexhe93- 0.02, or 1nS) 

    elif manipulation == "g_AMPA=0":
        
        g_AMPA=0
        g_GABAb=0.01
        g_GABAa=0.069 #mS/cm^2 (from 0.069)

    

    '''TRN to TC'''
    #initialise synaptic equation variables for currents from TRN to TC: GABAa and GABAb
    
    # GABABb
    g=np.zeros(len(T))
    r=np.zeros(len(T))

    g[0]=0
    r[0]=0

    #state varaible equation (r,g)
    k1= 0.5 #(mM^-1 ms^-1)
    k2=0.0012 #(ms^-1)
    k3=0.18 #(ms^-1)
    k4= 0.034 #(ms^-1)

    # GABAa
    s_GABAa=np.zeros(len(T))
    s_GABAa[0]=0 #find out the initial value, assuming that this denotes the initial concentration of the neurotramsitter (NT), this can be taken as 0

    #reversal potentials for currents from TRN --> TC 
    V_GABAb= -95 #mV, '' GABAb
    V_GABAa=-80 #mV, reveral potential of GABAa
    tau_GABAa=5 #ms, decay constant

    '''TC to TRN'''
   #initialise synaptic equation variables for currents from TC to TRN: AMPA

    s_AMPA=np.zeros(len(T))
    s_AMPA[0]=0 #find out the initial value, assuming that this denotes the initial concentration of the neurotramsitter (NT), this can be taken as 0

    V_AMPA=0 #mV, reveral potential of AMPA
    tau_AMPA=2 #ms, decay constant for AMPA


    return g_GABAa, g_GABAb, g_AMPA, s_GABAa, s_AMPA, g, r, k1, k2, k3, k4, V_GABAa, V_AMPA, V_GABAb, tau_GABAa, tau_AMPA

    
