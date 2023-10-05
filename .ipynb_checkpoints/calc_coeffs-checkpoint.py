
from numpy import *

def calc_coeffs(T,S):
#'''
#    function [k1,k2,k0,kb,kw,BT]=calc_coeffs(T,S)
#       T   = temperature (degrees C)
#       S   = salinity (PSU)
#'''
# total boron, BT
# estimated concentration of borate based on salinity
    scl=S/1.80655
    BT=0.000232 * scl/10.811
    
    # Meshgrid of temperature and salinity
    # [T,S] = np.meshgrid(T,S)
    
    # some definitions 
    S2 = S * S 
    sqrtS = S**0.5
    T = T + 273.15
    invT = 1.0/T 
    T1 = T/100.0 
    
    # K1, K2 Millero (1995) using Mehrbach data
    k1 = 10**(-1*( 3670.7*invT - \
        62.008 + 9.7944*(log(T)) - \
        0.0118*S + 0.000116*S2)) 
    
    k2 = 10**(-1*(1394.7*invT + 4.777 - \
        0.0184*S + 0.000118*S2)) 
    
    # Kb, Millero (1995) using data from Dickson
    kb = exp((-8966.90 - 2890.53*sqrtS - 77.942*S + \
        1.728*S**1.5 - 0.0996*S2)*invT + \
        (148.0248 + 137.1942*sqrtS + 1.62142*S) + \
        (-24.4344 - 25.085*sqrtS - 0.2474*S) * \
        log(T) + 0.053105*T*sqrtS) 
    
    # Kw, Millero (1995)
    kw = exp(-13847.26*invT + 148.9652 - \
        23.6521*log(T) + \
        (118.67*invT - 5.977 + 1.0495*log(T)) * \
        sqrtS - 0.01615*S) 
    
    # solubility of CO2, K0 (Weiss, 1974)
    # mol kg-1 atmos-1 
    k0 = exp(93.4517/T1 - 60.2409 + \
       23.3585 * log(T1) + \
       S * (0.023517 - 0.023656*T1 + \
       0.0047036*T1*T1))
    
    data = {
        'k0': k0,
        'k1': k1,
        'k2': k2,
        'kb': kb,
        'kw': kw,
        'BT': BT
        }
    
    return data
