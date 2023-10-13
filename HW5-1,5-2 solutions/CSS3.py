
from calc_coeffs import *
from pylab import *


def H_poly(H, AT, CO2_star, k0,k1,k2,kb,BT):
    value=H**3*AT+H**2*(AT*kb-kb*BT-CO2_star*k1)+ H*(-kb*k1*CO2_star-2*k1*k2*CO2_star) \
        - 2*k1*k2*kb*CO2_star
    return value
    
def CO2sys_ALKpCO2(pCO2, AT,T,S):
    
    num_iter=100
    diff_per=0.001
    # Calculate coefficients
    coeffs = calc_coeffs(T,S)
    
    # Calculate CO2_star
    CO2_star=coeffs['k0']*pCO2 # mol/kg
    
    # Use bracket bisection to find root
    
    # Initial guess
    h_low=10**-10
    h_high=10**-6
    
    for i in arange(num_iter):
        
        h_mid=(h_high+h_low)/2
        
        MidRoot=H_poly(h_mid, AT, CO2_star, coeffs['k0'],coeffs['k1'],coeffs['k2'],coeffs['kb'],coeffs['BT'])
        
        if MidRoot>0:
            # New low root is the mid root
            h_high=h_mid
        elif MidRoot<0:
            h_low=h_mid
        else:
            print('Error')
        
        if abs((h_high-h_low)/h_high*100)<diff_per:
            H_final=h_high
            break
        
    ## Solve system
    HCO3=(coeffs['k1']*CO2_star)/H_final
    CO3=(coeffs['k2']*HCO3)/H_final
    H_plus=H_final
    C_sat=CO2_star+HCO3+CO3
    
    return CO2_star, HCO3, CO3, C_sat, H_plus
    