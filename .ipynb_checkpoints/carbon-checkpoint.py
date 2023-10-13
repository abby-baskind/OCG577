import numpy as np

def correct_dataarray(data):
    """
    """
    
    import numpy 
    
    if not isinstance(data,numpy.ndarray):
        data = numpy.asarray(data)
    return data

def pHsolver_TA_DIC(TA = None, 
                    temperature = None, 
                    salinity = None, 
                    DIC = None, 
                    pHlo = None, pHhi = None):
    
    """
    This function, as it is currently written, takes TA [umol/kg], temperature [degC], 
    salinity [PSU], and DIC [umol/kg] to solve for pH [total] and HCO3, CO3, CO2*, and pCO2,
    all in [umol/kg].
    
    Arguments for TA, T, S, and DIC can be of any length as long as they are all the same
    length as each other and do not have NaN values. 
    
    Required dependencies:
        - calc_coeffs.py: calculates carbonate system coefficients
            function [k1,k2,k0,kb,kw,BT]=calc_coeffs(T,S)
            T   = temperature (degrees C)
            S   = salinity (PSU)
        - H_poly2.py: H polynomial root finder taking TA and DIC as inputs
            H_poly2(H, AT, DIC,K1,K2,KB,BT)
    
    Inputs:
        - total alkalinity in umol/kg
        - DIC in umol/kg
        - temperature in °C
        - salinity in PSU
        - OPTIONAL: pHlo and pHhi in total scale
            - if pHlo and pHhi not provided, default to 6 and 9
    
    Returns:
        - [H+] in umol/kg
        - pH in total scale
        - [HCO3] in umol/kg
        - [CO3] in umol/kg
        - [CO2*] in umol/kg
        - pCO2 in uatm
    """
    
    import calc_coeffs as co2
    import H_poly as hpoly
    import H_poly2 as hpoly2
    import numpy as np 
    import math
    
    # Check that all the inputs are provided
    if np.asarray(TA).any() == None:
        raise KeyError('Please provide TA [umol/kg] with kwarg "TA = *value*".')
    if np.asarray(temperature).any() == None:
        raise KeyError('Please provide temperature [degC] with kwarg "temperature = *value*".')
    if np.asarray(salinity).any() == None:
        raise KeyError('Please provide salinity [PSU] with kwarg "salinity = *value*".')
    if np.asarray(DIC).any() == None:
        raise KeyError('Please provide DIC [umol/kg] with kwarg "DIC = *value*".')
        
    # If args are not already numpy ndarrays,
    # Convert them to ndarrays
    TA = correct_dataarray(TA)
    temperature = correct_dataarray(temperature)
    salinity = correct_dataarray(salinity)
    DIC = correct_dataarray(DIC)
    
        
    # Check for NaN values
    if np.isnan(np.min(TA)) or np.isnan(np.min(temperature)) or np.isnan(np.min(salinity)) or np.isnan(np.min(DIC)):
        raise Exception("No NaN values allowed.") 
    
    # Check all input arrays are same size
    if not (TA.size == temperature.size and TA.size == salinity.size and TA.size == DIC.size):
        raise Exception("Input arrays are not the same size.") 
        
    # Check magnitude of TA and DIC to guess if units are right
    if not TA.min() > 1:
        raise Exception('This function takes TA units of umol/kg. Based on the magnitude of your input, we assume you entered TA in mol/kg.')
    if not DIC.min() > 1:
        raise Exception('This function takes DIC units of umol/kg. Based on the magnitude of your input, we assume you entered DIC in mol/kg.')
    TA = TA*1e-6
    DIC = DIC*1e-6
    
    # Calculate the coefficients
    coeffs = co2.calc_coeffs(temperature, salinity)
    k0 = coeffs['k0']
    k1 = coeffs['k1']
    k2 = coeffs['k2']
    kb = coeffs['kb']
    kw = coeffs['kw']
    BT = coeffs['BT']
    
    # Initialize arrays to store results
    pH = np.zeros(TA.size)
    H = np.zeros(TA.size)
    HCO3 = np.zeros(TA.size)
    CO3 = np.zeros(TA.size)
    pCO2 = np.zeros(TA.size)
    co2star = np.zeros(TA.size)
    
    # If inputs are size 1 are provided
    if TA.size == 1:
        Hhi = 10**(-pHlo)
        Hlo = 10**(-pHhi)
        Hmid = np.mean([Hhi,Hlo])
        threshold = Hhi-Hlo
        while np.abs(threshold) > 5e-23:
            Hmid = np.mean([Hhi,Hlo])
            # H, AT, DIC,K1,K2,KB,BT
            sol = hpoly2.H_poly2(Hmid, TA, DIC, k1, k2, kb, BT)
            if sol >= 0:
                Hhi = Hmid
            elif sol < 0:
                Hlo = Hmid
            if Hlo > Hhi:
                print('mix up')
                break
            threshold = Hhi-Hlo
        Hmid = np.mean([Hhi,Hlo])
        H = Hmid
        pH = -math.log10(H)
        
        # Calculate pCO2 using DIC and H
        # pCO2 = [DIC / k0] * [H^2 / (H^2 + k1 * H + k1k2)]
        # Sarmiento & Gruber (2006) Table 8.2.1 Eq 16
        pCO2 = (DIC / k0) * ((H**2)/((H**2) + (k1 * H) + (k1 * k2)))
        
        # Calculate CO2* from pCO2
        # K0 = [CO2*]/pCO2
        co2star = k0 * pCO2
        
        # K1 = [HCO3][H]/[CO2*]
        HCO3 = (k1 * co2star)/H
        
        # k2 = [CO3][H]/[HCO3]
        CO3 = (k2 * HCO3)/H

    
    # If longer inputs provided
    else:
        for i in np.arange(0,TA.size)-1: 
            Hhi = 10**(-pHlo)
            phlo = 6
            Hlo = 10**(-pHhi)
            Hmid = np.mean([Hhi,Hlo])
            threshold = Hhi-Hlo
            while np.abs(threshold) >= 3e-23:
                Hmid = np.mean([Hhi,Hlo])
                sol = hpoly2.H_poly2(Hmid, TA[i], DIC[i], k1[i], k2[i], kb[i], BT[i])
                # print(sol)
                if sol > 0:
                    Hhi = Hmid
                elif sol <= 0:
                    Hlo = Hmid
                if Hlo > Hhi:
                    print('mix up')
                    break
                threshold = Hhi-Hlo
            Hmid = np.mean([Hhi,Hlo])
            H[i] = Hmid
            pH[i] = -math.log10(H[i])
            
            # Calculate pCO2 using DIC and H
            # pCO2 = [DIC / k0] * [H^2 / (H^2 + k1 * H + k1k2)]
            # Sarmiento & Gruber (2006) Table 8.2.1 Eq 16
            pCO2[i] = (DIC[i] / k0[i]) * ((H[i]**2)/(H[i]**2 + k1[i] * H[i] + k1[i] * k2[i]))

            # Calculate CO2* from pCO2
            # K0 = [CO2*]/pCO2
            co2star[i] = k0[i] * pCO2[i]
        
            # K1 = [HCO3][H]/[CO2*]
            HCO3[i] = (k1[i] * co2star[i])/H[i]
        
            # k2 = [CO3][H]/[HCO3]
            CO3[i] = (k2[i] * HCO3[i])/H[i]
    
    # Return all the data in a dictionary 
    data = {
        '[H+]': H/1e-6,
        'pH': pH,
        '[HCO3]': HCO3/1e-6,
        '[CO3]': CO3/1e-6,
        '[CO2*]': co2star/1e-6,
        'pCO2': pCO2/1e-6,
        'DIC': DIC/1e-6,
        'TA': TA/1e-6,
        'k0': k0
        }
    return data

def pHsolver_TA_pCO2(TA = None, 
                     temperature = None, 
                     salinity = None, 
                     pCO2 = None, 
                     pHlo = None, pHhi = None):
    
    """
    This function, as it is currently written, takes TA [umol/kg], temperature [degC], 
    salinity [PSU], and pCO2 [uatm] to solve for pH [total] and HCO3, CO3, CO2*, and DIC,
    all in [umol/kg].
    
    Arguments for TA, T, S, and pCO2 can be of any length as long as they are all the same
    length as each other and do not have NaN values. 
    
    Required dependencies:
        - calc_coeffs.py: calculates carbonate system coefficients
            function [k1,k2,k0,kb,kw,BT]=calc_coeffs(T,S)
            T   = temperature (degrees C)
            S   = salinity (PSU)
        - H_poly.py: H polynomial root finder
    
    Inputs:
        - total alkalinity in umol/kg
        - pCO2 in uatm
        - temperature in °C
        - salinity in PSU
        - OPTIONAL: pHlo and pHhi in total scale
            - if pHlo and pHhi not provided, default to 6 and 9
    
    Returns:
        - [H+] in umol/kg
        - pH in total scale
        - [HCO3] in umol/kg
        - [CO3] in umol/kg
        - [CO2*] in umol/kg
        - Csat = [HCO3] + [CO3] + [CO2*] in umol/kg
    """
    
    import calc_coeffs as co2
    import H_poly as hpoly
    import H_poly2 as hpoly2
    import numpy as np 
    import math
    
    # Check that all the inputs are provided
    if np.asarray(TA).any() == None:
        raise KeyError('Please provide TA [umol/kg] with kwarg "TA = *value*".')
    if np.asarray(temperature).any() == None:
        raise KeyError('Please provide temperature [degC] with kwarg "temperature = *value*".')
    if np.asarray(salinity).any() == None:
        raise KeyError('Please provide salinity [PSU] with kwarg "salinity = *value*".')
    if np.asarray(pCO2).any() == None:
        raise KeyError('Please provide pCO2 [uatm] with kwarg "pCO2 = *value*".')
        
    # If args are not already numpy ndarrays,
    # Convert them to ndarrays
    TA = correct_dataarray(TA)
    temperature = correct_dataarray(temperature)
    salinity = correct_dataarray(salinity)
    pCO2 = correct_dataarray(pCO2)
    
        
    # Check for NaN values
    if np.isnan(np.min(TA)) or np.isnan(np.min(temperature)) or np.isnan(np.min(salinity)) or np.isnan(np.min(pCO2)):
        raise Exception("No NaN values allowed.") 
    
    # Check all input arrays are same size
    if not (TA.size == temperature.size and TA.size == salinity.size and TA.size == pCO2.size):
        raise Exception("Input arrays are not the same size.") 
    # If pH limits are not explicitly stated
    # Default to 6 and 9
    if pHlo == None or pHhi == None:
        pHhi = 9
        pHlo = 6
        
    # Check magnitude of TA and pCO2 to guess if units are right
    if not TA.min() > 1:
        raise Exception('This function takes TA units of umol/kg. Based on the magnitude of your input, we assume you entered TA in mol/kg.')
    if not pCO2.min() > 1:
        raise Exception('This function takes pCO2 units of uatm. Based on the magnitude of your input, we assume you entered pCO2 in atm.')
    TA = TA*1e-6
    pCO2 = pCO2*1e-6
    
    # Calculate the coefficients
    coeffs = co2.calc_coeffs(temperature, salinity)
    k0 = coeffs['k0']
    k1 = coeffs['k1']
    k2 = coeffs['k2']
    kb = coeffs['kb']
    kw = coeffs['kw']
    BT = coeffs['BT']
    
    # Calculate CO2* from pCO2
    # K0 = [CO2*]/pCO2
    co2star = k0 * pCO2
    
    # Initialize arrays to store results
    pH = np.zeros(TA.size)
    H = np.zeros(TA.size)
    HCO3 = np.zeros(TA.size)
    CO3 = np.zeros(TA.size)
    Csat = np.zeros(TA.size)
    
    # If inputs are size 1 are provided
    if TA.size == 1:
        Hhi = 10**(-pHlo)
        phlo = 6
        Hlo = 10**(-pHhi)
        Hmid = np.mean([Hhi,Hlo])
        threshold = Hhi-Hlo
        while np.abs(threshold) > 5e-23:
            Hmid = np.mean([Hhi,Hlo])
            sol = hpoly.H_poly(Hmid, TA, co2star, k0, k1, k2, kb, BT)
            if sol >= 0:
                Hhi = Hmid
            elif sol < 0:
                Hlo = Hmid
            if Hlo > Hhi:
                print('mix up')
                break
            threshold = Hhi-Hlo
        Hmid = np.mean([Hhi,Hlo])
        H = Hmid
        pH = -math.log10(H)
    
        # K1 = [HCO3][H]/[CO2*]
        HCO3 = (k1 * co2star)/H
        # k2 = [CO3][H]/[HCO3]
        CO3 = (k2 * HCO3)/H
        # Csat = CO3 + CO2* + HCO3
        Csat = CO3 + co2star + HCO3 
    
    # If longer inputs provided
    else:
        for i in np.arange(0,TA.size)-1: 
            Hhi = 10**(-pHlo)
            phlo = 6
            Hlo = 10**(-pHhi)
            Hmid = np.mean([Hhi,Hlo])
            threshold = Hhi-Hlo
            while np.abs(threshold) > 5e-23:
                Hmid = np.mean([Hhi,Hlo])
                sol = hpoly.H_poly(Hmid, TA[i], co2star[i], k0[i], k1[i], k2[i], kb[i], BT[i])
                # print(sol)
                if sol >= 0:
                    Hhi = Hmid
                elif sol < 0:
                    Hlo = Hmid
                if Hlo > Hhi:
                    print('mix up')
                    break
                threshold = Hhi-Hlo
            Hmid = np.mean([Hhi,Hlo])
            H[i] = Hmid
            pH[i] = -math.log10(H[i])
    
            # K1 = [HCO3][H]/[CO2*]
            HCO3[i] = (k1[i] * co2star[i])/H[i]
            # k2 = [CO3][H]/[HCO3]
            CO3[i] = (k2[i] * HCO3[i])/H[i]
            # Csat = CO3 + CO2* + HCO3
            Csat[i] = CO3[i] + co2star[i] + HCO3[i]
    
    # Return all the data in a dictionary 
    data = {
        '[H+]': H/1e-6,
        'pH': pH,
        '[HCO3]': HCO3/1e-6,
        '[CO3]': CO3/1e-6,
        '[CO2*]': co2star/1e-6,
        'DIC': Csat/1e-6,
        'TA': TA/1e-6,
        'pCO2': pCO2/1e-6
        }
    return data

# def pHsolver(TA = None, 
#              temperature = None, 
#              salinity = None, 
#              pCO2 = None, 
#              DIC = None,
#              pHlo = None, pHhi = None):
def pHsolver(temperature = None, salinity = None, pHhi = None, pHlo = None, **kwargs):
    
    """
    This function, as it is currently written, takes TA [umol/kg], temperature [degC], 
    salinity [PSU], and pCO2 [uatm] to solve for pH [total] and HCO3, CO3, CO2*, and DIC,
    all in [umol/kg].
    
    Arguments for TA, T, S, and pCO2 can be of any length as long as they are all the same
    length as each other and do not have NaN values. 
    
    Required dependencies:
        - calc_coeffs.py: calculates carbonate system coefficients
            function [k1,k2,k0,kb,kw,BT]=calc_coeffs(T,S)
            T   = temperature (degrees C)
            S   = salinity (PSU)
        - H_poly.py: H polynomial root finder
    
    Inputs:
        - total alkalinity in umol/kg
        - pCO2 in uatm
        - temperature in °C
        - salinity in PSU
        - OPTIONAL: pHlo and pHhi in total scale
            - if pHlo and pHhi not provided, default to 6 and 9
    
    Returns:
        - [H+] in umol/kg
        - pH in total scale
        - [HCO3] in umol/kg
        - [CO3] in umol/kg
        - [CO2*] in umol/kg
        - Csat = [HCO3] + [CO3] + [CO2*] in umol/kg
    """
    
    import calc_coeffs as co2
    import H_poly as hpoly
    import H_poly2 as hpoly2
    import numpy as np 
    import math
    
    # If pH limits are not explicitly stated
    # Default to 6 and 9
    if pHlo == None or pHhi == None:
        pHhi = 9
        pHlo = 6
        
    temperature = correct_dataarray(temperature)
    
    
    if temperature.size >= 1:
    
        # if TA.all() != None and pCO2.all() != None:
        if 'TA' in kwargs and 'pCO2' in kwargs:
            TA = kwargs.get('TA')
            pCO2 = kwargs.get('pCO2')
            data = pHsolver_TA_pCO2(TA = TA, 
                                    temperature = temperature, 
                                    salinity = salinity, 
                                    pCO2 = pCO2, 
                                    pHlo = pHlo, pHhi = pHhi)
        # elif TA.any() != None and DIC.any() != None:
        if 'TA' in kwargs and 'DIC' in kwargs:
            TA = kwargs.get('TA')
            DIC = kwargs.get('DIC')
            data = pHsolver_TA_DIC(TA = TA, 
                                   temperature = temperature, 
                                   salinity = salinity, 
                                   DIC = DIC, 
                                   pHlo = pHlo, pHhi = pHhi)
    elif temperature.size < 1:
    
        if 'TA' in kwargs and 'pCO2' in kwargs:
            TA = kwargs.get('TA')
            pCO2 = kwargs.get('pCO2')
            data = pHsolver_TA_pCO2(TA = TA, 
                                    temperature = temperature, 
                                    salinity = salinity, 
                                    pCO2 = pCO2, 
                                    pHlo = pHlo, pHhi = pHhi)
        if 'TA' in kwargs and 'DIC' in kwargs:
            TA = kwargs.get('TA')
            DIC = kwargs.get('DIC')
            data = pHsolver_TA_DIC(TA = TA, 
                                   temperature = temperature, 
                                   salinity = salinity, 
                                   DIC = DIC, 
                                   pHlo = pHlo, pHhi = pHhi)
    return data
        
    

        