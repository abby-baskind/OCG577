import numpy as np

def correct_dataarray(data):
    """
    """
    
    import numpy 
    
    if not isinstance(data,numpy.ndarray):
        data = numpy.asarray(data)
    return data

def pHsolver_pH_TA(pH = None, 
                    temperature = None, 
                    salinity = None, 
                    TA = None, 
                    pHlo = None, pHhi = None):
    
    
    import calc_coeffs as co2
    import H_poly as hpoly
    import H_poly2 as hpoly2
    import numpy as np 
    import math
    import numpy.ma as ma
    
    pH = correct_dataarray(pH)
    temperature = correct_dataarray(temperature)
    salinity = correct_dataarray(salinity)
    TA = correct_dataarray(TA)
    
    # Calculate the coefficients
    coeffs = co2.calc_coeffs(temperature, salinity)
    k0 = coeffs['k0']
    k1 = coeffs['k1']
    k2 = coeffs['k2']
    kb = coeffs['kb']
    kw = coeffs['kw']
    BT = coeffs['BT']
    
    # Solve for H
    H = 10**(-pH)
    
    # Solve for B(OH)4
    # # kb = [H][BOH4]/[H3BO3]
    # # TB = BOH4 + H3BO3
    # # H3BO3 = TB - BOH4
    # # kb = [H][BOH4]/(TB - BOH4)
    # # kb * TB - kb * BOH4 = [H][BOH4]
    # # kb  * TB = kb * BOH4 + H * BOH4
    # # kb * TB = BOH4 (kb + H)
    # # BOH4 = (kb * TB) / (H + kb)
    BOH4 = ((kb * BT)/(H * kb))*1e-6
    
    # Solve for OH
    # # kw = [H][OH]
    # # kw/[H] = [OH]
    OH = kw/H
 
    # Solve for HCO3
    # # Sarmiento & Gruber 2006, Table 8.2.1, Eq. 12
    # # TA = [HCO3] + 2[CO3] + [OH] - [H] + [BOH4]
    # # # Rewrite [CO3] in terms of [HCO3]
    # # # k2 = [H][CO3]/[HCO3]
    # # # [CO3] = k2[HCO3]/[H]
    # # TA = [HCO3] + 2 * k2[HCO3]/[H] + [OH] - H + [BOH4]
    # # TA + [H] - [OH] - [BOH4] = [HCO3] (1 + 2 * k2/[H])
    # # [HCO3] = (TA + [H] - [OH] - [BOH4])/((1 + 2 * k2/[H]))
    HCO3 = (TA + H - OH - BOH4)/(1 + 2*k2/H)
    
    # Solve for CO3
    # # k2 = [H][CO3]/[HCO3]
    # # [CO3] = k2[HCO3]/[H]
    CO3 = k2 * HCO3/H
    
    # Solve for [CO2*]
    # # k1 = [H][HCO3]/[CO2*]
    # # [CO2*] = [H][HCO3]/k1
    co2star = (H * HCO3)/k1
    
    # Solve for pCO2
    pCO2 = co2star / k0 
    
    # Solve for DIC
    # # DIC = [HCO3] + [CO3] + [CO2*]
    DIC = HCO3 + CO3 + co2star
    
    # Aragonite saturation Ωarag
    # # Sarmiento & Gruber 2006, Eq. 9.3.2
    # # Ωarag = [CO3][Ca]/KAr
    # Ca: Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967
    # in mol/kg
    Ca = (.0213/40.078)*(salinity/1.80655)
    # Temperature in kelvin
    TempK = np.asarray(temperature) + 273.15
    # CO3 in mol/kg
    CO3molkg = CO3 / 1e6
    # Solubility of aragonite (mol/kg)^2
    # Mucci 1983 qtd.
    # Sarmiento & Gruber 2006, Table 9.3.1, Eq. 2
    lnKAr_1 = -395.918 + (6685.079/TempK) + 71.595 * np.log(TempK) - 0.17959 * TempK 
    lnKAr_2 = (-0.157481 + 202.938/TempK + 0.003978 * TempK) * (np.power(salinity,0.5))
    lnKAr_3 = -0.23067 * salinity + 0.0136808 * (np.power(salinity,3/2))
    lnKAr = lnKAr_1 + lnKAr_2 + lnKAr_3
    KAr = np.exp(lnKAr) 
    OmegaAr = (CO3molkg * Ca)/KAr
        
    Revelle = (3 * TA * DIC - 2 * DIC**2)/((2*DIC - TA)*(TA - DIC))
    
        
    # Return all the data in a dictionary 
    data = {
        '[H+]': H,
        'pH': pH,
        '[HCO3]': HCO3,
        '[CO3]': CO3,
        '[CO2*]': co2star,
        'pCO2': pCO2,
        'DIC': DIC,
        'TA': TA,
        'k0': k0,
        'k1': k1,
        'k2': k2,
        'kw': kw,
        'kb': kb,
        'BT': BT,
        'OH': OH,
        'BOH4': BOH4,
        'OmegaAr': OmegaAr,
        'KAr': KAr,
        'Revelle Factor': Revelle
        }
    return data

def pHsolver_pH_DIC(pH = None, 
                    temperature = None, 
                    salinity = None, 
                    DIC = None, 
                    pHlo = None, pHhi = None):
    
    
    import calc_coeffs as co2
    import H_poly as hpoly
    import H_poly2 as hpoly2
    import numpy as np 
    import math
    import numpy.ma as ma
    
    pH = correct_dataarray(pH)
    temperature = correct_dataarray(temperature)
    salinity = correct_dataarray(salinity)
    DIC = correct_dataarray(DIC)
    
    # Calculate the coefficients
    coeffs = co2.calc_coeffs(temperature, salinity)
    k0 = coeffs['k0']
    k1 = coeffs['k1']
    k2 = coeffs['k2']
    kb = coeffs['kb']
    kw = coeffs['kw']
    BT = coeffs['BT']
    
    # Solve for H
    H = 10**(-pH)
 
    # Solve for HCO3
    # # Sarmiento & Gruber 2006, Table 8.2.1, Eq. 14
    # # DIC = ([H][HCO3]/k1) + ([HCO3]) + (k2[HCO3]/[H])
    # # DIC = [HCO3] ([H]/k1 + 1 + k2/[H])
    # # [HCO3] = DIC/([H]/k1 + 1 + k2/[H])
    HCO3 = (DIC)/((H/k1) + 1 + (k2/H))
    
    # Solve for CO3
    # # Sarmiento & Gruber 2006, Table 8.2.1, Eq. 15
    # # DIC = [H]^2[CO3]/k1k2 + [H][CO3]/k2 + [CO3]
    # # DIC = [CO3] ([H]^2/k1k2 + [H]/k2 + 1)
    # # [CO3] = DIC / ([H]^2/k1k2 + [H]/k2 + 1)
    CO3 = (DIC)/((H**2)/(k1*k2) + (H/k2) + 1)
    
    # Solve for pCO2
    # # Sarmiento & Gruber 2006, Table 8.2.1, Eq. 16
    # # pCO2 = (DIC/k0) * (([H]^2)/([H]^2 + k1[H] + k1k2))
    pCO2 = (DIC/k0) * ((H**2)/(H**2 + k1*H + k1*k2))
    
    # Solve for co2star
    co2star = k0 * pCO2
    
    # Solve for B(OH)4
    # # kb = [H][BOH4]/[H3BO3]
    # # TB = BOH4 + H3BO3
    # # H3BO3 = TB - BOH4
    # # kb = [H][BOH4]/(TB - BOH4)
    # # kb * TB - kb * BOH4 = [H][BOH4]
    # # kb  * TB = kb * BOH4 + H * BOH4
    # # kb * TB = BOH4 (kb + H)
    # # BOH4 = (kb * TB) / (H + kb)
    BOH4 = ((kb * BT)/(H * kb))*1e-6
    
    # Solve for OH
    # # kw = [H][OH]
    # # kw/[H] = [OH]
    OH = kw/H
    
    # Solve for TA
    # # Alk = [HCO3] + 2[CO3] + [OH] - [H] + [BOH4]
    TA = HCO3 + 2 * CO3 + OH + BOH4
    
    # Aragonite saturation Ωarag
    # # Sarmiento & Gruber 2006, Eq. 9.3.2
    # # Ωarag = [CO3][Ca]/KAr
    # Ca: Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967
    # in mol/kg
    Ca = (.0213/40.078)*(salinity/1.80655)
    # Temperature in kelvin
    TempK = np.asarray(temperature) + 273.15
    # CO3 in mol/kg
    CO3molkg = CO3 / 1e6
    # Solubility of aragonite (mol/kg)^2
    # Mucci 1983 qtd.
    # Sarmiento & Gruber 2006, Table 9.3.1, Eq. 2
    lnKAr_1 = -395.918 + (6685.079/TempK) + 71.595 * np.log(TempK) - 0.17959 * TempK 
    lnKAr_2 = (-0.157481 + 202.938/TempK + 0.003978 * TempK) * (np.power(salinity,0.5))
    lnKAr_3 = -0.23067 * salinity + 0.0136808 * (np.power(salinity,3/2))
    lnKAr = lnKAr_1 + lnKAr_2 + lnKAr_3
    KAr = np.exp(lnKAr) 
    OmegaAr = (CO3molkg * Ca)/KAr
        
    Revelle = (3 * TA * DIC - 2 * DIC**2)/((2*DIC - TA)*(TA - DIC))
    
        
    # Return all the data in a dictionary 
    data = {
        '[H+]': H,
        'pH': pH,
        '[HCO3]': HCO3,
        '[CO3]': CO3,
        '[CO2*]': co2star,
        'pCO2': pCO2,
        'DIC': DIC,
        'TA': TA,
        'k0': k0,
        'k1': k1,
        'k2': k2,
        'kw': kw,
        'kb': kb,
        'BT': BT,
        'OH': OH,
        'BOH4': BOH4,
        'OmegaAr': OmegaAr,
        'KAr': KAr,
        'Revelle Factor': Revelle
        }
    return data

def pHsolver_pH_pCO2(pH = None, 
                    temperature = None, 
                    salinity = None, 
                    pCO2 = None, 
                    pHlo = None, pHhi = None):
    
    
    import calc_coeffs as co2
    import H_poly as hpoly
    import H_poly2 as hpoly2
    import numpy as np 
    import math
    import numpy.ma as ma
    
    pH = correct_dataarray(pH)
    temperature = correct_dataarray(temperature)
    salinity = correct_dataarray(salinity)
    pCO2 = correct_dataarray(pCO2)
    
    # Calculate the coefficients
    coeffs = co2.calc_coeffs(temperature, salinity)
    k0 = coeffs['k0']
    k1 = coeffs['k1']
    k2 = coeffs['k2']
    kb = coeffs['kb']
    kw = coeffs['kw']
    BT = coeffs['BT']
    
    # # Initialize arrays to store results
    # TA = np.zeros(pH.shape)
    # H = np.zeros(pH.shape)
    # HCO3 = np.zeros(pH.shape)
    # CO3 = np.zeros(pH.shape)
    # pCO2 = np.zeros(pH.shape)
    # co2star = np.zeros(pH.shape)
    # DIC = np.zeros(pH.shape)
    
    # Solve for H
    H = 10**(-pH)
 
    # Solve for co2star
    co2star = k0 * pCO2
    
    # Solve for HCO3
    HCO3 = (k1 * co2star)/H
    
    # Solve for CO3
    CO3 = (k2 * HCO3)/H
    
    # Solve for B(OH)4
    # # kb = [H][BOH4]/[H3BO3]
    # # TB = BOH4 + H3BO3
    # # H3BO3 = TB - BOH4
    # # kb = [H][BOH4]/(TB - BOH4)
    # # kb * TB - kb * BOH4 = [H][BOH4]
    # # kb  * TB = kb * BOH4 + H * BOH4
    # # kb * TB = BOH4 (kb + H)
    # # BOH4 = (kb * TB) / (H + kb)
    BOH4 = ((kb * BT)/(H * kb))*1e-6
    
    # Solve for OH
    # # kw = [H][OH]
    # # kw/[H] = [OH]
    OH = kw/H
    
    # Solve for DIC
    DIC = HCO3 + co2star + CO3
    
    # Solve for TA
    # # Alk = [HCO3] + 2[CO3] + [OH] - [H] + [BOH4]
    TA = HCO3 + 2 * CO3 + OH + BOH4
    
    # Aragonite saturation Ωarag
    # # Sarmiento & Gruber 2006, Eq. 9.3.2
    # # Ωarag = [CO3][Ca]/KAr
    # Ca: Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967
    # in mol/kg
    Ca = (.0213/40.078)*(salinity/1.80655)
    # Temperature in kelvin
    TempK = np.asarray(temperature) + 273.15
    # CO3 in mol/kg
    CO3molkg = CO3 / 1e6
    # Solubility of aragonite (mol/kg)^2
    # Mucci 1983 qtd.
    # Sarmiento & Gruber 2006, Table 9.3.1, Eq. 2
    lnKAr_1 = -395.918 + (6685.079/TempK) + 71.595 * np.log(TempK) - 0.17959 * TempK 
    lnKAr_2 = (-0.157481 + 202.938/TempK + 0.003978 * TempK) * (np.power(salinity,0.5))
    lnKAr_3 = -0.23067 * salinity + 0.0136808 * (np.power(salinity,3/2))
    lnKAr = lnKAr_1 + lnKAr_2 + lnKAr_3
    KAr = np.exp(lnKAr) 
    OmegaAr = (CO3molkg * Ca)/KAr
        
    Revelle = (3 * TA * DIC - 2 * DIC**2)/((2*DIC - TA)*(TA - DIC))
    
        
    # Return all the data in a dictionary 
    data = {
        '[H+]': H,
        'pH': pH,
        '[HCO3]': HCO3,
        '[CO3]': CO3,
        '[CO2*]': co2star,
        'pCO2': pCO2,
        'DIC': DIC,
        'TA': TA,
        'k0': k0,
        'k1': k1,
        'k2': k2,
        'kw': kw,
        'kb': kb,
        'BT': BT,
        'OH': OH,
        'BOH4': BOH4,
        'OmegaAr': OmegaAr,
        'KAr': KAr,
        'Revelle Factor': Revelle
        }
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
    length as each other. 
    
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
    import numpy.ma as ma
    
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
    
        
    # # Check for NaN values
    # if np.isnan(np.min(TA)) or np.isnan(np.min(temperature)) or np.isnan(np.min(salinity)) or np.isnan(np.min(DIC)):
    #     raise Exception("No NaN values allowed.") 
    
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
    pH = np.zeros(TA.shape)
    H = np.zeros(TA.shape)
    HCO3 = np.zeros(TA.shape)
    CO3 = np.zeros(TA.shape)
    pCO2 = np.zeros(TA.shape)
    co2star = np.zeros(TA.shape)
    
    # If inputs are size 1 are provided
    if TA.size == 1:
        Hhi = 10**(-pHlo)
        Hlo = 10**(-pHhi)
        Hmid = np.mean([Hhi,Hlo])
        threshold = Hhi-Hlo
        while np.abs(threshold) > 1e-22:
            Hmid = np.mean([Hhi,Hlo])
            # H, AT, DIC,K1,K2,KB,BT
            sol = hpoly2.H_poly2(Hmid, TA, DIC, k1, k2, kb, BT)
            print(sol)
            if sol >= 0:
                Hhi = Hmid
            elif sol < 0:
                Hlo = Hmid
            if Hlo > Hhi:
                print('mix up')
                break
            threshold = Hhi-Hlo
            print(threshold)
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
    elif TA.ndim == 1:
        for i in np.arange(0,TA.size)-1: 
            Hhi = 10**(-pHlo)
            phlo = 6
            Hlo = 10**(-pHhi)
            Hmid = np.mean([Hhi,Hlo])
            threshold = Hhi-Hlo
            while np.abs(threshold) > 1e-19:
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
                # print(threshold)
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
            
    elif TA.ndim == 2:
        for i in np.arange(0,TA.shape[0]):
            for j in np.arange(0,TA.shape[1]):
                Hhi = 10**(-pHlo)
                Hlo = 10**(-pHhi)
                Hmid = np.mean([Hhi,Hlo])
                threshold = Hhi-Hlo
                while np.abs(threshold) > 1e-19:
                    Hmid = np.mean([Hhi,Hlo])
                    sol = hpoly2.H_poly2(Hmid, TA[i,j], DIC[i,j], k1[i,j], k2[i,j], kb[i,j], BT[i,j])
                    if not ma.is_masked(sol):
                        if sol > 0:
                            Hhi = Hmid
                        elif sol <= 0:
                            Hlo = Hmid
                        if Hlo > Hhi:
                            print('mix up')
                            break
                        threshold = Hhi-Hlo
                    else:
                        threshold = 0
                        Hhi = np.nan
                        Hlo = np.nan
                # print(threshold)
                Hmid = np.mean([Hhi,Hlo])
                H[i,j] = Hmid
                pH[i,j] = -math.log10(H[i,j])
                
                # Calculate pCO2 using DIC and H
                # pCO2 = [DIC / k0] * [H^2 / (H^2 + k1 * H + k1k2)]
                # Sarmiento & Gruber (2006) Table 8.2.1 Eq 16
                pCO2[i,j] = (DIC[i,j] / k0[i,j]) * ((H[i,j]**2)/(H[i,j]**2 + k1[i,j] * H[i,j] + k1[i,j] * k2[i,j]))

                # Calculate CO2* from pCO2
                # K0 = [CO2*]/pCO2
                co2star[i,j] = k0[i,j] * pCO2[i,j]
        
                # K1 = [HCO3][H]/[CO2*]
                HCO3[i,j] = (k1[i,j] * co2star[i,j])/H[i,j]
        
                # k2 = [CO3][H]/[HCO3]
                CO3[i,j] = (k2[i,j] * HCO3[i,j])/H[i,j]
    elif TA.ndim == 3:
        for i in np.arange(0,TA.shape[0]):
            for j in np.arange(0,TA.shape[1]):
                for k in np.arange(0,TA.shape[2]):
                    Hhi = 10**(-pHlo)
                    Hlo = 10**(-pHhi)
                    Hmid = np.mean([Hhi,Hlo])
                    threshold = Hhi-Hlo
                    while np.abs(threshold) > 1e-19:
                        Hmid = np.mean([Hhi,Hlo])
                        sol = hpoly2.H_poly2(Hmid, TA[i,j,k], DIC[i,j,k], k1[i,j,k], k2[i,j,k], kb[i,j,k], BT[i,j,k])
                        if not ma.is_masked(sol) or not math.isnan(sol):
                            if sol > 0:
                                Hhi = Hmid
                            elif sol <= 0:
                                Hlo = Hmid
                            if Hlo > Hhi:
                                print('mix up')
                                break
                            threshold = Hhi-Hlo
                        else:
                            threshold = 0
                            Hhi = np.nan
                            Hlo = np.nan
                    # print(threshold)
                    Hmid = np.mean([Hhi,Hlo])
                    H[i,j,k] = Hmid
                    pH[i,j,k] = -math.log10(H[i,j])
                    # Calculate pCO2 using DIC and H
                    # pCO2 = [DIC / k0] * [H^2 / (H^2 + k1 * H + k1k2)]
                    # Sarmiento & Gruber (2006) Table 8.2.1 Eq 16
                    pCO2[i,j,k] = (DIC[i,j,k] / k0[i,j,k]) * ((H[i,j,k]**2)/(H[i,j,k]**2 + k1[i,j,k] * H[i,j,k] + k1[i,j,k] * k2[i,j,k]))
                    # Calculate CO2* from pCO2
                    # K0 = [CO2*]/pCO2
                    co2star[i,j,k] = k0[i,j,k] * pCO2[i,j,k]
                    # K1 = [HCO3][H]/[CO2*]
                    HCO3[i,j,k] = (k1[i,j,k] * co2star[i,j,k])/H[i,j,k]
                    # k2 = [CO3][H]/[HCO3]
                    CO3[i,j,k] = (k2[i,j,k] * HCO3[i,j,k])/H[i,j,k]
    else:
        raise Exception("This function currently does not have the capability to process data higher than 3 dimensions.") 
        
    # Aragonite saturation Ωarag
    # # Sarmiento & Gruber 2006, Eq. 9.3.2
    # # Ωarag = [CO3][Ca]/KAr
    # Ca: Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967
    # in mol/kg
    Ca = (.0213/40.078)*(salinity/1.80655)
    # Temperature in kelvin
    TempK = np.asarray(temperature) + 273.15
    # CO3 in mol/kg
    CO3molkg = CO3 / 1e6
    # Solubility of aragonite (mol/kg)^2
    # Mucci 1983 qtd.
    # Sarmiento & Gruber 2006, Table 9.3.1, Eq. 2
    lnKAr_1 = -395.918 + (6685.079/TempK) + 71.595 * np.log(TempK) - 0.17959 * TempK 
    lnKAr_2 = (-0.157481 + 202.938/TempK + 0.003978 * TempK) * (np.power(salinity,0.5))
    lnKAr_3 = -0.23067 * salinity + 0.0136808 * (np.power(salinity,3/2))
    lnKAr = lnKAr_1 + lnKAr_2 + lnKAr_3
    KAr = np.exp(lnKAr) 
    OmegaAr = (CO3molkg * Ca)/KAr
        
    Revelle = (3 * TA * DIC - 2 * DIC**2)/((2*DIC - TA)*(TA - DIC))
    
    # Solve for OH
    # # kw = [H][OH]
    # # kw/[H] = [OH]
    OH = kw/H
    
    # Solve for B(OH)4
    # # kb = [H][BOH4]/[H3BO3]
    # # TB = BOH4 + H3BO3
    # # H3BO3 = TB - BOH4
    # # kb = [H][BOH4]/(TB - BOH4)
    # # kb * TB - kb * BOH4 = [H][BOH4]
    # # kb  * TB = kb * BOH4 + H * BOH4
    # # kb * TB = BOH4 (kb + H)
    # # BOH4 = (kb * TB) / (H + kb)
    BOH4 = ((kb * BT)/(H * kb))*1e-6
    
        
    # Return all the data in a dictionary 
    data = {
        '[H+]': H*1e6,
        'pH': pH,
        '[HCO3]': HCO3*1e6,
        '[CO3]': CO3*1e6,
        '[CO2*]': co2star*1e6,
        'pCO2': pCO2*1e6,
        'DIC': DIC*1e6,
        'TA': TA*1e6,
        'k0': k0,
        'k1': k1,
        'k2': k2,
        'kw': kw,
        'kb': kb,
        'BT': BT*1e6,
        'OH': OH*1e6,
        'BOH4': BOH4,
        'OmegaAr': OmegaAr,
        'KAr': KAr,
        'Revelle Factor': Revelle
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
    import numpy.ma as ma
    
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
    
        
    # # Check for NaN values
    # if np.isnan(np.min(TA)) or np.isnan(np.min(temperature)) or np.isnan(np.min(salinity)) or np.isnan(np.min(pCO2)):
    #     raise Exception("No NaN values allowed.") 
    
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
    pH = np.zeros(TA.shape)
    H = np.zeros(TA.shape)
    HCO3 = np.zeros(TA.shape)
    CO3 = np.zeros(TA.shape)
    Csat = np.zeros(TA.shape)
    
    # If inputs are size 1 are provided
    if TA.size == 1:
        Hhi = 10**(-pHlo)
        phlo = 6
        Hlo = 10**(-pHhi)
        Hmid = np.mean([Hhi,Hlo])
        threshold = Hhi-Hlo
        while np.abs(threshold) > 1e-22:
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
    elif TA.ndim == 1:
        for i in np.arange(0,TA.size)-1: 
            Hhi = 10**(-pHlo)
            phlo = 6
            Hlo = 10**(-pHhi)
            Hmid = np.mean([Hhi,Hlo])
            threshold = Hhi-Hlo
            while np.abs(threshold) > 1e-19:
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
            
    elif TA.ndim == 2:
        for i in np.arange(0,TA.shape[0]):
            for j in np.arange(0,TA.shape[1]):
                Hhi = 10**(-pHlo)
                Hlo = 10**(-pHhi)
                Hmid = np.mean([Hhi,Hlo])
                threshold = Hhi-Hlo
                while np.abs(threshold) > 1e-19:
                    Hmid = np.mean([Hhi,Hlo])
                    sol = hpoly.H_poly(Hmid, TA[i,j], co2star[i,j], k0[i,j], k1[i,j], k2[i,j], kb[i,j], BT[i,j])
                    if not ma.is_masked(sol) or not math.isnan(sol):
                        if sol >= 0:
                            Hhi = Hmid
                        elif sol < 0:
                            Hlo = Hmid
                        if Hlo > Hhi:
                            print('mix up')
                            break
                        threshold = Hhi-Hlo
                    else:
                        threshold = 0
                        Hhi = np.nan
                        Hlo = np.nan
                # print(threshold)
                Hmid = np.mean([Hhi,Hlo])
                H[i,j] = Hmid
                pH[i,j] = -math.log10(H[i,j])
                # K1 = [HCO3][H]/[CO2*]
                HCO3[i,j] = (k1[i,j] * co2star[i,j])/H[i,j]
                # k2 = [CO3][H]/[HCO3]
                CO3[i,j] = (k2[i,j] * HCO3[i,j])/H[i,j]
                # Csat = CO3 + CO2* + HCO3
                Csat[i,j] = CO3[i,j] + co2star[i,j] + HCO3[i,j]
    elif TA.ndim == 3:
        for i in np.arange(0,TA.shape[0]):
            for j in np.arange(0,TA.shape[1]):
                for k in np.arange(0,TA.shape[2]):
                    Hhi = 10**(-pHlo)
                    Hlo = 10**(-pHhi)
                    Hmid = np.mean([Hhi,Hlo])
                    threshold = Hhi-Hlo
                    while np.abs(threshold) > 1e-19:
                        Hmid = np.mean([Hhi,Hlo])
                        sol = hpoly.H_poly(Hmid, TA[i,j,k], co2star[i,j,k], k0[i,j,k], k1[i,j,k], k2[i,j,k], kb[i,j,k], BT[i,j,k])
                        if not ma.is_masked(sol) or not math.isnan(sol):
                            if sol >= 0:
                                Hhi = Hmid
                            elif sol < 0:
                                Hlo = Hmid
                            if Hlo > Hhi:
                                print('mix up')
                                break
                            threshold = Hhi-Hlo
                        else:
                            threshold = 0
                            Hhi = np.nan
                            Hlo = np.nan
                    # print(threshold)
                    Hmid = np.mean([Hhi,Hlo])
                    H[i,j,k] = Hmid
                    pH[i,j,k] = -math.log10(H[i,j])
                    # K1 = [HCO3][H]/[CO2*]
                    HCO3[i,j,k] = (k1[i,j,k] * co2star[i,j,k])/H[i,j,k]
                    # k2 = [CO3][H]/[HCO3]
                    CO3[i,j,k] = (k2[i,j,k] * HCO3[i,j,k])/H[i,j,k]
                    # Csat = CO3 + CO2* + HCO3
                    Csat[i,j,k] = CO3[i,j,k] + co2star[i,j,k] + HCO3[i,j,k]
    else:
        raise Exception("This function currently does not have the capability to process data higher than 3 dimensions.") 
    
    # Aragonite saturation Ωarag
    # # Sarmiento & Gruber 2006, Eq. 9.3.2
    # # Ωarag = [CO3][Ca]/KAr
    # Ca: Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967
    # in mol/kg
    Ca = (.0213/40.078)*(salinity/1.80655)
    # Temperature in kelvin
    TempK = np.asarray(temperature) + 273.15
    # CO3 in mol/kg
    CO3molkg = CO3 / 1e6
    # Solubility of aragonite (mol/kg)^2
    # Mucci 1983 qtd.
    # Sarmiento & Gruber 2006, Table 9.3.1, Eq. 2
    lnKAr_1 = -395.918 + (6685.079/TempK) + 71.595 * np.log(TempK) - 0.17959 * TempK 
    lnKAr_2 = (-0.157481 + 202.938/TempK + 0.003978 * TempK) * (np.power(salinity,0.5))
    lnKAr_3 = -0.23067 * salinity + 0.0136808 * (np.power(salinity,3/2))
    lnKAr = lnKAr_1 + lnKAr_2 + lnKAr_3
    KAr = np.exp(lnKAr) 
    OmegaAr = (CO3molkg * Ca)/KAr
    
    Revelle = (3 * TA * Csat - 2 * Csat**2)/((2*Csat - TA)*(TA - Csat))
    
    # Return all the data in a dictionary 
    data = {
        '[H+]': H,
        'pH': pH,
        '[HCO3]': HCO3,
        '[CO3]': CO3,
        '[CO2*]': co2star,
        'pCO2': pCO2,
        'DIC': DIC,
        'TA': TA,
        'k0': k0,
        'k1': k1,
        'k2': k2,
        'kw': kw,
        'kb': kb,
        'BT': BT,
        'OH': OH,
        'BOH4': BOH4,
        'OmegaAr': OmegaAr,
        'KAr': KAr,
        'Revelle Factor': Revelle
        }
    return data

# def pHsolver(TA = None, 
#              temperature = None, 
#              salinity = None, 
#              pCO2 = None, 
#              DIC = None,
#              pHlo = None, pHhi = None):
def solver(temperature = None, salinity = None, pHhi = None, pHlo = None, **kwargs):
    
    """
    This function is a wrapper function that eventually solves the marine carbonate system. Depending on the input arguments, it calls one of the following: [A] pHsolver_TA_pCO2(), [B] pHsolver_TA_DIC(), [C] pHsolver_pH_TA(), [D] pHsolver_pH_DIC(), or [E] pHsolver_pH_pCO2.
    
    Required dependencies:
        - calc_coeffs.py: calculates carbonate system coefficients
            function [k1,k2,k0,kb,kw,BT]=calc_coeffs(T,S)
            T   = temperature (degrees C)
            S   = salinity (PSU)
        - H_poly.py: H polynomial root finder for pCO2
        - H_poly2.py: H polynomial root finder for DIC
        - numpy
        - math
        - numpy.ma
    
    Inputs:
        - 2 of the following:
            - "TA": total alkalinity in umol/kg
            - "pCO2": pCO2 in uatm
            - "DIC": dissolved inorganic carbon in umol/kg
            - "pH" in total scale
            - EXCLUDING DIC and pCO2
        - "temperature" in °C
        - "salinity" in PSU
        - OPTIONAL: pHlo and pHhi in total scale
            - if pHlo and pHhi not provided, default to 6 and 9
    
    Returns:
        - [H+] in umol/kg
        - pH in total scale
        - [HCO3] in umol/kg
        - [CO3] in umol/kg
        - [CO2*] in umol/kg
        - DIC = [HCO3] + [CO3] + [CO2*] in umol/kg
        - Total Alkalinity (TA) in umol/kg
        - k0
        - k1
        - k2 
        - kw
        - kb
        - total borate
        - B(OH)4-
        - Aragonite saturation (OmegaAr)
        - Aragonite solubility (KAr)
        - approximate Revelle factor
        
    Sources:
        - Riley, J. P., & Tongudai, M. (1967). The major cation/chlorinity ratios in sea water. Chemical Geology, 2, 263-269.
        - Friis, K., Körtzinger, A., & Wallace, D. W. (2003). The salinity normalization of marine inorganic carbon chemistry data. Geophysical Research Letters, 30(2).
        - Mucci, A. (1983). The solubility of calcite and aragonite in seawater at various salinities, temperatures, and one atmosphere total pressure. Am. J. Sci, 283(7), 780-799.
        - Millero, F. J. (1995). Thermodynamics of the carbon dioxide system in the oceans. Geochimica et Cosmochimica Acta, 59(4), 661-677.
        - Weiss, R. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. Marine chemistry, 2(3), 203-215.
        - Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. Limnology and Oceanography: Methods, 12(6), 351-362.
    """
    
    import calc_coeffs as co2
    import H_poly as hpoly
    import H_poly2 as hpoly2
    import numpy as np 
    import math
    import numpy.ma as ma
    
    # If pH limits are not explicitly stated
    # Default to 6 and 9
    if pHlo == None:
        pHlo = 6
    if pHhi == None:
        pHhi = 9
        
        
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
        elif 'TA' in kwargs and 'DIC' in kwargs:
            TA = kwargs.get('TA')
            DIC = kwargs.get('DIC')
            data = pHsolver_TA_DIC(TA = TA, 
                                   temperature = temperature, 
                                   salinity = salinity, 
                                   DIC = DIC, 
                                   pHlo = pHlo, pHhi = pHhi)
            
        elif 'TA' in kwargs and 'pH' in kwargs:
            TA = kwargs.get('TA')
            pH = kwargs.get('pH')
            data = pHsolver_pH_TA(TA = TA, 
                                   temperature = temperature, 
                                   salinity = salinity, 
                                   pH = pH)
        elif 'DIC' in kwargs and 'pH' in kwargs:
            DIC = kwargs.get('DIC')
            pH = kwargs.get('pH')
            data = pHsolver_pH_DIC(DIC = DIC, 
                                   temperature = temperature, 
                                   salinity = salinity, 
                                   pH = pH)
        elif 'pCO2' in kwargs and 'pH' in kwargs:
            pCO2 = kwargs.get('pCO2')
            pH = kwargs.get('pH')
            data = pHsolver_pH_pCO2(pCO2 = pCO2, 
                                   temperature = temperature, 
                                   salinity = salinity, 
                                   pH = pH)
        else: 
            raise KeyError('Insufficient arguments provided. Please provide either \n(A) TA and pCO2, \n(B) TA and DIC, \n(C) TA and pH, \n(D) DIC and pH or \n(E) pCO2 and pH.')
            
        if 'O2' in kwargs:
            O2 = kwargs.get('O2')
            # Reference Salinity
            # Default = 35 PSU
            if 'Sref' in kwargs:
                Sref = kwargs.get('Sref')
            else:
                Sref = 35
    
            # Reference salinity normalized alkalinity
            # Default = 2298 umol/kg
            # from Chen et al (2022)
            if 'sAlkref' in kwargs:
                sAlkref = kwargs.get('sAlkref')
            else:
                sAlkref = 2298 #umol/kg, Chen et al, 2022
    
            # Reference salinity normalized DIC
            # Default = 1967 umol/kg
            # from Chen et al (2022)
            if 'sDICref' in kwargs:
                sDICref = kwargs.get('sDICref')
            else:
                sDICref = 1967 #umol/kg, Chen et al, 2022
    
            # Estimation of Cant in the surface ocean
            # https://bg.copernicus.org/articles/10/2169/2013/bg-10-2169-2013.pdf
            # S. Khatiwala et al (2013, Biogeosciences): See figure 4
            if 'Cant' in kwargs:
                Cant = kwargs.get('Cant')
            else:
                Cant = 40
                
            data2 = pumps(data['TA'], data['DIC'], temperature, salinity, O2, 
                          Sref = Sref, sAlkref = sAlkref, sDICref = sDICref, Cant = Cant)
        else:
            data2 = {}
            
        DATA = {**data, **data2}

    return DATA

def pumps(TA, DIC, T, S, O2, **kwargs):
    
    # Reference Salinity
    # Default = 35 PSU
    if 'Sref' in kwargs:
        Sref = kwargs.get('Sref')
    else:
        Sref = 35
    
    # Reference salinity normalized alkalinity
    # Default = 2298 umol/kg
    # from Chen et al (2022)
    if 'sAlkref' in kwargs:
        sAlkref = kwargs.get('sAlkref')
    else:
        sAlkref = 2298 #umol/kg, Chen et al, 2022
    
    # Reference salinity normalized DIC
    # Default = 1967 umol/kg
    # from Chen et al (2022)
    if 'sDICref' in kwargs:
        sDICref = kwargs.get('sDICref')
    else:
        sDICref = 1967 #umol/kg, Chen et al, 2022
    
    # Estimation of Cant in the surface ocean
    # https://bg.copernicus.org/articles/10/2169/2013/bg-10-2169-2013.pdf
    # S. Khatiwala et al (2013, Biogeosciences): See figure 4
    if 'Cant' in kwargs:
        Cant = kwargs.get('Cant')
    else:
        Cant = 40 
    
    TA = correct_dataarray(TA)
    T = correct_dataarray(T)
    S = correct_dataarray(S)
    DIC = correct_dataarray(DIC)
    O2 = correct_dataarray(O2)
    
    O2sat = o2sat(T,S)
    AOU = O2sat - O2
    
    # Redfield ratio 
    # C:N:P:O2
    # 117:16:1:150
    C = 117
    N = 16
    P = 1
    O = 150
    CtoO2 = C/O
    TAtoO2 = (N + P)/O2
    TAtoC = (N + P)/C
    
    # Soft tissue pump
    # Roughly based on Sarmiento & Gruber 2006, Eq. 8.4.1 but in terms of O2
    dCsoft = CtoO2 * AOU
    dTAsoft = -TAtoO2 * AOU
    
    # Salinity normalization
    # nX = (X/S) * Sref
    # Sarmiento & Gruber (2006) and Chen et al (2022) use Sref = 35
    sDIC = (DIC/S) * Sref
    sAlk = (TA/S) * Sref
    
    # Carbonate pump
    dCcarb = 0.5 * (sAlk - sAlkref)
    dTAcarb = sAlk - sAlkref - dTAsoft
    
    # Gas exchange
    dCgas = sDIC - sDICref - Cant - dCsoft - dCcarb
    dTAgas = -TAtoC * dCgas
    
    data2 = {
        'S_ref': Sref,
        'sDIC': sDIC,
        'sDIC_ref': sDICref,
        'sAlk': sAlk,
        'sAlk_ref': sAlkref,
        'C_ant': Cant,
        '∆C_soft': dCsoft,
        '∆TA_soft': dTAsoft,
        'O2_sat': O2sat,
        'AOU': AOU,
        '∆C_carb': dCcarb,
        '∆TA_carb': dTAcarb,
        '∆C_gasex': dCgas,
        '∆TA_gasex': dTAgas,
        'C:O2': CtoO2,
        'TA:O2': TAtoO2,
        'TA:C': TAtoC,
        }
    return data2

def o2sat(temp, S):
    
    import numpy as np
    
    """
    This function computes oxygen solubility in seawater in umol/kg according to Garcia and Gordon (1992). Results are valid for 0°C ≤ T ≤ 40°C and 0 PSU ≤ S ≤ 42 PSU.
    
    Usage is as follows:
        o2_sat = o2sat(temperature, salinity),
        where o2_sat is oxygen saturation in umol/kg,
        temperature is in °C, and salinity in PSU.
    
    Source:
        Garcia, H. E., & Gordon, L. I. (1992). Oxygen solubility in seawater: Better fitting equations. Limnology and oceanography, 37(6), 1307-1312.
    """
    
    A0 = 5.80818
    A1 = 3.20684
    A2 = 4.11890
    A3 = 4.93845
    A4 = 1.01567
    A5 = 1.41575
    B0 = -7.01211e-3
    B1 = -7.25958e-3
    B2 = -7.93334e-3
    B3 = -5.54491e-3
    C0 = -1.32412e-7
    
    T = np.log((298.15 - temp)*(273.15 + temp)**(-1))
    
    lnC = A0 + (A1 * T) + (A2 * T**2) + (A3 * T**3) + (A4 * T**4) + (A5 * T**5) + \
        S * (B0 + (B1 * T) + (B2 * T**2) + (B3 * T**3)) + C0 * S**2
    C = np.exp(lnC)
    
    return C

def density_1atm(T,S):
    
    """
    This function calculates the density (kg/m^3) of seawater at 1atm of pressure using temperature (°C) and salinity (PSU).
    
    The equation of state is of the form:
        rho = rho0 + AS + BS^(3/2) + CS
        rho0 = a + bT + cT^2 + dT^3 + eT^4 + fT^5
        A = g + hT + iT^2 + jT^3 + kT^4
        B = l + mT + nT^2
        C = constant
    
    Usage is as follows:
        rho = density_1atm(temperature, salinity),
        where rho is density at 1atm in kg/m^3,
        temperature is in °C, and salinity in PSU.
        
    Source:
        Millero, F. J., & Poisson, A. (1981). International one-atmosphere equation of state of seawater. Deep Sea Research Part A. Oceanographic Research Papers, 28(6), 625-629.
    """
    
    # Reference density calculation 
    # Bigg, 1967, British Journal of Applied Physics, 8, 521–537
    a = 999.842594
    b = 6.793952e-2
    c = -9.095290e-3
    d = 1.001685e-4
    e = -1.120083e-6
    f = 6.536339e-9
    rho0 = a + (b * T) + (c * T**2) + (d * T**3) + (e * T**4) + (f * T**5)
    
    # First coeff A
    g = 8.24493e-1
    h = -4.0899e-3
    i = 7.6438e-5
    j = -8.2467e-7
    k = 5.3875e-9
    A = g + (h * T) + (i * T**2) + (j * T**3) + (k * T**4)
    
    # Second coeff B
    l = -5.72466e-3
    m = 1.0227e-4
    n = -1.6546e-6
    B = l + (m * T) + (n * T**2)
    
    # Third coeff C
    C = 4.8313e-4
    
    # Equation of state
    rho = rho0 + (A * S) + (B * np.power(S, 1.5)) + (C * S)
    
    return rho
    
    
    
        
    

        