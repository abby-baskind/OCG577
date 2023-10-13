#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import *
from matplotlib.pyplot import *
from calc_coeffs import *


#First part=========
T = 20
S = 34.5
Alk = 2350 * 1e-6
K = calc_coeffs(T,S)
pCO2 = 380*1e-6
print(K['k0'],K['k1'],K['k2'],K['BT'])
a = K['k0']
print(a)
#Second part========

T = arange(0,35+1e-10)
S = arange(30,38+1e-10)

k0_box = zeros((len(T),len(S)))*nan
k1_box = copy(k0_box)
k2_box = copy(k0_box)

for i in arange(len(T)):
    for j in arange(len(S)):
        
        K = calc_coeffs(T[i],S[j])
        
        k0_box[i,j] = K['k0']
        k1_box[i,j] = K['k1']
        k2_box[i,j] = K['k2']

figure(1,figsize=(10,4))
subplot(1,3,1)
pcolormesh(T,S,k0_box.T)
xlabel('T ($\N{DEGREE SIGN}$C)')
ylabel('S (o/oo)')
title('k0 (mol/kg)')
colorbar()

subplot(1,3,2)
pcolormesh(T,S,k1_box.T)
xlabel('T ($\N{DEGREE SIGN}$C)')
ylabel('S (o/oo)')
title('k1 (mol/kg)')
colorbar()

subplot(1,3,3)
pcolormesh(T,S,k2_box.T)
xlabel('T ($\N{DEGREE SIGN}$C)')
ylabel('S (o/oo)')
title('k2 (mol/kg)')
colorbar()

subplots_adjust(wspace=.7,top=0.8 , bottom=0.2)

show()
