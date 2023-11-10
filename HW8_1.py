'''
Created on Aug 23, 2021

@author: keiin
'''

from numpy import *
from matplotlib.pyplot import *

dx = array([17e6, 17e6, 17e6])  #(m) longitude
dy = array([8e6, 8e6, 8e6 *2]) #(m) latitude
dz = array([100, 100, 5000])  #(m) depth

vol = dx*dy*dz #(m3)

phi = 20e6 #m3/s
k13 = 1e6 #m3/s
k23 = 1e6 #m3/s
k12 = 10e6 #m3/s

k31 = k13
k32 = k23
k21 = k12

years = 100*86400*365 #(s)
dt = 1*86400*365 #(s) timestep = 1year
t = arange(0,years+10e-10,dt)

C1 = zeros(size(t))*nan
C2 = zeros(size(t))*nan
C3 = zeros(size(t))*nan

C1[0] = 1
C2[0] = 1
C3[0] = 1

for i in arange(size(t)-1):

    dC1 = 1/vol[0]*(phi*(C3[i] - C1[i]) + k31*(C3[i] - C1[i]) + k21*(C2[i] - C1[i]))
    dC2 = 1/vol[1]*(phi*(C1[i] - C2[i]) + k12*(C1[i] - C2[i]) + k32*(C3[i] - C2[i]))
    dC3 = 1/vol[2]*(phi*(C2[i] - C3[i]) + k23*(C2[i] - C3[i]) + k13*(C1[i] - C3[i]))
    
    C1[i+1] = C1[i] + dC1*dt
    C2[i+1] = C2[i] + dC2*dt
    C3[i+1] = C3[i] + dC3*dt

ty = t/(86400*365)

figure('1')
plot(ty, C1,label='C1')
plot(ty, C2,label='C2')
plot(ty, C3,label='C3')
legend()
xlabel('year')

show()





