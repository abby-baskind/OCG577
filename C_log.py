'''
Created on Nov 3, 2021

@author: keiin
'''


from numpy import *
from matplotlib.pyplot import *

#Analytical solution
t_analytical = arange(1,100+1e-10,0.001)
C_analytical = log(t_analytical)  #    #computational solution of C

#Setting up time step and array frame
dt = 0.01 
t = arange(1,100+1e-10,dt) 
C = zeros(size(t))*nan 

#Initial condition
C[0] = log(1)

import time

time0 = time.time()
#Time stepping
for i in arange(size(t)-1):
    dCdt = 1/t[i]
    C[i+1] = C[i] + dt*dCdt

time1 = time.time()

print(time1-time0)


#Plotting
figure(1)
plot(t_analytical,C_analytical,label='analytical')
plot(t,C,label='computational')
xlabel('t')
ylabel('C')
legend()

show()