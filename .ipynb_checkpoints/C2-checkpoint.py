'''
Created on Nov 3, 2021

@author: keiin
'''


from pylab import *

#Analytical solution
t_analytical = arange(100)
C_analytical = t_analytical**2  #    #computational solution of C

#Setting up time step and array frame
dt = 0.00001
t = arange(0,100+1e-20,dt)
C = zeros(size(t))*nan

#Initial condition
C[0] = 0

#Time stepping
for i in arange(size(t)-1):
    dCdt = 2*t[i]
    C[i+1] = C[i] + dt*dCdt

#Plotting
figure(1)
plot(t_analytical,C_analytical,label='analytical')
plot(t,C,label='computational')
xlabel('t')
ylabel('C')
legend()

show()