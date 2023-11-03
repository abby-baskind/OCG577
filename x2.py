'''
Created on Nov 3, 2021

@author: keiin
'''


from numpy import *
from matplotlib.pylab import *


x1 = arange(100)
y1 = x1**2  #    #computational solution of y

dx = 0.1
x = arange(0,100+1e-10,dx)
dydx = 2*x

y = zeros(size(x))*nan

y[0] = 0

for i in arange(size(x)-1):
    y[i+1] = y[i]+dx*dydx[i] 

figure(1)
plot(x1,y1,label='analytical')
plot(x,y,label='computational')
xlabel('x')
ylabel('y')
legend()

show()