
from numpy import *
from matplotlib.pyplot import *

#surface data:
S = loadtxt('SurfS.txt')
T = loadtxt('SurfT.txt')
DIC = loadtxt('SurfDIC.txt')
Alk = loadtxt('SurfAlk.txt')

#corresponding lat and lon for grid corners:
lat = loadtxt('latitude.txt')
lon = loadtxt('longitude.txt')

#plot DIC:
figure()
pcolormesh(lon,lat,ma.masked_invalid(DIC))
colorbar()
xlabel('Longitude',fontsize=14)
ylabel('Latitude',fontsize=14)
title('Surface DIC',fontsize=16)
show()
