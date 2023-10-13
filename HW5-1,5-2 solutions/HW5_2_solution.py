
from numpy import *
from CSS3 import *

T = 20
S = 34.5
Alk = 2350 * 1e-6
pCO2 = 380*1e-6

CO2_star, HCO3, CO3, Csat, H = CO2sys_ALKpCO2(pCO2, Alk, T,S)
pH = -log10(H)



print(pH)
print(CO2_star*1e6)
print(HCO3*1e6)
print(CO3*1e6)
print(Csat*1e6)

show()
