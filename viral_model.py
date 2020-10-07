import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

dt = 0.01
time = np.arange(0, 7, dt)
pV = 210
beta = 5e-7
betap = 3e-8
V0 = 1e+4
T0 = 7e+7
gT = 0.8
deltaV = 5
deltaI = 2

ka = 11.04
ke = 2.64
EC50 = 30
emax = 0.98
omega = 4.63


V = np.zeros(1,len(time))
V[0] = V0
T = V;
T[0] = 7e+7;
I = T;
I[0] = 0;
De = I;
D = I;

Dadmin= 75 #75 mg tablet
Td0 = 28/24 #start of treatment day
Td = 12/24 #period of dosage (1 dose per 12 hours)
#Td= 24/24 (one dose per 24 hours)

def viral_model():
    dVdx = (1-emax*y(5)/(y(5)+ EC50))*pV*y(3)-deltaV*y(1)-beta*y(1)*y(2); #dV/dt
    dTdt = gT*y(2)*(1-y(2)+y(3))/T0-betap*y(1)*y(2); #dT/dt
    dIdt = betap*y(1)*y(2)-deltaI*y(3); #dI/dt
    #ynew(4) = ka*y(4);
#ynew(5) = omega*ka*y(4)-ke*y(5); %dD/dt