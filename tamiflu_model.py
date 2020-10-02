import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

ka = 1.011  # OP blood absorption rate
kf = 0.684  # OP to OC conversion rate
ke = 0.136  # OC elimintation rate

Go = 75  # initial dose of 75mg
OPo = Go * ka
t = np.arange(0, 24, .1)  # 24 hours, iterate by tenth of hour

tdose = 12 # time between doses in hours



def dose(t):
    if t % tdose <= .1:
        return Go
    else:
        return 0

def OP_metabolism(Y, t):
    G = Y[0]
    OP = Y[1]
    OC = Y[2]

    dGdt = -1 * ka * G
    dOPdt = ka * G - kf * OP
    dOCdt = kf * OP - ke * OC

    return [dGdt, dOPdt, dOCdt]


def plot(G, OP, OC):
    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot(1, 1, 1)
    # Plot using red circles
    ax.plot(t, G, 'b-', label='Oral OP Cocentration', markevery=10)
    ax.plot(t, OP, 'r-', label='Plasma OP Concentration', markevery=10)
    ax.plot(t, OC, 'g-', label='Plasma OC Concentration', markevery=10)

    # ax.plot(t, go, 'g-', label='Other Conductance', markevery=10)

    # Set labels and turn grid on
    ax.set(xlabel='Time $t$, hrs', ylabel=r'Concentration', title='Tamiflu Uptake')
    ax.grid(True)
    ax.legend(loc='best')
    # Use space most effectively
    fig.tight_layout()
    # Save as a PNG file
    fig.savefig('Oseltamivir_Metabolism.png')


Yo = [Go, 0, 0]  # initial conditions of G, dose, is 75mg and OP and OC are both 0

out = odeint(OP_metabolism, Yo, t)
G = out[:, 0]
OP = out[:, 1]
OC = out[:, 2]

plot(G, OP, OC)
