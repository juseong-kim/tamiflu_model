import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

ka = 1.011  # OP blood absorption rate
kf = 0.684  # OP to OC conversion rate
ke = 0.136  # OC elimintation rate

Go = 75  # initial dose of 75mg
OPo = Go * ka

tdose = 12  # time between doses in hours
num_doses = 4
dt = .1  # dose time


def OP_metabolism_multidose(Y, t):
    G = Y[0]
    OP = Y[1]
    OC = Y[2]

    dGdt = -1 * ka * G
    dOPdt = ka * G - kf * OP
    dOCdt = kf * OP - ke * OC

    return [dGdt, dOPdt, dOCdt]


def plot(G, OP, OC, t):
    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot(1, 1, 1)
    # Plot using red circles
    ax.plot(t, G, 'b-', label='Oral OP Cocentration', markevery=10)
    ax.plot(t, OP, 'r-', label='Plasma OP Concentration', markevery=10)
    ax.plot(t, OC, 'g-', label='Plasma OC Concentration', markevery=10)

    # ax.plot(t, go, 'g-', label='Other Conductance', markevery=10)

    # Set labels and turn grid on
    ax.set(xlabel='Time $t$, hrs', ylabel=r'Concentration', title='Tamiflu Single Dose')
    ax.grid(True)
    ax.legend(loc='best')
    # Use space most effectively
    fig.tight_layout()
    # Save as a PNG file
    fig.savefig('Oseltamivir_Metabolism.png')


def dosing(tdose, num_doses):
    dose = 1  # Indexing from 1 feels dirty, but its the way that I have to store time data
    G = []  # Outer variable storage for overall concentrations
    OP = []
    OC = []
    T = []
    Yo = [75, 0, 0]  # initial conditions of G, dose, is 75mg and OP and OC are both 0
    t = np.arange(0, tdose, .01)  # 12 hours, iterate by the minute
    while dose <= num_doses:
        # This function runs a new ODE for each dose, saving the final condition of one dose
        # and feeding it in as the initial condition of the next dose
        # lower case letters represent concentration after a given dose, upper case are total concentrations

        out = odeint(OP_metabolism_multidose, Yo, t)
        g = out[:, 0]
        op = out[:, 1]
        oc = out[:, 2]
        G.extend(g)
        OP.extend(op)
        OC.extend(oc)

        T.extend(t + dose * tdose)
        dose += 1
        added_dose = g[len(g) - 1] + 75
        Yo = [added_dose, op[len(op) - 1], oc[len(oc) - 1]]

    print(G)
    print(OP)
    print(OC)
    plot(G, OP, OC, T)


dosing(12, 7)  # Function call for the dosing equation, frequency of dosing, number of doses
