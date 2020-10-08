import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

ka = 1.011  # OP blood absorption rate
kf = 0.684  # OP to OC conversion rate
ke = 0.136  # OC elimintation rate

Go = 75  # initial dose of 75mg
OPo = Go * ka

tdose = 12  # time between doses in hours
dt = .1  # dose time

To = 4e8  # initial number of epithelial cells in upper respiratory
Vo = 1e4  # initial free virus
Io = 0  # initial infected cells

gt = .03333  # Basal growth rate of healthy cells, h^-1
InfectionRate = 1.253e-9  # Rate of infection of target cells by virus, (EID50/mL)-1 h^-1
InfectedDeathRate = .0833333  # h^-1
MaxAntiviralEffect = .98

EC50 = 30  # Concentration of drug reaching half-maximal effect,
ViralProdRate = 8.75  # *EID50/mL) cell^-1 h^-1
ViralClearanceRate = .2167  # Nonspecific viral clearance h^-1
ViralConsumption = 2.08333e-8  # Rate of viral consumption by binding to target cells, cell^-1 h^-1


# T is target cells
# I is infected cells
# V is free virus


def OP_metabolism_multidose(Y, t):
    G = Y[0]
    OP = Y[1]
    OC = Y[2]

    dGdt = -1 * ka * G
    dOPdt = ka * G - kf * OP
    dOCdt = kf * OP - ke * OC

    T = Y[3]
    I = Y[4]
    V = Y[5]

    urtOC = .95 * OC  # for the below calculations, OC is taken to be in ng/ml, which is equal to ug/L
    AntiviralEffect = (MaxAntiviralEffect * urtOC) / (urtOC + EC50)
    dTdt = gt * T * (1 - (T + I) / To) - InfectionRate * V * T
    dIdt = InfectionRate * V * T - InfectedDeathRate * I
    dVdT = (1 - AntiviralEffect) * ViralProdRate * I - ViralClearanceRate * V - ViralConsumption * V * T

    return [dGdt, dOPdt, dOCdt, dTdt, dIdt, dVdT]


def plot(G, OP, OC, t, T, I, V, plot_name, forgotten_dose):
    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot(1, 1, 1)
    # Plot using red circles
    # ax.plot(t, G, 'b-', label='Oral OP Cocentration', markevery=10)
    ax.plot(t, OP, 'r-', label='Plasma OP Concentration', markevery=10)
    ax.plot(t, OC, 'g-', label='Plasma OC Concentration', markevery=10)

    if forgotten_dose:
        plt.axvline(forgotten_dose, color='r', linestyle='--')

    # Set labels and turn grid on
    ax.set(xlabel='Time $t$, hrs', ylabel=r'Concentration', title='Tamiflu Multiple Doses')
    ax.grid(True)
    ax.legend(loc='best')
    # Use space most effectively
    fig.tight_layout()
    # Save as a PNG file
    fig.savefig('Oseltamivir_Concentration_{}.png'.format(plot_name))

    # Create Figure for Viral Data

    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot(2, 1, 1)

    # Plot using red circles
    ax.plot(t, T, 'b-', label='Target Cells', markevery=10)
    ax.plot(t, I, 'r-', label='Infected Cells', markevery=10)
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(t, np.log(V), 'g-', label='Free Virus', markevery=10)

    # Set labels and turn grid on
    ax.set(xlabel='Time $t$, hrs', ylabel=r'Concentration', title='Viral Model')
    ax.grid(True)
    ax.legend(loc='best')

    if forgotten_dose:
        plt.axvline(forgotten_dose, color='r', linestyle='--')

    ax2.set(xlabel='Time $t$, hrs', ylabel=r'Concentration Log(V)')
    ax2.grid(True)
    ax2.legend(loc='best')
    # Use space most effectively
    plt.ylim([0, 25])  # set bounds on Y axis for free virus concentration
    fig.tight_layout()
    # Save as a PNG file
    fig.savefig('Viral_Model_{}.png'.format(plot_name))


def dosing(dose_size, tdose, num_doses, dose_delay, plot_name, forgotten_dose=0):
    dose = 1  # Indexing from 1 feels dirty, but its the way that I have to store time data
    G = []  # Outer variable storage for overall concentrations
    OP = []
    OC = []
    T = []
    I = []
    V = []

    time = []
    Yo = [0, 0, 0, To, Io, Vo]  # initial conditions of G, dose, is 75mg and OP and OC are both 0
    t = np.arange(0, tdose, .01)  # 12 hours, iterate by the minute
    while dose <= num_doses:

        # This function runs a new ODE for each dose, saving the final condition of one dose
        # and feeding it in as the initial condition of the next dose
        # lower case letters represent concentration after a given dose, upper case are total concentrations

        if not dose_delay:
            Yo[0] = dose_size  # if there's no delay for treatment, set the initial drug dose to dose size

        out = odeint(OP_metabolism_multidose, Yo, t)
        g = out[:, 0]
        op = out[:, 1]
        oc = out[:, 2]
        target_cells = out[:, 3]
        i = out[:, 4]
        v = out[:, 5]

        G.extend(g)
        OP.extend(op)
        OC.extend(oc)
        T.extend(target_cells)
        I.extend(i)
        V.extend(v)  # adding the individual run values to the overall lists

        time.extend(t + dose * tdose)
        dose += 1

        if dose == forgotten_dose or dose < dose_delay:
            added_dose = g[len(g) - 1]
        elif dose >= dose_delay:
            added_dose = g[len(g) - 1] + dose_size

        Yo = [added_dose, op[len(op) - 1], oc[len(oc) - 1], target_cells[len(target_cells) - 1],
              i[len(i) - 1], v[len(v) - 1]]
        # set initial conditions for the next loop equal to the last values from the previous loop

    """print(G)
    print(OP)
    print(OC)"""
    print(np.argmin(T))
    print(np.min(T))
    plot(G, OP, OC, time, T, I, V, plot_name, forgotten_dose * tdose)


def viral_growth(Y, t):
    T = Y[0]
    I = Y[1]
    V = Y[2]

    AntiviralEffect = 0
    dTdt = gt * T * (1 - (T + I) / To) - InfectionRate * V * T
    dIdt = InfectionRate * V * T - InfectedDeathRate * I
    dVdT = (1 - AntiviralEffect) * ViralProdRate * I - ViralClearanceRate * V - ViralConsumption * V * T

    return [dTdt, dIdt, dVdT]


'''
t = np.arange(0, 200, .01)

out = odeint(viral_growth, [To, Io, Vo], t)
print(out)
plot(t, t , t, t, out[:,0], out[:,1], out[:,2], 0)'''

dosing(dose_size=Go, tdose=tdose, num_doses=14, dose_delay=1, plot_name="14Dose_1DoseDelay_Dose5Forget_75mg", forgotten_dose=5)
# Function call for the dosing equation, dose size in mg, frequency of dosing, number of doses, number of doses
# missed after onset of infection

