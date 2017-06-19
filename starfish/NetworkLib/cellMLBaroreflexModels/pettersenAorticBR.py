# Prevent style checking of autogenerated code
# pylint: disable-all
from __future__ import print_function, absolute_import
from math import *
from numpy.core import *
from matplotlib.cbook import iterable
from scipy.optimize import fsolve
# Size of variable arrays:
sizeAlgebraic = 14
sizeStates = 6
sizeConstants = 32

inputID = 29 # ID of the input quantity, which is the wall strain calculated in STARFiSh
outputArray = "algebraic" # array in which the output (heart period) of the BR model is found
outputID = 13 # id of heart period in output array

#
n_ID = {'type':'alg', 'val': 3}
Tsym_ID = {'type':'alg', 'val' : 7}
Tparasym_ID = {'type':'alg', 'val': 8}
c_nor_ID = {'type':'state', 'val': 3}
c_ach_ID = {'type':'state', 'val': 4}


def createLegends():
    legend_y = [""] * sizeStates
    legend_ydot = [""] * sizeStates
    legend_alg = [""] * sizeAlgebraic
    legend_t = ""
    legend_p = [""] * sizeConstants
    legend_t = "Time in component Environmental (second)"
    legend_p[0] = "Tsmax in component Parameters (AU)"
    legend_p[1] = "Tpmax in component Parameters (AU)"
    legend_p[2] = "Tsmin in component Parameters (AU)"
    legend_p[3] = "Tpmin in component Parameters (AU)"
    legend_p[4] = "Gcns in component Parameters (dimensionless)"
    legend_p[5] = "Gs in component Parameters (per_Hertz)"
    legend_p[6] = "Gp in component Parameters (per_Hertz)"
    legend_p[7] = "tau_nor in component Parameters (second)"
    legend_p[8] = "tau_ach in component Parameters (second)"
    legend_p[9] = "tau_HR_nor in component Parameters (second)"
    legend_p[10] = "tau_HR_ach in component Parameters (second)"
    legend_p[11] = "HRo in component Parameters (Beats_per_min)"
    legend_p[12] = "HRmax in component Parameters (Beats_per_min)"
    legend_p[13] = "HRmin in component Parameters (Beats_per_min)"
    legend_p[14] = "Beta in component Parameters (dimensionless)"
    legend_p[15] = "delta_th in component Parameters (dimensionless)"
    legend_p[16] = "q_nor in component Parameters (per_s)"
    legend_p[17] = "q_ach in component Parameters (per_s)"
    legend_p[18] = "K_nor in component Parameters (AU)"
    legend_p[19] = "K_ach in component Parameters (AU)"
    legend_p[20] = "Gamma in component Parameters (dimensionless)"
    legend_p[21] = "alpha_rf in component Parameters (dimensionless)"
    legend_p[22] = "tau1 in component Parameters (second)"
    legend_p[23] = "tau2 in component Parameters (second)"
    legend_p[24] = "L0 in component Parameters (dimensionless)"
    legend_p[25] = "n0 in component Parameters (hertz)"
    legend_p[26] = "g in component Parameters (hertz)"
    legend_alg[7] = "alpha_cns in component Nervous_System (hertz)"
    legend_alg[5] = "n in component Receptive_Field (hertz)"
    legend_p[27] = "alpha_s0 in component Nervous_System (hertz)"
    legend_p[28] = "alpha_p0 in component Nervous_System (hertz)"
    legend_p[29] = "Eps_wall in component Receptive_Field (dimensionless)"
    legend_y[0] = "L1 in component Receptive_Field (dimensionless)"
    legend_y[1] = "L2 in component Receptive_Field (dimensionless)"
    legend_alg[0] = "L in component Receptive_Field (dimensionless)"
    legend_alg[3] = "F in component Receptive_Field (hertz)"
    legend_alg[9] = "Ts in component PNS_tones (AU)"
    legend_alg[10] = "Tp in component PNS_tones (AU)"
    legend_y[2] = "c_nor in component Norepinephrine (AU)"
    legend_y[3] = "C_ach in component Acetylcholine (AU)"
    legend_alg[1] = "delta_HR_ss in component Heart_Response_Nor (Beats_per_min)"
    legend_p[30] = "delta_HR_smax in component Heart_Response_Nor (Beats_per_min)"
    legend_y[4] = "delta_HR_s in component Heart_Response_Nor (Beats_per_min)"
    legend_alg[2] = "delta_HR_ps in component HR_ach (Beats_per_min)"
    legend_p[31] = "delta_HR_pmax in component HR_ach (Beats_per_min)"
    legend_alg[4] = "delta_HR_pfast in component HR_ach (Beats_per_min)"
    legend_y[5] = "delta_HR_pslow in component HR_ach (Beats_per_min)"
    legend_alg[6] = "delta_HR_p in component HR_ach (Beats_per_min)"
    legend_alg[12] = "HR in component HR_Combined (Beats_per_min)"
    legend_alg[11] = "HR_p in component HR_Combined (Beats_per_min)"
    legend_alg[8] = "HR_s in component HR_Combined (Beats_per_min)"
    legend_alg[13] = "Period in component HR_Combined (Sec_per_Beat)"
    legend_ydot[0] = "d/dt L1 in component Receptive_Field (dimensionless)"
    legend_ydot[1] = "d/dt L2 in component Receptive_Field (dimensionless)"
    legend_ydot[2] = "d/dt c_nor in component Norepinephrine (AU)"
    legend_ydot[3] = "d/dt C_ach in component Acetylcholine (AU)"
    legend_ydot[4] = "d/dt delta_HR_s in component Heart_Response_Nor (Beats_per_min)"
    legend_ydot[5] = "d/dt delta_HR_pslow in component HR_ach (Beats_per_min)"
    return legend_y, legend_alg, legend_t, legend_p

def initConsts(t=0, init_strain=0.01):
    p = zeros(sizeConstants)

    y = zeros(sizeStates)
    alg = zeros(sizeAlgebraic)
    p[0] = 4.12 #Tsmax
    p[1] = 4.994 #Tpmax
    p[2] = 0.5 #Tsmin
    p[3] = 1.6 # Tpmin
    p[4] = 1 # Gcns
    p[5] = 0.178 # Gs
    p[6] = 0.492 # Gp
    p[7] = 9.1 # tau_nor
    p[8] = 0.2 # tau_ach
    p[9] = 2.1 # tau_HR_nor
    p[10] = 2.5 # tau_HR_ach
    p[11] = 107 # intrinsic heart rate
    p[12] = 194 # max heart rate / could be adjusted for different ages
    p[13] = 50 #minimal heart rate / might be very low, might have to be increased
    p[14] = 0.175 # beta
    p[15] = 0 # delta_th
    p[16] = 0.1099 # q_nor
    p[17] = 5 # q_ach
    p[18] = 1.12 # K_nor
    p[19] = 0.65 # K_ach
    p[20] = 0.75 # Gamma
    p[21] = 2.5 # alfa_rf
    p[22] = 0.1 # tau1
    p[23] = 0.5 # tau2
    p[24] = 0.52 #0.32 #threshold value L0 values for young (39 years) and old group (75 years)
    p[25] = 30 # background firing rate, estimated with CellML simulation of Bugenhagen
    p[26] = 87.5 #100 #gain # both age groups again
    p[27] = 58.6 # alpha_s0
    p[28] = 76.019 # alpha_p0
    p[29] = init_strain # Epsilon_wall --> the input
    y[0] = 0.608 # inital value for L1
    y[1] = 0.626 # initial value for L2
    y[2] = 2.407 #  1.5672 # norepinephrine release initial state
    y[3] = 1.95 # 1.60 # acetylcholine release initial state
    y[4] = 71.205 # 57.587 # delta_HR_s initial value
    y[5] = 13.24  # 12.2467 # delta_HR_pslow initial value
    p[30] = p[12]-p[11]
    p[31] = p[11]-p[13]
    # y = steadyStates(y, p, t)
    return y, p

def steadyStates(y,p,t):
    alg = computeAlgebraic(p, y, t)
    y[0] = p[29] # inital value for L1
    y[1] = p[29] # initial value for L2
    y[2] = p[16]*alg[9]*p[7]
    y[3] = p[17]*alg[10]*p[8]
    y[4] = alg[1]
    y[5] = (1-p[20])*alg[2]
    return y
    
def computeRates(t, y, p):
    ydot = zeros(sizeStates)
    alg = zeros(sizeAlgebraic)
    ydot[0] = (p[29]-y[0])/p[22]
    ydot[1] = (p[29]-y[1])/p[23]
    alg[1] = (p[30]*(power(y[2], 2.00000)))/(power(p[18], 2.00000)+power(y[2], 2.00000))
    ydot[4] = (-y[4]+alg[1])/p[9]
    alg[2] = (p[31]*(power(y[3], 2.00000)))/(power(p[19], 2.00000)+power(y[3], 2.00000))
    ydot[5] = (-y[5]+(1.00000-p[20])*alg[2])/p[10]
    alg[0] = (p[21]*y[0]-y[1])/(p[21]-1.00000)
    alg[3] =  (p[26]*(alg[0]-p[24]) if (alg[0] >= p[24]) else  0.00000 if (alg[0] < p[24]) else nan)
    alg[5] = p[25]+alg[3]
    alg[7] = p[4]*alg[5]
    alg[9] = p[2]+(p[0]-p[2])/(exp(p[5]*(alg[7]-p[27]))+1.00000)
    ydot[2] = -(y[2]/p[7])+p[16]*alg[9]
    alg[10] = p[3]+(p[1]-p[3])/(exp(-p[6]*(alg[7]-p[28]))+1.00000)
    ydot[3] = -(y[3]/p[8])+p[17]*alg[10]
    return ydot

def computeAlgebraic(p, y, t):
    alg = zeros(sizeAlgebraic)
    alg[1] = (p[30]*(power(y[2], 2.00000)))/(power(p[18], 2.00000)+power(y[2], 2.00000))
    alg[2] = (p[31]*(power(y[3], 2.00000)))/(power(p[19], 2.00000)+power(y[3], 2.00000))
    alg[0] = (p[21]*y[0]-y[1])/(p[21]-1.00000)
    alg[3] =  (p[26]*(alg[0]-p[24]) if (alg[0] >= p[24]) else  0.00000 if (alg[0] < p[24]) else nan)
    alg[5] = p[25]+alg[3]
    alg[7] = p[4]*alg[5]
    alg[9] = p[2]+(p[0]-p[2])/(exp(p[5]*(alg[7]-p[27]))+1.00000)
    alg[10] = p[3]+(p[1]-p[3])/(exp(-p[6]*(alg[7]-p[28]))+1.00000)
    alg[4] = p[20]*alg[2]
    alg[6] = alg[4]+y[5]
    alg[8] = p[11]+y[4]
    alg[11] = p[11]-alg[6]
    alg[12] = alg[11]+((alg[8]-p[11])*(alg[11]-p[14]*p[13]))/(p[11]-p[14]*p[13])
    alg[13] = 60.0000/alg[12]
    return alg

def solve_model():
    """Solve model with ODE solver"""
    from scipy.integrate import ode
    # Initialise constants and state variables
    (init_states, constants) = initConsts()

    # Set timespan to solve over
    voi = linspace(0, 100, 1000)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = zeros((len(voi), sizeStates))
    algebraic = zeros((len(voi), sizeAlgebraic))
    states[0] = init_states
    for i, t in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[i + 1] = r.y
        else:
            break

    # Compute algebraic variables
    for i, (s, v) in enumerate(zip(states, voi)):
        algebraic[i] = computeAlgebraic(constants, s, v)
    return (voi, states, algebraic)

def plot_model(voi, states, algebraic):
    """Plot variables against variable of integration"""
    import matplotlib.pyplot as plt
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()
    plt.figure()
    plt.plot(voi, states)
    plt.plot(voi, algebraic)
    plt.xlabel(legend_voi)
    plt.legend(legend_states + legend_algebraic, loc='best')
    plt.show()

### BEGIN added for cgptoolbox

# @todo: The following module-level variables are shared across instances.
#        It might be better to wrap them in a class, allowing each instance of 
#        the same model to have its own parameter vector.

import sys
import numpy as np

ftype = np.float64 # explicit type declaration, can be used with cython
y0 = np.zeros(sizeStates, dtype=ftype)
ydot = np.zeros(sizeStates, dtype=ftype)
p = np.zeros(sizeConstants, dtype=ftype)
algebraic = alg = np.zeros(sizeAlgebraic, dtype=ftype)

y0[:], p[:] = initConsts()

# Sundials calling convention: https://computation.llnl.gov/casc/sundials/documentation/cv_guide/node6.html#SECTION00661000000000000000

def ode(t, y, ydot, f_data):
    """
    Compute rates of change for differential equation model.
    
    Rates are written into ydot[:]. 
    f_data is ignored, but required by the CVODE interface.
    
    The function returns 0 on success and -1 on failure.
    
    >>> ode(None, None, None, None)
    -1
    
    For debugging in case of failure, exception info is stored in the 
    module-level variable exc_info. (The message ends in "unsubscriptable" 
    under Python 2.6 but "not subscriptable" under Python 2.7, hence the 
    ellipsis.) Unfortunately, this is currently not implemented in a compiled 
    ODE. It will check the type of arguments before executing, but I am not 
    sure what happens in case of run-time errors inside the ODE.
    
    >>> exc_info
    (<type 'exceptions.TypeError'>,
    TypeError("'NoneType' object is ...subscriptable",),
    <traceback object at 0x...>)
    """
    global traceback
    traceback = None
    try:
        ydot[:] = computeRates(t, y, p)
        return 0
    except StandardError:
        import traceback
        ode.traceback = traceback.format_exc()
        return -1

def rates_and_algebraic(t, y):
    """
    Compute rates and algebraic variables for a given state trajectory.
    
    Unfortunately, the CVODE machinery does not offer a way to return rates and 
    algebraic variables during integration. This function re-computes the rates 
    and algebraics at each time step for the given state.
    
    This returns a simple float array; 
    :meth:`cgp.physmod.cellmlmodel.Cellmlmodel.rates_and_algebraic`
    will cast them to structured arrays with named fields.
    
    This version is pure Python; 
    :func:`~cgp.physmod.cythonize.cythonize`
    will generate a faster version.
    """
    t = np.atleast_1d(t)
    imax = len(t)
    # y can be NVector, unstructured or structured Numpy array.
    # If y is NVector, its data will get copied into a Numpy array.
    y = np.array(y).view(float)
    ydot = np.zeros_like(y)
    a = np.zeros((imax, len(alg)))
    for i in range(imax):
        ydot[i] = computeRates(t[i], y[i], p)
        if len(alg):
            # need np.atleast_1d() because computeAlgebraic() uses len(t)
            a[i] = computeAlgebraic(p, y[i], np.atleast_1d(t[i])).squeeze()
    return ydot, a

### END added for cgptoolbox

def solver2(timeArray,initial_states,constants):
    """
    Solve model with ODE solver
    solver function defined for the use with STARFiSh
    """
    
    from scipy.integrate import ode
    # Initialise constants and state variables
    #(init_states, constants) = initConsts()

    # Set timespan to solve over
    voi = timeArray
    max_step = timeArray[1]-timeArray[0]
    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=max_step)
    r.set_initial_value(initial_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = zeros((len(voi), sizeStates))
    algebraic = zeros((len(voi), sizeAlgebraic))
    states[0] = initial_states
    for i, t in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[i + 1] = r.y
        else:
            break

    # Compute algebraic variables
    for i, (s, v) in enumerate(zip(states, voi)):
        algebraic[i] = computeAlgebraic(constants, s, v)
    return (voi, states, algebraic)


if __name__ == "__main__":
    
    (states, constants) = initConsts()
    timeArray = linspace(0, 1, 2500)
    strain = 0.3419557616+0.0671788693*sin((78.9000/60.0000*6.28)*timeArray)
    timeArray2 = linspace(0,0.0004,2)
    constants[29] = strain[0]
    voi, states, algebraic = solver2(timeArray2,states,constants)
    
    HR = np.zeros(np.size(strain))
    HR[0] = algebraic[-1][12]
    
    for t in xrange((size(strain)-1)):
        
        constants[29] = strain[t+1]
        voi, states, algebraic = solver2(timeArray2,states[-1],constants)
        HR[t+1] = algebraic[-1][12]
        
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(timeArray, HR)
    #plt.plot(timeArray, strain)
    plt.show()
    
    file = open("solv","w")
    np.save(file,HR)
    file = open("solvTime","w")
    np.save(file,timeArray)
    file = open("solvStrain","w")
    np.save(file,strain)
    
