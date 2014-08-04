import numpy as np 
#from numpy.linalg import solve
#from scipy.optimize import fsolve
#from numpy.linalg import inv

import copy

import pprint

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')

class BaroReceptor(object):
    def __init__(self, dt, n, Tsteps, vesselNi, A):
        '''
        Baroreceptor model
        
        '''
        #System and Vessel Variables
        self.dt = dt
        self.n = n
        self.Tsteps = Tsteps
        
        self.vesselNi = vesselNi
        self.A = A
        
        self.A0 = copy.deepcopy(A[vesselNi])
        
        # constants of baroreceptor
        self.phi0  = 0
        self.F0    = 0 #called L0 in article
        self.g     = 0
        self.alpha = 0
        self.tau1  = 0
        self.tau2  = 0
        
        self.F1 = F0
        self.F2 = F0
        
    def __call__(self):
        
        # create local variables for this timestep
        dt = self.dt
        n = self.n[0]
        
        # 1. calculate strain
        An = self.A[n][self.vesselNi]
        A0 = self.A0
        epsilon = (An-A0)/A0
        
        ## 2. solve F1,F2  (== L1 L2) functions
        # explicit euler
        self.F1 = self.F1 + dt/self.tau1 (epsilon-self.F1)
        self.F2 = self.F2 + dt/self.tau2 (epsilon-self.F2)
        
        # 3. calculate L(L1,L2)
        F = (self.alpha * self.F1 - self.F2)*(alpha-1)
        
        # 4. calculate firing rate phi
        phi = self.phi0 + self.g*(F-self.F0)* 0.5 * (np.sign((F-self.F0)) + 1)
        
        # 5. calcualte phi_sn phi_pn
        
        ## should go into the classBoundaryConditions.py
        # 6. (adrelanin ODE) phi_sn -> calculate heart rate // apply to boundaryConditions
        
        