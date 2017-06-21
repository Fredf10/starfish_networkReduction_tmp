from __future__ import print_function, absolute_import
from future.utils import iteritems, iterkeys, viewkeys, viewitems, itervalues, viewvalues
from builtins import input as input3
import sys,os
import numpy as np
import logging
logger = logging.getLogger('starfish')

class Field():
    
    def __init__(self, vessel, currentMemoryIndex, dt, rigidArea, solvingSchemeField = 'MacCormack_Flux'):
        """
        Constructor of Field object
        
        calculates the interior field of a vessel
        with a MackKormack Predictor-Corrector Schmea
        """
        
        self.name = ' '.join(['Field',str(vessel.Id)])
        
        self.vessel = vessel
        #System and Vessel Variables
        self.dz = vessel.dz
        self.dt = dt
        # current time index in memory
        self.currentMemoryIndex = currentMemoryIndex
        
        self.AFunction = vessel.A
        self.rigidArea = rigidArea
        
        #SolutionVariables
        self.P = vessel.Psol
        self.P_pre = vessel.Psol[0].copy()
        self.Q = vessel.Qsol
        self.Q_pre = vessel.Qsol[0].copy()
        self.A = vessel.Asol
        self.A_pre = np.ones_like(vessel.Asol[0])
        
        gamma = self.vessel.gamma
        self.alpha = 1 #(gamma+2.0)/(gamma+1.0) # Velocity profile correction, uncomment to activate
        self.step = "predictor"
                       
        if solvingSchemeField == "MacCormack_Flux":
            self._callfct = self.MacCormackFlux
        elif solvingSchemeField == "MacCormack_TwoStep":
            self._callfct = self.MacCormackTwoStep
            self.predict = self.MacCormackPredictor
            self.correct = self.MacCormackCorrector
        else:
            raise ValueError('Classfields51: error, scheme for solving field not correct')
    
    def __call__(self):
        return self._callfct()
        
    def F(self, u, A, C, Aconst, Cconst):
    
        """Flux based on conservative form of governing mass and momentum equations"""

        flux=np.zeros_like(u)
        p=u[0,:]
        q=u[1,:]
        rho = self.vessel.rho
        flux[0,:] = q/Cconst
        flux[1,:] = Aconst*p/rho + self.alpha*q**2/A
        return flux 
    
    def MacCormackFlux(self):
        """ This is an implementation of the MacCormack scheme as proposed in Fredrik Eikeland Fossans Master Thesis"""
        
        dt = self.dt
        
        # the current position in solution memory
        n = self.currentMemoryIndex[0]
        
        P = self.P[n]
        Q = self.Q[n]
        A = self.A[n]
        C = self.vessel.C(P)
        netGravity = self.vessel.netGravity[n]
                
        
        dx = self.dz[0] #currently only equidistant spacing is implemented and thus dz can be taken as float instead of vector
        
        rho = self.vessel.rho
        gamma = self.vessel.gamma
        
        my = self.vessel.my
        
        u = np.zeros((2, len(P))) #matrix for storing pressure and flow to caalculate fluxes
        
        u[0,:] = P
        u[1,:] = Q
        up=u.copy()
        
        Utemp1=u[1, :-1]/A[:-1]
        Aconst1 = A[:-1]
        Cconst1 = C[:-1]
        A1 = A[1:]
        A2 = A[:-1]
        C1 = C[1:]
        C2 = C[:-1]
        up[:,:-1] = u[:, :-1] - dt*(self.F(u[:, 1:], A1, C1, Aconst1, Cconst1) - self.F(u[:, :-1], A2, C2, Aconst1, Cconst1))/dx 
        
        A_grav = A[:-1]
        up[1,:-1] = up[1, :-1]-dt*2*(gamma + 2)*my*Utemp1*np.pi/rho +  dt*A_grav*netGravity
        
        if self.rigidArea == True:
            A_p = A
            C_p = C
        else:
            try:
                A_p = self.vessel.A(up[0, :])
                C_p = self.vessel.C(up[0, :])
            except FloatingPointError as E:
                print("Floating Point error in Field {}".format(self.name))
                raise E
        
        A_p1 = A_p[1:]
        A_p2 = A_p[:-1]
        C_p1 = C_p[1:]
        C_p2 = C_p[:-1]
        Utemp2=up[1, 1:]/A_p1
        Aconst2 = A_p[1:]
        Cconst2 = C_p[1:]
        u[:,1:] = .5*(u[:, 1:] + up[:, 1:] -  dt/dx*(self.F(up[:, 1:], A_p1, C_p1, Aconst2, Cconst2) - self.F(up[:, :-1], A_p2, C_p2, Aconst2, Cconst2)))

        u[1,1:] = u[1, 1:] - 0.5*dt*2*(gamma + 2)*my*Utemp2*np.pi/rho + 0.5*dt*A_p1*netGravity
        
        Pnewinterior = u[0, :]
        Qnewinterior = u[1, :]
        
        if self.rigidArea == True:
            Anewinterior = A
        else:
            Anewinterior = self.vessel.A(Pnewinterior)
        
        Pnewinterior = Pnewinterior[1:-1]
        Qnewinterior = Qnewinterior[1:-1]
        Anewinterior=Anewinterior[1:-1]
        
        self.P[n + 1][1:-1] = Pnewinterior
        self.Q[n + 1][1:-1] = Qnewinterior
        self.A[n + 1][1:-1] = Anewinterior

        #TODO: Please explain this if statement in a comment.
        if (self.P[n+1] < 0).any():
            raise ValueError("ERROR: {} calculated negative pressure in corrector step at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt))

    def MacCormackTwoStep(self):
        self.MacCormackPredictor()
        self.MacCormackCorrector()
    
    def MacCormackPredictor(self):
        """ TPredictor Step of MacCorMack scheme as proposed in Fredrik Eikeland Fossans Master Thesis"""
        
        dt = self.dt
        
        # the current position in solution memory
        n = self.currentMemoryIndex[0]
        
        P = self.P[n]
        Q = self.Q[n]
        A = self.A[n]
        C = self.vessel.C(P)
        netGravity = self.vessel.netGravity[n]
                
        
        dx = self.dz[0] #currently only equidistant spacing is implemented and thus dz can be taken as float instead of vector
        
        rho = self.vessel.rho
        gamma = self.vessel.gamma
        
        my = self.vessel.my
        
        u = np.zeros((2, len(P))) #matrix for storing pressure and flow to caalculate fluxes
        
        u[0,:] = P
        u[1,:] = Q
        
        self.u = u
        
        up = u.copy()
        
        Utemp1=u[1, :-1]/A[:-1]
        Aconst1 = A[:-1]
        Cconst1 = C[:-1]
        A1 = A[1:]
        A2 = A[:-1]
        C1 = C[1:]
        C2 = C[:-1]
        up[:, :-1] = u[:, :-1] - dt*(self.F(u[:, 1:], A1, C1, Aconst1, Cconst1) - self.F(u[:, :-1], A2, C2, Aconst1, Cconst1))/dx 
         
        A_grav = A[:-1]
        up[1,:-1] = up[1, :-1]-dt*2*(gamma + 2)*my*Utemp1*np.pi/rho +  dt*A_grav*netGravity
        
        
        self.up = up
        self.P_pre = up[0, :]
        self.Q_pre =up[1, :]
        
    def MacCormackCorrector(self):
        
        n = self.currentMemoryIndex[0]
        
        netGravity = self.vessel.netGravity[n]
        
        dt = self.dt
        dx = self.dz[0]
        
        rho = self.vessel.rho
        my = self.vessel.my
        gamma = self.vessel.gamma
        
        up = self.up
        u = self.u
        A_p = self.vessel.A(up[0, :])
        C_p = self.vessel.C(up[0, :])
        
        A_p1 = A_p[1:]
        A_p2 = A_p[:-1]
        C_p1 = C_p[1:]
        C_p2 = C_p[:-1]
        Utemp2=up[1, 1:]/A_p1
        Aconst2 = A_p[1:]
        Cconst2 = C_p[1:]
        u[:,1:] = .5*(u[:, 1:] + up[:, 1:] -  dt/dx*(self.F(up[:, 1:], A_p1, C_p1, Aconst2, Cconst2) - self.F(up[:, :-1], A_p2, C_p2, Aconst2, Cconst2)))

        u[1,1:] = u[1, 1:] - 0.5*dt*2*(gamma + 2)*my*Utemp2*np.pi/rho + 0.5*dt*A_p1*netGravity
        
        Pnewinterior = u[0, :]
        Qnewinterior = u[1, :]
        
        if self.rigidArea == True:
            Anewinterior = self.A[n]
        else:
            Anewinterior = self.vessel.A(Pnewinterior)
        
        Pnewinterior = Pnewinterior[1:-1]
        Qnewinterior = Qnewinterior[1:-1]
        Anewinterior=Anewinterior[1:-1]
        
        self.P[n + 1][1:-1] = Pnewinterior
        self.Q[n + 1][1:-1] = Qnewinterior
        self.A[n + 1][1:-1] = Anewinterior

        #TODO: Please explain this if statement in a comment.
        if (self.P[n+1] < 0).any():
            raise ValueError("ERROR: {} calculated negative pressure in corrector step at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt))
