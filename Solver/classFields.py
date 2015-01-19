import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')

import numpy as np

import time

from copy import copy as copy 

class Field(object):
    
    def __init__(self, vessel, currentMemoryIndex, dt, systemEquation, rigidArea):
        '''
        Constructor of Field object
        
        calculates the interior field of a vessel
        with a MackKormack Predictor-Corrector Schmea
        '''
        
        self.name = ' '.join(['Field',str(vessel.Id)])
        
        #System and Vessel Variables
        self.dz = vessel.dz
        self.dt = dt
        # current time index in memory
        self.currentMemoryIndex = currentMemoryIndex
        
        self.systemEquation = systemEquation
        self.AFunction = vessel.A
        self.rigidArea = rigidArea
        
        #SolutionVariables
        self.P = vessel.Psol
        self.P_pre = np.ones_like(vessel.Psol[0])
        self.Q = vessel.Qsol
        self.Q_pre = np.ones_like(vessel.Qsol[0])
        self.A = vessel.Asol
        self.A_pre = np.ones_like(vessel.Asol[0])
        
        self.step = "predictor"
        
    def __call__(self):
        '''
        Mack Kormac Predictor-Corrector
        '''
        # solve vessel objects
        dt = self.dt
        
        # the current position in solution memory
        n = self.currentMemoryIndex[0]
        
        P = self.P[n]
        Q = self.Q[n]
        A = self.A[n]
                
        # set bc values of predictor step arrays
        self.P_pre[0]   = P[0]
        self.P_pre[-1]  = P[-1] 
        self.Q_pre[0]   = Q[0]
        self.Q_pre[-1]  = Q[-1] 
        self.A_pre[0]   = A[0]
        self.A_pre[-1]  = A[-1]
        
        #''' Predictor Step '''            
        # update matrices               
        m12,m21,m22,b2 = self.systemEquation.updateSystem(P,Q,A)
         
        # create lokal variables for predictor varibales with
        P_pre = self.P_pre
        Q_pre = self.Q_pre 
        A_pre = self.A_pre
          
        dzPre = self.dz[1::]
        dzCor = self.dz[0:-1]
        
        # calculate derivatives forward 
        dPdz  = (P[2::] - P[1:-1])/dzPre
        dQdz  = (Q[2::] - Q[1:-1])/dzPre
        dQ2A  = pow(Q,2.)/A
        dQ2dz = (dQ2A[2::] - dQ2A[1:-1])/dzPre
        
        # solve vessel
        P_pre[1:-1] = (P[1:-1] - (m12*dQdz)*dt)
        Q_pre[1:-1] = (Q[1:-1] - (m21*dPdz + m22*dQ2dz - b2 )*dt)
         
        # check pressure solution
        if P_pre.any()<0:
            print "ERROR: {} calculated negativ pressure in predictor step at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print P_pre
            exit()
             
        # solve area
        if self.rigidArea == True:
            A_pre = A
        else:
            A_pre[1:-1] = self.AFunction(P_pre)[1:-1]        
                         
        #'''Corrector Step'''    
        # update matrices  
        m12,m21,m22,b2 = self.systemEquation.updateSystem(P_pre,Q_pre,A_pre)
         
        # calculate derivatives backward 
        dPdz  = (P_pre[1:-1] - P_pre[0:-2])/dzCor
        dQdz  = (Q_pre[1:-1] - Q_pre[0:-2])/dzCor
        dQ2A  = pow(Q_pre,2.)/A_pre
        dQ2dz = (dQ2A[1:-1] - dQ2A[0:-2])/dzCor
                 
        #solve vessel
        self.P[n+1][1:-1] = (P[1:-1] + P_pre[1:-1] - (m12*dQdz)*dt)/2.0
        self.Q[n+1][1:-1] = (Q[1:-1] + Q_pre[1:-1] - (m21*dPdz + m22*dQ2dz - b2 )*dt )/2.0
         
        # check pressure solution
        if self.P[n+1].any() < 0:
            print "ERROR: {} calculated negative pressure in corrector step at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print self.P[n+1]
            exit()
         
        # solve area
        if self.rigidArea == True:
            self.A[n+1] = A_pre
        else:
            self.A[n+1][1:-1] = self.AFunction(self.P[n+1])[1:-1]
