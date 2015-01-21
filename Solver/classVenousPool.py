from classBoundaryConditions import *
import numpy as np
import math


class venousPool(object):
    """
    Very simple model of the venous side, considering the veins as one big compliant reservoir,
    and assuming a pure pressure gain between CVP and LAP
    The Baroreflex regulates the unstretched volume of the venous side, through which the CVP and ultimately
    the LAP are changed
    
    self.V is the blood volume in the veins
    self.Vusv: unstretched Volume of Veins
    self.P0: constant
    self.k: constant
    self.pressureGain: pressure gain from CVP to LAP --> a pure gain is used, value according to Bell
    self.P: CVP
    self.P_LA: pressure in left atrium (LAP)
    self.Qin: inflow
    self.Qout: outflow
    """
    
    __init__(self):
        
        self.dt = 0 #will be updated with update method
        self.n = 0 # current time step
        
        self.boundarys = 0 # make it a dictionary/needs to be initialized in FlowSolver
        
        self.update()
        
        self.V = 5600.0 * 10e-6 * 0.61 # estimated blood volume on venous side under normal conditions
        self.Vusv0 = 3213 * 10e-6 # unstretched volume at reference state
        self.P0 = 2.0 * 133.322368 # pressure constant for calculation of P venous
        
        self.k = 0.1124*10e-9 # constant
        
        self.pressureGain = 1.0/0.228 # pressure gain between CVP and LAP - Bell paper
        
        self.Vusv = 2378.0 *10e-3 # initial states for Vus, CVP and LAP
        self.P = self.P0*(math.exp(self.k*math.pow((self.V-self.Vusv),1.5)/self.V))
        self.P_LA = self.pressureGain*self.P
        
        self.Qin = 0.0 # in and outflow to the venous pool
        self.Qout = 0.0
        
            
    def estimateInflow(self):
        """
        calculate the inflow to the venous side, from terminal boundary conditions
        """
        n = self.n[0]
        Qin = 0
        
        
        for boundaryID, boundary in self.boundarys:
            
            if boundary.position == -1: #distal boundaries    
                Qin = Qin + boundary.Q[n]
                
        self.Qin = Qin
        
        
    
    def estimateOutflow
        """
        calculate outflow
        """
        n = self.n[0]
        Qout = 0
        
        for boundaryID,boundary in self.boundarys:
            
            if boundary.position == 0: #proximal boundaries
                
                if boundary.bcType2[0].name == 'VaryingElastanceHeart':     
                    Qout = Qout + boundary.bcType2[0].mitralQ[n]
                
                else:
                    Qout = Qout + boundary.Q[n]
                           
        self.Qout = Qout
        
        
    
    def updateVenousPool(self):
        """
        update the state of the venous pool (volume and pressure)
        the model is from Ursino_1999
        """
        
        self.V = self.V + self.dt*(-self.Qin + self.Qout)
        self.P = self.P0*(math.exp(self.k*math.pow((self.V-self.Vusv),1.5)/self.V))
        self.P_LA = self.pressureGain*self.P
        
    
        
    def updateBoundaryConditions(self):
        """
        update all the boundary conditions
        """
        
        for boundaryID,boundary in self.boundarys:
            
            for bc in boundary.bcType2:
                if bc.name == 'VaryingElastanceHeart':
                    bc.atriumPressure = self.P_LA
                    
                elif bc.name == 'Windkessel-2Elements':
                    bc.venousPressure = self.P
                    
                elif bc.name == 'Windkessel-3Elements':
                    bc.venousPressure = self.P
                    
                elif bc.name == 'Resistance':
                    bc.venousPressure = self.P
                
                
             
    def __call__(self):
        """
        call function for Venous pool, updates the venous pool volume and pressure
        and updates the boundary conditions with new pressure values
        """
        self.estimateInflow()
        self.estimateOutflow()
        self.updateVenousPool()
        self.updateBoundaryConditions()
        
    
        
    def update(self,dataDict):
        '''
        updates the data
        Dict = {'variableName': value}
        '''
        for key,value in dataDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except: 
                print 'ERROR venousPool.update(): wrong key: %s, could not set up venousPool' %key    
            
            
            
            
    
    