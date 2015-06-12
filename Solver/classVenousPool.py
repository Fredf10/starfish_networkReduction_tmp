from NetworkLib.classBoundaryConditions import *
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
    
    def __init__(self,dataDict):
        
        self.dt = 0 #will be updated with update method
        self.currentTimeStep = 0 # current time step
        self.currentMemoryIndex = 0
        self.nTsteps = 0
        
        self.boundarys = 0 # make it a dictionary/needs to be initialized in FlowSolver
        
        self.update(dataDict)
        
        self.V = 5600.0e-6 * 0.61 # estimated blood volume on venous side under normal conditions
        self.Vusv0 = 3213e-6 # unstretched volume at reference state
        self.P0 = 2.0 * 133.322368 # pressure constant for calculation of P venous
        
        self.k = 0.1124 #0.1124e-9 # constant
        
        self.pressureGain = 1.0/0.228 # pressure gain between CVP and LAP - Bell paper
        
        self.Vusv = 2335.46e-6 #3037.0e-6 # initial states for Vus, CVP and LAP
        self.P = self.P0*(math.exp(self.k*math.pow((self.V*1e6-self.Vusv*1e6),1.5)/(self.V*1e6)))
        self.P_LA = self.pressureGain*self.P
        
        self.Qin = 0.0 # in and outflow to the venous pool
        self.Qout = 0.0
        
        ### vectors for export
        self.Vvector = np.ones(self.nTsteps+1)*self.V
        self.Pvector = np.ones(self.nTsteps+1)*self.P
        self.P_LAvector = np.ones(self.nTsteps+1)*self.P_LA
        self.Vusvvector = np.ones(self.nTsteps+1)*self.Vusv
        
            
    def estimateInflow(self):
        """
        calculate the inflow to the venous side, from terminal boundaries
        """
        n = self.currentTimeStep[0]
        Qin = 0
        
        for key in self.boundarys:
            for x in range(len(self.boundarys[key])):
                if self.boundarys[key][x].position == -1: #distal boundaries    
                    Qin = Qin + self.boundarys[key][x].Q[self.currentMemoryIndex,-1]

        self.Qin = Qin
        
        
    
    def estimateOutflow(self):
        """
        calculate outflow
        """
        n = self.currentTimeStep[0]
        Qout = 0
        
        for key in self.boundarys:
            for x in range(len(self.boundarys[key])):
                if self.boundarys[key][x].position == 0: #proximal boundaries
                    
                    if self.boundarys[key][x].bcType2[0].name == 'VaryingElastanceHeart':     
                        Qout = Qout + self.boundarys[key][x].bcType2[0].mitralQ[n]
                        #print "VP 78"
                        #print self.boundarys[key][x].bcType2[0].mitralQ[n]
                        
                    else:
                        Qout = Qout + self.boundarys[key][x].Q[self.currentMemoryIndex,0]
                        
                        
        self.Qout = Qout
        
    
    def updateVenousPool(self):
        """
        update the state of the venous pool (volume and pressure)
        the model is from Ursino_1999
        """
        self.V = self.V + self.dt*(+self.Qin - self.Qout)
        self.P = self.P0*(math.exp(self.k*math.pow((self.V[0]*1e6-self.Vusv*1e6),1.5)/(self.V[0]*1e6)))
        self.P_LA = self.pressureGain*self.P
        
        n = self.currentTimeStep[0]
        
        self.Vvector[n+1] = self.V
        self.Pvector[n+1] = self.P
        self.P_LAvector[n+1] = self.P_LA
        self.Vusvvector[n+1] = self.Vusv
        
        
        
    def updateBoundaryConditions(self):
        """
        update all the boundary conditions
        """
        n = self.currentTimeStep[0]
        
        for key in self.boundarys:
            for x in range(len(self.boundarys[key])):
                
                if self.boundarys[key][x].position == 0:
                    if self.boundarys[key][x].type == 'VaryingElastanceHeart':
                        self.boundarys[key][x].atriumPressure = self.P_LA
                        
                    else: pass

                if self.boundarys[key][x].position == -1:
                    if self.boundarys[key][x].type == 'Windkessel-2Elements':
                        self.boundarys[key][x].bcType2[0].venousPressure[n+1] = self.boundarys[key][x].bcType2[0].venousPressure[n+1]-self.boundarys[key][x].bcType2[0].venousPressure[0]+self.P
                            
                    if self.boundarys[key][x].type == 'Windkessel-3Elements':
                        if n < (np.size(self.boundarys[key][x].bcType2[0].venousPressure)-1):
                            self.boundarys[key][x].bcType2[0].venousPressure[n+1] = self.boundarys[key][x].bcType2[0].venousPressure[n+1]-self.boundarys[key][x].bcType2[0].venousPressure[0]+self.P
                    
                    if self.boundarys[key][x].type == 'Resistance':
                        self.boundarys[key][x].bcType2[0].venousPressure[n+1] = self.boundarys[key][x].bcType2[0].venousPressure[n+1]-self.boundarys[key][x].bcType2[0].venousPressure[0]+self.P
                
                
             
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
            
            
            
            
    
    
