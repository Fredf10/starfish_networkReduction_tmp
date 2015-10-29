import sys, os
#from NetworkLib.classBoundaryConditions import *
### should work without this, but won't delete yet just in case.
import numpy as np
import math

cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')
import UtilityLib.classStarfishBaseObject as cSBO

class StaticVenousPool(cSBO.StarfishBaseObject):
    """
    A venous pool model with only fixed values. Mimics the venousPool class
    Very simple model of the venous side, considering the veins as one big compliant reservoir,
    and assuming a pure pressure gain between CVP and LAP
    The Baroreflex regulates the unstretched volume of the venous side, through which the CVP and ultimately
    the LAP are changed
    
    self.V is the blood volume in the veins
    self.Vusv: unstretched Volume of Veins with zero external pressure
    self.P0: constant
    self.k: constant
    self.pressureGain: pressure gain from CVP (right atrial) to LAP (left atrial) --> a pure gain is used, value according to Bell
    self.P: Central Venouse Pressure i.e. right atrial pressure
    self.P_LA: pressure in left atrium (LAP)
    self.Qin: inflow
    self.Qout: outflow
    """
    
    def __init__(self,dataDict):
        self.dt = 0 #will be updated with update method
        self.currentTimeStep = 0 # current time step
        self.currentMemoryIndex = 0
        self.nTsteps = 0
        
        self.boundarys = {} # make it a dictionary/needs to be initialized in FlowSolver
        
        self.update(dataDict)
        self.veinId  = 0
        self.V = 3770.4477485970647e-6 # estimated blood volume on venous side under normal conditions
        self.Vusv0 = 3400e-6#  ? 3213e-6 # unstretched volume at reference state
        self.P0 = 2.0 * 133.322368 # pressure constant for calculation of P venous
        
        self.k = 0.1124 #0.1124e-9 # constant
        
        self.pressureGain = 3. #1.0/0.228 # pressure gain between CVP and LAP - Bell paper
        
        self.Vusv = 3400e-6 #3037.0e-6 # initial states for Vus, CVP and LAP
        self.P = self.P0*(math.exp(self.k*math.pow((self.V*1e6-self.Vusv*1e6),1.5)/(self.V*1e6)))
        self.P_LA = self.pressureGain*self.P
        
        self.Qin = 0.0 # in and outflow to the venous pool
        self.Qout = 0.0
        
        ### vectors for export
        self.Vvector = np.ones(self.nTsteps+1)*self.V
        self.Pvector = np.ones(self.nTsteps+1)*self.P
        self.P_LAvector = np.ones(self.nTsteps+1)*self.P_LA
        self.Vusvvector = np.ones(self.nTsteps+1)*self.Vusv
    
    def initializeForSimulation(self,flowSolver, vascularNetwork):
        
        self.dsetGroup               = vascularNetwork.solutionDataFile.create_group('Venous')
        self.currentTimeStep         = flowSolver.currentTimeStep
        self.currentMemoryIndex      = flowSolver.currentMemoryIndex
        self.dt                      = flowSolver.dt
        self.nTsteps                 = flowSolver.nTsteps
        self.boundarys               = flowSolver.boundarys

        self.P
        self.P_LA
        self.Qin
        self.Qout
                ### vectors for export
        self.Vusvvector = np.ones(self.nTsteps+1)*self.Vusv
        self.dsetGroup.create_dataset("Vus", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Vvector = np.ones(self.nTsteps+1)*self.V    
        self.dsetGroup.create_dataset("V", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Pvector = np.ones(self.nTsteps+1)*self.P
        self.dsetGroup.create_dataset("CVP", (vascularNetwork.savedArraySize,), dtype='float64')
        self.P_LAvector = np.ones(self.nTsteps+1)*self.P_LA
        self.dsetGroup.create_dataset("LAP", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Qin_vector = np.zeros(self.nTsteps+1)
        self.dsetGroup.create_dataset("Qin", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Qout_vector = np.zeros(self.nTsteps+1)
        self.dsetGroup.create_dataset("Qout", (vascularNetwork.savedArraySize,), dtype='float64')
        
    def flushSolutionData(self,saving, nDB,nDE,nSB,nSE,nSkip):

        if saving:
            self.dsetGroup["V"][nDB:nDE]  = self.Vvector[nSB:nSE:nSkip]
            self.dsetGroup["Vus"][nDB:nDE]  = self.Vusvvector[nSB:nSE:nSkip]
            self.dsetGroup["CVP"][nDB:nDE]  = self.Pvector[nSB:nSE:nSkip]
            self.dsetGroup["LAP"][nDB:nDE] = self.P_LAvector[nSB:nSE:nSkip]
            self.dsetGroup["Qin"][nDB:nDE] = self.Qin_vector[nSB:nSE:nSkip]
            self.dsetGroup["Qout"][nDB:nDE] = self.Qout_vector[nSB:nSE:nSkip]


             
    def __call__(self):
        """
        call function for Venous pool, updates the venous pool volume and pressure
        and updates the boundary conditions with new pressure values
        """
        pass
    
        
    def update(self,dataDict):
        """
        updates the data
        Dict = {'variableName': value}
        """
        for key,value in dataDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except Exception:
                self.warning("StaticVenousPool.update(): wrong key: %s, could not set up venousPool" %key)

class venousPool(cSBO.StarfishBaseObject):
    """
    Very simple model of the venous side, considering the veins as one big compliant reservoir,
    and assuming a pure pressure gain between CVP and LAP
    The Baroreflex regulates the unstretched volume of the venous side, through which the CVP and ultimately
    the LAP are changed
    
    self.V is the blood volume in the veins
    self.Vusv: unstretched Volume of Veins with zero external pressure
    self.P0: constant
    self.k: constant
    self.pressureGain: pressure gain from CVP (right atrial) to LAP (left atrial) --> a pure gain is used, value according to Bell
    self.P: Central Venouse Pressure i.e. right atrial pressure
    self.P_LA: pressure in left atrium (LAP)
    self.Qin: inflow
    self.Qout: outflow
    """
    
    def __init__(self,dataDict):
        
        self.dt = 0 #will be updated with update method
        self.currentTimeStep = 0 # current time step
        self.currentMemoryIndex = 0
        self.nTsteps = 0
        
        self.boundarys = {} # make it a dictionary/needs to be initialized in FlowSolver
        
        self.update(dataDict)
        self.veinId  = 0
       

        '''
        ### FINDING INITIAL VALUES FOR VOLUME
        V = 3892. # 5600.0e-6 *0.61# * 0.61 # estimated blood volume on venous side under normal conditions
        Vusv0 = 2378. #  ? 3213e-6 # unstretched volume at reference state
        P0 = 2.0 # pressure constant for calculation of P venous    
        k = 0.1124 
              
        Vusv = Vusv0  
        def pFct(V):
            P = P0*(np.exp(k*(V-Vusv)**1.5/(V)))
            return P
            
        from scipy import optimize as opt
        fct = lambda V: pFct(V) - 3.0
        opt.brentq(fct, 2400,4000)
        # 2850.912397321067
        fct = lambda V: pFct(V) - 4.0
        opt.brentq(fct, 2400,4000)
        # 3091.6779241832583
        fct = lambda V: pFct(V) - 5.0
        opt.brentq(fct, 2400,4000)
        # 3270.4477485970647
        fct = lambda V: pFct(V) - 6.0
        opt.brentq(fct, 2400,4000)
        # 3414.6023352947177
        fct = lambda V: pFct(V) - 7.0
        opt.brentq(fct, 2400,4000)
        # 3536.118289840244
        fct = lambda V: pFct(V) - 8.0
        opt.brentq(fct, 2400,4000)
        # 3641.5166289669837
        '''

        # TODO figure out how to have BRX not need these as maximal values
        self.V = 3770.4477485970647e-6 # 3892e-6 # 5600.0e-6 *0.61# * 0.61 # estimated blood volume on venous side under normal conditions
        self.Vusv0 =  3400.e-6 # 2378e-6 #  ? 3213e-6 # unstretched volume at reference state
        
        
        self.P0 = 2.0 * 133.322368 # pressure constant for calculation of P venous
        
        self.k = 0.1124 #0.1124e-9 # constant
        
        self.pressureGain = 1.0/0.228 # pressure gain between CVP and LAP - Bell paper
        
        self.Vusv = self.Vusv0 # 2335.46e-6 #3037.0e-6 # initial states for Vus, CVP and LAP
        self.P = self.pressureFromVolume(self.V)
        # self.P0*(math.exp(self.k*math.pow((self.V*1e6-self.Vusv*1e6),1.5)/(self.V*1e6)))
        self.P_LA = self.pressureGain*self.P
        
        self.Qin = 0.0 # in and outflow to the venous pool
        self.Qout = 0.0
        
        ### vectors for export
        self.Vvector = np.ones(self.nTsteps+1)*self.V
        self.Pvector = np.ones(self.nTsteps+1)*self.P
        self.P_LAvector = np.ones(self.nTsteps+1)*self.P_LA
        self.Vusvvector = np.ones(self.nTsteps+1)*self.Vusv
    
    def initializeForSimulation(self,flowSolver, vascularNetwork):
        
        self.dsetGroup               = vascularNetwork.solutionDataFile.create_group('Venous')
        self.currentTimeStep         = flowSolver.currentTimeStep
        self.currentMemoryIndex      = flowSolver.currentMemoryIndex
        self.dt                      = flowSolver.dt
        self.nTsteps                 = flowSolver.nTsteps
        self.boundarys               = flowSolver.boundarys
        self.boundaryCondtions = vascularNetwork.boundaryConditions

        ### vectors for export
        self.Vusvvector = np.ones(self.nTsteps+1)*self.Vusv
        self.dsetGroup.create_dataset("Vus", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Vvector = np.ones(self.nTsteps+1)*self.V    
        self.dsetGroup.create_dataset("V", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Pvector = np.ones(self.nTsteps+1)*self.P
        self.dsetGroup.create_dataset("CVP", (vascularNetwork.savedArraySize,), dtype='float64')
        self.P_LAvector = np.ones(self.nTsteps+1)*self.P_LA
        self.dsetGroup.create_dataset("LAP", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Qin_vector = np.zeros(self.nTsteps+1)
        self.dsetGroup.create_dataset("Qin", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Qout_vector = np.zeros(self.nTsteps+1)
        self.dsetGroup.create_dataset("Qout", (vascularNetwork.savedArraySize,), dtype='float64')
        
#         print "DEBUG"
#         print "self.Vusv0", self.Vusv0
#         print "self.Vusv", self.Vusv
#         print "self.V", self.V
#         print "self.P", self.P
#         
#         exit()
        
    def flushSolutionData(self,saving, nDB,nDE,nSB,nSE,nSkip):

        if saving:
            self.dsetGroup["V"][nDB:nDE]  = self.Vvector[nSB:nSE:nSkip]
            self.dsetGroup["Vus"][nDB:nDE]  = self.Vusvvector[nSB:nSE:nSkip]
            self.dsetGroup["CVP"][nDB:nDE]  = self.Pvector[nSB:nSE:nSkip]
            self.dsetGroup["LAP"][nDB:nDE] = self.P_LAvector[nSB:nSE:nSkip]
            self.dsetGroup["Qin"][nDB:nDE] = self.Qin_vector[nSB:nSE:nSkip]
            self.dsetGroup["Qout"][nDB:nDE] = self.Qout_vector[nSB:nSE:nSkip]
            
    def estimateInflow(self):
        """
        calculate the inflow to the venous side, from terminal boundaries
        """
        
        Qin = 0
        
        for key in self.boundarys:
            for x in range(len(self.boundarys[key])):
                if self.boundarys[key][x].position == -1: #distal boundaries    
                    bcCondition = self.boundarys[key][x].bcType2[0]
                    # TODO: Add other types
                    if bcCondition.name == 'Windkessel-3Elements':                    
                        deltaP = self.boundarys[key][x].P[self.currentMemoryIndex,-1] - bcCondition.venousPressure[self.currentMemoryIndex]
                        Qin = Qin + deltaP/bcCondition.Rtotal
                        
                    

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
        
    def pressureFromVolume(self,V):
        return self.P0*(np.exp(self.k*(V*1e6-self.Vusv*1e6)**1.5/(V*1e6)))
        
    def updateVenousPool(self):
        """
        update the state of the venous pool (volume and pressure)
        the model is from Ursino_1999
        """
        self.V = self.V + self.dt*(+self.Qin - self.Qout)
        try:
            # P = P0*(np.exp(k*(V*1e6-Vusv*1e6)**1.5/(V*1e6))-1.0)
            self.P = self.pressureFromVolume(self.V[0])
        except ValueError:
            if self.V[0] - self.Vusv < 0:
                self.exception(
                   "Venous volume {} is lower the unstressed venous volume {} at time {} and time step {}.".format(
                        self.V[0], self.Vusv, self.dt*self.currentTimeStep[0], self.currentTimeStep[0]))
            else:
                raise
                
        self.P_LA = self.pressureGain*self.P
        
        n = self.currentTimeStep[0]
        
        self.Vvector[n+1] = self.V
        self.Vusvvector[n+1] = self.Vusv
        self.Pvector[n+1] = self.P
        self.P_LAvector[n+1] = self.P_LA
        self.Vusvvector[n+1] = self.Vusv
        self.Qout_vector[n+1] = self.Qout
        self.Qin_vector[n+1] = self.Qin

        
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
                        #
                        # Replaces precalculated Pv(t) = Pcv_init + pgh(t) + -Pext by 
                        # Pv[n+1] = Pv[n+1] - Pv[0] + CVP[n+1],
                        # where Pv[n+1] - Pv[0] = Pcv_init + pgh[n+1] + -Pext - (Pcv_init + pgh[0] + -Pext)
                        # which leaves Pv[n+1] - Pv[0] = pgh[n+1] - pgh[0], and thus
                        # Pv[n+1] = CVP[n+1] + rho gh[n+1] - rho gh[0], which is correct so long as rho gh[0] = 0 and Pext = 0
                        # 
                        
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
        """
        updates the data
        Dict = {'variableName': value}
        """
        for key,value in dataDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except Exception:
                self.warning("venousPool.update(): wrong key: %s, could not set up venousPool" %key)
            
