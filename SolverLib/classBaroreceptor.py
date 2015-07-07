import numpy as np 
import math

# import cellMLBaroreflexModels # Used later on in classes
import cellMLBaroreflexModels.pettersenAorticBR
import cellMLBaroreflexModels.bugenhagenAorticBR

class Baroreceptor(object):
    """
    Mother class for all baroreceptor models
    """
    
    def __init__(self, BaroDict):
        '''
        Baroreceptor model initialisation
        '''
        ## Solver related variables
        #System and Vessel Variables
        self.dt = 0
        self.currentTimeStep = 0
        self.currentMemoryIndex = 0
        self.nTsteps = 0
        
        # Configuration and solution data variables
        self.modelName = ''
        self.baroId = None
            
        self.boundaryCondition = 0
        self.boundaryConditionII = 0
        self.boundaryConditionIIout = {}
        
        self.newCycles = None
        self.Pheart = None
        self.Vheart = None
        self.venousPool = None
        
        self.dsetGroup = None
        self.Venous_dsetGroup = None

        # update with BaroDict --> see also the initialize method in FlowSolver
        self.update(BaroDict)
                
    def updateConstants(self,constants):
        """
        Method for extending to specific child classes to allow modification of cellML parameters and constants
        """
        self.constants = constants
                    
    def initializeForSimulation(self,flowSolver, vascularNetwork): 
        """
        Configures class members for simulation, using data from the network 
        which may not have been available at construction.
        """
        self.currentTimeStep         = flowSolver.currentTimeStep
        self.currentMemoryIndex      = flowSolver.currentMemoryIndex
        self.dt                      = flowSolver.dt
        self.nTsteps                 = flowSolver.nTsteps
        
        self.dsetGroup = vascularNetwork.BrxDataGroup.create_group('Baroreflex - '+str(self.baroId))
      
                     
        bc2out = {}
        terminalBoundaries = 0
        
        for bcId,bcs in vascularNetwork.boundaryConditions.iteritems():
            
            if bcId == vascularNetwork.root:
                
                for bc in bcs:
                    if bc.type == 1:
#                                 baroData['boundaryCondition'] = bc
                    
                        self.boundaryCondition = bc
                    
                    elif bc.type == 2:
#                                 baroData['boundaryConditionII'] = bc
                        self.boundaryConditionII = bc

                    else:
                        print "Warning Baroreceptor: vascularNetwork.root has the wrong type of boundary condition"   
        
            elif bcId != vascularNetwork.root:
                for bc in bcs:
                    if bc.type == 2: # type 2 BC, outflow or Varying Elastance heart
                        bc2out[bcId] = bc
                        terminalBoundaries = terminalBoundaries + 1
        
#                 baroData['boundaryConditionIIout'] = bc2out
        self.boundaryConditionIIout = bc2out
        self.terminalBoundaries = terminalBoundaries # number of terminal boundaries used to calculate delta_R for each WK at the distal end of a network
        self.venousPool = flowSolver.venousPool # venous pool object for the update of Vusv
            
        ### quantities in the left ventricle --> to extract for postprocessing
        self.newCycles = np.zeros(1)
        self.dsetGroup.create_dataset("newCycles", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Pheart = np.zeros(self.nTsteps+1)
        self.dsetGroup.create_dataset("Pheart", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Vheart = np.zeros(self.nTsteps+1)
        self.dsetGroup.create_dataset("Vheart", (vascularNetwork.savedArraySize,), dtype='float64')
   
        self.Venous_dsetGroup = self.dsetGroup.create_group('venousPool')
        self.Venous_dsetGroup.create_dataset("V_vein", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Venous_dsetGroup.create_dataset("CVP", (vascularNetwork.savedArraySize,), dtype='float64')
        self.Venous_dsetGroup.create_dataset("LAP", (vascularNetwork.savedArraySize,), dtype='float64')
        
        
    def update(self,baroDict):
            '''
            updates the Baroreceptor using a dictionary in form of 
            baroDict = {'variableName': value}
            '''
            for key,value in baroDict.iteritems():
                try:
                    self.__getattribute__(key)
                    self.__setattr__(key,value)
                except: 
                    print "ERROR 139 baroreceptor.update(): wrong key: %s, could not set up baroreceptor" %key
    
    
    def getVariableValue(self,variableName):
        '''
        Returns value of variable with name : variableName
        States Error if not such variable
        '''
        try:
            return self.__getattribute__(variableName)
        except: 
            print "ERROR CarotidBaroreceptor.getVariable() : CarotidBaroreceptor has no variable {}".format(variableName)
                       
                       
    def getVariableDict(self):
        '''
        Returns a deep copy of the class variable dict
        '''
        return self.__dict__
    
    
    ### Begin solver type functions
    ### solving of imported CellML model
    def solveCellML(self, strain):
        """
        solve CellML model, with the solver2 method that has been added to the CellML export
        """
        # import the CellML Module
        baroreceptorCellML = self.baroreceptorCellML
        
        nbElements = np.shape(strain)[0]            
        #timeArray = np.linspace(0,(nbElements-2)*self.dt,(nbElements-1))
        
        # set the input values to the model, which are defined in input
        for i in xrange(0,(nbElements)):
                
            self.constants[self.cellMLinputID] = strain[i]
            self.voi, self.states, self.algebraic = baroreceptorCellML.solver2([0,self.dt],self.states[-1][:],self.constants)
        
        
        return self.voi, self.states, self.algebraic
    
    
    def calcAndupdatePeriodTypeIcellML(self,n,strain):
        """
        Function to calculate new period for BC of type 1
        Calls self.solveCellML()
        is used by self.__call__()
        """
        if n == (round(self.newUpdateTime-2)): 
            # solve the cellML system using the method defined above
            self.voi, self.states, self.algebraic = self.solveCellML(strain)
            
            Tperiod  = self.algebraic[-1][self.cellMLoutputID]
            
            # update period in runtime using the method defined in the BC type 1 class
            self.boundaryCondition.updatePeriodRuntime(Tperiod,(self.newUpdateTime-2)*self.dt)
            
            # set new update time
            self.oldUpdateTime = self.newUpdateTime
            self.newUpdateTime = self.oldUpdateTime + round(self.boundaryCondition.Tperiod/self.dt)
        
        
    def calcAndupdatePeriodTypeIIcellML(self,strain):
        """
        Function to calculate new period of BC type 2 - Varying Elastance Heart (only implemented in the simple heart)
        Calls self.solveCellML()
        is used by self.__call__()
        """
        if self.boundaryConditionII.name == 'VaryingElastanceSimple':
            self.voi, self.states, self.algebraic = self.solveCellML(strain)
                
            if self.boundaryConditionII.newCycle == True:
                self.boundaryConditionII.T = self.algebraic[-1][self.cellMLoutputID] # update period
                self.boundaryConditionII.Tpeak = 0.43*self.boundaryConditionII.T # update Tpeak
                self.newCycles = np.append(self.newCycles,self.currentTimeStep) # to save the time step where new cycles start to solution data
                
        else: 
            print "BR line 124: not a varying elastance heart"
            
         
class AorticBaroreceptor(Baroreceptor):
    """
    for models of the AorticBaroreceptors
    Aortic Baroreceptor models with strain input and period of the heart cycle as output
    """
    
    def __init__(self,BaroDict):
        """
        constructor method of an AorticBaroreceptor object
        """
        # intialize with mother class constructor
        super(AorticBaroreceptor, self).__init__(BaroDict)
        
        # CellML expected to exist initialized in initializeForSimulation
        self.cellMLBaroreceptorModel = True
        self.cellMLimport = '' # for the import of the CellML file
        self.baroreceptorCellML = None
        self.voi = None
        self.states = None 
        self.algebraic = None 
        
        self.vesselIds = None
        self.Area1 = 0
        self.Area2 = 0
        self.Pressure1 = 0
        self.Pressure2 = 0
        self.Strain = 0
        self.MStrain = 0
        self.T = 0
        
        self.n = None
        self.Tsym = None
        self.Tparasym = None
        self.c_nor = None
        self.c_ach = None
              
        # update with dictionary 
        self.update(BaroDict)
        
    def initializeForSimulation(self,flowSolver, vascularNetwork): 
        """
        Configures class members for simulation, using data from the network 
        which may not have been available at construction.
        """
        # Do basic stuff
        super(AorticBaroreceptor,self).initializeForSimulation(flowSolver, vascularNetwork)
        
        #initialize the CellML Baroreceptor model if given
        if self.cellMLBaroreceptorModel == True:
            baroreceptorCellML = self.baroreceptorCellML
            # input to the CellML model and output from it --> defined in the header of the Python CellML export
            self.cellMLinputID = baroreceptorCellML.inputID
            self.cellMLoutputArray = baroreceptorCellML.outputArray
            self.cellMLoutputID = baroreceptorCellML.outputID
         
           
            # initialize the model (constant parameters and initial values of states) 
            (iniStates, constants) = baroreceptorCellML.initConsts()
            self.updateConstants(constants)
            
            timeArray = np.linspace(0,self.dt,2)
            self.voi, self.states, self.algebraic = baroreceptorCellML.solver2(timeArray,iniStates,self.constants)
            
                     
        self.Area1 = vascularNetwork.vessels[self.vesselIds[0]].Asol
        self.Area2 = vascularNetwork.vessels[self.vesselIds[1]].Asol
        self.Pressure1 = vascularNetwork.vessels[self.vesselIds[0]].Psol
        self.Pressure2 = vascularNetwork.vessels[self.vesselIds[1]].Psol
        
        self.Strain = np.zeros([self.nTsteps + 1,np.shape(self.Area1)[1]+np.shape(self.Area2)[1]])
        self.MStrain = np.zeros(self.nTsteps + 1)
        self.T       = np.zeros(self.nTsteps + 1)
        
        # to save the timesteps where new cycles start in the solution data
        self.newCycles = np.zeros(1)
        
        # intial areas from the inlet node of the vessels
        self.Ao_1 = self.Area1[0] 
        self.Ao_2 = self.Area2[0]
        
        #initial pressure --> used for estimation of unstretched radius
        self.Po_1 = self.Pressure1[0]
        self.Po_2 = self.Pressure2[0]
            
        # unstretched radius
        self.Ro_1 = None # np.zeros(np.shape(self.Area1)[0])
        self.Ro_2 = None # np.zeros(np.shape(self.Area2)[0])
        
        # estimate the unstretched Radii R0 for the calculation of the strain
        self.estimateUnstretchedRadius()
        
        # update Time for the period of the inflow function (type 1 BC's)
        self.oldUpdateTime = 0
        if self.boundaryCondition != 0:
            self.newUpdateTime = round(self.boundaryCondition.Tperiod/self.dt)


                                       
        
    def flushSolutionData(self,saving, nDB,nDE,nSB,nSE):
        
        if saving:
            
            ### quantities of the Baroreflex loop
            self.dsetGroup['T'][nDB:nDE] =  self.T[nSB:nSE]
            self.dsetGroup["n"][nDB:nDE] = self.n[nSB:nSE]
            self.dsetGroup["Tsym"][nDB:nDE] = self.Tsym[nSB:nSE]
            self.dsetGroup["Tparasym"][nDB:nDE] = self.Tparasym[nSB:nSE]
            self.dsetGroup["c_nor"][nDB:nDE] = self.c_nor[nSB:nSE]
            self.dsetGroup["c_ach"][nDB:nDE] = self.c_ach[nSB:nSE]
            
                                                
            
            self.dsetGroup["newCycles"][nDB:nDE] =    self.newCycles[nSB:nSE]
            self.dsetGroup["Pheart"][nDB:nDE]    =self.Pheart[nSB:nSE]
            self.dsetGroup["Vheart"][nDB:nDE]    = self.Vheart[nSB:nSE]
    
            self.Venous_dsetGroup["V_vein"][nDB:nDE] = self.venousPool.Vvector[nSB:nSE]
            self.Venous_dsetGroup["CVP"][nDB:nDE] = self.venousPool.Pvector[nSB:nSE]
            self.Venous_dsetGroup["LAP"][nDB:nDE] = self.venousPool.P_LAvector[nSB:nSE]
            
            
    def estimateUnstretchedRadius(self):
        """
        Function to estimate the unstretched radius of the Vessels of the Aortic Arch for the calculation of strain.
        Takes information about compliance and initial vessel geometry.
        Takes the initial compliance and assumes it as constant, uses the initial Area to calculate 
        the Area at zero pressure, where the vessel is in its unstressed state, so that the true
        strain can be calculated.
        """
        
        ### these following lines did not work with the given initial 
        
        #compliance1 = self.initialCompliance1.C(self.Po_1)
        #compliance2 = self.initialCompliance2.C(self.Po_2)
        
        #deltaA1 = compliance1 * self.Po_1
        #deltaA2 = compliance2 * self.Po_2

        #self.Ro_1 = np.power((self.Ao_1 - deltaA1)/math.pi,0.5)
        #self.Ro_2 = np.power((self.Ao_2 - deltaA2)/math.pi,0.5)
        
        
        f = 1.3 # estimated value, probably somewhere between 1.35 and 1.65
        
        # initial radii of the two vessels given
        self.Ro_1 = np.power((self.Ao_1)/math.pi,0.5)/f
        self.Ro_2 = np.power((self.Ao_2)/math.pi,0.5)/f
        
    
       
    def __call__(self):
        """
        call function of AorticBaroreceptor
        """
        
        n = self.currentTimeStep[0]
        n_mem = self.currentMemoryIndex[0]
        
        ## read area, calculate radius and strain
        A1 = self.Area1[n_mem]
        R1 = np.power(A1/math.pi,0.5)
       
        A2 = self.Area2[n_mem]
        R2 = np.power(A2/math.pi,0.5)
        
        epsilon1 = (R1 - self.Ro_1)/self.Ro_1 # strains in the two vessels
        epsilon2 = (R2 - self.Ro_2)/self.Ro_2
        
        # concatenate the strain arrays of the two vessels of the Aortic Arch
        epsilon = np.concatenate((epsilon1,epsilon2),axis=0)
        epsMean = np.mean(epsilon) # calculate one mean strain for the input to the Baroreceptor model
        
        
        # update the strain and the mean strain
        if n < self.nTsteps-1:
            self.Strain[n+1]  = epsilon
            self.MStrain[n+1] = epsMean
            
            
        # use CellML model to calculate new period, the method calcAndupdatePeriodcellML
        # does this on a beat-to-beat basis
        if self.cellMLBaroreceptorModel == True:
            
            if self.boundaryCondition != 0:
                self.calcAndupdatePeriodTypeIcellML(n,self.MStrain[self.oldUpdateTime:(self.newUpdateTime-2)])
                print "BR255"
                print self.algebraic[-1][self.cellMLoutputID]
                print self.boundaryCondition.Tperiod
                
            else:
                self.calcAndupdatePeriodTypeIIcellML(self.MStrain[n:(n+1)])
            
            # update the current period, update the period depending on the array in which it is defined in the CellML export
            if self.cellMLoutputArray == 'algebraic':
                
                self.T[n+1] = self.algebraic[-1][self.cellMLoutputID]
            
            
            elif self.cellMLoutputArray == 'states':
                
                self.T[n+1] = self.states[-1][self.cellMLoutputID]

        

class bugenhagenAorticBR(AorticBaroreceptor):
    '''
    for models of the AorticBaroreceptors
    Aortic Baroreceptor models with strain input and period of the heart cycle as output
    '''
    
    def __init__(self,BaroDict):
        """
        constructor method of an AorticBaroreceptor object
        """
        # intialize with mother class constructor
        super(bugenhagenAorticBR, self).__init__(BaroDict)
        # Configuration and solution data variables
        self.modelName = 'bugenhagenAorticBR'
        self.baroreceptorCellML = cellMLBaroreflexModels.bugenhagenAorticBR
        
        self.update(BaroDict)
        
    def initializeForSimulation(self, flowSolver, vascularNetwork):
        AorticBaroreceptor.initializeForSimulation(self, flowSolver, vascularNetwork)

        # arrays  used to save BR quantities to solution data
        # the saving is done at the end of the solver method in classFlowSolver
        self.n = np.ones(self.nTsteps+1)*self.algebraic[0][3]
        self.Tsym = np.ones(self.nTsteps+1)*self.algebraic[0][7]
        self.Tparasym = np.ones(self.nTsteps+1)*self.algebraic[0][8]
        self.c_nor = np.ones(self.nTsteps+1)*self.states[0][3]
        self.c_ach = np.ones(self.nTsteps+1)*self.states[0][4]
            
            
        ### quantities of the Baroreflex loop
        self.dsetGroup.create_dataset("T", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("n", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("Tsym", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("Tparasym", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("c_nor", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("c_ach", (vascularNetwork.savedArraySize,), dtype='float64')
        

    def __call__(self):
        '''
        Implements bugenhagen specific actions for the numerical object call
        '''
        super(bugenhagenAorticBR,self).__call__()
        
        n = self.currentTimeStep[0]
        n_mem = self.currentMemoryIndex[0]
        self.n[n+1] = self.algebraic[-1][3]
        self.Tsym[n+1] = self.algebraic[-1][7]
        self.Tparasym[n+1] = self.algebraic[-1][8]
        self.c_nor[n+1] = self.states[-1][3]
        self.c_ach[n+1] = self.states[-1][4]
        

class pettersenAorticBR(AorticBaroreceptor):
    '''
    for models of the AorticBaroreceptors
    Aortic Baroreceptor models with strain input and period of the heart cycle as output
    '''
    
    def __init__(self,BaroDict):
        """
        constructor method of an AorticBaroreceptor object
        """
        # intialize with mother class constructor
        super(pettersenAorticBR, self).__init__(BaroDict)
        self.modelName = 'pettersenAorticBR'
        self.baroreceptorCellML = cellMLBaroreflexModels.pettersenAorticBR
        
        self.L0 = None
        self.n0 = None
        self.g = None
        self.tau1 = None
        self.tau2 = None
        self.Gp = None
        self.Gs = None
        self.delta_HR_smax = None
        self.delta_HR_pmax = None
    
        self.update(BaroDict)
    
    def updateConstants(self, constants):
        """
        Assign parameter values to those used by the cellML solver
        """
        self.constants = constants
        # Pettersen
        self.constants[30] = self.delta_HR_smax
        self.constants[30] = self.delta_HR_smax
        self.constants[31] = self.delta_HR_pmax 
        self.constants[24] = self.L0
        self.constants[25] = self.n0
        self.constants[26] = self.g
        self.constants[22] = self.tau1
        self.constants[23] = self.tau2
        self.constants[6]  = self.Gp
        self.constants[5]  = self.Gs
        
    def initializeForSimulation(self, flowSolver, vascularNetwork):
        super(pettersenAorticBR,self).initializeForSimulation(flowSolver, vascularNetwork)
        # arrays  used to save BR quantities to solution data
        # the saving is done at the end of the solver method in classFlowSolver
        self.n = np.ones(self.nTsteps+1)*self.algebraic[0][5]
        self.Tsym = np.ones(self.nTsteps+1)*self.algebraic[0][9]
        self.Tparasym = np.ones(self.nTsteps+1)*self.algebraic[0][10]
        self.c_nor = np.ones(self.nTsteps+1)*self.states[0][2]
        self.c_ach = np.ones(self.nTsteps+1)*self.states[0][3]
            
        ### quantities of the Baroreflex loop
        self.dsetGroup.create_dataset("T", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("n", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("Tsym", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("Tparasym", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("c_nor", (vascularNetwork.savedArraySize,), dtype='float64')
        self.dsetGroup.create_dataset("c_ach", (vascularNetwork.savedArraySize,), dtype='float64')
        
    def __call__(self):
        """
        Implements pettersen specific actions for the numerical object call
        """      
        super(pettersenAorticBR,self).__call__()
        
        n = self.currentTimeStep[0]
        n_mem = self.currentMemoryIndex[0]
        
        self.n[n+1] = self.algebraic[-1][5]
        self.Tsym[n+1] = self.algebraic[-1][9]
        self.Tparasym[n+1] = self.algebraic[-1][10]
        self.c_nor[n+1] = self.states[-1][2]
        self.c_ach[n+1] = self.states[-1][3]

class CarotidBaroreceptor(Baroreceptor):
    '''
    for models of the Carotid Baroreceptors
    '''
    def __init__(self,BaroDict):
        
        """
        constructor for the CarotidBaroreceptor
        """
        #System and Vessel Variables
        super(CarotidBaroreceptor, self).__init__()
        
        
        self.vesselIdLeft = 0
        self.vesselIdRight = 0
        self.terminalBoundaries = 0
        self.venousPool = 0
        
        # pressure in left and right Carotid sinus
        self.pressureLeft = 0
        self.pressureRight = 0
        
        ### total resistance
        self.Res0 = {} # Windkessel resistances of all boundaries --> for updating the resistances
        self.ResTot0 = 0.0 # total resistance of the windkessels
        
#         ### for update of period in type 1 boundary condition
#         self.oldUpdateTime = 0
#         if self.boundaryCondition != 0:
#             self.newUpdateTime = round(self.boundaryCondition.Tperiod/self.dt)
#         
#         
        #######################################
        """
        Parameters for the Ursino model of the carotid baroreflex
        """
        #######################################
        
        ### model parameters afferent part of the Ursino model - Ursino 1999
        
        self.pn = 135*133.32 # adjusted from 92 mmHg as given in paper, chose a value close to the initialisation
        self.ka = 1567. # parameter associated to the static pressure - firing rate function
        self.fmin = 2.51 # lower saturation of afferent firing rate
        self.tau_z= 6.37 # time constant of zero of first order dynamical block of afferent part
        self.fmax = 47.78 # upper saturation of afferent firing rate
        self.tau_p= 2.076 # time constant of pole of first order dynamical block of afferent part
        
        # model parameters efferent part of the Ursino model  - Ursino 1999
        
        self.fe_inf = 2.10 # 
        self.fe_0 = 16.11 #
        self.ke = 0.0675 # 
        self.fe_min = 2.66 # threshold value for efferent firing used in the effector parts
        
        # model parameters for the TPR effector part of the Ursino model  - Ursino 1999
        
        self.cR = 0.317*133.32e6 # gain
        self.tauR = 6 # time constant
        self.DR = 2.0 # delay
        self.R0 = 0.61*133.32e6 # value in absence of sympathetic drive
        
        
        
        # model parameters for the EmaxLV effector part of the Ursino model  - Ursino 1999   
        self.cE = 0.475 * 133.322368e6 # gain
        self.tauE = 8. # time constant
        self.DE = 2.0 # pure delay
        self.E0 = 2.392*133.322368e6 # value in absence of sympathetic drive
        
        # model parameters for the T effector part of the Ursino model  - Ursino 1999
        
        self.cT = -0.08 # gain
        self.tauT = 2.0 # time constant
        self.DT = 2.0 # pure delay
        self.T0 = 0.58 # value in absence of sympathetic drive 
        
        
        # model parameters for the Vusv effector part of the Ursino model  - Ursino 1999
        
        self.cVusv = -625.0e-6 # gain
        self.tauVusv = 20.0 # time constant
        self.DVusv = 3.0 # pure delay
        self.Vusv0 = 3213e-6 # value in absence of sympathetic drive
        
        
        ### initial values for affected quantities
        ### values set close to steady state values which were found by open-loop simulation
        ### and readjusted after closing the loop, in order to minimize the transients
        
        self.dTPRin = 5.85e7
        self.dTin = -0.1114
        self.dEmaxin = 87.65e6
        self.dVusvin = -877.54e-6
        self.Ptildin = 10600.
        # ratio between total WK resistance of Network and value given by Ursino 
        self.ratio = 1.0 # Updated when initializing for simulation 
        # self.ratio = (self.ResTot0)/(self.R0 + self.dTPRin)
        
#         ### states of the Ursino model###
#         # afferent part: left and right side
#         self.PtildLeft = np.zeros(self.nTsteps+1)
#         self.PtildLeft[0] = self.Ptildin # intial value for first oprder dynamical block
#         
#         self.PtildRight = np.zeros(self.nTsteps+1)
#         self.PtildRight[0] = self.Ptildin
#         
#         self.F_cs_left = np.zeros(self.nTsteps+1)
#         self.F_cs_right = np.zeros(self.nTsteps+1)
#         
#         self.F_cs = np.zeros(self.nTsteps+1) # mean value of afferent firing (mean of left and right carotid sinus)
#         
#         
#         # efferent part
#         self.F_efferent = np.zeros(self.nTsteps+1)
#         
#         
#         # effector parts: initilized with the initial value of the respective quantities
#         self.delta_TPR = np.ones(self.nTsteps+1)*self.dTPRin
#         self.delta_Emax = np.ones(self.nTsteps+1)*self.dEmaxin
#         self.delta_Vusv = np.ones(self.nTsteps+1)*self.dVusvin
#         self.delta_T = np.ones(self.nTsteps+1)*self.dTin
#         
#         
#         ### quantities in the left ventricle --> to extract for postprocessing
#         self.newCycles = np.zeros(1)
#         self.Pheart = np.zeros(self.nTsteps+1)
#         self.Vheart = np.zeros(self.nTsteps+1)
    
            ### states of the Ursino model###
        # afferent part: left and right side
        self.PtildLeft = None
        self.PtildRight = None
        
        self.F_cs_left = None
        self.F_cs_right = None
        self.F_cs = None # mean value of afferent firing (mean of left and right carotid sinus)
        
        
        # efferent part
        self.F_efferent = None
                
        # effector parts: initilized with the initial value of the respective quantities
        self.delta_TPR = None
        self.delta_Emax = None
        self.delta_Vusv = None
        self.delta_T = None
        
        ### quantities in the left ventricle --> to extract for postprocessing
        self.newCycles = None
        self.Pheart = None
        self.Vheart = None

        # update with dictionary
        self.update(BaroDict)

        
        
    def initializeForSimulation(self,flowSolver, vascularNetwork): 
        """
        Configures class members for simulation, using data from the network 
        which may not have been available at construction.
        """
        # Do basic stuff
        Baroreceptor.initializeForSimulation(self, flowSolver, vascularNetwork)
                               
        # New Solver Initialization
        self.pressureLeft   = vascularNetwork.vessels[self.vesselIdLeft].Psol
        self.pressureRight  = vascularNetwork.vessels[self.vesselIdRight].Psol
        
        reciprocal_sum = 0.0     
        for key in self.boundaryConditionIIout:
            self.Res0[key] = self.boundaryConditionIIout[key].Rtotal
            reciprocal_sum = reciprocal_sum + 1/self.boundaryConditionIIout[key].Rtotal # summation of parallel resisitances
                   
        # total resisistance of the Windkessels (reciprocal of the sum above - R in parallel!) 
        self.ResTot0 = 1.0/reciprocal_sum
        
        self.ratio = (self.ResTot0)/(self.R0 + self.dTPRin)
        
        ### for update of period in type 1 boundary condition
        self.oldUpdateTime = 0
        if self.boundaryCondition != 0:
            self.newUpdateTime = round(self.boundaryCondition.Tperiod/self.dt)
        
        
        ### states of the Ursino model###
        # afferent part: left and right side
        self.PtildLeft = np.zeros(self.nTsteps+1)
        self.PtildLeft[0] = self.Ptildin # intial value for first oprder dynamical block
        
        self.PtildRight = np.zeros(self.nTsteps+1)
        self.PtildRight[0] = self.Ptildin
        
        self.F_cs_left = np.zeros(self.nTsteps+1)
        self.F_cs_right = np.zeros(self.nTsteps+1)
        
        self.F_cs = np.zeros(self.nTsteps+1) # mean value of afferent firing (mean of left and right carotid sinus)
        self.dsetGroup.create_dataset("F_cs", (vascularNetwork.savedArraySize,), dtype='float64')
        
        # efferent part
        self.F_efferent = np.zeros(self.nTsteps+1)
        self.dsetGroup.create_dataset("F_efferent", (vascularNetwork.savedArraySize,), dtype='float64')
        
        # effector parts: initilized with the initial value of the respective quantities
        self.delta_TPR = np.ones(self.nTsteps+1)*self.dTPRin
        self.dsetGroup.create_dataset("delta_TPR", (vascularNetwork.savedArraySize,), dtype='float64')
        
        self.delta_Emax = np.ones(self.nTsteps+1)*self.dEmaxin
        self.dsetGroup.create_dataset("delta_Emax", (vascularNetwork.savedArraySize,), dtype='float64')
        self.delta_Vusv = np.ones(self.nTsteps+1)*self.dVusvin
        self.dsetGroup.create_dataset("delta_Vusv", (vascularNetwork.savedArraySize,), dtype='float64')
        self.delta_T = np.ones(self.nTsteps+1)*self.dTin
        self.dsetGroup.create_dataset("delta_T", (vascularNetwork.savedArraySize,), dtype='float64')
        
        

    
    def flushSolutionData(self,saving, nDB,nDE,nSB,nSE):
        
        if saving:
            # save solution for carotid baroreceptor type
            self.dsetGroup["F_cs"][nDB:nDE] = self.F_cs[nSB:nSE]
         
            # efferent part
            self.dsetGroup["F_efferent"][nDB:nDE] =self.F_efferent[nSB:nSE]
            
    
            self.dsetGroup["delta_TPR"][nDB:nDE] = self.delta_TPR[nSB:nSE]
            self.dsetGroup["delta_Emax"][nDB:nDE] = self.delta_Emax[nSB:nSE]
            self.dsetGroup["delta_Vusv"][nDB:nDE] = self.delta_Vusv[nSB:nSE]
            self.dsetGroup["delta_T"][nDB:nDE] = self.delta_T[nSB:nSE]
                                                
            
            self.dsetGroup["newCycles"][nDB:nDE] =    self.newCycles[nSB:nSE]
            self.dsetGroup["Pheart"][nDB:nDE]    =self.Pheart[nSB:nSE]
            self.dsetGroup["Vheart"][nDB:nDE]    = self.Vheart[nSB:nSE]
    
            self.Venous_dsetGroup["V_vein"][nDB:nDE] = self.venousPool.Vvector[nSB:nSE]
            self.Venous_dsetGroup["CVP"][nDB:nDE] = self.venousPool.Pvector[nSB:nSE]
            self.Venous_dsetGroup["LAP"][nDB:nDE] = self.venousPool.P_LAvector[nSB:nSE]
        
        
    #############################################################################
    
    """
    methods defining the different parts of the Ursino Baroreceptor model
    - Afferent
    - Efferent
    - Effector
    """
    
    def UrsinoAfferent(self,n,n_mem):
        """
        Afferent part of the Baroreceptor model as defined by Ursino 1999
        The firing rate of the left and right carotid sinus are averaged to give a global firing rate
        the firing rate of the carotid sinus will be processed by the efferent part
        """
        
        pLeft = np.mean(self.pressureLeft[n_mem]) # mean pressure value for left side vessel
        pRight = np.mean(self.pressureRight[n_mem]) # mean pressure value for right side vessel
        
        if n == 0:
            
            dpLeft =  0.
            dpRight = 0.
            
        else:
            
            dpLeft = (pLeft - np.mean(self.pressureLeft[n_mem-1]))/self.dt
            dpRight = (pRight - np.mean(self.pressureRight[n_mem-1]))/self.dt
            
        ### forward Euler method to integrate differential equation
        self.PtildLeft[n+1] = (pLeft + self.tau_z*dpLeft-self.PtildLeft[n])*self.dt/self.tau_p+self.PtildLeft[n]
        self.PtildRight[n+1] = (pRight + self.tau_z*dpRight-self.PtildRight[n])*self.dt/self.tau_p+self.PtildRight[n]
        
        self.F_cs_left[n+1] = (self.fmin + self.fmax*math.exp((self.PtildLeft[n+1]-self.pn)/self.ka))/(1+math.exp((self.PtildLeft[n+1]-self.pn)/self.ka))
        self.F_cs_right[n+1] = (self.fmin + self.fmax*math.exp((self.PtildRight[n+1]-self.pn)/self.ka))/(1+math.exp((self.PtildRight[n+1]-self.pn)/self.ka))
        
        ### averaging the left and right firing rate
        self.F_cs[n+1] = 0.5*(self.F_cs_left[n+1]+self.F_cs_right[n+1])
 
              
    def UrsinoEfferent(self,n):
        """
        Efferent part of the Baroreceptor model as defined by Ursino 1999
        Calculation of Efferent firing rate resp. sympathetic activity
        """
        self.F_efferent[n+1] = self.fe_inf + (self.fe_0-self.fe_inf)*math.exp(-self.ke*self.F_cs[n+1])

            
    def UrsinoResistanceEffector(self,n):
        """
        Effector part for TPR of the Baroreceptor model as defined by Ursino 1999
        Calculation of effect of efferent firing rate (sympathetic activity) on system quantity 
        """
        delay = round(self.DR/self.dt)
        vR = 0.0 # "forcing term" of differential equation in effector part
        
        # determine "forcing term" with respect to fe_min and delay of the effector
        if n <= delay:
            vR = self.dTPRin # vR initialized to cancel RHS of effector ODE -> output constant before the pure delay is over -> reduce transients
            
        elif n > delay:
            deltaF = self.F_efferent[n+1-delay] - self.fe_min
            
            if deltaF >= 0.0:
                vR = self.cR * math.log(deltaF + 1)
            
            elif deltaF < 0.0:
                vR = 0.0
        
        # Forward Euler to integrate differential equation
        self.delta_TPR[n+1] = self.dt/self.tauR *(-self.delta_TPR[n] + vR) + self.delta_TPR[n]
        
        
    
    def UrsinoEmaxLVEffector(self,n):
        """
        Effector of maximal contractility of left ventricle
        """    
        delay = round(self.DE/self.dt)
        vE = 0
        
        if n <= delay:
            vE = self.dEmaxin
            
        elif n > delay:
            
            deltaF = self.F_efferent[n+1-delay] - self.fe_min
            
            if deltaF >= 0.0:
                vE = self.cE * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                vE = 0.0
        
        # Forward Euler to integrate differential equation
        self.delta_Emax[n+1] = self.dt/self.tauE *(-self.delta_Emax[n] + vE) + self.delta_Emax[n]
    
    
    
    def UrsinoTEffector(self,n):
        """
        Effector of maximal contractility of left ventricle
        """    
        delay = round(self.DT/self.dt)
        vT = 0
        
        if n <= delay:
            vT = self.dTin
            
        elif n > delay:
            
            deltaF = self.F_efferent[n+1-delay] - self.fe_min
            
            if deltaF >= 0.0:
                vT = self.cT * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                vT = 0.0
        
        # Forward Euler to integrate differential equation
        self.delta_T[n+1] = self.dt/self.tauT *(-self.delta_T[n] + vT) + self.delta_T[n]
    
    
    
    def UrsinoVusvEffector(self,n):
        """
        Effector for venous unstressed volume
        """
        delay = round(self.DVusv/self.dt)
        vVusv = 0.0
        
        if n <= delay:
            vVusv = self.dVusvin
            
        elif n > delay:
            deltaF = self.F_efferent[n+1-delay] - self.fe_min
            
            if deltaF >= 0.0:
                vVusv = self.cVusv * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                vVusv = 0.0
        
        # Forward Euler to integrate differential equation
        self.delta_Vusv[n+1] = self.dt/self.tauVusv *(-self.delta_Vusv[n] + vVusv) + self.delta_Vusv[n]
        
    
    
    def updateBC1(self,n):
        """
        update period in BC of type 1 at inflow boundaries
        uses the method updatedPeriodRuntime of the BoundaryCondition class
        uses newUpdateTime
        """
        
        if n == (round(self.newUpdateTime-2)):
             
            newPeriod = self.T0 + self.delta_T[n+1]
            
            self.boundaryCondition.updatePeriodRuntime(newPeriod,(self.newUpdateTime-2)*self.dt)
            
            # update the updateTime for the update of the next period
            self.oldUpdateTime = self.newUpdateTime
            self.newUpdateTime = self.oldUpdateTime + round(self.boundaryCondition.Tperiod/self.dt)
    
    
    def updateBC2(self,n):
        """
        update Boundary condition of type 2
        TPR and Emax for Varying Elastance heart model
        loop through all boundary conditions of type 2 and update the resistance at every timestep
        --> the single resistances are updated with respect to their initial value dRtot/Rtot0 = dRi/Ri0
        update the heart properties at completion of the beat
        """
        
        newTotalResistance = self.ratio*(self.R0 + self.delta_TPR[n+1])
         
        for key in self.boundaryConditionIIout:
            if self.boundaryConditionIIout[key].name == 'Windkessel-3Elements':
                self.boundaryConditionIIout[key].Rtotal = (newTotalResistance)/self.ResTot0*self.Res0[key]  
                   
            ## for Resistance outlets - not tested      
            elif self.boundaryConditionIIout[key].name == 'Resistance':
                self.boundaryConditionIIout[key].Rc = (newTotalResistance)/self.ResTot0*self.Res0[key]
                   
            ## for WK 2 outlets - not tested      
            elif self.boundaryConditionIIout[key].name == 'Windkessel-2Elements':    
                self.boundaryConditionIIout[key].Rc = (newTotalResistance)/self.ResTot0*self.Res0[key]
            
        ## update inlet BC's of type 2 --> Varying Elastance heart           
        if self.boundaryConditionII != 0:
            if self.boundaryConditionII.newCycle == True:
                self.boundaryConditionII.Emax = self.E0 + self.delta_Emax[n+1] # update the VaryingElastance heart (only implememnted in simple heart)
                self.boundaryConditionII.T = self.T0 + self.delta_T[n+1]
                self.newCycles = np.append(self.newCycles,n) # to save to solution data
                self.Pheart[n] = self.boundaryConditionII.pressure[n] # to save to solution data
                self.Vheart[n] = self.boundaryConditionII.volume[n] # to save to solution data
                
                
            else:
                self.Pheart[n] = self.boundaryConditionII.pressure[n] # to save to solution data
                self.Vheart[n] = self.boundaryConditionII.volume[n]
                

    def updateVenousSide(self,n):
        """
        to update the Venous unstretched Volume
        """
        ### first update when the delay in the effector of Vusv is complete
        if self.currentTimeStep*self.dt > self.DVusv:
            self.venousPool.Vusv = self.Vusv0 + self.delta_Vusv[n]
     
    def UrsinoBRmodel(self,n,n_mem):
        """
        Calculate Baroreflex as described by Ursino_1999
        Effected quantities are TPR, Emax of the left heart, Vusv
        Update the affected quantities in their classes
        """
          
        self.UrsinoAfferent(n,n_mem) # calculate the afferent signals
        self.UrsinoEfferent(n) # calculate the efferent signal
        
        self.UrsinoResistanceEffector(n) # the different effector parts
        self.UrsinoEmaxLVEffector(n) 
        self.UrsinoVusvEffector(n)
        self.UrsinoTEffector(n)
        
        #self.updateBC1(n) # for update of BC type 1 at inlet
        self.updateBC2(n) # update BC type 2
        self.updateVenousSide(n) # update the venous side with Vusv

     
    def __call__(self):
        
        """
        call function
        takes the Baroreceptor model (Ursino_1999) and evaluates it
        """
        
        n = self.currentTimeStep[0]
        n_mem = self.currentMemoryIndex[0]
        
        if self.modelName == 'Ursino':
            
            self.UrsinoBRmodel(n,n_mem)
#             print "BR646"
#             print self.F_efferent[n]
#             print self.delta_TPR[n]
#             print self.delta_T[n]
#             print self.delta_Emax[n]
#             print self.delta_Vusv[n]