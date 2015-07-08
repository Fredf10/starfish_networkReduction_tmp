import sys, os
from reportlab.lib.validators import isNumber

#from duplicity.tarfile import TUREAD
# from UtilityLib.saveSimulationDataToCSV import vesselId
# set the path relative to THIS file not the executing file!
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur + '/../')

#sys.path.append(cur + '/../NetworkLib')
#sys.path.append(cur+'/../SolverLib')

from classVessel import Vessel
from classBoundaryConditions import *
from SolverLib import classBaroreceptor as cBRX
#sys.path.append(cur + '/../UtilityLib')
import UtilityLib.moduleFilePathHandler as mFPH

#sys.path.append(cur + '/../VascularPolynomialChaosLib')
from VascularPolynomialChaosLib.classRandomInputManager import RandomInputManager

import numpy as np
from scipy import interpolate
from math import pi, cos, sin
import pprint

import h5py

class VascularNetwork(object):
    '''
    Class representing a vascular Network
    The vascular network consits out of vessels defined in classVessel::Vessel()
    Additional Topologie, BoundaryConditions and the SimulationContext are saved.
    '''
    def __init__(self, quiet=False):

        # # vascularNetwork variables to set via XML
        self.name = 'vascularNetwork'  # name of the network
        self.description = ''  # description of the current case
        self.dataNumber = 'xxx'  # data number of the network
        self.quiet = quiet  # bool to suppress output
        
        # saving options
        self.pathSolutionDataFilename = None
        self.timeSaveBegin = 0.0  # time when to start saving
        self.timeSaveEnd = 2.0  # time when to end saving
        self.maxMemory = 20  # maximum memory in MB 
        self.saveInitialisationPhase = False  # bool to enable saving of the initPhase
                
        self.vesselsToSave = {}
        self.nSaveBegin = None
        self.nSaveEnd  = None
        self.savedArraySize = None
        self.nDCurrent = None
        self.memoryArraySizeTime = None  # memory array size for the time arrays
        self.solutionDataFile = None  # file name of the solution data
        
        # keep track of time points loaded in memory
        self.tsol =None
        
        # running options
        self.cycleMode = False
        
        # simulation Context
        self.totalTime = 1.0  # simulation time in seconds
        self.CFL = 0.85  # maximal initial CFL number
        self.dt = None  # time step of the simulation determined by the solver
        self.nTsteps = None  # number of timesteps of the simulation case determined by the solver
        self.simulationTime = None  # array with simulation Time
        
                
        # self.motion         = {'keyframe': [0, 0.1, 1.0],
        #                        'X1'     : [0, 45, 90]}
        # dict defining the movment by change of angle unsing keyframes
        # {'keyframe': [t0, t1, tend],
        #  'X1: [0, 45, 90]} ## <- correspond to 90 degree change of angleXtoMother of vessel 1        
        self.motionAngles = {}
        
        # gravity controls
        self.gravitationalField = False  # bool, turn gravity on or off
        self.gravityConstant = -9.81  # earth gravity
                
        # venous system
        self.centralVenousPressure = 0.0  # central venous pressure
        self.minimumVenousPressure = 0.0  # minimum allowed venous pressure
        
        # the solver calibration 
        self.rigidAreas = False  # # 'True' 'False' to change
        self.simplifyEigenvalues = False  #
        self.riemannInvariantUnitBase = 'Pressure'  # 'Pressure' or 'Flow'
        self.automaticGridAdaptation = True  # False True        
        # self.solvingSchemeField       = 'MacCormack' # MacCormack
        self.solvingSchemeConnections = 'Linear'  # 'Linear'
                
        # initialization controls
        self.initialsationMethod = 'Auto'  # 'Auto', 'MeanFlow', 'MeanPressure', 'ConstantPressure'
        self.initMeanFlow = 0.0  # initial mean flow value (at inflow point)
        self.initMeanPressure = 0.0  # initial pressure value (at inflow point)
        self.initialisationPhaseExist = True  # bool is False only for 'ConstantPressure'
        self.initPhaseTimeSpan = 0.0  # time span of the init phase
        self.nTstepsInitPhase = 0  # number of timesteps of the initPhase
        
        self.estimateWindkesselCompliance = 'Tree'  # 'Tree', 'Sys', 'Wk3', 'None'
        self.compPercentageWK3 = 0.3  # Cwk3 percentage on total Csys
        self.compPercentageTree = 0.8  # Ctree percentage on total Csys
        self.compTotalSys = 5.0  # total Csys
                
        self.optimizeTree = False  # optmizie areas of vessels to minimize reflections in root direction 
                
        # # dictionaries for network components 
        self.vessels = {}  # Dictionary with containing all vessel data,  key = vessel id; value = vessel::Vessel()
                
        self.boundaryConditions = {}
                
        self.globalFluid = {'my': 1e-6, 'rho': 1050., 'gamma': 2.0}  # dictionary containing the global fluid data if defined
        
        self.baroreceptors = {}  # dictionary with baroreceptors
        
        self.communicators = {}  # dictionary with communicators, key = communicator id; values = {communicator data} 
        
        # # interal calculated variables
        self.root = None  # the root vessel (mother of the mothers)
        self.boundaryVessels = []  # includes all vessels with terminal boundaryConditions (except start of root)
        self.treeTraverseList = []  # tree traverse list 
        self.treeTraverseConnections = []  # tree travers list including connections [ LM, RM , LD, RD ]
     
        self.initialValues = {}
        self.Rcum = {}  # Dictionary with all cumultative resistances
        self.Cends = {}  # Dictionary with the area compliances of all terminating vessels (at ends)
        self.totalTerminalAreaCompliance = None  # the sum of all Cends
        self.TotalVolumeComplianceTree = None  # total volume compliance of all vessels
        
#         ### random variables TODO: move out of here to global class        
        self.randomInputManager = None        
        
    # all classes concerning vessel
    def addVessel(self, vesselId=None, dataDict=False):
        '''
        adds vessel to the Network
        if no id, a random id is choosen
        if no DataDict, no values are assigned
        '''
        # set id to 1 + highest id of existing vessels
        if vesselId == None: 
            try: vesselId = max(self.vessels.keys()) + 1
            except: vesselId = 0
            
        # check Id
        if vesselId not in self.vessels:
            vessel = Vessel(Id=vesselId , name=('vessel_' + str(vesselId)))  # create vessel with given variables
            if dataDict:
                vessel.update(dataDict)  # set vesselData if available
            self.vessels[vessel.Id] = vessel  # add vessel to network
        else:  
            print "Error vascularNetwork.addVessel: vessel with Id {} exists already! Could not add vessel".format(vesselId)  # raise error if Id is set doubled
 
    def deleteVessel(self, inputId):
        '''
        Remove vessel from network and delete it
        '''
        try:
            del self.vessels[inputId]
        except:
            print "ERROR vascularNetwork.deleteVessel(): vessel with Id {} does not exist! Could not remove vessel".format(inputId)
    
    def addBaroreceptor(self, baroId=None, dataDict=False):
        '''
        adds vessel to the Network
        if no id, a random id is choosen
        if no DataDict, no values are assigned
        '''
        # set id to 1 + highest id of existing vessels
        if baroId == None: 
            try: baroId = max(self.baroreceptors.keys()) + 1
            except: baroId = 0
             
        # check Id
        if baroId not in self.baroreceptors:
            baroType = dataDict['modelName']
            instance = getattr(cBRX, baroType)(dataDict)
            self.baroreceptors[baroId] = instance
            
        else:  
            print "Error vascularNetwork.addBaroreceptor: baroreceptor with Id {} exists already! Could not add baroreceptor".format(baroId)  # raise error if Id is set doubled
  
    def update(self, vascularNetworkData):
        '''
        updates the vascularNetwork data using a dictionary in form of 
        vascularNetworkData = {'variableName': value}
        '''
        for key, value in vascularNetworkData.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key, value)
            except: 
                print 'WARNING vascularNetwork.update(): wrong key: %s, could not update vascularNetwork' % key 
                
    def getVariableValue(self, variableName):
        '''
        Returns value of variable with name : variableName
        States Error if not such variable
        '''
        try:
            return self.__getattribute__(variableName)
        except: 
            print "ERROR vascularNetwork.getVariable() : VascularNetwork has no variable {}".format(variableName)
    
    def updateNetwork(self, updateDict):
        '''
        Update vascular Network with an Dictionary: updateDict
        
        updateDict = {'vascularNetworkData': {},
                      'globalFluid': {},
                      'globalFluidPolyChaos': {},
                      'communicators': {},
                      'vesselData': {},
                      'baroreceptors': {}}

            'vascularNetworkData'  := dict with all vascularNetwork variables to update
            'globalFluid'          := dict with all global fluid properties
            'communicators'        := netCommunicators}
            'vesselData'           := { vessel.id : DataDict}
            'baroreceptors'        := { baroreceptor.id : DataDict}
        '''
        
        for dictName in ['vascularNetworkData']:
            try: self.update(updateDict[dictName])
            except: pass
            
#         for dictName in ['globalFluid', 'communicators', 'baroreceptors']:
        for dictName in ['globalFluid', 'communicators']: 
            try: self.getVariableValue(dictName).update(updateDict[dictName])
            except: pass
                        
        if 'vesselData' in updateDict:
            for vesselId, vesselData in (updateDict['vesselData']).iteritems():
                try:
                    self.vessels[vesselId].update(vesselData)
                except:
                    self.addVessel(vesselId, vesselData)    
                     
        if 'baroreceptors' in updateDict:
            for baroId, baroData in (updateDict['baroreceptors']).iteritems():
                try:
                    self.baroreceptors[baroId].update(baroData)
                except:
                    self.addBaroreceptor(baroId, baroData)
            print self.baroreceptors
                      
    def showVessels(self):
        '''
        writes the Vesseldata for each vessel to console (calls printToConsole() from each vessel)
        '''
        print " Vessels in Network:"
        for vessel in self.vessels.itervalues():
            vessel.printToConsole()
    
    def showNetwork(self):
        '''
        writes Network properties (without vesselData) to console
        '''
        print "-------------------"
        print " vascularNetwork ", self.name, "\n"
        for variable, value in self.__dict__.iteritems():
            try: 
                print " {:<20} {:>8}".format(variable, value)
            except: print " {:<20} {:>8}".format(variable, 'None')
        
    def initialize(self, initializeForSimulation=False):
        '''
        Initializes vascular network: the compliance of the vessels and the position of the call function of boundary type 2
        Check if boundaryConditions and globalFluid properties are defined in a rigth manner;
        '''   
        # # refresh all connections and generate traversing lists
        self.evaluateConnections()
        
        # # checks if gravity is turned on
        if self.gravitationalField == False: self.gravityConstant = 0.
                        
        # ## check global fluid properties
        for fluidItem, value in self.globalFluid.iteritems():
            if value == None:
                if fluidItem == 'dlt':
                    try:
                        gamma = self.globalFluid['gamma']
                        self.globalFluid['dlt'] = (gamma + 2.0) / (gamma + 1)
                    except:
                        print "ERROR: VascularNetwork.initialize(): global fluid properties are not proper defined! Check:"
                        pprint.pprint(self.globalFluid)
                        exit()
                else:
                    print "ERROR: VascularNetwork.initialize(): global fluid properties are not proper defined! Check:"
                    pprint.pprint(self.globalFluid)
                    exit()
        # ## initialize vessels
        for vessel in self.vessels.itervalues():
            vessel.initialize(self.globalFluid)
            vessel.update({'gravityConstant': self.gravityConstant})
        
            
        # ## check and initialize boundary conditions
        if self.boundaryConditions != {}:
            # # Set position of boundary conditions
            # check position if one Vessel
            if len(self.vessels) == 1:
                vesselId = self.boundaryConditions.keys()[0]
                if vesselId != self.root: print "Error Wrong Root found"
                for bc in self.boundaryConditions[vesselId]:
                    if '_' not in bc.name[0]: bc.setPosition(0)
                    else: bc.setPosition(-1)                    
            else:
                for vesselId, bcs in self.boundaryConditions.iteritems():
                    if vesselId == self.root:
                        for bc in bcs:
                            bc.setPosition(0)
                            
                    elif vesselId in self.boundaryVessels:
                        for bc in bcs:
                            bc.setPosition(-1)       
             
        definedButNotatBC = set(self.boundaryConditions.keys()).difference([self.root] + self.boundaryVessels)
        atBCButNotDefined = set([self.root] + self.boundaryVessels).difference(self.boundaryConditions.keys())
        if len(definedButNotatBC.union(atBCButNotDefined)) > 0:
            print "ERROR: VascularNetwork.initialize(): BoundaryConditions are not proper defined:"
            if len(definedButNotatBC) > 0:print "       for Vessel(s) {} boundaryConditions are defined but \n        Vessel(s) is(are) not at the Boundary!!!".format(list(definedButNotatBC))
            if len(atBCButNotDefined) > 0:print "       for Vessel(s) {} no BoundaryConditions are defined!!!".format(list(atBCButNotDefined))
            exit()
        if len(self.vessels) == 1:
            bcPositions = []
            for Id, bcs in self.boundaryConditions.iteritems():
                for bc in bcs:
                    bcPositions.append(bc.position)
            if 1 not in bcPositions and -1 not in bcPositions:
                print "ERROR: VascularNetwork.initialize(): BoundaryConditions are not proper defined Vessel {} at least one boundaryCondition at both ends! system exit".format(self.vessels[0].name)
                exit() 
        
        # initialize boundary conditions of type 1
        for Id, bcs in self.boundaryConditions.iteritems():
            for bc in bcs:
                try: bc.initialize({})
                except: pass
               
        windkesselExist = False
        for Id, bcs in self.boundaryConditions.iteritems():
            
            for bc in bcs:
                if bc.name in ['_Velocity-Gaussian', 'Velocity-Gaussian']: bc.update({'area':self.vessels[Id].A0})
                # relink the positionFunction
                if bc.type == 2: bc.setPosition(bc.position)
                # initialise windkessel # venousPressure
                if bc.name in ['_Windkessel-2Elements', '_Windkessel-2Elements', '_Windkessel-3Elements', 'Windkessel-3Elements']:
                    windkesselExist = True
                # initialise 
                if bc.name in ['VaryingElastanceHeart']:
                    try:
                        bc.mitral.rho = self.globalFluid['rho']
                    except: print "WARNING: VascularNetwork.initialize(): could not set blood density for mitral valve!"
                    try:
                        bc.aortic.rho = self.globalFluid['rho']
                    except: print "WARNING: VascularNetwork.initialize(): could not set blood density for aortic valve!"
        
        
        # # initialize 3d positions of the vascularNetwork
        self.calculate3DpositionsAndGravity(nSet=0) 
        
        # ## initialize for simulation
        if initializeForSimulation == True:
              
            # # initialize venous pressure and checks central venous pressure
            self.initializeVenousGravityPressure() 
            
            # # print 3D positions
            if self.quiet == False: self.print3D()
                         
            # calculate the cumulative network resistances and vessel resistances of the network
            if self.initialsationMethod != 'ConstantPressure':
                self.calculateNetworkResistance()    
                
            # calculate the initial values of the network
            self.calculateInitialValues()
            
            # show wave speed of network
            if self.quiet == False: self.showWaveSpeedOfNetwork()
            
            # optimize tree reflection coefficients BADDDDD
            if self.optimizeTree: self.optimizeTreeRefelctionCoefficients()
            
            if self.quiet == False: self.showReflectionCoefficientsConnectionInitialValues()
            
            if self.estimateWindkesselCompliance != 'No' and windkesselExist:
                # calculate terminal vessel compliance
                self.evaluateWindkesselCompliance()
                
    def initializeNetworkForSimulation(self):
        '''
        Method to initialize the network for a simulation.
        Creates hdf5 File and groups for the vessels
        Enforces memory allocation.
        Set initial values for the simulations.
        '''
            
        # create solution file
        if self.pathSolutionDataFilename == None:
            self.pathSolutionDataFilename = mFPH.getFilePath('solutionFile', self.name, self.dataNumber, 'write')
        self.solutionDataFile = h5py.File(self.pathSolutionDataFilename, "w")
        
        self.vesselDataGroup = self.solutionDataFile.create_group('vessels')
        
        # initialize saving indices
        if  self.timeSaveEnd < 0 or self.timeSaveEnd > self.totalTime:
            print "ERROR: VascularNetwork.initializeSolutionMatrices(): timeSaveEnd not in [0, totalTime], exit()"
            exit()
        
        if self.timeSaveBegin < 0 or self.timeSaveBegin > self.timeSaveEnd:
            print "WARNING: VascularNetwork.initializeSolutionMatrices(): timeSaveBegin not in [0, timeSaveEnd], exit()"
            exit()
        
        self.nSaveBegin = int(np.floor(self.timeSaveBegin / self.dt))
        self.nSaveEnd = int(np.ceil(self.timeSaveEnd / self.dt))
        # set save counter to the correct parts
        if self.initialisationPhaseExist:
            self.nSaveEnd += self.nTstepsInitPhase
            if self.saveInitialisationPhase:
                self.nSaveBegin = 0
            else:
                self.nSaveBegin += self.nTstepsInitPhase
        self.savedArraySize = self.nSaveEnd-self.nSaveBegin+1   
        self.nDCurrent = 0
        
        # # -> derive  number int(maxMemory / (vessels*3 arrays per vessel*vessel.N)) = memoryArraySizeTime    
        estimatedMemorySolutionDataSpace = 0
        for vessel in self.vessels.itervalues():
            estimatedMemorySolutionDataSpace += vessel.N * 8 * 3  # byte
        
        self.memoryArraySizeTime = int(np.floor(self.maxMemory * 1024.*1024. / estimatedMemorySolutionDataSpace))       
        # Don't allocate more memory than needed
        if self.memoryArraySizeTime > (self.nTsteps + 1):
            self.memoryArraySizeTime = self.nTsteps + 1
                
        # initialize for simulation
        for vesselId, vessel in self.vessels.iteritems():
            # initialize the vessel for simulation
            vessel.initializeForSimulation(self.initialValues[vesselId],
                                           self.memoryArraySizeTime,
                                           self.nTsteps)       
            # Put a reference to the dsetGroup into the saving dictionary if needed
            if vessel.save == True:
                # create a new group in the data file
                dsetGroup = self.vesselDataGroup.create_group(' '.join([vessel.name, ' - ', str(vessel.Id)]))
                nGridPoints =  vessel.N
                dsetP = dsetGroup.create_dataset("Pressure", (self.savedArraySize,nGridPoints), dtype='float64')
                dsetQ = dsetGroup.create_dataset("Flow", (self.savedArraySize,nGridPoints), dtype='float64')
                dsetA = dsetGroup.create_dataset("Area", (self.savedArraySize,nGridPoints), dtype='float64')
     
                if self.nSaveBegin==0:
                    dsetP[0] = vessel.Psol[0]
                    dsetQ[0] = vessel.Qsol[0]
                    dsetA[0] = vessel.Asol[0]
                    self.nDCurrent = 1
                    
                self.vesselsToSave[vesselId] =  dsetGroup      

        
        # # initialize varying elastance model
        # # initialize boundary condition type 1: initial phase
        for vesselId, boundaryConditions in self.boundaryConditions.iteritems():
            for bC in boundaryConditions:
                if bC.name in ['VaryingElastanceHeart', 'VaryingElastanceSimple']:
                    Qm = self.initialValues[vesselId]['Flow']
                    bC.update({'aorticFlowPreviousTimestep':Qm})
                    bC.initializeSolutionVectors(self.nTsteps)
                if bC.type == 1:
                    if self.initialisationPhaseExist:
                        bC.update({'initialisationPhaseExist': True,
                                   'nTstepsInitPhase': self.nTstepsInitPhase})
                    
        # # initialize gravity and 3d positions over time
        # create motion description out of motion dict of vascularNetwork
        # self.motion = [] # [ { vesselId : { angleXMother: ax, angleYMother: ay, angleZMotheraz }_n ] for all n in range (0,Tsteps-1)
        
        # define motion
        motionDict = {}
        headUpTilt = False
        # # head up tilt
        if headUpTilt == True:
            tSteps4 = int(self.nTsteps / 40.0)
            start = self.vessels[1].angleXMother
            end = start - 70 * np.pi / 180
            startAngle = np.ones(tSteps4 * 20.0) * start
            endAngle = np.ones(18.0*tSteps4) * end
            tiltAngle = np.linspace(start, end, self.nTsteps - 38. * tSteps4)
             
            angleXSystem = np.append(startAngle, np.append(tiltAngle, endAngle))
                     
            motionDict = {1:{'angleXMotherTime': angleXSystem}}
         
        for vesselId, angleDict in motionDict.iteritems():
            self.vessels[vesselId].update(angleDict)
            
        # # calculate gravity and positions   
        self.calculate3DpositionsAndGravity(nTsteps=self.nTsteps)
            
        # # calculate venous pressure for windkessel
        self.initializeVenousGravityPressureTime(self.nTsteps)
        
        
        ##
        for vesselId, vessel in self.vessels.iteritems():
            dsetGroup = self.vesselsToSave[vesselId]
            dsetPos = dsetGroup.create_dataset("PositionStart", (self.savedArraySize,3), dtype='float64')
            dsetRot = dsetGroup.create_dataset("RotationToGlobal", (self.savedArraySize,3,3), dtype='float64')
            dsetGravity = dsetGroup.create_dataset("NetGravity", (self.savedArraySize,1), dtype='float64')
            
            # TODO: Verify that numpy's index range protocal... it seems to cutoff the final value in the range selected.
            dsetPos[:] = vessel.positionStart[self.nSaveBegin:self.nSaveEnd+1]
            dsetRot[:] = vessel.rotToGlobalSys[self.nSaveBegin:self.nSaveEnd+1]
            del vessel.positionStart, vessel.rotToGlobalSys # free memory not used during simulation
            #TODO: do this nicer way
            vessel.positionStart  = np.zeros((1,3))         # instanteanoues position of vessel start point in the global system
            vessel.positionEnd    = np.zeros((1,3))         # instanteanoues position of vessel start point in the global system
            vessel.rotToGlobalSys = np.array([np.eye(3)])
            dsetGravity[:] = vessel.netGravity[self.nSaveBegin:self.nSaveEnd+1]
        
        self.BrxDataGroup = self.solutionDataFile.create_group('Baroreflex')
        print self.baroreceptors
        
        
        
    def flushSolutionMemory(self, currentTimeStep, currentMemoryIndex, chunkCount):
        """
        saving utility function to determine if solution data needs to be sent to the outputfile,
        and to calculate the correct indices between solution memory and the data output file.
        """
        
        """
        Explanation of index variables
        nCB,nCE where the beginning and end of the current solution data in memory would lie
         in a full time history of the solution.
        nMB,nME what indices of the current memory mark the beginning and end of what should be saved 
         nME is a slice index, i.e. position + 1
        nSB,nSE where the beginning and end of the current data to save would lie in the whole of the
         simulation. nSE is a slice index, i.e. position + 1
        nDB,nDE where does the current selection of data belong in the whole of the saved data
        """
        
        memoryArraySize = self.memoryArraySizeTime;
        offset = (memoryArraySize - 1) * chunkCount
        # indices mapping beginning and end of memory to the absolute number time steps in solution
        nCB = offset+1
        nCE = offset+memoryArraySize-1 # nCE == currentTimeStep+1
    
        # check if we need to save
        saving = not(nCE < self.nSaveBegin or nCB > self.nSaveEnd) # not(not saving)
        
        if saving:
            ## memory indices
            nMB = 1
            nSB = nCB
            if (self.nSaveBegin-nCB)>0:
                nMB = 1+(self.nSaveBegin-nCB)
                nSB = self.nSaveBegin
            
            #determine length to write
            # assume we write out through the end of memory
            lengthToWrite = self.memoryArraySizeTime - nMB 
                
            nME = memoryArraySize
            nSE = nCE + 1
            # correct this if save index is less than the current time step
            if (self.nSaveEnd-nCE)<0:
                # set the index to end saving
                nME -= (nCE - self.nSaveEnd)
                nSE -= (nCE - self.nSaveEnd)
                lengthToWrite -= (nCE - self.nSaveEnd) # -(-nME) as nME is negative

            nDB = self.nDCurrent
            nDE = self.nDCurrent+lengthToWrite
                                        
            self.nDCurrent += lengthToWrite
            
        # For vessels in the saving dictionary 
        for vesselId,dsetGroup in self.vesselsToSave.iteritems():
            # access each variable to save.
            # TODO: Is there a better way to define these in the vessel class
            vessel = self.vessels[vesselId]
            if saving:
                dsetGroup["Pressure"][nDB:nDE] = vessel.Psol[nMB:nME]
                dsetGroup["Flow"][nDB:nDE] = vessel.Qsol[nMB:nME]
                dsetGroup["Area"][nDB:nDE] = vessel.Asol[nMB:nME]
            # roll the end of the buffer
            vessel.Psol[0] = vessel.Psol[-1]
            vessel.Qsol[0] = vessel.Qsol[-1]
            vessel.Asol[0] = vessel.Asol[-1]
            
        for baro in self.baroreceptors.itervalues():
            baro.flushSolutionData(saving,nDB,nDE,nSB,nSE)
        
                    
    def saveSolutionData(self):
        '''
        # solution of the system over time 
        # {vesselID: { 'Psol' : [ [solution at N nodes]<-one array for each timePoint , ...  ], ..  }
        '''        
        globalData = self.solutionDataFile.create_group('VascularNetwork')
        
        globalData.attrs['dt'] = self.dt
        globalData.attrs['nTsteps'] = self.nTsteps
        globalData.attrs['nTstepsInitPhase'] = self.nTstepsInitPhase
        globalData.attrs['simulationDescription'] = self.description        
        
        savedArraySize = self.nSaveEnd - self.nSaveBegin + 1
        dsetTime = globalData.create_dataset('Time', (savedArraySize,), dtype='float64')
                
        # find start and end time of the time vector of the solution data
        if self.initialisationPhaseExist:
            if self.saveInitialisationPhase:
                startTime = -self.dt * self.nTstepsInitPhase
                endTime = self.dt * (self.nSaveEnd - self.nTstepsInitPhase)
            else:
                startTime = self.dt * (self.nSaveBegin - self.nTstepsInitPhase)
                endTime = self.dt * (self.nSaveEnd - self.nTstepsInitPhase)
        else:
            startTime = self.nSaveBegin * self.dt
            endTime = self.dt * self.nSaveEnd
            
        dsetTime[:] = np.linspace(startTime, endTime, savedArraySize).reshape(savedArraySize,)
        
        self.solutionDataFile.close()
        
        
    def linkSolutionData(self):
        '''
        This function prepares the solution data when the network is loaded
        assigning the appropriate information to allow the user to call 
        classVascularNetwork::loadSolutionDataRange to get specific values 
        loaded into memory.
        
        '''

        if self.pathSolutionDataFilename == None:
            self.pathSolutionDataFilename = mFPH.getFilePath('solutionFile', self.name, self.dataNumber, 'read')
        # TODO, what if this fails? do we know?
        self.solutionDataFile = h5py.File(self.pathSolutionDataFilename, "r")
        
        vesselId = None
        for groupName, group in self.solutionDataFile.iteritems():            
            if groupName == 'VascularNetwork':
                self.dt = group.attrs['dt']
                self.nTsteps = group.attrs['nTsteps']
                self.simulationTime = group['Time'][:]
                
            elif groupName == 'Baroreflex':
                # This works perfectly as long as the variables are the same in the group as in the class __init__
                for subGroupName, subGroup in group.iteritems():
                    baroId = int(subGroupName.split(' - ')[-1])
                    self.baroreceptors[baroId].update(subGroup)
            
            elif groupName == 'Heart':
                pass
            
            elif groupName == 'Vein':
                pass
            
            elif groupName == 'vessels' or '-' in groupName: # '-' is loads older hdf5 data files
                for subGroupName, subGroup in group.iteritems():
                    vesselId = int(subGroupName.split(' - ')[-1])
                    # try:
                    # link data
                    self.vesselsToSave[vesselId] = subGroup
                        # except: 
                            # print "WARNING: vascularNetwork.loadSolutionData() could not link solution data of vessel {}".format(vesselId)
                    # except: print "WARNING: could not read in solution data for vessel {}".format(groupName)
            else:
                print "classVascularNetwork::linkSolutionData() Unable to identify data group", groupName
        
        self.initialize()
        
               
    def _checkAccessInputs(self,t1,t2, mindt):
        
        # Check if the time span is valid
        startTime = self.simulationTime[0]
        endTime = self.simulationTime[-1]
        
        
        # Assume inputs are valid, otherwise flag invalid inputs
        inputsAreValid = True
        if t1>t2 :
            print 'ERROR:Invalid time range t1=', t1, '> t2=', t2
            inputsAreValid = False
                 
        if t1 < startTime :
            print 'ERROR:Invalid start time t1=', t1, 'before beginning of saved data t=', startTime
            inputsAreValid = False
        
        if t2 > endTime:
            print 'ERROR:Invalid end time t2=', t2, 'after end of saved data t=', endTime
            inputsAreValid = False
            
        if isNumber(mindt) and mindt > endTime - startTime:
            inputsAreValid = False
            print 'ERROR: Invalid minimum time step ', mindt , ' larger than solution time span.'
         
        if self.dt > t2-t1:
            inputsAreValid = False
            print 'ERROR: Invalid time range t2-t1=', t2-t1, 'is smaller than the solution time step dt'
          
        return inputsAreValid
    
    
    def getSolutionData(self,vesselId, variables, tvals, xvals):
        '''
        Get interpolated solution data
        Inputs:
        vesselId - the vessel from which the data is wanted
        variables - a list of strings with desired variables
            "Pressure",
            "Flow", 
            "Area", 
            "WaveSpeed", 
            "MeanVelocity",
            "ForwardFlow",
            "BackwardFlow",
            "ForwardPressure",
            "BackwardPressure"
        tvals - a numpy array (or python list) of times at which the values are desired
        xvals - a numpy array (or python list) of positions at which the values are desired
        
        Returns: A dictionary with keys corresponding to the input variables, and values are
            numpy arrays with rows corresponding to times(tvals) and columns corresponding to position(xvals) 
        '''
        
        tspan = [np.min(tvals),np.max(tvals)]
        mindt=None

        if "ForwardPressure" in variables or "BackwardPressure" in variables or "ForwardFlow" in variables or  "BackwardFlow" in variables:
            variables.append('linearWavesplit')
        
        self.loadSolutionDataRange([vesselId], tspan, mindt, variables) 
        
        data_dict = {}
        # Create Interpolating Function
        # interpolate.interp2d(self.tsol,self.vessels[vesselId].z,self.vessels,kind='linear',copy=False)
        if 'Pressure' in variables:
            interpfct= interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].Psol,kind='linear',copy=False)
            data_dict['Pressure'] = interpfct(xvals,tvals)
        if 'Flow' in variables:
            interpfct = interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].Qsol,kind='linear',copy=False)
            data_dict['Flow'] = interpfct(xvals,tvals)
        if  'Area' in variables:
            interpfct= interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].Asol,kind='linear',copy=False)
            data_dict['Area'] = interpfct(xvals,tvals) 
        if 'WaveSpeed' in variables:
            interpfct = interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].csol,kind='linear',copy=False)
            data_dict['WaveSpeed'] = interpfct(xvals,tvals) 
        if 'MeanVelocity' in variables:
            interpfct = interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].vsol,kind='linear',copy=False)
            data_dict['MeanVelocity'] = interpfct(xvals,tvals) 
        if 'ForwardPressure' in variables:
            interpfct  = interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].PsolF,kind='linear',copy=False) 
            data_dict['ForwardPressure'] = interpfct(xvals,tvals) 
        if 'BackwardPressure' in variables:
            interpfct = interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].PsolB,kind='linear',copy=False) 
            data_dict['BackwardPressure'] = interpfct(xvals,tvals) 
        if 'ForwardFlow' in variables:
            interpfct = interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].QsolF,kind='linear',copy=False)
            data_dict['ForwardFlow']  = interpfct(xvals,tvals) 
        if 'BackwardFlow' in variables:
            interpfct = interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].QsolB,kind='linear',copy=False)
            data_dict['BackwardFlow'] = interpfct(xvals,tvals) 
        return data_dict
    
    def loadSolutionDataRange(self, vesselIds = None, tspan=None, mindt=None, 
                                  values=["All",
                                  "Pressure",
                                  "Flow", 
                                  "Area", 
                                  "WaveSpeed", 
                                  "MeanVelocity",
                                  "Gravity",
                                  "Position",
                                  "Rotation"]):
        '''
        loads the solution data of the vessels specified into memory for the times 
            specified and drops any other previously loaded data.
        Inputs:
            vesselIds - a list of vessel Ids to load
                if vesselIds = None, data of all vessels is loaded
            tspan=[t1,t2] - a time range to load into memory t2 must be greater than t1.
                if tspan=None, all times are loaded
            values = a dictionary specifying which quantities to load entries keys are booleans and may be 'loadAll', 
                'loadPressure', 'loadArea', 'loadFlow', 'loadWaveSpeed', and 'loadMeanVelocity'. If 'All' 
                is in the list all quantities are loaded. Inputs are case insensitive.
            mindt := the minimum spacing in time between successively loaded points if 
                none is specified, the solution time step is used.
        Effects and Usage:
            loads the specified values into memory such that they may be accessed as
            vascularNetwork.vessels[vesselId].Pressure, etc, returning a matrix of 
            solution values corresponding to the time points in vascularNetwork.tsol.
            Accessing vessels and values not set to be loaded will produce errors.
        '''
        # Update loaded data tracking if inputs are valid
        # We could do this value = d.get(key, False) returns the value or False if it doesn't exist

        values = set(values) 
        if 'All' in values:
            values.update(["All",
                        "Pressure",
                        "Flow", 
                        "Area", 
                        "WaveSpeed", 
                        "MeanVelocity",
                        "linearWavesplit",
                        "Gravity",
                        "Position",
                        "Rotation"])
        elif 'WaveSpeed' in values:
            values.update(['Pressure', 'Area'])
        elif 'MeanVelocity' in values:
            values.update(['Pressure','Flow'])
        elif "linearWavesplit" in values:
            values.update(['Pressure','Flow','Area',"WaveSpeed"])
            
        if tspan is not None:
            t1 = tspan[0]
            t2 = tspan[1]
        else:
            t1 = self.simulationTime[0]
            t2 = self.simulationTime[-1]
        
        
        if self._checkAccessInputs(t1, t2, mindt):
            nSelectedBegin, nSelectedEnd =self.getFileAccessIndices(t1, t2)
            
            if isNumber(mindt):
                nTStepSpaces = int(np.ceil(mindt / self.dt))
            else:
                nTStepSpaces = 1
                
            self.tsol = self.simulationTime[nSelectedBegin:nSelectedEnd:nTStepSpaces]
            # check if all vessels should be loaded
            if vesselIds == None: vesselIds = self.vessels.keys()
            # Update selected vessels
            for vesselId in vesselIds:
                if vesselId in self.vesselsToSave:
                    vessel = self.vessels[vesselId]
                    dsetGroup = self.vesselsToSave[vesselId]
                    if dsetGroup['Pressure'].shape[1] == vessel.N:
                        del vessel.Psol, vessel.Qsol, vessel.Asol
                        # TODO Implement h5py direct_read method to improve speed
                        if 'Pressure' in values: 
                            vessel.Psol = dsetGroup['Pressure'][nSelectedBegin:nSelectedEnd:nTStepSpaces] 
                        if 'Flow' in values:
                            vessel.Qsol = dsetGroup['Flow'][nSelectedBegin:nSelectedEnd:nTStepSpaces]    
                        if  'Area' in values:
                            vessel.Asol = dsetGroup['Area'][nSelectedBegin:nSelectedEnd:nTStepSpaces] 
                        if 'WaveSpeed' in values:
                            #vessel.csol = vessel.waveSpeed(vessel.Asol,vessel.C(vessel.Psol))
                            vessel.postProcessing(['WaveSpeed'])
                        if 'MeanVelocity' in values:
                            #vessel.vsol = vessel.Qsol/vessel.Asol
                            vessel.postProcessing(["MeanVelocity"])
                        if "linearWavesplit" in values:
                            vessel.postProcessing(["linearWavesplit"])
                        if 'Gravity' in values:
                            try: vessel.netGravity = dsetGroup['NetGravity'][nSelectedBegin:nSelectedEnd:nTStepSpaces]
                            except: print "WARNING vascularNetwork.loadSolutionDataRange():  no netGravity stored in solutiondata file"
                        if 'Rotation' in values:
                            try: vessel.rotToGlobalSys = dsetGroup['RotationToGlobal'][nSelectedBegin:nSelectedEnd:nTStepSpaces]
                            except: print "WARNING vascularNetwork.loadSolutionDataRange():  no rotation matrices stored in solutiondata file"
                        if 'Position' in values:
                            try: vessel.positionStart = dsetGroup['PositionStart'][nSelectedBegin:nSelectedEnd:nTStepSpaces]
                            except: print "WARNING vascularNetwork.loadSolutionDataRange():  no positionStart stored in solutiondata file"
                    else:
                        print 'classVascularNetwork::loadSolutionDataRangeVessel Warning: vessel ', vesselId, 'not in saved data'
                        print "this is a very bad exception text, as it is raised if the saved number of gridpoints is different from the xml file"
            
        else:
            print 'classVascularNetwork::loadSolutionDataRangeVessel Error: Inputs were not valid you should not get here'
            exit()


        
    def getFileAccessIndices(self,t1,t2):
        '''
        Helper method to convert times to indices in the saved data.
        Input:
        t1,t2 the beginning and ending times to access
        Output:
        nSelectedBegin, nSelectedEnd - the indices corresponding to t1 and t2 in the file
        '''
        startTime = self.simulationTime[0]
        nSelectedBegin = int(np.floor((t1 - startTime) / self.dt))
        nSelectedEnd = int(np.ceil((t2 - startTime) / self.dt))+1
        return nSelectedBegin, nSelectedEnd
              
        
    def findRootVessel(self):
        '''
        Finds the root of a network, i.e. the vessel which is not a daughter of any vessel
        Evaluates a startRank for the evaulation of the network
        '''
        daughters = []
        approximatedBif = 0
        for vessel in self.vessels.itervalues():
            try:
                if vessel.leftDaughter != None: 
                    daughters.append(vessel.leftDaughter)
                try:
                    if vessel.rightDaughter != None: 
                        daughters.append(vessel.rightDaughter)
                        approximatedBif += 1
                except: pass
            except: pass
                
        # find startRank by approximation of numbers of generations
        approxGen = len(set(daughters)) - 2 * approximatedBif + int(np.sqrt(approximatedBif))
        self.startRank = 2.0 ** (approxGen - 1)
        # find root with difference between daughters and all vessels as root is never daughter
        roots = list(set(self.vessels.keys()).difference(daughters))
        try:
            self.root = roots[0]
        except: print "ERROR: vascularNetwork.searchRoot(): could not find a root node, system exit", exit()
        if len(roots) > 1:
            print "ERROR: vascularNetwork.searchRoot(): found several roots: {}, check network again system exit!".format(roots), exit()
        
    def checkDaughterDefinition(self):
        '''
        Method to check if all daughters are defined in the correct way, i.e. if a vessel has only 1 daughter 
        it should be defined as a leftDaughter, if it is defined as a rightDaughter, this method will rename it!
        additional check if there is a vessel with this id, if not remove daughter
        '''
        for vessel in self.vessels.itervalues():
            if vessel.leftDaughter == None and vessel.rightDaughter != None:
                print "WARNING vascularNetwork.checkDaughterDefiniton(): Vessel {} has no leftDaughter but a rightDaughter {}, this daughter is now assumed to be leftDaughter".format(vessel.Id, vessel.rightDaughter)
                vessel.leftDaughter = vessel.rightDaughter
                vessel.rightDaughter = None
            # check if daughter vessel exists
            if vessel.leftDaughter != None:
                try:
                    self.vessels[vessel.leftDaughter]
                except: 
                    print "WARNING: vascularNetwork.checkDaugtherDefinition():\n      leftDaughter with Id {} of vessel {} does not exist".format(vessel.leftDaughter, vessel.Id,)
                    vessel.leftDaughter = None 
                    
            if vessel.rightDaughter != None:
                try:
                    self.vessels[vessel.rightDaughter]
                except: 
                    print "WARNING: vascularNetwork.checkDaugtherDefinition():\n       rightDaughter with Id {} of vessel {} does not exist".format(vessel.rightDaughter, vessel.Id,)
                    vessel.rightDaughter = None
        
    def evaluateConnections(self):
        '''
        Method to evaluate all connections:
        
        - check for right daughter definition (call)
        - find root of the network (call)
        - evalualte all connections link, bifurcation, anastomosis
        - apply mothers to all vessels (call)
        
        Method traverses tree with defined daughters,
        - finds mothers and connections
        
        -> creates treeTraverseList breadth first traversing list
        -> creates treeTraverseConnections list of connections [ [LeftMother, rightMother, leftDaughter, rightDaughter], ..]
        '''
        
        # check for proper definition: if one daughter := leftDaughter ..
        self.checkDaughterDefinition()
        # find the current root
        self.findRootVessel()
        
        self.treeTraverseList = []
        self.treeTraverseConnections = []
        self.boundaryVessels = []
    
        root = self.root
        toVisit = []
        generation = 0
        rankGeneration = self.startRank
        ranking = {}
        mothers = {}
                
        if self.vessels[root].leftDaughter != None:  
            toVisit.append(root)  # Add root to the 'toVisit'-vessels if root has daughters:
            toVisit.append('nextGeneration')  # add nextGeneration marker                       
        else:   
            self.boundaryVessels.append(root)  # append root to ends as it has no left and right daughters  
                
        self.treeTraverseList.append(root)
        
        ranking[root] = rankGeneration
        rankGeneration = rankGeneration / 2.0
        
        # loop through tree until all daughters are conected
        while len(toVisit) != 0:                
                        
            # check if next generation has come
            motherVessel = toVisit.pop(0)
            if motherVessel == 'nextGeneration':
                try: motherVessel = toVisit.pop(0)
                except: break
                # set new brakepoint after the current generation
                toVisit.append('nextGeneration')
                generation += 1
                rankGeneration = rankGeneration / 2.0
            
            # current connection List reads [leftMother, rightMother, leftDaughter, rightDaughter]
            currentConnectionList = [motherVessel, None]  # at first each mother is assumed to be leftMother
            
            # Grab left daughter
            leftDaughter = self.vessels[motherVessel].leftDaughter   
            
            if leftDaughter != None:
                # adjust ranking
                rankingLeftDaughter = ranking[motherVessel] - rankGeneration
                
                # # check if exists in list (if so -> anastomsis!!)
                if leftDaughter not in self.treeTraverseList:
                    # # normal daughter: link or connection
                    # apply values to treeTraverseList, ranking, mothers, currentConnectionList
                    self.treeTraverseList.append(leftDaughter)
                    ranking[leftDaughter] = rankingLeftDaughter
                    mothers[leftDaughter] = [motherVessel]
                    currentConnectionList.append(leftDaughter)
                else:
                    # # either anastomosis or vessel has to moved to its real generation 
                    # 1.remove leftDaughter from treeTraversingList
                    self.treeTraverseList.remove(leftDaughter)
                    
                    existingMothers = mothers[leftDaughter]
                    existingRanking = ranking[leftDaughter]
                    if len(existingMothers) == 1:
                        
                        if existingMothers[0] == motherVessel:
                            # 2a.if the same mothers, just move it to its real generation and add it again
                            self.treeTraverseList.append(leftDaughter)
                            ranking[leftDaughter] = rankingLeftDaughter
                            currentConnectionList.append(leftDaughter)
                        else:
                            # 2b.  different mothers --> anastomosis!!!
                            #      check ranking: lower rank -> left mother; 
                            if existingRanking < rankingLeftDaughter:
                                # 2.1 existing is left mother, new ranking
                                self.treeTraverseList.append(leftDaughter)
                                mothers[leftDaughter] = [existingMothers[0], motherVessel]
                                self.treeTraverseConnections.remove([existingMothers[0], None, leftDaughter, None])
                                currentConnectionList = [existingMothers[0], motherVessel, leftDaughter, None]
                                ranking[leftDaughter] = rankingLeftDaughter
                                
                            elif existingRanking > rankingLeftDaughter:
                                # 2.2 existing is right mother, new ranking
                                self.treeTraverseList.append(leftDaughter)
                                mothers[leftDaughter] = [motherVessel, existingMothers[0]]
                                self.treeTraverseConnections.remove([existingMothers[0], None, leftDaughter, None])
                                currentConnectionList = [motherVessel, existingMothers[0], leftDaughter, None]
                                ranking[leftDaughter] = rankingLeftDaughter
                                
                            else:  # existingRanking == rankingLeftDaughter
                                # 2.3 existing is left mother, mean ranking 
                                self.treeTraverseList.append(leftDaughter)
                                mothers[leftDaughter] = [existingMothers[0], motherVessel]
                                self.treeTraverseConnections.remove([existingMothers[0], None, leftDaughter, None])
                                currentConnectionList = [existingMothers[0], motherVessel, leftDaughter, None]
                                ranking[leftDaughter] = (rankingLeftDaughter + existingRanking) / 2.0
                                                
                    elif len(existingMothers) == 2:
                        self.treeTraverseList.append(leftDaughter)
                        ranking[leftDaughter] = rankingLeftDaughter
                        currentConnectionList = [existingMothers[0], existingMothers[1], leftDaughter, None]
                
                # check if leftDaughter has also daughters which should be visualized
                if self.vessels[leftDaughter].leftDaughter != None:
                    toVisit.append(leftDaughter)
                else:
                    # append vessel to ends as it has no left and right daughters
                    if leftDaughter not in self.boundaryVessels: self.boundaryVessels.append(leftDaughter)
                
                rightDaughter = self.vessels[motherVessel].rightDaughter
                
                if rightDaughter != None:
                    # adjust ranking
                    rankingRightDaughter = ranking[motherVessel] + rankGeneration
                    
                    # # check if exists in list (if so -> anastomsis!!)
                    if rightDaughter not in self.treeTraverseList:
                        # # normal daughter: link or connection
                        # apply values to treeTraverseList, ranking, mothers, currentConnectionList
                        self.treeTraverseList.append(rightDaughter)
                        ranking[rightDaughter] = rankingRightDaughter
                        mothers[rightDaughter] = [motherVessel]
                        currentConnectionList.append(rightDaughter)
                    else:
                        # # either anastomosis or vessel has to moved to its real generation 
                        # 1.remove leftDaughter from treeTraversingList
                        self.treeTraverseList.remove(rightDaughter)
                        
                        existingMothers = mothers[rightDaughter]
                        existingRanking = ranking[rightDaughter]
                        if len(existingMothers) == 1:
                            
                            if existingMothers[0] == motherVessel:
                                # 2a.if the same mothers, just move it to its real generation and add it again
                                self.treeTraverseList.append(rightDaughter)
                                ranking[rightDaughter] = rankingRightDaughter
                                currentConnectionList.append(rightDaughter)
                            else:
                                print "ERROR right daughter forced to anastomosis, not possible"
                            
                        elif len(existingMothers) == 2:
                            self.treeTraverseList.append(rightDaughter)
                            ranking[rightDaughter] = rankingRightDaughter
                            currentConnectionList = [existingMothers[0], existingMothers[1], rightDaughter, None]
                    
                    # check if rightDaughter has also daughters which should be visualized
                    if self.vessels[rightDaughter].leftDaughter != None: toVisit.append(rightDaughter)
                    else: 
                        if rightDaughter not in self.boundaryVessels: self.boundaryVessels.append(rightDaughter)
                        # append vessel to ends as it has no left and right daughters
                        
                else:
                    if len(currentConnectionList) == 3:
                        currentConnectionList.append(None)
                                           
            if len(currentConnectionList) == 4:
                # check if already in list -> remove it
                if currentConnectionList in self.treeTraverseConnections : self.treeTraverseConnections.remove(currentConnectionList)
                # add current list
                self.treeTraverseConnections.append(currentConnectionList)
        
        self.applyMothersToVessel()
        
    def applyMothersToVessel(self):
        '''
        Functions traverses the self.treeTraverseConnections and saves the id of the
        left and right mother of the vessel
        '''        
        for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:  
            
            self.vessels[leftDaughter].leftMother = leftMother
            self.vessels[leftDaughter].rightMother = rightMother
            try:
                self.vessels[rightDaughter].leftMother = leftMother
                self.vessels[rightDaughter].rightMother = rightMother
            except: pass
    
    def findStartAndEndNodes(self):
        '''
        Function traverses self.treeTraverseConnections and creates start- and
        end-nodes for all vessels in the network 
        '''
        nodeCount = 0
        self.vessels[self.root].startNode = nodeCount
        # add end node for root vessel
        nodeCount += 1
        self.vessels[self.root].endNode = nodeCount
        
        # # add rest of the vessels by traversing the connections
        for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:  
            # # link
            if rightMother == None and rightDaughter == None:
                # set start of LD
                self.vessels[leftDaughter].startNode = self.vessels[leftMother].endNode
                # set end of LD
                nodeCount += 1
                self.vessels[leftDaughter].endNode = nodeCount
                
            # # bifurcation
            elif rightMother == None:
                # set start of LD & RD
                self.vessels[leftDaughter].startNode = self.vessels[leftMother].endNode
                self.vessels[rightDaughter].startNode = self.vessels[leftMother].endNode
                # set end of LD
                nodeCount += 1
                self.vessels[leftDaughter].endNode = nodeCount
                # set end of RD
                nodeCount += 1
                self.vessels[rightDaughter].endNode = nodeCount
                
            # # anastomosis
            elif rightDaughter == None:
                # end node of right is changed to the one of the left
                self.vessels[rightMother].endNode = self.vessels[leftMother].endNode
                # set start if LD
                self.vessels[leftDaughter].startNode = self.vessels[leftMother].endNode
                # set end of LD
                nodeCount += 1
                self.vessels[leftDaughter].endNode = nodeCount

            
    def calculateNetworkResistance(self):
        '''
        This function travers the network tree and calculates the 
        cumultative system resistances Rcum for each vessel in the Network.
        '''
                
        # # travers tree and create the R_cum list including the cumltative values
        for vesselId in self.treeTraverseList:
            if vesselId in self.boundaryVessels:
                boundaryResistance = 0
                for bc in self.boundaryConditions[vesselId]:
                    # # if Rtotal is not given evaluate Rtotal := Rc + Zc_vessel
                    try:
                        # # windkessel 3 elements
                        if bc.Rtotal == None:
                            if bc.Z == 'VesselImpedance':
                                P = np.ones(self.vessels[vesselId].N) * self.vessels[vesselId].Ps  # 158.747121018*133.32 #97.4608013004*133.32#
                                compliance = self.vessels[vesselId].C(P)
                                area = self.vessels[vesselId].A(P)
                                waveSpeed = self.vessels[vesselId].c(area, compliance)
                                Z = 1.0 / (compliance * waveSpeed)[-1]
                            else: Z = bc.Z
                            Rtotal = bc.Rc + Z
                            bc.update({'Rtotal':Rtotal})
                            print "vessel {} : estimated peripheral windkessel resistance (Rtotal) {}".format(vesselId, Rtotal / 133.32 * 1.e-6)
                    except: pass
                    # # add resistance to the value
                    try: boundaryResistance = boundaryResistance + bc.Rtotal
                    except:
                        # # winkessel 2 elements and single resistance
                        try:
                            if bc.Rc == 'VesselImpedance':
                                P = np.ones(self.vessels[vesselId].N) * self.vessels[vesselId].Ps  # 158.747121018*133.32 #97.4608013004*133.32#
                                compliance = self.vessels[vesselId].C(P)
                                area = self.vessels[vesselId].A(P)
                                waveSpeed = self.vessels[vesselId].c(area, compliance)
                                Z = 1.0 / (compliance * waveSpeed)[-1]
                                boundaryResistance = boundaryResistance + Z
                        except: pass
                        try:
                            # # winkessel 2 elements and single resistance
                            boundaryResistance = boundaryResistance + bc.Rc
                        except: pass
                    
                # print 'boundaryResistance',boundaryResistance/133.32*1.e-6
                if boundaryResistance == 0: 
                    print "\n Boundary Condition at end of vessel {} has no resistance".format(vesselId)
                    # # set boundaryresistance to 1/133.32*1.e6
                    print "The resistance is set to 1*133.32*1.e6 \n"
                    boundaryResistance = 1.*133.32 * 1.e6
                    
                self.Rcum[vesselId] = self.vessels[vesselId].resistance + boundaryResistance
            else:
                self.Rcum[vesselId] = None
        
        # # travers trhough the connections backwards to evaluate the cumulative resistance
        for leftMother, rightMother, leftDaughter, rightDaughter in reversed(self.treeTraverseConnections):  
            # # link
            if rightMother == None and rightDaughter == None:
                self.Rcum[leftMother] = self.vessels[leftMother].resistance + self.Rcum[leftDaughter]
            # # bifurcation
            elif rightMother == None:
                self.Rcum[leftMother] = self.vessels[leftMother].resistance + 1.0 / (1.0 / self.Rcum[leftDaughter] + 1.0 / self.Rcum[rightDaughter])
            # # anastomosis
            elif rightDaughter == None:
                print "\n WARNING: no method for resistance calculation for anastomosis is implemented!!! \n"
                
        
    
    def calculateInitialValues(self):
        '''
        This function travers the network tree and calculates the 
        esitimates the initial flow and pressure values for each vessel in the Network
        based on the meanflow/pressure value at the root node using the cumultative resistance
        '''
        
        initialValues = {}
                    
        root = self.root
        meanInflow = None
        meanInPressure = None
        
        
        ## find root inflow boundary condition, ie. bc condition with type 1:
        # varying elastance is type 2 and is only initialized with constant pressure
        inflowBoundaryCondition = None
        for bc in self.boundaryConditions[root]:
            if bc.type == 1:
                inflowBoundaryCondition = bc
        
        if self.venousSystemCollaps == True and self.initialsationMethod != 'ConstantPressure':
            print '\nERROR: Auto, MeanFlow, Mean Pressure: initialization not implemented for collapsing venous system! \n'
            exit()
        
        if self.initialsationMethod == 'Auto':
            try:
                meanInflow, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(quiet=self.quiet)
                self.initialisationPhaseExist = True
            except:
                print "Error: classVascularNetwork: Unable to calculate mean flow at inflow point"
                exit()
                
        elif self.initialsationMethod == 'MeanFlow':
            try:
                meanInflow = self.initMeanFlow
                # # addjust bc condition
                xxx, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(meanInflow, quiet=self.quiet)
                self.initialisationPhaseExist = True
            except:
                print "Error: classVascularNetwork: Unable to set given meanFlow at inflow point"
                exit()
                
        elif self.initialsationMethod == 'MeanPressure':
            try:
                meanInPressure = self.initMeanPressure
                self.initialisationPhaseExist = True
            except:
                print "Error: classVascularNetwork: Unable to set given meanFlow at inflow point"
                exit()
        
        elif self.initialsationMethod == 'ConstantPressure':
            
            constantPressure = self.initMeanPressure
            try:
                constantPressure = self.initMeanPressure
                if inflowBoundaryCondition != None:
                    xxx, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(0.0, quiet=self.quiet)
                
                self.initialisationPhaseExist = False 
                if self.initPhaseTimeSpan > 0:
                    self.initialisationPhaseExist = True 
                     
            except:
                print "Error: classVascularNetwork: Unable to evaluate time shift to 0 at inflow point"
                exit()
            #############################Inititalisation Method constant pressure #############
            initialValues[root] = {}
            initialValues[root]['Pressure'] = [constantPressure, constantPressure]
            initialValues[root]['Flow'] = 0
            
            # # set initial values of the vessels by traversing the connections
            for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:  
                calcDaughters = [leftDaughter]
                if rightDaughter != None: calcDaughters.append(rightDaughter)
                for daughter in calcDaughters:
                    initialValues[daughter] = {}
                    initialValues[daughter]['Pressure'] = [constantPressure, constantPressure]
                    initialValues[daughter]['Flow'] = 0
            
            # # adjust pressure with venous pressure
            for initialArray in initialValues.itervalues():
                initialArray['Pressure'][0] = initialArray['Pressure'][0] + self.centralVenousPressure
                initialArray['Pressure'][1] = initialArray['Pressure'][1] + self.centralVenousPressure          
            
            # # Check if gravity is on and if user would correct for hydrostatic pressure
            if self.gravitationalField == True:
                
                input = ['K']
                while input not in ['y', 'n']:
                    input = str(raw_input("\n Adjust for hydrostatic pressure(y/n): "))
                   
                if input == 'y':  # 'y' Adjust ConstantPressure to correct for hydrostatic pressure
                    initialValuesWithGravity = self.initializeGravityHydrostaticPressure(initialValues, root)                   
                    self.initialValues = initialValuesWithGravity
                    
                else:  # # if input is 'n'
                    self.initialValues = initialValues
            else:  # with no gravity       
                self.initialValues = initialValues
            return
            
         
        #############################Inititalisation Method Tree travers######################
        
        ###### initialize refelctionCoefficientTimeVarying --> move to boundary ? condition ?
        bcdict = {}
        for boundaryCondition in self.boundaryConditions[root]:
            if boundaryCondition.type == 1:
                bcdict = boundaryCondition.__dict__
                
        for boundaryCondition in self.boundaryConditions[root]:
            if boundaryCondition.name == 'ReflectionCoefficientTimeVarying':
                boundaryCondition.update(bcdict)
        ######
        
        if self.venousSystemCollaps == True:
            "WARNING: no method  for venous collapsing system is implemented to initialize network with method 'Tree' !!!"
            pass
        
        if meanInflow != None:
            p0 = self.Rcum[root] * meanInflow 
            p1 = p0 - self.vessels[root].resistance * meanInflow
            
        elif meanInPressure != None:
            meanInflow = meanInPressure / self.Rcum[root]  # calculate mean flow
            p0 = meanInPressure  
            # # addjust bc condition
            try:    xxx, self.initPhaseTimeSpan = self.boundaryConditions[root][0].findMeanFlowAndMeanTime(meanInflow, quiet=self.quiet)
            except:
                print "Error: VascularNetwork: Unable to adjust calculated meanFlow at inflow point boundary condition !"
                exit()
            p1 = p0 - self.vessels[root].resistance * meanInflow
        else: 
            print "Error: Neither flow or pressure value given at inflow point!"
            exit()
        
        initialValues[root] = {}
        initialValues[root]['Pressure'] = [p0, p1]
        initialValues[root]['Flow'] = meanInflow
                
        # # calculate initial values of the vessels by traversing the connections
        for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:  
            
            # # link & bifurcation
            if rightMother == None:
                p0 = initialValues[leftMother]['Pressure'][1]
                calcDaughters = [leftDaughter]
                # check for bifurcation
                if rightDaughter != None:
                    calcDaughters.append(rightDaughter)
                for daughter in calcDaughters:
                    qm = p0 / self.Rcum[daughter]
                    p1 = p0 - self.vessels[daughter].resistance * qm
                    
                    initialValues[daughter] = {}
                    initialValues[daughter]['Pressure'] = [p0, p1]
                    initialValues[daughter]['Flow'] = qm  
                                        
            # # anastomosis
            elif rightDaughter == None:
                print "WARNING: VascularNetwork: no method for anastomosis is implemented to initialize network !!!"
                pass

        
        # # adjust pressure with venous pressure
        for initialArray in initialValues.itervalues():
            initialArray['Pressure'][0] = initialArray['Pressure'][0] + self.centralVenousPressure
            initialArray['Pressure'][1] = initialArray['Pressure'][1] + self.centralVenousPressure
        
        # # adjust pressure for gravity pressure 
        initialValuesWithGravity = self.initializeGravityHydrostaticPressure(initialValues, root)
                      
        self.initialValues = initialValuesWithGravity
        
    def evaluateWindkesselCompliance(self):
        
        self.TotalVolumeComplianceTree = 0.0
        self.totalTerminalAreaCompliance = 0.0
        for vesselId, vessel_i in self.vessels.iteritems():
            # vessel_i = self.vessels[vesselId]
            
            p0, p1 = self.initialValues[vesselId]['Pressure']
            initialPressure = np.linspace(p0, p1, int(vessel_i.N))
            C = vessel_i.C(initialPressure)
            if vesselId in self.boundaryVessels:
                self.totalTerminalAreaCompliance = self.totalTerminalAreaCompliance + C[-1]
            
            self.Cends[vesselId] = C[-1]
            
            Cvol = sum((C[1::] + C[0:-1]) / 2.0) * vessel_i.dz[0]  # ## works only if equidistant grid
            # Cvol = C[-1]*vessel_i.length
            
            # print sum(C[1:-1])*vessel_i.dz[0], Cvol2, C[0]*vessel_i.length
            # print C[-1]*vessel_i.length - Cvol2
            # print sum(C[1:-1])*vessel_i.dz[0]- C[-1]*vessel_i.length
            # print 'Cvol ',vesselId, ' ',Cvol,' ',C[-1]
            self.TotalVolumeComplianceTree = self.TotalVolumeComplianceTree + Cvol
                    
        # # calculate C_wkTotal according to choosen method
        if self.estimateWindkesselCompliance == 'System':
            C_wkTotal = self.compTotalSys - self.TotalVolumeComplianceTree
        
        elif self.estimateWindkesselCompliance == 'Tree':
            a = self.compPercentageTree
            C_wkTotal = (1. - a) / a * self.TotalVolumeComplianceTree
        
        elif self.estimateWindkesselCompliance == 'Wk3':
            b = self.compPercentageWK3
            C_wkTotal = b / (1 - b) * self.TotalVolumeComplianceTree
        else:
            print "Error: VascularNetwork in calculating C_wkTotal!"
            exit()
        
        if self.quiet == False:
            print '====================================='
            print '__________total compliances________'
            print '               Compliance'
            print "TerminalArea     {:5.3}".format(self.totalTerminalAreaCompliance * 133.32 * 1.e6)
            print "TreeVolume       {:5.3}".format(self.TotalVolumeComplianceTree * 133.32 * 1.e6)
            print "Total System     {:5.3}".format((self.TotalVolumeComplianceTree + C_wkTotal) * 133.32 * 1.e6)
            print "Total WK's       {:5.3}".format(C_wkTotal * 133.32 * 1.e6)
        
        wk3CompPrintList = {}
        # calculate wk3Compliance and apply it to boundaryCondition
        for vesselId in self.boundaryVessels:
            wk3Compliance = C_wkTotal * self.Cends[vesselId] / self.totalTerminalAreaCompliance
            if self.boundaryConditions[vesselId][-1].name in ['_Windkessel-3Elements', 'Windkessel-3Elements']:
                Cdef = self.boundaryConditions[vesselId][-1].C
                self.boundaryConditions[vesselId][-1].C = wk3Compliance
                Rt = self.boundaryConditions[vesselId][-1].Rtotal
                wk3CompPrintList[vesselId] = [Rt / 133.32 * 1.e-6, wk3Compliance * 133.32 * 1.e6 * 1e5, Cdef * 133.32 * 1.e6 * 1e5]
                          
                #### set Z to Z   = 'VesselImpedance'
                # self.boundaryConditions[vesselId][-1].Z   = 'VesselImpedance'
                                
                if wk3Compliance < 0:
                    print "ERROR: Windkessel Compliance at vessel {}:  {} < 0!".format(vesselId, wk3Compliance)
                    exit()
        if self.quiet == False:
            print '________estimated compliances________'
            print ' vesselId       Rt       C     Cdef'
            for vesselId in self.vessels.keys():
                try: print "{:3} {:10.3f} {:10.3f} {:10.3f}".format(vesselId, wk3CompPrintList[vesselId][0], wk3CompPrintList[vesselId][1], wk3CompPrintList[vesselId][2])
                except: print "{:3}".format(vesselId)
        
    def calculateReflectionCoefficientConnection(self, mothers, daughters):
        '''
        Function calculates reflection coefficient of a vessel connection
        
        Input:
            motherVessels   = [ [Id mother1, pressure mother1] ...  ]
            daughterVessels = [ [Id daughter1, pressure daughter1] ...  ]
        
        Return: reflectionCoefficient
        '''
        
        admittanceM = 0
        admittanceD = 0
        
        for motherId, motherPressure in mothers:
            # calculate addmintance of current mother
            impedanceM = self.vessels[motherId].Impedance(motherPressure)
            admittanceM = admittanceM + 1.0 / impedanceM[-1]
        
        for daughterId, daughterPressure in daughters:
            # calculate addmintance of current daughter
            impedanceLD = self.vessels[daughterId].Impedance(daughterPressure)
            admittanceD = admittanceD + 1.0 / impedanceLD[0]
            
        # calculate reflection coefficient
        reflectionCoefficient = (admittanceM - admittanceD) / (admittanceM + admittanceD)
        # TransmissionCoeffLeftDaughter = (- AdmittanceM + AdmittanceD) / (AdmittanceM+AdmittanceD)        
        
        return reflectionCoefficient
        
    def optimizeTreeRefelctionCoefficients(self):
        '''
        Calculates the optimal reflection coeffiecients for the network
        
        addapted from article Reymond et al.2009
        (very poor and instable method)
        '''  
        
        # # add rest of the vessels by traversing the connections
        for leftMother, rightMother, leftDaughter, rightDaughter  in self.treeTraverseConnections:  
            #### to be changed 
            
            maxReflectionCoeff = 0.005  # values in reymonds code
            toleranceReflectionCoeff = 0.0  # values in reymonds code
            reflectionCoefficient = 10.0  # start value to get while running 
                  
            radiusLeftDaughterInit = self.vessels[leftDaughter].radiusProximal
            
            # print "connection:",leftMother,rightMother, leftDaughter, rightDaughter
            # while (abs(reflectionCoefficient)-maxReflectionCoeff) > toleranceReflectionCoeff:
            while abs(reflectionCoefficient) > maxReflectionCoeff or reflectionCoefficient < 0:    
                # # setup initial pressure for left mother
                p0, p1 = self.initialValues[leftMother]['Pressure']
                initialPressureLM = np.linspace(p0, p1, int(self.vessels[leftMother].N))
                try:
                    # # setup initial pressure for right daughter used if anastomosis     
                    p0, p1 = self.initialValues[rightMother]['Pressure']
                    initialPressureRM = np.linspace(p0, p1, int(self.vessels[rightMother].N))
                except:pass
                # # setup initial pressure for left daughter         
                p0, p1 = self.initialValues[leftDaughter]['Pressure']
                initialPressureLD = np.linspace(p0, p1, int(self.vessels[leftDaughter].N))
                # # setup initial pressure for right daughter used if bifurcation
                try: 
                    p0, p1 = self.initialValues[rightDaughter]['Pressure']
                    initialPressureRD = np.linspace(p0, p1, int(self.vessels[rightDaughter].N))
                except: pass
                # # calculate reflection coefficient
                if rightMother == None and rightDaughter == None:
                    reflectionCoefficient = self.calculateReflectionCoefficientConnection([[leftMother, initialPressureLM, ]],
                                                                                    [[leftDaughter, initialPressureLD]])
                elif  rightMother == None:
                    reflectionCoefficient = self.calculateReflectionCoefficientConnection([[leftMother, initialPressureLM, ]],
                                                                                    [[leftDaughter, initialPressureLD],
                                                                                     [rightDaughter, initialPressureRD]])
                elif  rightDaughter == None:
                    reflectionCoefficient = self.calculateReflectionCoefficientConnection([[leftMother, initialPressureLM],
                                                                                          [rightMother, initialPressureRM]],
                                                                                         [[leftDaughter, initialPressureLD]])
                # adjust daughter radii
                if reflectionCoefficient > maxReflectionCoeff:
                    for vesselId in [leftDaughter, rightDaughter]:
                        try:
                            self.vessels[vesselId].radiusProximal = self.vessels[vesselId].radiusProximal * 1.005
                            try: self.vessels[vesselId].radiusDistal = self.vessels[vesselId].radiusDistal * 1.005
                            except: pass
                            self.vessels[vesselId].initialize({})
                        except: pass
                else:
                    for vesselId in [leftDaughter, rightDaughter]:
                        try:
                            self.vessels[vesselId].radiusProximal = self.vessels[vesselId].radiusProximal * 0.995
                            try: self.vessels[vesselId].radiusDistal = self.vessels[vesselId].radiusDistal * 0.995
                            except: pass
                            self.vessels[vesselId].initialize({})
                        except: pass
            print " new Reflection Coeff area ratio", radiusLeftDaughterInit, self.vessels[leftDaughter].radiusProximal, 1 - (radiusLeftDaughterInit) / self.vessels[leftDaughter].radiusProximal
                # print "      new Reflection coefficient {}, areas".format(reflectionCoefficient), self.vessels[leftDaughter].radiusProximal #, self.vessels[rightDaughter].radiusProximal 
            # print
        
    def showReflectionCoefficientsConnectionInitialValues(self): 
        if self.quiet == False:
            print '====================================='
            print '________Reflection Coefficients______'
            print ' LM RM LD RD   Reflection coefficient'
        # # add rest of the vessels by traversing the connections
        for leftMother, rightMother, leftDaughter, rightDaughter  in self.treeTraverseConnections:  
                p0, p1 = self.initialValues[leftMother]['Pressure']
                initialPressureLM = np.linspace(p0, p1, int(self.vessels[leftMother].N))
                try:
                    # # setup initial pressure for right daughter used if anastomosis     
                    p0, p1 = self.initialValues[rightMother]['Pressure']
                    initialPressureRM = np.linspace(p0, p1, int(self.vessels[rightMother].N))
                except:pass
                # # setup initial pressure for left daughter         
                p0, p1 = self.initialValues[leftDaughter]['Pressure']
                initialPressureLD = np.linspace(p0, p1, int(self.vessels[leftDaughter].N))
                # # setup initial pressure for right daughter used if bifurcation
                try: 
                    p0, p1 = self.initialValues[rightDaughter]['Pressure']
                    initialPressureRD = np.linspace(p0, p1, int(self.vessels[rightDaughter].N))
                except: pass
                # # calculate reflection coefficient
                if rightMother == None and rightDaughter == None:
                    reflectionCoefficient = self.calculateReflectionCoefficientConnection([[leftMother, initialPressureLM, ]],
                                                                                    [[leftDaughter, initialPressureLD]])
                elif  rightMother == None:
                    reflectionCoefficient = self.calculateReflectionCoefficientConnection([[leftMother, initialPressureLM, ]],
                                                                                    [[leftDaughter, initialPressureLD],
                                                                                     [rightDaughter, initialPressureRD]])
                elif  rightDaughter == None:
                    reflectionCoefficient = self.calculateReflectionCoefficientConnection([[leftMother, initialPressureLM],
                                                                                          [rightMother, initialPressureRM]],
                                                                                         [[leftDaughter, initialPressureLD]])
                if rightMother == None: rightMother = '-'
                if rightDaughter == None: rightDaughter = '-'
                print "{:3} {:3} {:3} {:3}      {:.4}".format(leftMother, rightMother, leftDaughter, rightDaughter, reflectionCoefficient)
        
        
        
    def showWaveSpeedOfNetwork(self, Pressure=None, Flow=None):
        print '====================================='
        print '__________initial wave speed_________'
        print ' vessel    wave speed c(Pinit)   A(Pinit)    As(Pinit)      Dw(Pinit)      Re(Pinit)'
        for vesselId, vessel in self.vessels.iteritems():
            if Pressure == None:
                # calc initial pressure
                p0, p1 = self.initialValues[vesselId]['Pressure']
                pressureVessel = np.linspace(p0, p1, int(vessel.N))
            else:
                pressureVessel = Pressure[vesselId]
                                       
            A = vessel.A(pressureVessel)
            C = vessel.C(pressureVessel)
            c = np.max(vessel.c(A, C))
            Dw = np.max(C / A)
            As = np.max(vessel.compliance.As)
            
            if Flow == None:
                v = self.initialValues[vesselId]['Flow'] / A
            else:
                v = Flow / A
                
            Re = np.max(np.sqrt(A / np.pi) * 2.0 * v / self.globalFluid['my'] * self.globalFluid['rho'])    
            print ' {:3}            {:5.4}            {:5.4}     {:5.4}     {:4.4}    {:5.0f}'.format(vesselId, c, np.max(A), As, Dw, Re)

                
    def initializeGravityHydrostaticPressure(self, initialValues, root):
        '''
        Traverse the tree and initialize the nodes with the steady state hydrostatic pressure distribution
        ''' 
        negativePressure = False
        
        # # root vessel   
        p0, p1 = initialValues[root]['Pressure']
        p1 = p1 + self.vessels[root].netGravity[0] * self.vessels[root].length
        
        if p1 < 0. : print """ERROR: classVascularNetwork.initializeGravityHydrostaticPressure(), \n 
                            calculated negative pressure in initialization of vessel {} with inital values {}""".format(root, [p0, p1]) ; exit()
                    
        initialValues[root]['Pressure'] = [p0, p1]

        # # traverse tree to calculate the pressure influence of gravity
        for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:  
            
            # # link & anastomosis
            calcDaughters = [leftDaughter]
            
            # # add for bifucation
            if rightDaughter != None: 
                calcDaughters.append(rightDaughter)
                
            for daughter in calcDaughters:
                # initial pressure gradiant due to viscos effects without gravity 
                initialPressureDiff = initialValues[daughter]['Pressure'][1] - initialValues[daughter]['Pressure'][0]
                # update p0 with new p1 from mother including gravity
                p0 = initialValues[leftMother]['Pressure'][1]
                
                p1 = p0 + initialPressureDiff + self.vessels[daughter].netGravity[0] * self.vessels[daughter].length
                
                initialValues[daughter]['Pressure'] = [p0, p1]
                
                if p1 < 0. : print """ERROR: classVascularNetwork.initializeGravityHydrostaticPressure(), \n 
                                    calculated negative pressure in initialization of vessel {} with inital values {}""".format(daughter, [p0, p1]) ; exit()
                              
        return initialValues
                    
            
    def initializeVenousGravityPressure(self):
        '''
        Calculate and initialze the venous pressure depending on gravity for the 2 and 3 element windkessel models
        '''
        self.venousSystemCollaps = False
        
        # calculate absolute and relative venous pressure at boundary nodes
        for vesselId in self.boundaryVessels:
            
            relativeVenousPressure = self.centralVenousPressure + self.globalFluid['rho'] * self.vessels[vesselId].positionEnd[0][2] * self.gravityConstant - self.vessels[vesselId].externalPressure
            
            if self.minimumVenousPressure != None:
                if round(relativeVenousPressure, 2) < round(self.minimumVenousPressure, 2):  # round off everthing after 2 decimal points x.xx
                    relativeVenousPressure = self.minimumVenousPressure
                    self.venousSystemCollaps = True
                    print 'Warning: Venous system showing collapsing dynamics! \n'
                        
            for bc in self.boundaryConditions[vesselId]:
                # update venous pressure at boundary nodes
                if bc.name in ['_Windkessel-2Elements', 'Windkessel-2Elements', '_Windkessel-3Elements', 'Windkessel-3Elements']:
                    bc.update({'venousPressure':relativeVenousPressure})
        
        # # print out of method
        if self.quiet == False:
            print '\n============================================================='
            print '_______________Venous Pressures _____________________________'
            print '%s %36.1f' % ('Central venous pressure:', round(self.centralVenousPressure, 2))
       
            if self.gravitationalField == True:
                # for vesselId in sorted(self.boundaryVessels):
                #    print '%s %2i %15s %20.1f' % ('Boundary vessel',vesselId,',relative pressure  :', venousPressure[vesselId]/133.32)
                # checks if venous system is collapsing(having less than minimum allowed negative pressure)
                if self.venousSystemCollaps == True:
                    print '\n'
                    print 'Warning: Venous system showing collapsing dynamics! \n'
        
    def print3D(self):
    
        # # print
        print '==========================================================================================  \n'
        print '__________________________Vessel Id: position, net Gravity________________________________'
     
        # traverse vascular network
        for vesselId in sorted(self.treeTraverseList):           
            # positionStart = self.vessels[vesselId].positionStart
            # positionEnd = self.vessels[vesselId].positionEnd
            # print 'Start position  : vessel  {} {:19.3f} {:20.3f} {:21.3f}'.format(vesselId, positionStart[0],   positionStart[1],   positionStart[2])
            # print 'End position    : vessel  {} {:19.3f} {:20.3f} {:21.3f}'.format(vesselId, positionEnd[0],     positionEnd[1],     positionEnd[2])
            if self.gravitationalField == True:
                print '%s %2i %19.3f' % ('Net gravity     : vessel ', vesselId, self.vessels[vesselId].netGravity[0])
        
    def calculate3DpositionsAndGravity(self, nTsteps=None, nSet=None):
        '''
        Initializing the position and rotation of each vessel in 3D space
        Initializing netGravity of the vessels.
        '''
        if nSet != None:
            nTsteps = 0
            
        for n in xrange(nTsteps+1):
        
            if nSet != None: n = nSet
        
            if n == 0:
                self.vessels[self.root].angleXMother = 90.*np.pi / 180.
                self.vessels[self.root].angleYMother = 0  # 45*np.pi/180.
                self.vessels[self.root].angleZMother = 0  # 45*np.pi/180.
                                    
            positionEndMother = np.zeros(3)
            rotToGlobalSysMother = np.eye(3)
            self.vessels[self.root].caculatePositionAndGravity(n, positionEndMother, rotToGlobalSysMother)
            
            for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:
                # initialize left daughter
                positionEndMother = self.vessels[leftMother].positionEnd[n]
                rotToGlobalSysMother = self.vessels[leftMother].rotToGlobalSys[n]
                self.vessels[leftDaughter].caculatePositionAndGravity(n, positionEndMother, rotToGlobalSysMother) 
                # initiaize right daughter
                if rightDaughter != None:                               
                    self.vessels[rightDaughter].caculatePositionAndGravity(n, positionEndMother, rotToGlobalSysMother) 
                
                if rightMother != None:
                    if np.sum(self.vessels[rightMother].positionEnd - self.vessels[leftMother].positionEnd) < 3.e-15:
                        print 'ERROR: 3d positions of anastomosis {} {} {} is not correct!'.format(leftMother, rightMother, leftDaughter)
        
    def initializeVenousGravityPressureTime(self, nTsteps):
        '''
        Calculate and initialze the venous pressure depending on gravity for the 2 and 3 element windkessel models
        '''
        
        self.venousSystemCollaps = False
        
        # calculate absolute and relative venous pressure at boundary nodes
        for vesselId in self.boundaryVessels:
            relativeVenousPressure = np.empty(nTsteps+1)
            for n in xrange(nTsteps+1):
                           
                relativeVP = self.centralVenousPressure + self.globalFluid['rho'] * self.vessels[vesselId].positionEnd[n][2] * self.gravityConstant - self.vessels[vesselId].externalPressure
                
                if self.minimumVenousPressure != None:
                    if round(relativeVP, 2) < round(self.minimumVenousPressure, 2):  # round off everthing after 2 decimal points x.xx
                        relativeVP = self.minimumVenousPressure
                        print 'Warning: Venous system showing collapsing dynamics! \n'
                        
                relativeVenousPressure[n] = relativeVP
                        
            # update bc
            for bc in self.boundaryConditions[vesselId]:
                # update venous pressure at boundary nodes
                if bc.name in ['_Windkessel-2Elements', 'Windkessel-2Elements', '_Windkessel-3Elements', 'Windkessel-3Elements']:
                    bc.update({'venousPressure':relativeVenousPressure})
                
