from __future__ import print_function, absolute_import
from builtins import range
from future.utils import iteritems, iterkeys, viewkeys, viewitems, itervalues, viewvalues
from builtins import input as input3
import sys
import os
from sympy.series.tests.test_order import test_eval
# set the path relative to THIS file not the executing file!

#import classWKinitializers as cWKinit

#cur = os.path.dirname(os.path.realpath(__file__))
#sys.path.append(cur + '/../')

from starfish.UtilityLib import classStarfishBaseObject as cSBO
from starfish.NetworkLib import classVessel as cVes
from starfish.NetworkLib import classBaroreceptor as cBRX
from starfish.NetworkLib import classVenousPool as classVenousPool
from starfish.UtilityLib import moduleFilePathHandler as mFPH

import numpy as np
import math
from scipy import interpolate
from copy import deepcopy
import pprint
import h5py
from starfish.NetworkLib.classBoundaryConditions import *

from starfish.UtilityLib import classRuntimeMemoryManager
import logging
logger = logging.getLogger('starfish')

class VascularNetwork(cSBO.StarfishBaseObject):
    """
    Class representing a vascular Network
    The vascular network consists out of vessels defined in classVessel::Vessel()
    Additional Topology, BoundaryConditions and the SimulationContext are saved.
    """
    solutionMemoryFields    = ["simulationTime", "arterialVolume"]
    solutionMemoryFieldsToSave = ["simulationTime", "arterialVolume"]

    def __init__(self, quiet=True):

        # # vascularNetwork variables to set via XML
        self.name = 'vascularNetwork'  # name of the network
        self.description = ''  # description of the current case
        self.dataNumber = 'xxx'  # data number of the network
        self.quiet = quiet  # bool to suppress output


        self.dsetGroup = None
        self.tiltAngle = None # Angle the network is tilted relative to supine position
        # keep track of time points loaded in memory
        self.tsol = np.zeros(0)
        self.arterialVolume = np.zeros(0)


        # running options
        self.cycleMode = False


        # simulation Context
        self.totalTime = 1.0  # simulation time in seconds
        self.dt = -1.0
        self.CFL = 0.85  # maximal initial CFL number
        
        self.nTSteps = None  # number of timesteps of the simulation case determined by the solver
        self.simulationTime = np.zeros(0)  # array with simulation Time
        self.currentMemoryIndex = None
        # TODO: Remove when refactored
        self.initDataManagement()

        # self.motion         = {'keyframe': [0, 0.1, 1.0],
        #                        'X1'     : [0, 45, 90]}
        # dict defining the movement by change of angle using keyframes
        # {'keyframe': [t0, t1, tend],
        #  'X1: [0, 45, 90]} ## <- correspond to 90 degree change of angleXtoMother of vessel 1
        self.motionAngles = {}

        # gravity controls
        self.gravitationalField = False  # bool, turn gravity on or off
        self.gravityConstant = -9.81  # earth gravity

        # the solver calibration
        self.rigidAreas = False  # # 'True' 'False' to change
        self.simplifyEigenvalues = False  #
        self.riemannInvariantUnitBase = 'Pressure'  # 'Pressure' or 'Flow'
        self.automaticGridAdaptation = True  # False True
        self.solvingSchemeField       = 'MacCormack_Flux' # MacCormack_Flux or MacCormack_Matrix
        self.solvingSchemeConnections = 'NonLinear'  # 'Linear'

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

        self.optimizeTree = False  # optimize areas of vessels to minimize reflections in root direction

        #dictionaries for network components
        self.vessels = {}  # Dictionary with containing all vessel data,  key = vessel id; value = vessel::Vessel()

        self.venousPool = None
        # Todo: figure out better handling of venous pressure. How to set this in xml?
        self.centralVenousPressure = 5*133.32
        
        self.heart = None

        self.boundaryConditions = {}

        self.globalFluid = {'my': 1e-6, 'rho': 1050., 'gamma': 2.0}  # dictionary containing the global fluid data if defined

        self.externalStimuli = {}

        self.baroreceptors = {}  # dictionary with baroreceptors

        self.communicators = {}  # dictionary with communicators, key = communicator id; values = {communicator data}

        # internally calculated variables
        self.root = None  # the root vessel (mother of the mothers)
        self.anastomosisExists = False
        self.boundaryVessels = []  # includes all vessels with terminal boundaryConditions (except start of root)
        self.treeTraverseList = []  # tree traverse list
        self.treeTraverseList_sorted = []
        self.treeTraverseConnections = []  # tree traversal list including connections [ LM, RM , LD, RD ]
        
        self.nodes = []
        self.connectionNodes = []

        self.initialValues = {}
        self.lumpedValues = {}
        self.Rcum = {}  # Dictionary with all cumulative resistances
        self.Cends = {}  # Dictionary with the area compliances of all terminating vessels (at ends)
        self.totalTerminalAreaCompliance = None  # the sum of all Cends
        self.TotalVolumeComplianceTree = None  # total volume compliance of all vessels

        # random variables TODO: move out of here to global class
        self.randomInputManager = None
        self.measurementRoutine = None

    def initDataManagement(self):
        """
        Refactoring all "data management" code and variables into "independent functions"
        """
        self.pathSolutionDataFilename = None
        self.timeSaveBegin = 0.0  # time when to start saving
        self.timeSaveEnd = None # time when to end saving
        self.maxMemory = 20  # maximum memory in MB
        self.saveInitialisationPhase = True  # bool to enable saving of the initPhase

        self.vesselsToSave = {}
        self.nSaveSkip = 1
        self.nSkipShift = 0
        self.minSaveDt = -1
        self.saveDt = None
        self.nSaveBegin = None
        self.nSaveEnd  = None
        self.savedArraySize = None
        self.nDCurrent = None
        self.memoryArraySizeTime = None  # memory array size for the time arrays
        self.solutionDataFile = None  # file name of the solution data

        self.runtimeMemoryManager = None

    # all classes concerning vessel
    def addVessel(self, vesselId=None, dataDict=False):
        """
        adds vessel to the Network
        if no id, a random id is chosen
        if no DataDict, no values are assigned
        """
        # set id to 1 + highest id of existing vessels
        if vesselId == None:
            try: vesselId = max(self.vessels.keys()) + 1
            except: vesselId = 0

        # check Id
        if vesselId not in self.vessels:
            vessel = cVes.Vessel(Id=vesselId , name=('vessel_' + str(vesselId)))  # create vessel with given variables
            if dataDict:
                vessel.update(dataDict)  # set vesselData if available
            self.vessels[vessel.Id] = vessel  # add vessel to network
        else:
            self.warning("vascularNetwork.addVessel: vessel with Id {} exists already! Could not add vessel".format(vesselId), noException= True)

    def deleteVessel(self, inputId):
        """
        Remove vessel from network and delete it
        """
        try:
            del self.vessels[inputId]
        except Exception:
            logger.warn("vascularNetwork.deleteVessel(): vessel with Id {} does not exist! Could not remove vessel".format(inputId))
            self.warning("vascularNetwork.deleteVessel(): vessel with Id {} does not exist! Could not remove vessel".format(inputId))

    def addBaroreceptor(self, baroId=None, dataDict=False):
        """
        adds vessel to the Network
        if no id, a random id is chosen
        if no DataDict, no values are assigned
        """
        # set id to 1 + highest id of existing vessels
        if baroId == None:
            try:
                baroId = max(self.baroreceptors.keys()) + 1
            except ValueError: # Empty dict
                baroId = 0

        if baroId not in self.baroreceptors:
            baroType = dataDict['modelName']
            instance = getattr(cBRX, baroType)(dataDict)
            self.baroreceptors[baroId] = instance

        else:
            self.warning("vascularNetwork.addBaroreceptor: baroreceptor with Id {} exists already! Could not add baroreceptor".format(baroId),noException=True)

    def update(self, vascularNetworkData):
        """
        updates the vascularNetwork data using a dictionary in form of
        vascularNetworkData = {'variableName': value}
        """
        for key, value in iteritems(vascularNetworkData):
            try:
                self.__getattribute__(key)
                self.__setattr__(key, value)
            except Exception:
                self.warning("vascularNetwork.update(): wrong key: %s, could not update vascularNetwork" % key)

    def getVariableValue(self, variableName):
        """
        Returns value of variable with name : variableName
        States Error if not such variable
        """
        try:
            return self.__getattribute__(variableName)
        except Exception:
            self.warning("vascularNetwork.getVariable() : VascularNetwork has no variable {}".format(variableName))

    def updateNetwork(self, updateDict):
        """
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
        """

        for dictName in ['vascularNetworkData']:
            try:
                self.update(updateDict[dictName])
            except Exception:
                pass #self.warning("old except: pass clause #1 in classVascularNetwork.updateNetwork", oldExceptPass= True)

        for dictName in ['globalFluid', 'communicators', 'externalStimuli']:
            try:
                self.getVariableValue(dictName).update(updateDict[dictName])
            except Exception:
                pass #self.warning("old except: pass clause #2 in classVascularNetwork.updateNetwork", oldExceptPass= True)

        if 'vesselData' in updateDict:
            for vesselId, vesselData in iteritems((updateDict['vesselData'])):
                try:
                    self.vessels[vesselId].update(vesselData)
                except KeyError:
                    self.addVessel(vesselId, vesselData)

        if 'baroreceptors' in updateDict:
            for baroId, baroData in iteritems((updateDict['baroreceptors'])):
                try:
                    self.baroreceptors[baroId].update(baroData)
                except KeyError:
                    self.addBaroreceptor(baroId, baroData)

    def showVessels(self):
        """
        writes the Vesseldata for each vessel to console (calls printToConsole() from each vessel)
        """
        print(" Vessels in Network:")
        for vessel in itervalues(self.vessels):
            vessel.printToConsole()

    def showNetwork(self):
        """
        writes Network properties (without vesselData) to console
        """
        print("-------------------")
        print(" vascularNetwork ", self.name, "\n")
        for variable, value in iteritems(self.__dict__):
            try:
                print(" {:<20} {:>8}".format(variable, value))
            except Exception: print(" {:<20} {:>8}".format(variable, 'None'))

    def initialize(self, initializeForSimulation=False):
        """
        Initializes vascular network: the compliance of the vessels and the position of the call function of boundary type 2
        Check if boundaryConditions and globalFluid properties are defined in a right manner;
        """
        # # refresh all connections and generate traversing lists
        self.evaluateConnections()

        # # checks if gravity is turned on
        if self.gravitationalField == False: self.gravityConstant = 0.

        # ## check global fluid properties
        for fluidItem, value in iteritems(self.globalFluid):
            if value == None:
                if fluidItem == 'dlt':
                    try:
                        gamma = self.globalFluid['gamma']
                        self.globalFluid['dlt'] = (gamma + 2.0) / (gamma + 1)
                    except Exception: #TODO: Should htis exception be propagated?
                        logger.error("ERROR: VascularNetwork.initialize(): global fluid properties are not properly defined! Check:" 
                                + '\n' 
                                + pprint.pformat(self.globalFluid) + '\n'
                                + 'Please fix network file and try again')
                        exit()
                else:
                    logger.error("ERROR: VascularNetwork.initialize(): global fluid properties are not properly defined! Check:" 
                                + '\n' 
                                + pprint.pformat(self.globalFluid) + '\n'
                                + 'Please fix network file and try again')
                    exit()

        # ## initialize vessels
        for vessel in itervalues(self.vessels):
            vessel.initialize(self.globalFluid)
            vessel.update({'gravityConstant': self.gravityConstant})

        # ## update wall models from measurment data
        if self.measurementRoutine is not None:
            self.measurementRoutine.adaptationToPatientSpecificCondition(self)

        # ## check and initialize boundary conditions
        if self.boundaryConditions != {}:
            # # Set position of boundary conditions
            # check position if one Vessel
            if len(self.vessels) == 1:
                vesselId = list(self.boundaryConditions.keys())[0]
                if vesselId != self.root: logger.error("Error Wrong Root found") #TODO: should this stop something?
                for bc in self.boundaryConditions[vesselId]:
                    if '_' not in bc.name[0]: bc.setPosition(0)
                    else: bc.setPosition(-1)
            else:
                for vesselId, bcs in iteritems(self.boundaryConditions):
                    if vesselId == self.root:
                        for bc in bcs:
                            bc.setPosition(0)
                            if "Elastance" in bc.name:
                                self.heart = bc

                    elif vesselId in self.boundaryVessels:
                        for bc in bcs:
                            bc.setPosition(-1)

        definedButNotatBC = set(self.boundaryConditions.keys()).difference([self.root] + self.boundaryVessels)
        atBCButNotDefined = set([self.root] + self.boundaryVessels).difference(self.boundaryConditions.keys())
        if len(definedButNotatBC.union(atBCButNotDefined)) > 0:
            tmpstring = "VascularNetwork.initialize(): BoundaryConditions are not properly defined:"
            if len(definedButNotatBC) > 0:tmpstring = tmpstring + "for Vessel(s) {} boundaryConditions are defined but \n   Vessel(s) is(are) not at the Boundary!".format(list(definedButNotatBC))
            if len(atBCButNotDefined) > 0:tmpstring = tmpstring + "for Vessel(s) {} no BoundaryConditions are defined!".format(list(atBCButNotDefined))
            logger.error(tmpstring)
            raise RuntimeError(tmpstring)
        if len(self.vessels) == 1:
            bcPositions = []
            for Id, bcs in iteritems(self.boundaryConditions):
                for bc in bcs:
                    bcPositions.append(bc.position)
            if 1 not in bcPositions and -1 not in bcPositions:
                error_msg = "VascularNetwork.initialize(): BoundaryConditions are not properly defined Vessel {} at least one boundaryCondition at both ends! system exit".format(self.vessels[0].name)
                logger.error(error_msg)
                raise RuntimeError(error_msg)

        # initialize boundary conditions of type 1
        for Id, bcs in iteritems(self.boundaryConditions):
            for bc in bcs:
                try: bc.initialize({})
                except Exception: pass #self.warning("old except: pass clause in VascularNetwork.initialize", oldExceptPass= True)

        windkesselExist = False
        for Id, bcs in iteritems(self.boundaryConditions):

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
                    except Exception:
                        self.warning("VascularNetwork.initialize(): could not set blood density for mitral valve!")
                    try:
                        bc.aortic.rho = self.globalFluid['rho']
                    except Exception:
                        self.warning("VascularNetwork.initialize(): could not set blood density for aortic valve!")

        # # initialize 3d positions of the vascularNetwork
        if self.anastomosisExists:
            logger.debug("WARNING: The network contain one or more anastomosis; 3DpositionsAndGravity will not be calculated. line 410 classVascularNetwork")
        else:
            self.calculate3DpositionsAndGravity(nSet=0)

        # ## initialize for simulation
        # TODO: Can this be moved?
        if initializeForSimulation == True:
            
            # # initialize venous pressure and checks central venous pressure
            self.initializeVenousGravityPressure()

            # # print(3D positions)
            if self.quiet == False:
                self.print3D()

            # calculate the cumulative network resistances and vessel resistances of the network
            if self.initialsationMethod not in ['ConstantPressure', 'AutoLinearSystem', 'FromSolution']:
                self.calculateNetworkResistance()

            # calculate the initial values of the network
            self.calculateInitialValues()

            if self.quiet == False:
                # evaluate the total arterial compiance and resistacne
                self.evaluateNetworkResistanceAndCompliance()

                # show wave speed of network
                #self.showWaveSpeedOfNetwork()

            # optimize tree reflection coefficients BADDDDD
            if self.optimizeTree:
                self.optimizeTreeRefelctionCoefficients()

            if self.quiet == False:
                self.showReflectionCoefficientsConnectionInitialValues()

            if self.estimateWindkesselCompliance != 'No' and windkesselExist:
                # calculate terminal vessel compliance
                self.evaluateWindkesselCompliance()

    def initializeNetworkForSimulation(self):
        """
        Method to initialize the network for a simulation.
        Creates hdf5 File and groups for the vessels
        Enforces memory allocation.
        Set initial values for the simulations.
        """

        # initialize saving indices
        self.timeSaveEnd = self.totalTime

        if self.timeSaveBegin < 0 or self.timeSaveBegin > self.timeSaveEnd:
            raise ValueError("VascularNetwork.initializeNetworkForSimulation(): timeSaveBegin not in [0, timeSaveEnd]")

        self.nSaveSkip = max(int(np.ceil(self.minSaveDt/self.dt)),1)
        self.saveDt = self.nSaveSkip*self.dt
        self.nSaveBegin = int(np.floor(self.timeSaveBegin / self.dt))
        self.nSaveEnd = int(np.ceil(self.timeSaveEnd / self.dt))

        # set save counter to the correct parts
        if self.initialisationPhaseExist:
            self.nSaveEnd += self.nTstepsInitPhase
            if self.timeSaveBegin > 0:
                self.nSaveBegin += self.nTstepsInitPhase
                self.saveInitialisationPhase = False
            else:
                self.saveInitialisationPhase = True
                self.nSaveBegin += self.nTstepsInitPhase

        self.savedArraySize = (self.nSaveEnd-self.nSaveBegin)//self.nSaveSkip + 1

        self.runtimeMemoryManager = classRuntimeMemoryManager.RuntimeMemoryManager(self.nSaveBegin,
                                                                                   self.nSaveEnd,
                                                                                   self.nSaveSkip,
                                                                                   self.nTSteps,
                                                                                   self.maxMemory)


        # Register all objects with the memory manager
        sizes = self.getSolutionMemorySizes()
        self.runtimeMemoryManager.registerDataSize(sizes)

        for vessel in itervalues(self.vessels):
            sizes = vessel.getSolutionMemorySizes()
            self.runtimeMemoryManager.registerDataSize(sizes)

        for bcList in itervalues(self.boundaryConditions):
            for bc in bcList:
                self.runtimeMemoryManager.registerDataSize(bc.getSolutionMemorySizes())

        for baroData in itervalues(self.baroreceptors):
            sizes = baroData.getSolutionMemorySizes()
            self.runtimeMemoryManager.registerDataSize(sizes)

        if self.venousPool is not None:
            sizes = self.venousPool.getSolutionMemorySizes()
            self.runtimeMemoryManager.registerDataSize(sizes)

        self.memoryArraySizeTime = self.runtimeMemoryManager.memoryArraySizeTime


        # Initialize solution file and data set groups
        # create solution file
        if self.pathSolutionDataFilename == None:
            self.pathSolutionDataFilename = mFPH.getFilePath('solutionFile', self.name, self.dataNumber, 'write')

        self.solutionDataFile = h5py.File(self.pathSolutionDataFilename, "w")
        self.dsetGroup = self.solutionDataFile.create_group('VascularNetwork')
        self.allocate(self.runtimeMemoryManager)

        # TODO: Integrate precalculated data into data saving framework
        self.dsetGroup.create_dataset('TiltAngle', (self.savedArraySize,),dtype='float64')
        self.tiltAngle = np.zeros(self.nTSteps)

        self.simulationTime[0] = -self.nTstepsInitPhase*self.dt


        logger.debug("cVN::InitializeNetworkForSimulation")
        logger.debug("nTSteps {}".format(self.nTSteps))
        logger.debug("Saving ={}:{}:{}".format(self.nSaveBegin,self.nSaveEnd,self.nSaveSkip))


        self.vesselDataGroup = self.solutionDataFile.create_group('vessels')

        # initialize objects for simulation
        for vesselId, vessel in iteritems(self.vessels):
            # initialize the vessel for simulation
            vessel.initializeForSimulation(self.initialValues[vesselId],
                                           self.runtimeMemoryManager,
                                           self.nTSteps,
                                           self.vesselDataGroup)



        for vesselId, boundaryConditions in iteritems(self.boundaryConditions):
            for bC in boundaryConditions:
                try:
                    bC.initializeSolutionVectors(self.runtimeMemoryManager, self.solutionDataFile)
                except AttributeError:
                    pass # bC doesn't have solution vector data
                bC.update({'initialisationPhaseExist': self.initialisationPhaseExist,
                                     'nTstepsInitPhase': self.nTstepsInitPhase})


        self.BrxDataGroup = self.solutionDataFile.create_group('Baroreflex')
        for baroData in itervalues(self.baroreceptors):
            baroData.initializeForSimulation(self)
        
        try:
            # Not all venous classes use this
            self.venousPool.initializeForSimulation(self)
        except AttributeError:
            logger.debug("Using static venous system")

        # # initialize gravity and 3d positions over time
        for stimulus in itervalues(self.externalStimuli):
            if stimulus['type'] == "headUpTilt":
                self.initializeHeadUpTilt(stimulus)


        # calculate gravity and positions
        if self.anastomosisExists:
            self.initializeVenousGravityPressureTime(self.nTSteps)
            logger.debug("WARNING: The network contain one or more anastomosis; lines 573-593 in classVascularNetwork will not be run")
        else:
            self.calculate3DpositionsAndGravity(nTsteps=self.nTSteps)
    
            # calculate venous pressure for windkessel
            self.initializeVenousGravityPressureTime(self.nTSteps)
    
            # Save gravity data if appropriate
            for vesselId, vessel in iteritems(self.vessels):
                dsetGroup = vessel.dsetGroup
                if dsetGroup:
                    dsetPos = dsetGroup.create_dataset("PositionStart", (self.savedArraySize,3), dtype='float64')
                    dsetRot = dsetGroup.create_dataset("RotationToGlobal", (self.savedArraySize,3,3), dtype='float64')
                    dsetGravity = dsetGroup.create_dataset("NetGravity", (self.savedArraySize,1), dtype='float64')
    
                    dsetPos[:] = vessel.positionStart[self.nSaveBegin:self.nSaveEnd+1:self.nSaveSkip]
                    dsetRot[:] = vessel.rotToGlobalSys[self.nSaveBegin:self.nSaveEnd+1:self.nSaveSkip]
                    dsetGravity[:] = vessel.netGravity[self.nSaveBegin:self.nSaveEnd+1:self.nSaveSkip]
    
                # TODO: Better way to return this to normal, while clearing the data?
                vessel.positionStart  = np.zeros((1,3))         # instantaneous position of vessel start point in the global system
                vessel.positionEnd    = np.zeros((1,3))         # instantaneous position of vessel end point in the global system
                vessel.rotToGlobalSys = np.array([np.eye(3)])

    class WholeBodyTilt(cSBO.StarfishBaseObject):
        """Encapsulates data related to the tilting motion of the network.

        A WholeBodyTilt object specifies the action of tilting the entire
        network about the root vessel of the network.

        Attributes:
            startTime  (float): The time in seconds when the tilt begins
            duration  (float): The length in second of the tilt
            stopAngle (float): The angle swept out by the tilt, relative to
             a supine position with a positive angle meaning an elevation
             of the feet above the head.
        """
        def __init__(self):
            self.startTime
            self.duration
            self.stopAngle

    def initializeHeadUpTilt(self, headUpTilt):
        """
        Takes a head upt tilt stimulus object and applies the specification to
        generate vessel positions over the simulation
        """

        tstart = headUpTilt['startTime']
        duration = headUpTilt['duration']
        tiltAngle = headUpTilt['stopAngle']

        tstop = tstart + duration

        # TODO: Is vessels[1] the root?
        start = self.vessels[1].angleXMother
        end = start + tiltAngle
        nStepsStart = int(math.floor(tstart/self.dt))
        nStepsTilt = int(math.ceil(tstop/self.dt)) - nStepsStart
        # TODO determine appropriate behaviour if simulation time is shorter that head up tilt time
        assert tstop < self.totalTime, 'tstop > totalTime'
        # nStepsEnd = ceil((self.nTSteps*self.dt - tstop)/self.dt)
        nStepsEnd = int(math.ceil(self.totalTime/self.dt))- nStepsTilt - nStepsStart

        startAngle = np.ones(nStepsStart)*start
        endAngle = np.ones(nStepsEnd)*end
        tiltAngle = np.linspace(start, end, nStepsTilt+1) #Account for time points not time steps
        angleXSystem = np.append(startAngle, np.append(tiltAngle, endAngle))
        # TODO: Why is the key "1" here?
        motionDict = {1:{'angleXMotherTime': angleXSystem}}

        self.tiltAngle = angleXSystem

        # TODO: Do these belong here? and do they need to happen every simulation?
        for vesselId, angleDict in iteritems(motionDict):
            self.vessels[vesselId].update(angleDict)

    def __call__(self):
        # Global Compliance
        # Global Impedance?
        nmem = self.currentMemoryIndex[0]
        self.simulationTime[nmem+1] = self.simulationTime[nmem] + self.dt
        
        # TODO: Pressure update assumes happening last
        if self.heart:
            if self.venousPool is not None:
                if len(self.venousPool.P_LA)>1:
                    self.heart.atriumPressure[nmem+1] = self.venousPool.P_LA[nmem+1]
                else:
                    self.heart.atriumPressure[nmem+1] = self.venousPool.pressureGain*self.venousPool.P[0]
            else:
                logger.debug("Using static atrial pressure")

        # TODO: Volume calculation assumes all other objects have been updated for the current time step!!!!
        self.arterialVolume[nmem+1] = self.calculateNetworkVolume(nmem+1)


    def calculateNetworkVolume(self, n):
        # Adds the volume of all compartments in the network
        cumVolume = 0.0
        for vesselId,vessel in iteritems(self.vessels):
            
            # access each variable to save.
            # TODO: Is there a better way to define these in the vessel class
            # vessel = self.vessels[vesselId]

            # calculate vessel volume
            A1 = self.vessels[vesselId].Asol[n:n+1,0:-1]
            A2 = self.vessels[vesselId].Asol[n:n+1,1:]
            volume = np.sum(vessel.dz*(A1+A2+np.sqrt(A1*A2))/3.0,axis=1)
            cumVolume += volume
        try:
            cumVolume += self.venousPool.V[n]
        except AttributeError:
            logger.debug("venous pool has no volume")

        # TODO: add heart handling
        # cumVolume += self.heart.V[n]
        return cumVolume


    def saveSolutionData(self):
        """
        # solution of the system over time
        # {vesselID: { 'Psol' : [ [solution at N nodes]<-one array for each timePoint , ...  ], ..  }
        """
        globalData = self.dsetGroup
        globalData.attrs['dt'] = self.dt
        globalData.attrs['nTSteps'] = self.nTSteps
        globalData.attrs['nTstepsInitPhase'] = self.nTstepsInitPhase
        globalData.attrs['simulationDescription'] = str(self.description)

        # dsetTime = globalData.create_dataset('Time', (savedArraySize,), dtype='float64')

        # TODO: Integrate this better with the chunking mechanism

        # dsetTime[:] = startTime + self.saveDt*np.arange(savedArraySize).reshape(savedArraySize,)

        self.solutionDataFile.close()



    def linkSolutionData(self):
        """
        This function prepares the solution data when the network is loaded
        assigning the appropriate information to allow the user to call
        classVascularNetwork::loadSolutionDataRange to get specific values
        loaded into memory.

        """

        if self.pathSolutionDataFilename == None:
            self.pathSolutionDataFilename = mFPH.getFilePath('solutionFile', self.name, self.dataNumber, 'read')
        # TODO, what if this fails? do we know?
        self.solutionDataFile = h5py.File(self.pathSolutionDataFilename, "r")

        vesselId = None
        for groupName, group in iteritems(self.solutionDataFile):
            if groupName == 'VascularNetwork':
                self.dt = group.attrs['dt']
                self.nTSteps = group.attrs['nTSteps']
                self.simulationTime = group['simulationTime'][:]

            elif groupName == 'Baroreflex':
                # This works perfectly as long as the variables are the same in the group as in the class __init__
                for subGroupName, subGroup in iteritems(group):
                    baroId = int(subGroupName.split(' - ')[-1])
                    try:
                        self.baroreceptors[baroId].update(subGroup)
                    except KeyError: # will fail for nested data
                        pass

            elif groupName == 'Heart':
                pass

            elif groupName == 'Vein':
                pass

            elif groupName == 'vessels': # or '-' in groupName: # '-' is loads older hdf5 data files
                for subGroupName, subGroup in iteritems(group):
                    vesselId = int(subGroupName.split(' - ')[-1])
                    self.vesselsToSave[vesselId] = subGroup
                    self.vessels[vesselId].dsetGroup = subGroup
            else:
                logger.warning("classVascularNetwork::linkSolutionData() Unable to identify data group {}".format(groupName))

        self.initialize()


    def _checkAccessInputs(self,t1,t2, mindt):
        """
        Checks to ensure the data requested actually exists.

        Args:
            t1 (float): initial time of data requested
            t2 (float): final time of data requested
            mindt (float): the minimum time separating data points requested
        Raises:
            ValueError: If t1 or t2 lie outside the range of simulationTime,
             mindt is larger than the range of simulationTime, or the intrinsic
             time step, dt, is larger than t2-t1.
        """

        # Check if the time span is valid
        startTime = self.simulationTime[0];
        endTime = self.simulationTime[-1]


        # TODO: should these be errors?
        # Assume inputs are valid, otherwise flag invalid inputs
        inputsAreValid = True
        if t1>t2 :
            raise ValueError("ERROR:Invalid time range t1=%f > t2=%f" % (t1,t2))
            inputsAreValid = False

        if t1 < startTime :
            raise ValueError("ERROR:Invalid start time t1=%f before beginning of saved data t=%f" % (t1,startTime))
            inputsAreValid = False

        if t2 > endTime:
            raise ValueError("ERROR:Invalid end time t2=%f after end of saved data t=%f" % (t2, endTime))
            inputsAreValid = False

        if mindt is not None and mindt > endTime - startTime:
            inputsAreValid = False
            raise ValueError("ERROR: Invalid minimum time step %f larger than solution time span." % (mindt))

        if self.dt > t2-t1:
            inputsAreValid = False
            raise ValueError("ERROR: Invalid time range t2-t1=%f is smaller than the solution time step dt" %(t2-t1))

        return inputsAreValid


    def getSolutionData(self,vesselId, variables, tvals, xvals):
        """
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
            'Compliance'
        tvals - a numpy array (or python list) of times at which the values are desired
        xvals - a numpy array (or python list) of positions at which the values are desired

        Returns: A dictionary with keys corresponding to the input variables, and values are
            numpy arrays with rows corresponding to times(tvals) and columns corresponding to position(xvals)
        """
        #TODO: return full non interpolated solution
        
        tspan = [np.min(tvals),np.max(tvals)]
        mindt=None
        waveSplittingVariables =  ["ForwardPressure","BackwardPressure", "ForwardFlow","BackwardFlow"]
        if any(i in variables for i in waveSplittingVariables):
        # if "ForwardPressure" in variables or "BackwardPressure" in variables or "ForwardFlow" in variables or  "BackwardFlow" in variables:
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
        if 'Compliance' in variables:
            interpfct = interpolate.interp2d(self.vessels[vesselId].z,self.tsol,self.vessels[vesselId].Csol,kind='linear',copy=False)
            data_dict['Compliance'] = interpfct(xvals,tvals)
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
                                  values=["All"]):
        """
        loads the solution data of the vessels specified into memory for the times
            specified and drops any other previously loaded data.
        Inputs:
            vesselIds - a list of vessel Ids to load
                if vesselIds = None, data of all vessels is loaded
            tspan=[t1,t2] - a time range to load into memory t2 must be greater than t1.
                if tspan=None, all times are loaded
            values = a list specifying which quantities to load entries keys are booleans and may be 'loadAll',
                'loadPressure', 'loadArea', 'loadFlow', 'loadWaveSpeed', and 'loadMeanVelocity'. If 'All'
                is in the list all quantities are loaded. Inputs are case insensitive.
            mindt := the minimum spacing in time between successively loaded points if
                none is specified, the solution time step is used.
        Effects and Usage:
            loads the specified values into memory such that they may be accessed as
            vascularNetwork.vessels[vesselId].Pressure, etc, returning a matrix of
            solution values corresponding to the time points in vascularNetwork.tsol.
            Accessing vessels and values not set to be loaded will produce errors.
        """
        if tspan is not None:
            t1 = tspan[0]
            t2 = tspan[1]
        else:
            t1 = self.simulationTime[0]
            t2 = self.simulationTime[-1]


        if self._checkAccessInputs(t1, t2, mindt):
            nSelectedBegin, nSelectedEnd = self.getFileAccessIndices(t1, t2)

            if mindt is not None:
                nTStepSpaces = int(np.ceil(mindt / self.dt))
            else:
                nTStepSpaces = 1

            tSlice = np.s_[nSelectedBegin:nSelectedEnd:nTStepSpaces]

            self.tsol = self.simulationTime[tSlice]
            # check if all vessels should be loaded
            if vesselIds == None: vesselIds = self.vessels.keys()
            # Update selected vessels
            for vesselId in vesselIds:
                if vesselId in self.vesselsToSave:
                    vessel = self.vessels[vesselId]
                    vessel.loadSolutionDataRange(tSlice,values)

        else:
            raise ValueError("classVascularNetwork::loadSolutionDataRangeVessel Error: Inputs were not valid you should not get here")



    def getFileAccessIndices(self,t1,t2):
        """
        Helper method to convert times to indices in the saved data.
        Input:
        t1,t2 the beginning and ending times to access
        Output:
        nSelectedBegin, nSelectedEnd - the indices corresponding to t1 and t2 in the file
        """
        startTime = self.simulationTime[0]
        nSelectedBegin = int(np.floor((t1 - startTime) / self.dt))
        nSelectedEnd = int(np.ceil((t2 - startTime) / self.dt))+1
        return nSelectedBegin, nSelectedEnd


    def findRootVessel(self):
        """
        Finds the root of a network, i.e. the vessel which is not a daughter of any vessel
        Evaluates a startRank for the evaulation of the network
        """
        daughters = []
        approximatedBif = 0
        for vessel in itervalues(self.vessels):
            try:
                if vessel.leftDaughter is not None:
                    daughters.append(vessel.leftDaughter)
                try:
                    if vessel.rightDaughter is not None:
                        daughters.append(vessel.rightDaughter)
                        approximatedBif += 1
                except Exception: self.warning("old except: pass clause #1 in VascularNetwork.findRootVessel", oldExceptPass= True)
            except Exception: self.warning("old except: pass clause #2 in VascularNetwork.findRootVessel", oldExceptPass= True)

        # find startRank by approximation of numbers of generations
        approxGen = len(set(daughters)) - 2 * approximatedBif + int(np.sqrt(approximatedBif))
        self.startRank = 2.0 ** (approxGen - 1)
        # find root with difference between daughters and all vessels as root is never daughter
        roots = list(set(self.vessels.keys()).difference(daughters))
        try:
            self.root = roots[0]
        except Exception:
            self.exception("vascularNetwork.findRootVessel(): could not find a root node")
        if len(roots) > 1:
            raise ValueError("vascularNetwork.searchRoot(): found several roots: {}, check network again!".format(roots))

    def checkDaughterDefinition(self):
        """
        Method to check if all daughters are defined in the correct way, i.e. if a vessel has only 1 daughter
        it should be defined as a leftDaughter, if it is defined as a rightDaughter, this method will rename it!
        additional check if there is a vessel with this id, if not remove daughter
        """
        for vessel in itervalues(self.vessels):
            if vessel.leftDaughter == None and vessel.rightDaughter is not None:
                self.warning("vascularNetwork.checkDaughterDefiniton(): Vessel {} has no leftDaughter but a rightDaughter {}, this daughter is now assumed to be leftDaughter".format(vessel.Id, vessel.rightDaughter), noException= True)
                vessel.leftDaughter = vessel.rightDaughter
                vessel.rightDaughter = None
            # check if daughter vessel exists
            if vessel.leftDaughter is not None:
                try:
                    self.vessels[vessel.leftDaughter]
                except Exception:
                    self.warning("vascularNetwork.checkDaugtherDefinition():\n      leftDaughter with Id {} of vessel {} does not exist".format(vessel.leftDaughter, vessel.Id,))
                    vessel.leftDaughter = None

            if vessel.rightDaughter is not None:
                try:
                    self.vessels[vessel.rightDaughter]
                except Exception:
                    self.warning("vascularNetwork.checkDaugtherDefinition():\n       rightDaughter with Id {} of vessel {} does not exist".format(vessel.rightDaughter, vessel.Id,))
                    vessel.rightDaughter = None

    def evaluateConnections(self):
        """
        Method to evaluate all connections:

        - check for right daughter definition (call)
        - find root of the network (call)
        - evalualte all connections link, bifurcation, anastomosis
        - apply mothers to all vessels (call)

        Method traverses tree with defined daughters,
        - finds mothers and connections

        -> creates treeTraverseList breadth first traversing list
        -> creates treeTraverseConnections list of connections [ [LeftMother, rightMother, leftDaughter, rightDaughter], ..]
        """

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

        if self.vessels[root].leftDaughter is not None:
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

            if leftDaughter is not None:
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
                if self.vessels[leftDaughter].leftDaughter is not None:
                    toVisit.append(leftDaughter)
                else:
                    # append vessel to ends as it has no left and right daughters
                    if leftDaughter not in self.boundaryVessels: self.boundaryVessels.append(leftDaughter)

                rightDaughter = self.vessels[motherVessel].rightDaughter

                if rightDaughter is not None:
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
                                self.warning("right daughter forced to anastomosis, not possible", noException=True)

                        elif len(existingMothers) == 2:
                            self.treeTraverseList.append(rightDaughter)
                            ranking[rightDaughter] = rankingRightDaughter
                            currentConnectionList = [existingMothers[0], existingMothers[1], rightDaughter, None]

                    # check if rightDaughter has also daughters which should be visualized
                    if self.vessels[rightDaughter].leftDaughter is not None: toVisit.append(rightDaughter)
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
        """
        Functions traverses the self.treeTraverseConnections and saves the id of the
        left and right mother of the vessel. Also check if there are any anastomosis 
        in the network 
        """
        self.anastomosisExists = False
        
        for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:
            
            if leftMother is not None and rightMother is not None:
                self.anastomosisExists = True
                
            self.vessels[leftDaughter].leftMother = leftMother
            self.vessels[leftDaughter].rightMother = rightMother
            try:
                self.vessels[rightDaughter].leftMother = leftMother
                self.vessels[rightDaughter].rightMother = rightMother
            except Exception: self.warning("old except: pass clause in VascularNetwork.applyMothersToVessel", oldExceptPass= True)

    def findStartAndEndNodes(self):
        """
        Function traverses self.treeTraverseConnections and creates start- and
        end-nodes for all vessels in the network
        """
        self.treeTraverseList_sorted = self.treeTraverseList[:]
        self.treeTraverseList_sorted.sort()
        
        nodes = []
        nodeCount = 0
        nodes.append(nodeCount)
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

        connectionNodes = [0]
        for vesselID in self.treeTraverseList_sorted:
            nodes.append(self.vessels[vesselID].endNode)
            endNodeTmp = self.vessels[vesselID].endNode
            connection = False
            for vesselIDtmp in self.treeTraverseList_sorted:
                startNodetmp = self.vessels[vesselIDtmp].startNode
                if startNodetmp == endNodeTmp:
                    connection = True
            
            if connection:
                connectionNodes.append(endNodeTmp)
            #print "vessel{0}: startNode={1}, endNode={2}".format(vesselID, self.vessels[vesselID].startNode, self.vessels[vesselID].endNode)
        
        nodes.sort()
        self.nodes =  list(set(nodes))
        connectionNodes.sort()
        self.connectionNodes = list(set(connectionNodes))


    def calculateNetworkResistance(self):
        """
        This function travers the network tree and calculates the
        cumultative system resistances Rcum for each vessel in the Network.
        """

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
                            print("vessel {} : estimated peripheral windkessel resistance (Rtotal) {}".format(vesselId, Rtotal / 133.32 * 1.e-6))
                    except Exception: self.warning("Old except:pass clause #1 in VascularNetwork.calculateNetworkResistance", oldExceptPass= True)
                    # # add resistance to the value
                    try: boundaryResistance = boundaryResistance + bc.Rtotal
                    except Exception:
                        # # winkessel 2 elements and single resistance
                        try:
                            if bc.Rc == 'VesselImpedance':
                                P = np.ones(self.vessels[vesselId].N) * self.vessels[vesselId].Ps  # 158.747121018*133.32 #97.4608013004*133.32#
                                compliance = self.vessels[vesselId].C(P)
                                area = self.vessels[vesselId].A(P)
                                waveSpeed = self.vessels[vesselId].c(area, compliance)
                                Z = 1.0 / (compliance * waveSpeed)[-1]
                                boundaryResistance = boundaryResistance + Z
                        except Exception: self.warning("Old except: pass clause #2 in VascularNetwork.calculateNetworkResistance", oldExceptPass= True)
                        try:
                            # # winkessel 2 elements and single resistance
                            boundaryResistance = boundaryResistance + bc.Rc
                        except Exception: self.warning("Old except: pass clause #3 in VascularNetwork.calculateNetworkResistance", oldExceptPass= True)

                # print('boundaryResistance',boundaryResistance/133.32*1.e-6)
                if boundaryResistance == 0:
                    print("\n Boundary Condition at end of vessel {} has no resistance".format(vesselId))
                    # # set boundaryresistance to 1/133.32*1.e6
                    print("The resistance is set to 1*133.32*1.e6 \n")
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
                self.warning("no method for resistance calculation for anastomosis is implemented!!! \n")


    def calculateInitialValuesFromSolution(self):
        
        """ This function calculates initialValues based on values saved from previous simulations
            The initialvalues are stored in subfolder InitialValues of the network. In this folder 
            there are files for each vesselId e.g. 1.txt with format:
            
                    Qin, Qout, Pin, Pout
                    value, value, value, value
        """
        
        
        #initialValuePath = mFPH.getDirectory('initialValueFileDirectory', self.name, self.dataNumber, 'write')
        initialValuePath = "/home/fredrik/starfish_working_directory/Full96Model/InitialValues"
        initialValues =  {}
        for vesselId in self.treeTraverseList:
            
            strvesselId = str(vesselId)
            initialValueFile = ''.join([initialValuePath,'/',strvesselId,'.txt'])
            try:
                f = open(initialValueFile, 'r')
            except:
                print("ERROR: could not open file '{0}' for reading initialValue. Line 1350 classVascularNetwork".format(initialValueFile))
            
            try:
                for n, line in enumerate(f):
                    if n == 1:
                        Qin, Qout, Pin, Pout = line.split(',')
                        Qin, Qout = float(Qin), float(Qout)
                        Pin, Pout = float(Pin), float(Pout)
            
            except:
                print("ERROR: format of file '{0}' is wrong! Not possible to set initialValues.  Line 1360 classVascularNetwork".format(initialValueFile))
            
            
            initialValues[vesselId] = {}
            initialValues[vesselId]['Pressure'] = [Pin, Pout]
            initialValues[vesselId]['Flow'] = [Qin, Qout]
        
        self.initialValues = initialValues
        
    def calculateInitialValuesLinearSystem(self, Qmean, it=10, includeDynamic=False):
        """
        This function convert the system to a lumped model of resistors in series and paralell and calculate average pressure and flow values
        to be used as initial conditions. The system is reduced to a set of linear equations. For every vessel there is an equation for
        the pressure drop over the vessel [Pstart - Pend + Q*Rv = 0], and for every junction there is an equation for the conservation of mass
        [Qin - Qout = 0]
        """
        
        Qmean = Qmean*10**6
        
        
        
        Nunknowns = len(self.connectionNodes) + len(self.treeTraverseList_sorted) - 1
        
        vesselRangeList = range(1, len(self.treeTraverseList_sorted) + 1) # list of same length as number of vessels, but starting from 1 and ending at N (number of vessels)
        
        M = np.zeros((Nunknowns, Nunknowns)) #system Matrix
        
        RHS = np.zeros(Nunknowns) #system right hand side

        itMax = it
        
        initialValues_prev_it = None
        
        for it in range(itMax):
            initialValues = {}
            
            M = np.zeros((Nunknowns, Nunknowns)) #system Matrix
        
            RHS = np.zeros(Nunknowns) #system right hand side
            
            for n, vesselId in enumerate(self.treeTraverseList_sorted):
                """ iter through the list of vesseldicts and add the pressure equation:
                    p_start - p_end - Q*Rv = 0, and add the indices and floats in M and RHS """
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
                                #print "vessel {} : estimated peripheral windkessel resistance (Rtotal) {}".format(vesselId, Rtotal / 133.32 * 1.e-6)
                        except Exception: self.warning("Old except:pass clause #1 in VascularNetwork.calculateNetworkResistance", oldExceptPass= True)
                        # # add resistance to the value
                        try: boundaryResistance = boundaryResistance + bc.Rtotal
                        except Exception:
                            # # winkessel 2 elements and single resistance
                            try:
                                if bc.Rc == 'VesselImpedance':
                                    P = np.ones(self.vessels[vesselId].N) * self.vessels[vesselId].Ps  # 158.747121018*133.32 #97.4608013004*133.32#
                                    compliance = self.vessels[vesselId].C(P)
                                    area = self.vessels[vesselId].A(P)
                                    waveSpeed = self.vessels[vesselId].c(area, compliance)
                                    Z = 1.0 / (compliance * waveSpeed)[-1]
                                    boundaryResistance = boundaryResistance + Z
                            except Exception: self.warning("Old except: pass clause #2 in VascularNetwork.calculateNetworkResistance", oldExceptPass= True)
                            try:
                                # # winkessel 2 elements and single resistance
                                boundaryResistance = boundaryResistance + bc.Rc
                            except Exception: self.warning("Old except: pass clause #3 in VascularNetwork.calculateNetworkResistance", oldExceptPass= True)
    
                    # print 'boundaryResistance',boundaryResistance/133.32*1.e-6
                    if boundaryResistance == 0:
                        print("\n Boundary Condition at end of vessel {} has no resistance".format(vesselId))
                        # # set boundaryresistance to 1/133.32*1.e6
                        print("The resistance is set to 1*133.32*1.e6 \n")
                        boundaryResistance = 1.*133.32 * 1.e6
                    
                    if it == 0:
                        nodeToNodeResistance = self.vessels[vesselId].calcVesselResistance() + boundaryResistance
                    else:
                        
                        P_prev_it = initialValues_prev_it[vesselId]['Pressure']
                        nodeToNodeResistance = self.vessels[vesselId].calcVesselResistance(P=P_prev_it) + boundaryResistance
                        
                    boundaryVessel = True
                        
                else:
                    
                    if it == 0:
                        nodeToNodeResistance = self.vessels[vesselId].calcVesselResistance()
                    else:
                        
                        P_prev_it = initialValues_prev_it[vesselId]['Pressure']
                        nodeToNodeResistance = self.vessels[vesselId].calcVesselResistance(P=P_prev_it)
                    boundaryVessel = False
                
    
                Pstart_index = None # index in Matrix for startnode
                Pend_index = None # index in Matrix for endnode
                
                startNode = self.vessels[vesselId].startNode
                endNode = self.vessels[vesselId].endNode
                
                for pos, item in enumerate(self.connectionNodes):
                    # iterate through connectionNodes (unknown P's) and find their correct index in Matrix M

            
                    if item == startNode:
                        Pstart_index = pos
                    
                    elif item == endNode:
                        Pend_index = pos
                        
                vesselindex = self.treeTraverseList_sorted.index(vesselId) # find index of vesselID in rangeList
                Q_index = vesselindex + len(self.connectionNodes) - 1 # -1 due to Q1 is known and python index start at 0

                # TODO: 1) both flow and pressure BC at inlet
                if vesselId != self.root and boundaryVessel == False:
                    M[n, Pstart_index] = 1
                    M[n, Pend_index] = -1
                    M[n, Q_index] = - 1.e-6*nodeToNodeResistance/133.32
                elif vesselId == self.root and len(self.treeTraverseList_sorted) > 1:
                    M[n, Pstart_index] = 1
                    M[n, Pend_index] = -1
                    RHS[n] = Qmean*1.e-6*nodeToNodeResistance/133.32

                elif vesselId == self.root:
                    M[n, Pstart_index] = 1
                    #M[n, Pend_index] = -1
                    RHS[n] = Qmean*1.e-6*nodeToNodeResistance/133.32
                 
                elif boundaryVessel:
                    
                    M[n, Pstart_index] = 1
                    M[n, Q_index] = -1.e-6*nodeToNodeResistance/133.32
                    RHS[n] = 0 # could set to venous pressure

            for leftMother, rightMother, leftDaughter, rightDaughter in (self.treeTraverseConnections):
        
                """iter through the list of junctions and add the equation about conservation of mass:
                   Qin - Qout = 0, and add the indices and floats in M and RHS """
                
                n += 1
                
                
                #find Q in and Qout of junction
                if rightMother == None:
                    Qin = [leftMother]
                else:
                    Qin = [leftMother, rightMother]
                    
                if rightDaughter == None:
                    Qout = [leftDaughter]
                else:
                    Qout = [leftDaughter, rightDaughter]
                
                
                for vesselId in Qin:
                    vesselindex = self.treeTraverseList_sorted.index(vesselId)
                    Q_index = vesselindex + len(self.connectionNodes) - 1 #int(vesselId)
                    
                    if vesselId == self.root:
                        RHS[n] = - Qmean
                    else:
                        M[n, Q_index] = 1
            
                for vesselId in Qout:
                    vesselindex = self.treeTraverseList_sorted.index(vesselId)
                    Q_index = vesselindex + len(self.connectionNodes) - 1 #int(vesselId)
                    
                    M[n, Q_index] = - 1

            meanPandQ = np.linalg.solve(M, RHS)
                    
            for n, vesselId in enumerate(self.treeTraverseList_sorted):
                # post processing to assign initialvalues
                
                Pstart_index = None 
                Pend_index = None
                
                startNode = self.vessels[vesselId].startNode
                endNode = self.vessels[vesselId].endNode
                
                for pos, item in enumerate(self.connectionNodes):
                    # iterate through connectionNodes (unknown P's) and find their correct index in solution array meanPandQ
            
                    if item == startNode:
                        Pstart_index = pos
                    
                    elif item == endNode:
                        Pend_index = pos
                
                if vesselId != self.root:
                    vesselindex = self.treeTraverseList_sorted.index(vesselId)
                    Q_index = vesselindex + len(self.connectionNodes) - 1 #int(vesselId)
                    qm = meanPandQ[Q_index]*1e-6
                else:
                    qm = Qmean*1e-6
                
                p0 = meanPandQ[Pstart_index]*133.32
                
                
                if Pend_index != None:
                    p1 = meanPandQ[Pend_index]*133.32
                else:
                    
                    if it == 0:
                        Rv = self.vessels[vesselId].calcVesselResistance()
                    else:
                        P_prev_it = initialValues_prev_it[vesselId]['Pressure']
                        Rv = self.vessels[vesselId].calcVesselResistance(P=P_prev_it)
                    p1 = p0 - qm*Rv
                
                initialValues[vesselId] = {}
                initialValues[vesselId]['Pressure'] = [p0, p1]
                initialValues[vesselId]['Flow'] = [qm, qm]
                initialValues[vesselId]['R'] = [p0/qm, p1/qm]
                
                
                if it == 0:
                    radiusProximal = self.vessels[vesselId].radiusProximal
                    radiusDistal = self.vessels[vesselId].radiusDistal
                    areaProximal = np.pi*radiusProximal**2
                    areaDistal = np.pi*radiusDistal**2
                else:
                    radiusProximal = np.sqrt(self.vessels[vesselId].A_nID([p0, p1], 0)/np.pi)
                    radiusDistal = np.sqrt(self.vessels[vesselId].A_nID([p0, p1], -1)/np.pi)
                    areaProximal = np.pi*radiusProximal**2
                    areaDistal = np.pi*radiusDistal**2
                rho = self.vessels[vesselId].rho
                
                initialValues[vesselId]['Velocity'] = [qm/areaProximal, qm/areaDistal]
                initialValues[vesselId]['DynamicPressure'] = [0.5*rho*(qm/areaProximal)**2, 0.5*rho*(qm/areaDistal)**2]
                initialValues[vesselId]['radius'] = [radiusProximal, radiusDistal]
                
                self.Rcum[vesselId] = (p0)/qm
    
            
            # # adjust pressure with venous pressure and difference between mean and diastolic pressure
            for initialArray in initialValues.itervalues():
                initialArray['Pressure'][0] = initialArray['Pressure'][0] + self.centralVenousPressure
                initialArray['Pressure'][1] = initialArray['Pressure'][1] + self.centralVenousPressure
                
                    
            # # adjust pressure for gravity pressure
            initialValuesWithGravity = self.initializeGravityHydrostaticPressure(initialValues, self.root)
            
            initialValues_prev_it = initialValuesWithGravity.copy()
        

 
        if includeDynamic:
            initialValues = self.adjustForTotalPressure(initialValues_prev_it)
        else:
            initialValues = initialValues_prev_it
        return initialValues
    
    
    def adjustForTotalPressure(self, initialValues):
        """ 1) Traverese from root to periphery and adjust for total pressure. 2) Calculate the flow out that match the pressure at terminal sites.
        3 ) Scale the flow out so that it matches the flow in to the system, and traverese from periphery towards the root. 4) Do this process itMax times."""
        idToprint = 1
#         print initialValues[self.root]['Flow']
        #print initialValues[idToprint]['Pressure'][0]
        initialValues_new = initialValues.copy()
        
        Qin = initialValues[self.root]['Flow'][0]
        
        itMax = 100
        
        for it in range(itMax):

            toVisit = []
            leftDaughter = self.vessels[self.root].leftDaughter
            rightDaughter = self.vessels[self.root].rightDaughter
            
            if leftDaughter != None:
                toVisit.append(leftDaughter)
            if rightDaughter != None:
                toVisit.append(rightDaughter)
            # traverese from root to periphery to ajust for total pressure
            while len(toVisit) > 0:
                #print "first loop"
                for vesselId in toVisit:
                    
                    leftMother = self.vessels[vesselId].leftMother
                    P_mother = initialValues_new[leftMother]['Pressure'][-1]
                    P_dynamic_mother = initialValues_new[leftMother]['DynamicPressure'][-1]
                    
                    P_dynamic = initialValues_new[vesselId]['DynamicPressure'][0]
                    
                    p0 = P_mother + P_dynamic_mother - P_dynamic
                    
                    P_prev_it = initialValues_new[vesselId]['Pressure']
                    Qm_prev_it = initialValues_new[vesselId]['Flow'][0]
                    Rv = self.vessels[vesselId].calcVesselResistance(P=P_prev_it)
                    #print(p0, Qm_prev_it, Rv, vesselId)
                    p1 = p0 - Qm_prev_it*Rv
                    
                    initialValues_new[vesselId]['Pressure'] = [p0, p1]
    
                    radiusProximal = np.sqrt(self.vessels[vesselId].A_nID([p0, p1], 0)/np.pi)
                    radiusDistal = np.sqrt(self.vessels[vesselId].A_nID([p0, p1], -1)/np.pi)
                    areaProximal = np.pi*radiusProximal**2
                    areaDistal = np.pi*radiusDistal**2
                    rho = self.vessels[vesselId].rho
                    
                    initialValues_new[vesselId]['Velocity'] = [Qm_prev_it/areaProximal, Qm_prev_it/areaDistal]
                    initialValues_new[vesselId]['DynamicPressure'] = [0.5*rho*(Qm_prev_it/areaProximal)**2, 0.5*rho*(Qm_prev_it/areaDistal)**2]
                    initialValues_new[vesselId]['radius'] = [radiusProximal, radiusDistal]
                    #initialValues_new[vesselId]['Flow'] = [None, None]
                    
                    leftDaughter = self.vessels[vesselId].leftDaughter
                    rightDaughter = self.vessels[vesselId].rightDaughter
                    
                    if leftDaughter != None and leftDaughter not in toVisit:
                        toVisit.append(leftDaughter)
                    if rightDaughter != None and rightDaughter not in toVisit:
                        toVisit.append(rightDaughter)
                        
                    toVisit.remove(vesselId)
            
            

            Pv = self.centralVenousPressure
            
            Qout = 0
            
            # calculate flow out 
            for vesselId in self.boundaryVessels:

                bc = self.boundaryConditions[vesselId][-1]
                
                R = bc.Rc + bc.Z
                
                p1 = initialValues_new[vesselId]['Pressure'][-1]
                
                Qout += (p1 - Pv)/R
                
            Qscale = Qin/Qout
            print("Qscale:", Qscale)
            #print "Qscale: ", Qscale
                
            # traverse backwards to match flow with pressuredifference
            toVisit = self.boundaryVessels[:]
            
            while len(toVisit) > 0:
                #print "second loop"
                for vesselId in toVisit:

                    if vesselId in self.boundaryVessels:
                        
                        bc = self.boundaryConditions[vesselId][-1]
                        
                        R = bc.Rc + bc.Z
                        P_prev_it = initialValues_new[vesselId]['Pressure']
                        
                        p1 = P_prev_it[-1]
                        
                        qm = Qscale*(p1 - Pv)/R
                        Rv = self.vessels[vesselId].calcVesselResistance(P=P_prev_it)
                        p1 = qm*R + Pv
                        p0 = p1 + qm*Rv
                        
                        initialValues_new[vesselId]['Flow'] = [qm, qm]
                        initialValues_new[vesselId]['Pressure'] = [p0, p1]
                        initialValues_new[vesselId]['R'] = [(p0 - Pv)/qm, (p1 - Pv)/qm]
                        leftMother, rightMother = self.vessels[vesselId].leftMother, self.vessels[vesselId].rightMother
                        
                        if leftMother != None and leftMother not in toVisit:
                            toVisit.append(leftMother)
                        if rightMother != None and rightMother not in toVisit:
                            toVisit.append(rightMother)
                        
                        toVisit.remove(vesselId)
                        
                    else:
                        
                        leftDaughter = self.vessels[vesselId].leftDaughter
                        rightDaughter = self.vessels[vesselId].rightDaughter
                        
                        #bifurcation
                        if rightDaughter != None:
                            q1, q2 = initialValues_new[leftDaughter]['Flow'][0], initialValues_new[rightDaughter]['Flow'][0]
                            
                            if q1 != None and q2 != None:
                                
                                qm = q1 + q2
    
                                
                                P_prev_it = initialValues_new[vesselId]['Pressure']
                                Rv = self.vessels[vesselId].calcVesselResistance(P=P_prev_it)
                                
                                P_daughter = initialValues_new[leftDaughter]['Pressure'][0]
                                P_dynamic_daughter = initialValues_new[leftDaughter]['DynamicPressure'][0]
                    
                                P_dynamic = initialValues_new[vesselId]['DynamicPressure'][-1]
                    
                                p1 = P_daughter + P_dynamic_daughter - P_dynamic
                                p0 = p1 + qm*Rv
                                
                                initialValues_new[vesselId]['Flow'] = [qm, qm]
                                initialValues_new[vesselId]['Pressure'] = [p0, p1]
                                
                                leftMother, rightMother = self.vessels[vesselId].leftMother, self.vessels[vesselId].rightMother
                                
                                if leftMother != None and leftMother not in toVisit:
                                    toVisit.append(leftMother)
                                if rightMother != None and rightMother not in toVisit:
                                    toVisit.append(rightMother)
                                
                                toVisit.remove(vesselId)
                        # anastomosis
                        elif self.vessels[leftDaughter].rightMother != None:
                                #print('hello')
                                leftMotherAnastomosis = self.vessels[leftDaughter].leftMother
                                rightMotherAnastomosis = self.vessels[leftDaughter].rightMother
                                leftDaughterAnastomosis = leftDaughter
                                
                                P_daughter = initialValues_new[leftDaughterAnastomosis]['Pressure'][0]
                                P_dynamic_daughter = initialValues_new[leftDaughterAnastomosis]['DynamicPressure'][0]
                                Ptot_daughter = P_daughter + P_dynamic_daughter
                                Q_D = initialValues_new[leftDaughterAnastomosis]['Flow'][0]
                                
                                P_LM = initialValues_new[leftMotherAnastomosis]['Pressure'][-1]
                                P_RM = initialValues_new[rightMotherAnastomosis]['Pressure'][-1]
                                
                                A_LM = np.pi*initialValues_new[leftMotherAnastomosis]['radius'][-1]**2
                                A_RM = np.pi*initialValues_new[rightMotherAnastomosis]['radius'][-1]**2
                                
                                if Ptot_daughter - P_RM < 0:
                                    Q_RM = initialValues_new[rightMotherAnastomosis]['Flow'][-1]
                                else:
                                    Q_RM = np.sqrt(2.*(Ptot_daughter - P_RM)/rho)*A_RM
                                if Ptot_daughter - P_LM < 0:
                                    Q_LM = initialValues_new[leftMotherAnastomosis]['Flow'][-1]
                                else:
                                    Q_LM = np.sqrt(2.*(Ptot_daughter - P_LM)/rho)*A_LM
                                    
                                Q_scaleAnastomosis = Q_D/(Q_RM + Q_LM)
                                
                                if vesselId == leftMotherAnastomosis:
                                    qm = Q_LM*Q_scaleAnastomosis
                                elif vesselId == rightMotherAnastomosis:
                                    qm = Q_RM*Q_scaleAnastomosis
                                else:
                                    print("Error in adjustForTotalPressure. Wrong handling of anastomosis")
                                    
                                P_prev_it = initialValues_new[vesselId]['Pressure']
                                Rv = self.vessels[vesselId].calcVesselResistance(P=P_prev_it)
                                
                    
                                P_dynamic = initialValues_new[vesselId]['DynamicPressure'][-1]
                    
                                p1 = P_daughter + P_dynamic_daughter - P_dynamic
                                p0 = p1 + qm*Rv
                                
                                initialValues_new[vesselId]['Flow'] = [qm, qm]
                                initialValues_new[vesselId]['Pressure'] = [p0, p1]
                                
                                leftMother, rightMother = self.vessels[vesselId].leftMother, self.vessels[vesselId].rightMother

                                if leftMother != None and leftMother not in toVisit:
                                    toVisit.append(leftMother)
                                if rightMother != None and rightMother not in toVisit:
                                    toVisit.append(rightMother)
                                
                                toVisit.remove(vesselId)
                        # link
                        else:
    
                            q1 = initialValues_new[leftDaughter]['Flow'][0]
                            
                            if q1 != None:
                                
                                qm = q1
    
                                P_prev_it = initialValues_new[vesselId]['Pressure']
                                Rv = self.vessels[vesselId].calcVesselResistance(P=P_prev_it)
                                
                                P_daughter = initialValues_new[leftDaughter]['Pressure'][0]
                                P_dynamic_daughter = initialValues_new[leftDaughter]['DynamicPressure'][0]
                    
                                P_dynamic = initialValues_new[vesselId]['DynamicPressure'][-1]
                    
                                p1 = P_daughter + P_dynamic_daughter - P_dynamic
                                p0 = p1 + qm*Rv
                                
                                initialValues_new[vesselId]['Flow'] = [qm, qm]
                                initialValues_new[vesselId]['Pressure'] = [p0, p1]
                                
                                leftMother, rightMother = self.vessels[vesselId].leftMother, self.vessels[vesselId].rightMother
                                
                                if leftMother != None and leftMother not in toVisit:
                                    toVisit.append(leftMother)
                                if rightMother != None and rightMother not in toVisit:
                                    toVisit.append(rightMother)
                                
                                toVisit.remove(vesselId)
            
        return initialValues_new
                
    
    def calcReflectionJunction(self, lumpedValues=None):
        
        for vesselId in self.treeTraverseList_sorted:
            
            leftDaughter = self.vessels[vesselId].leftDaughter
            rightDaughter = self.vessels[vesselId].rightDaughter

            
            if leftDaughter != None and rightDaughter != None:
                
                if lumpedValues==None:
                    Am = np.pi*self.vessels[vesselId].radiusDistal**2
                    Ald = np.pi*self.vessels[leftDaughter].radiusProximal**2
                    Ard = np.pi*self.vessels[rightDaughter].radiusProximal**2
                    
                    cm = self.vessels[vesselId].calcVesselWavespeed_in_out()[-1]
                    cld = self.vessels[leftDaughter].calcVesselWavespeed_in_out()[0]
                    crd = self.vessels[rightDaughter].calcVesselWavespeed_in_out()[0]
                    
                    Am_prox = np.pi*self.vessels[vesselId].radiusProximal**2
                    cm_prox = self.vessels[vesselId].calcVesselWavespeed_in_out()[0]
                else:
                    Am = np.pi*lumpedValues[vesselId]['radius'][-1]**2
                    Ald = np.pi*lumpedValues[leftDaughter]['radius'][0]**2
                    Ard = np.pi*lumpedValues[rightDaughter]['radius'][0]**2
                    
                    Pm = lumpedValues[vesselId]['Pressure']
                    Pld = lumpedValues[leftDaughter]['Pressure']
                    Prd = lumpedValues[rightDaughter]['Pressure']
                    
                    cm = self.vessels[vesselId].calcVesselWavespeed_in_out(P=Pm)[-1]
                    cld = self.vessels[leftDaughter].calcVesselWavespeed_in_out(P=Pld)[0]
                    crd = self.vessels[rightDaughter].calcVesselWavespeed_in_out(P=Prd)[0]
                    
                    Am_prox = np.pi*lumpedValues[vesselId]['radius'][0]**2
                    cm_prox = self.vessels[vesselId].calcVesselWavespeed_in_out(P=Pm)[0]
                    
                rho = self.vessels[vesselId].rho
                
                Ym = Am/(rho*cm)
                Yld = Ald/(rho*cld)
                Yrd = Ard/(rho*crd)
                
                R_bif = (Ym - Yld -Yrd)/(Ym + Yld + Yrd)
                R_bif_ld = (Yld - Ym)/(Ym + Yld)
                R_bif_rd = (Yrd - Ym)/(Ym + Yrd)
                
                Ym_prox = Am_prox/(rho*cm_prox)
                
                R_inOut = (Ym_prox - Ym)/(Ym + Ym_prox)
                if abs(R_bif) > 0.05:
                    print("aiai \n")
                print("vessel: {0}, R_inOut: {1}, R_bif: {2}, R_bif_ld: {3}, R_bif_rd: {4}".format(vesselId, round(R_inOut, 2), round(R_bif, 2), round(R_bif_ld, 2), round(R_bif_rd, 2)))
    
    def sumInertance(self, L1, L2, connection):
        
        if connection == "link":
            L = L1
        elif connection == "bif":
            L = 1./(1./L1 + 1./L2)
        elif connection == "anastomosis":
            L = 2*L1
        
        return L

    def sumResistance(self, R1, R2, connection):
        
        if connection == "link":
            R = R1
        elif connection == "bif":
            R = 1./(1./R1 + 1./R2)
        elif connection == "anastomosis":
            R = 2*R1
        
        return R
    
    
    def sumCompliance(self, C1, C2, connection):
        
        if connection == "link":
            C = C1
        elif connection == "bif":
            C = C1 + C2
        elif connection == "anastomosis":
            C = 0.5*C1
        
        return C

    def impedanceWeightedCompliance(self, Cv, C, Rv, R1, R2):
        
        C_weight = (Cv*R2 + Cv*R1 + C*R2 + Cv*Rv)/(Rv + R1 + R2)
        
        return C_weight
        
    
    def calcComplianceAndInertance(self, lumpedValues, state='standard', weightOnlyPeripheral=False):
        
        """ Method that traverses the tree from the periphery to the root and calculate total peripheral values of
            R, C, Cw and L"""
        
        toVisit = self.boundaryVessels[:]
        # TODO: Fix so that we can use with anastomosis and base on qm from self.lumpedValues
        #print toVisit
        
        
        for vesselId in self.treeTraverseList:
            if vesselId in self.boundaryVessels:
                vessel = self.vessels[vesselId]
                
                if state == "standard":
                    L = self.vessels[vesselId].Lv
                    Cv = self.vessels[vesselId].Cv
                    Rv = self.vessels[vesselId].resistance
                    
                    radiusProximal = self.vessels[vesselId].radiusProximal
                    radiusDistal = self.vessels[vesselId].radiusDistal

                    c_prox = self.vessels[vesselId].cd_in
                    c_dist = self.vessels[vesselId].cd_out
                else:
                    P_lump = lumpedValues[vesselId]["Pressure"]
                    
                    L = self.vessels[vesselId].calcVesselInertance(P=P_lump)
                    Cv = self.vessels[vesselId].calcVesselCompliance(P=P_lump)
                    Rv = self.vessels[vesselId].calcVesselResistance(P=P_lump)
                    
                    radiusProximal = np.sqrt(self.vessels[vesselId].A_nID(P_lump, 0)/np.pi)
                    radiusDistal = np.sqrt(self.vessels[vesselId].A_nID(P_lump, -1)/np.pi)
                    #[radiusProximal, radiusDistal] = lumpedValues[vesselId]["radius"]
                    
                    c_prox, c_dist = self.vessels[vesselId].calcVesselWavespeed_in_out(P=P_lump)
                    
                C = self.boundaryConditions[vesselId][-1].C
            
                try:
                    R = self.boundaryConditions[vesselId][-1].Z + self.boundaryConditions[vesselId][-1].Rc
                except:
                    R = self.boundaryConditions[vesselId][-1].Rtotal
                
                rho = self.vessels[vesselId].rho
                    
                    
                Ad_prox = np.pi*radiusProximal**2
                Ad_dist = np.pi*radiusDistal**2

                Z_prox, Z_dist = rho*c_prox/Ad_prox, rho*c_dist/Ad_dist
                try: 
                    R1 = self.boundaryConditions[vesselId][-1].Z
                except:
                    R1 = Z_dist
                    
                R2 = R - R1
                
                C_w = self.impedanceWeightedCompliance(Cv, C, Rv, R1, R2)
                
                lumpedValues[vesselId]["R1"] = [Z_prox, R1]
                lumpedValues[vesselId]["Cw"] = [C_w, C]
                lumpedValues[vesselId]["R_new"] = [R + Rv, R]
                lumpedValues[vesselId]["C"] = [C + Cv, C]
                lumpedValues[vesselId]["L"] = [L, 0]
                
                leftMother = vessel.leftMother
                rightMother = vessel.rightMother
                if rightMother != None:
                    toVisit.append(leftMother)
                    toVisit.append(rightMother)
                    
                if leftMother not in toVisit:
                    toVisit.append(leftMother)
                toVisit.remove(vesselId)
                
            else:
                lumpedValues[vesselId]["C"] = None
                lumpedValues[vesselId]["L"] = None
                
                
        while len(toVisit)>0:
            #print toVisit
            for vesselId in toVisit:
                vessel = self.vessels[vesselId]
                
                leftDaughter = vessel.leftDaughter
                rightDaughter = vessel.rightDaughter
                
                if leftDaughter !=None and rightDaughter != None:
                    connection = "bif"
                    
                    if lumpedValues[leftDaughter]["C"] != None and lumpedValues[rightDaughter]["C"] != None:
                        calcNext = True
                    else:
                        calcNext = False
                
                else:
                    
                    leftMotherOfDaughter = self.vessels[leftDaughter].leftMother
                    rightMotherOfDaughter = self.vessels[leftDaughter].rightMother
                    
                    if leftMotherOfDaughter !=None and rightMotherOfDaughter != None:
                        connection = "anastomosis"
                        
                        if lumpedValues[leftDaughter]["C"] != None:
                            calcNext = True
                        
                        else:
                            calcNext = False
                    
                    else:
                        connection = "link"
                    
                        if lumpedValues[leftDaughter]["C"] != None:
                            calcNext = True
                        
                        else:
                            calcNext = False
                
                if calcNext:
                        
                    L = self.vessels[vesselId].Lv
                    Cv = self.vessels[vesselId].Cv
                    if connection == "bif":
                        R1_new = lumpedValues[leftDaughter]["R_new"][0]
                        R2_new = lumpedValues[rightDaughter]["R_new"][0]
                        C1 = lumpedValues[leftDaughter]["C"][0]
                        C2 = lumpedValues[rightDaughter]["C"][0]
                        Cw1 = lumpedValues[leftDaughter]["Cw"][0]
                        Cw2 = lumpedValues[rightDaughter]["Cw"][0]
                        L1 = lumpedValues[leftDaughter]["L"][0]
                        L2 = lumpedValues[rightDaughter]["L"][0]
                    
                    elif connection == "link":
                        R1_new = lumpedValues[leftDaughter]["R_new"][0]
                        R2_new = None
                        C1 = lumpedValues[leftDaughter]["C"][0]
                        C2 = None
                        Cw1 = lumpedValues[leftDaughter]["Cw"][0]
                        Cw2 = None
                        L1 = lumpedValues[leftDaughter]["L"][0]
                        L2 = None

                    elif connection == "anastomosis":
                        R1_new = lumpedValues[leftDaughter]["R_new"][0]
                        R2_new = None
                        C1 = lumpedValues[leftDaughter]["C"][0]
                        C2 = None
                        Cw1 = lumpedValues[leftDaughter]["Cw"][0]
                        Cw2 = None
                        L1 = lumpedValues[leftDaughter]["L"][0]
                        L2 = None
                    
                    if state == "standard":
                        L = self.vessels[vesselId].Lv
                        Cv = self.vessels[vesselId].Cv
                        Rv = self.vessels[vesselId].resistance
                        
                        radiusProximal = self.vessels[vesselId].radiusProximal
                        radiusDistal = self.vessels[vesselId].radiusDistal
    
                        c_prox = self.vessels[vesselId].cd_in
                        c_dist = self.vessels[vesselId].cd_out
                    else:
                        P_lump = lumpedValues[vesselId]["Pressure"]
                        
                        L = self.vessels[vesselId].calcVesselInertance(P=P_lump)
                        Cv = self.vessels[vesselId].calcVesselCompliance(P=P_lump)
                        Rv = self.vessels[vesselId].calcVesselResistance(P=P_lump)

                        radiusProximal = np.sqrt(self.vessels[vesselId].A_nID(P_lump, 0)/np.pi)
                        radiusDistal = np.sqrt(self.vessels[vesselId].A_nID(P_lump, -1)/np.pi)
                        #[radiusProximal, radiusDistal] = lumpedValues[vesselId]["radius"]
                        
                        c_prox, c_dist = self.vessels[vesselId].calcVesselWavespeed_in_out(P=P_lump)
                    
                    R_newDistal = self.sumResistance(R1_new, R2_new, connection)    
                    Cdistal = self.sumCompliance(C1, C2, connection)
                    CwDistal = self.sumCompliance(Cw1, Cw2, connection)
                    Ldistal = self.sumInertance(L1, L2, connection)
                    lumpedValues[vesselId]["C"] = [Cdistal + Cv, Cdistal]
                    lumpedValues[vesselId]["L"] = [Ldistal + L, Ldistal]
                    
                    rho = self.vessels[vesselId].rho
                    Ad_prox = np.pi*radiusProximal**2
                    Ad_dist = np.pi*radiusDistal**2
                    

                    Z_prox, Z_dist = rho*c_prox/Ad_prox, rho*c_dist/Ad_dist
                    R1 = Z_dist
                    R2 = R_newDistal - R1
                    if weightOnlyPeripheral:
                        CwProximal = CwDistal + Cv
                    else:
                        CwProximal = self.impedanceWeightedCompliance(Cv, CwDistal, Rv, R1, R2)
                    
                    lumpedValues[vesselId]["R_new"] = [R_newDistal + Rv, R_newDistal]
                    lumpedValues[vesselId]["R1"] = [Z_prox, Z_dist]
                    lumpedValues[vesselId]["C"] = [Cdistal + Cv, Cdistal]
                    lumpedValues[vesselId]["Cw"] = [CwProximal, CwDistal]
                    lumpedValues[vesselId]["L"] = [Ldistal + L, Ldistal]
                    leftMother = vessel.leftMother
                    rightMother = vessel.rightMother
                    if rightMother != None:
                        if rightMother not in toVisit and vesselId != self.root:
                            toVisit.append(rightMother)
                    if leftMother not in toVisit and vesselId != self.root:
                        toVisit.append(leftMother)
                    toVisit.remove(vesselId)
        


    def calculateInitialValues(self):
        """
        This function travers the network tree and calculates the
        esitimates the initial flow and pressure values for each vessel in the Network
        based on the meanflow/pressure value at the root node using the cumultative resistance
        """

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
                
        if self.venousSystemCollapse == True and self.initialsationMethod != 'ConstantPressure':
            raise NotImplementedError("Auto, MeanFlow, Mean Pressure: initialization not implemented for collapsing venous system! \n")
            #exit()

        if self.initialsationMethod == 'Auto':
            try:
                meanInflow, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(quiet=self.quiet)
                self.initialisationPhaseExist = True
            except Exception:
                self.exception("classVascularNetwork: Unable to calculate mean flow at inflow point")
                #exit()

        elif self.initialsationMethod == 'MeanFlow':
            try:
                meanInflow = self.initMeanFlow
                # # addjust bc condition
                xxx, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(meanInflow, quiet=self.quiet)
                self.initialisationPhaseExist = True
            except Exception:
                self.exception("classVascularNetwork: Unable to set given meanFlow at inflow point")
                #exit()

        elif self.initialsationMethod == 'MeanPressure':
            try:
                meanInPressure = self.initMeanPressure
                self.initialisationPhaseExist = True
            except Exception:
                self.exception("classVascularNetwork: Unable to set given meanFlow at inflow point")
                #exit()

        elif self.initialsationMethod == 'AutoLinearSystem':

            try:
                print(self.dt)
                #inflowBoundaryCondition.findMeanFlow()
                meanInflow, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(quiet=self.quiet)
                #exit()
                self.initialisationPhaseExist = True

            except Exception:
                self.exception("classVascularNetwork: Unable to evaluate time shift to 0 at inflow point")
            
            self.findStartAndEndNodes() # allocate start and end nodes to all vessels in the network
            self.initialValues = self.calculateInitialValuesLinearSystem(meanInflow)
            
            self.lumpedValues = deepcopy(self.initialValues)
            
            return

        elif self.initialsationMethod == 'FromSolution':

            try:
                meanInflow, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlow()
                #meanInflow, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(0.0, quiet=self.quiet)
                self.initialisationPhaseExist = False
                if self.initPhaseTimeSpan > 0:
                    self.initialisationPhaseExist

            except Exception:
                self.exception("classVascularNetwork: Unable to evaluate time shift to 0 at inflow point")
            
            self.findStartAndEndNodes() # allocate start and end nodes to all vessels in the network
            
            self.lumpedValues = self.calculateInitialValuesLinearSystem(meanInflow)
            self.calcComplianceAndInertance(self.lumpedValues, state='MeanValues', weightOnlyPeripheral=True)
            import pickle
            pickle.dump(self.lumpedValues, open('/home/fredrik/Documents/git/mypapers/networkReduction/tex/src/lumpedValues96_weighted_correct.p', 'wb'))
            self.calculateInitialValuesFromSolution()

            
            return

                
        elif self.initialsationMethod == 'ConstantPressure':
            meanInflowlumped, self.initPhaseTimeSpanlumped = inflowBoundaryCondition.findMeanFlow()

            #try:
            self.findStartAndEndNodes() # allocate start and end nodes to all vessels in the network
                
            self.lumpedValues = self.calculateInitialValuesLinearSystem(meanInflowlumped)
            meanInflow = self.initMeanFlow

            #except Exception:
            #    self.exception("classVascularNetwork (ConstantPressure init): Unable to estimate lumpedValues, meanInflow: "+ str(meanInflowlumped))
            constantPressure = self.initMeanPressure
            try:
                constantPressure = self.initMeanPressure
                ## TODO: uncomment again after MC simulations are done
                #if inflowBoundaryCondition is not None:
                    #xxx, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(0.0, quiet=self.quiet)

                self.initialisationPhaseExist = False
                if self.initPhaseTimeSpan > 0:
                    self.initialisationPhaseExist = True

            except Exception:
                self.exception("classVascularNetwork: Unable to evaluate time shift to 0 at inflow point")
                #exit()
            #############################Inititalisation Method constant pressure #############
            initialValues[root] = {}
            initialValues[root]['Pressure'] = [constantPressure, constantPressure]
            initialValues[root]['Flow'] = [meanInflow, 0]

            # # set initial values of the vessels by traversing the connections
            for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:
                calcDaughters = [leftDaughter]
                if rightDaughter is not None: calcDaughters.append(rightDaughter)
                for daughter in calcDaughters:
                    initialValues[daughter] = {}
                    initialValues[daughter]['Pressure'] = [constantPressure, constantPressure]
                    initialValues[daughter]['Flow'] = [0, 0]

            # # adjust pressure with venous pressure
            if self.venousPool is not None: 
                for initialArray in itervalues(initialValues):
                    initialArray['Pressure'][0] = initialArray['Pressure'][0] + self.venousPool.P[0]
                    initialArray['Pressure'][1] = initialArray['Pressure'][1] + self.venousPool.P[0]

            # # Check if gravity is on and if user wants to correct for hydrostatic pressure
            if self.gravitationalField == True:

                input_ = 'K'
                while input_ not in ['y', 'Y', 'yes', 'Yes', 'n', 'no', 'No', 'NO']:
                    input_ =input3("\n Adjust for hydrostatic pressure(y/n): ")

                if input_ in ['y', 'Y', 'yes', 'Yes']:  # 'y' Adjust ConstantPressure to correct for hydrostatic pressure
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

        if self.venousSystemCollapse == True:
            self.warning("no method for venous collapsing system is implemented to initialize network with method 'Tree' !!!", noException=True)

        if meanInflow is not None:
            p0 = self.Rcum[root] * meanInflow
            p1 = p0 - self.vessels[root].resistance * meanInflow
            print("DB cVN 1519: root init pressure", p0,p1, self.Rcum[root], meanInflow)

        elif meanInPressure is not None:
            meanInflow = meanInPressure / self.Rcum[root]  # calculate mean flow
            p0 = meanInPressure
            # # addjust bc condition
            try:    xxx, self.initPhaseTimeSpan = self.boundaryConditions[root][0].findMeanFlowAndMeanTime(meanInflow, quiet=self.quiet)
            except Exception:
                self.warning("VascularNetwork: Unable to adjust calculated meanFlow at inflow point boundary condition !")

                self.initialisationPhaseExist = False
                self.initPhaseTimeSpan = 0.

            p1 = p0 - self.vessels[root].resistance * meanInflow
        else:
            raise ValueError("Neither flow or pressure value given at inflow point!")
            #exit()

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
                if rightDaughter is not None:
                    calcDaughters.append(rightDaughter)
                for daughter in calcDaughters:
                    qm = p0 / self.Rcum[daughter]
                    p1 = p0 - self.vessels[daughter].resistance * qm

                    initialValues[daughter] = {}
                    initialValues[daughter]['Pressure'] = [p0, p1]
                    initialValues[daughter]['Flow'] = [qm, qm]

            # # anastomosis
            elif rightDaughter == None:
                self.warning("VascularNetwork: no method for anastomosis is implemented to initialize network !!!", noException= True)


        # # adjust pressure with venous pressure
        if self.venousPool is not None:
            for initialArray in itervalues(initialValues):
                initialArray['Pressure'][0] = initialArray['Pressure'][0] + self.venousPool.P[0]
                initialArray['Pressure'][1] = initialArray['Pressure'][1] + self.venousPool.P[0]

        # # adjust pressure for gravity pressure
        initialValuesWithGravity = self.initializeGravityHydrostaticPressure(initialValues, root)

        self.initialValues = initialValuesWithGravity
    
    def calcLumpedValues(self, it=1, includeDynamic=False):

        inflowBoundaryCondition = None
        for bc in self.boundaryConditions[self.root]:
            if bc.type == 1:
                inflowBoundaryCondition = bc
        try:
            meanInflow, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlow()
            #meanInflow, self.initPhaseTimeSpan = inflowBoundaryCondition.findMeanFlowAndMeanTime(0.0, quiet=self.quiet)
            self.initialisationPhaseExist = False
            if self.initPhaseTimeSpan > 0:
                self.initialisationPhaseExist

        except Exception:
            self.exception("classVascularNetwork: Unable to evaluate time shift to 0 at inflow point")
        
        self.findStartAndEndNodes() # allocate start and end nodes to all vessels in the network
        
        lumpedValues = self.calculateInitialValuesLinearSystem(meanInflow, it=it, includeDynamic=includeDynamic)
        
        self.lumpedValues = deepcopy(lumpedValues)
        
        return lumpedValues

        
        

    def evaluateNetworkResistanceAndCompliance(self):

        arterialCompliance = 0
        arterialCompliance120 = 0
        arterialCompliance80  = 0
        arterialCompliancePmean = 0
        arterialVolume = 0
        arterialVolume120 = 0
        arterialVolume80  = 0
        arterialVolumePmean = 0

        for vesselId, vessel_i in iteritems(self.vessels):

            p0, p1 = self.initialValues[vesselId]['Pressure']
            initialPressure = np.linspace(p0, p1, int(vessel_i.N))
            C = vessel_i.C(initialPressure)
            Cvol = sum((C[1::] + C[0:-1]) / 2.0) * vessel_i.dz[0]  # ## works only if equidistant grid
            A = vessel_i.A(np.linspace(p0, p1, int(vessel_i.N)))
            Avol = sum((A[1::] + A[0:-1]) / 2.0) * vessel_i.dz[0]  # ## works only if equidistant grid
            arterialVolume = arterialVolume + Avol
            arterialCompliance = arterialCompliance + Cvol

            p0 = 120*133.32
            p1 = p0
            C = vessel_i.C(np.linspace(p0, p1, int(vessel_i.N)))
            Cvol120 = sum((C[1::] + C[0:-1]) / 2.0) * vessel_i.dz[0]
            A = vessel_i.A(np.linspace(p0, p1, int(vessel_i.N)))
            Avol120 = sum((A[1::] + A[0:-1]) / 2.0) * vessel_i.dz[0]  # ## works only if equidistant grid
            arterialVolume120 = arterialVolume120 + Avol120
            arterialCompliance120 = arterialCompliance120 + Cvol120

            p0 = 75*133.32
            p1 = p0
            C = vessel_i.C(np.linspace(p0, p1, int(vessel_i.N)))
            Cvol80 = sum((C[1::] + C[0:-1]) / 2.0) * vessel_i.dz[0]
            A = vessel_i.A(np.linspace(p0, p1, int(vessel_i.N)))
            Avol80 = sum((A[1::] + A[0:-1]) / 2.0) * vessel_i.dz[0]  # ## works only if equidistant grid
            arterialVolume80 = arterialVolume80 + Avol80
            arterialCompliance80 = arterialCompliance80 + Cvol80
            

            numberEstimates = 20
            complianceEstimates = np.empty(numberEstimates)
            for index,p in zip(range(numberEstimates),np.linspace(65.,110.,numberEstimates)):
                pressure = np.linspace(p*133.32, p*133.32, int(vessel_i.N))
                C = vessel_i.C(pressure)
                complianceEstimates[index] = sum((C[1::] + C[0:-1]) / 2.0) * vessel_i.dz[0]

            arterialCompliancePmean = arterialCompliancePmean + np.mean(complianceEstimates)

        windkesselCompliance = 0
        for bcs in itervalues(self.boundaryConditions):
            for bc in bcs:
                if "Windkessel" in bc.name:
                    windkesselCompliance = windkesselCompliance + bc.C

        logger.info("{:6} - arterial Volume initPressure".format(arterialVolume*1e6))
        logger.info("{:6} - arterial Volume 120".format(arterialVolume120*1e6))
        logger.info("{:6} - arterial Volume 80".format(arterialVolume80*1e6))
        logger.info("")
        logger.info("--------------------------")
        logger.info("{:6} - arterial compliance initPressure".format(arterialCompliance*133.32*1e6))
        logger.info("{:6} - arterial compliance 120".format(arterialCompliance120*133.32*1e6))
        logger.info("{:6} - arterial compliance 80".format(arterialCompliance80*133.32*1e6))
        logger.info("{:6} - arterial compliance physiological MPA range 65-120 mmHg".format(arterialCompliancePmean*133.32*1e6))
        logger.info("")
        logger.info("{:6} - windkessel compliance".format(windkesselCompliance*133.32*1e6))
        logger.info("--------------------------")
        totalArterialCompliance = (arterialCompliancePmean+windkesselCompliance)
        logger.info("{:6} - total arterial compliance".format(totalArterialCompliance*133.32*1e6))
        logger.info("{:6} - ratio between arterial/total compliance".format(arterialCompliancePmean/totalArterialCompliance))
        #self.calculateNetworkResistance()
        if self.initialsationMethod not in ['ConstantPressure', 'FromSolution']:
            
            rootVesselResistance = self.vessels[self.root].resistance
            logger.info("{:6} - total arterial resistance".format(self.Rcum[self.root]/133.32*1e-6))
            logger.info("{:6} - root vessel resistance".format(rootVesselResistance/133.32*1e-6))
            logger.info("{:6} - total-root vessel resistance".format((self.Rcum[self.root]-rootVesselResistance)/133.32*1e-6))
            logger.info("{:6} - ratio root vessel / total".format(rootVesselResistance/self.Rcum[self.root]))


    def evaluateWindkesselCompliance(self):

        self.TotalVolumeComplianceTree = 0.0
        self.totalTerminalAreaCompliance = 0.0
        for vesselId, vessel_i in iteritems(self.vessels):
            # vessel_i = self.vessels[vesselId]

            p0, p1 = self.initialValues[vesselId]['Pressure']
            initialPressure = np.linspace(p0, p1, int(vessel_i.N))
            C = vessel_i.C(initialPressure)
            if vesselId in self.boundaryVessels:
                self.totalTerminalAreaCompliance = self.totalTerminalAreaCompliance + C[-1]

            self.Cends[vesselId] = C[-1]

            Cvol = sum((C[1::] + C[0:-1]) / 2.0) * vessel_i.dz[0] # ## works only if equidistant grid

            # Cvol = C[-1]*vessel_i.length

            # print(sum(C[1:-1])*vessel_i.dz[0], Cvol2, C[0]*vessel_i.length)
            # print(C[-1]*vessel_i.length - Cvol2)
            # print(sum(C[1:-1])*vessel_i.dz[0]- C[-1]*vessel_i.length)
            # print('Cvol ',vesselId, ' ',Cvol,' ',C[-1])
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
            raise ValueError("VascularNetwork in calculating C_wkTotal!")

        if self.quiet == False:
            print('=====================================')
            print('__________total compliances________')
            print('               Compliance')
            print("TerminalArea     {:5.3}".format(self.totalTerminalAreaCompliance * 133.32 * 1.e6))
            print("TreeVolume       {:5.3}".format(self.TotalVolumeComplianceTree * 133.32 * 1.e6))
            print("Total System     {:5.3}".format((self.TotalVolumeComplianceTree + C_wkTotal) * 133.32 * 1.e6))
            print("Total WK's       {:5.3}".format(C_wkTotal * 133.32 * 1.e6))

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
                    raise ValueError("Windkessel Compliance at vessel {}:  {} < 0!".format(vesselId, wk3Compliance))
                    #exit()
        if self.quiet == False:
            print('________estimated compliances________')
            print(' vesselId       Rt       C     Cdef')
            for vesselId in self.vessels.keys():
                try: print("{:3} {:10.3f} {:10.3f} {:10.3f}".format(vesselId, wk3CompPrintList[vesselId][0], wk3CompPrintList[vesselId][1], wk3CompPrintList[vesselId][2]))
                except Exception: print("{:3}".format(vesselId))

    def calculateReflectionCoefficientConnection(self, mothers, daughters):
        """
        Function calculates reflection coefficient of a vessel connection

        Input:
            motherVessels   = [ [Id mother1, pressure mother1] ...  ]
            daughterVessels = [ [Id daughter1, pressure daughter1] ...  ]

        Return: reflectionCoefficient
        """

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
        """
        Calculates the optimal reflection coeffiecients for the network

        addapted from article Reymond et al.2009
        (very poor and instable method)
        """

        # # add rest of the vessels by traversing the connections
        for leftMother, rightMother, leftDaughter, rightDaughter  in self.treeTraverseConnections:
            #### to be changed

            maxReflectionCoeff = 0.005  # values in reymonds code
            toleranceReflectionCoeff = 0.0  # values in reymonds code
            reflectionCoefficient = 10.0  # start value to get while running

            radiusLeftDaughterInit = self.vessels[leftDaughter].radiusProximal

            # print("connection:",leftMother,rightMother, leftDaughter, rightDaughter)
            # while (abs(reflectionCoefficient)-maxReflectionCoeff) > toleranceReflectionCoeff:
            while abs(reflectionCoefficient) > maxReflectionCoeff or reflectionCoefficient < 0:
                # # setup initial pressure for left mother
                p0, p1 = self.initialValues[leftMother]['Pressure']
                initialPressureLM = np.linspace(p0, p1, int(self.vessels[leftMother].N))
                try:
                    # # setup initial pressure for right daughter used if anastomosis
                    p0, p1 = self.initialValues[rightMother]['Pressure']
                    initialPressureRM = np.linspace(p0, p1, int(self.vessels[rightMother].N))
                except Exception: self.warning("Old except: pass clause #1 in VascularNetwork.optimizeTreeRefelctionCoefficients", oldExceptPass= True)
                # # setup initial pressure for left daughter
                p0, p1 = self.initialValues[leftDaughter]['Pressure']
                initialPressureLD = np.linspace(p0, p1, int(self.vessels[leftDaughter].N))
                # # setup initial pressure for right daughter used if bifurcation
                try:
                    p0, p1 = self.initialValues[rightDaughter]['Pressure']
                    initialPressureRD = np.linspace(p0, p1, int(self.vessels[rightDaughter].N))
                except Exception: self.warning("Old except: pass clause #2 in VascularNetwork.optimizeTreeRefelctionCoefficients", oldExceptPass= True)
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
                            except Exception: self.warning("Old except: pass clause #3 in VascularNetwork.optimizeTreeRefelctionCoefficients", oldExceptPass= True)
                            self.vessels[vesselId].initialize({})
                        except Exception: self.warning("Old except: pass clause #4 in VascularNetwork.optimizeTreeRefelctionCoefficients", oldExceptPass= True)
                else:
                    for vesselId in [leftDaughter, rightDaughter]:
                        try:
                            self.vessels[vesselId].radiusProximal = self.vessels[vesselId].radiusProximal * 0.995
                            try: self.vessels[vesselId].radiusDistal = self.vessels[vesselId].radiusDistal * 0.995
                            except Exception: self.warning("Old except: pass clause #5 in VascularNetwork.optimizeTreeRefelctionCoefficients", oldExceptPass= True)
                            self.vessels[vesselId].initialize({})
                        except Exception: self.warning("Old except: pass clause #6 in VascularNetwork.optimizeTreeRefelctionCoefficients", oldExceptPass= True)
            print(" new Reflection Coeff area ratio", radiusLeftDaughterInit, self.vessels[leftDaughter].radiusProximal, 1 - (radiusLeftDaughterInit) / self.vessels[leftDaughter].radiusProximal)
            # print("      new Reflection coefficient {}, areas".format(reflectionCoefficient), self.vessels[leftDaughter].radiusProximal #, self.vessels[rightDaughter].radiusProximal)
            # print

    def showReflectionCoefficientsConnectionInitialValues(self):
        if self.quiet == False:
            print('=====================================')
            print('________Reflection Coefficients______')
            print(' LM RM LD RD   Reflection coefficient')
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
                print("{:3} {:3} {:3} {:3}      {:.4}".format(leftMother, rightMother, leftDaughter, rightDaughter, reflectionCoefficient))



    def showWaveSpeedOfNetwork(self, Pressure=None, Flow=None):
        print('=====================================')
        print('__________initial wave speed_________')
        print(' vessel    wave speed c(P_init)   A(P_init)    As_init      Dw(P_init)      Re(P_init)')
        for vesselId, vessel in iteritems(self.vessels):
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
                v = self.initialValues[vesselId]['Flow'][0] / A
            else:
                v = Flow / A

            Re = np.max(np.sqrt(A / np.pi) * 2.0 * v / self.globalFluid['my'] * self.globalFluid['rho'])
            print(' {:3}            {:5.4}            {:5.4}     {:5.4}     {:4.4}    {:5.0f}'.format(vesselId, c, np.max(A), As, Dw, Re))


    def initializeGravityHydrostaticPressure(self, initialValues, root):
        """
        Traverse the tree and initialize the nodes with the steady state hydrostatic pressure distribution
        """

        # root vessel
        p0, p1 = initialValues[root]['Pressure']
        p1 = p1 + self.vessels[root].netGravity[0] * self.vessels[root].length

        if p1 < 0. :
            raise ValueError("classVascularNetwork.initializeGravityHydrostaticPressure(), \n calculated negative pressure in initialization of vessel {} with inital values {}".format(root, [p0, p1]))
            #exit()

        initialValues[root]['Pressure'] = [p0, p1]

        # traverse tree to calculate the pressure influence of gravity
        for leftMother, rightMother, leftDaughter, rightDaughter in self.treeTraverseConnections:

            # link & anastomosis
            calcDaughters = [leftDaughter]

            # add for bifucation
            if rightDaughter is not None:
                calcDaughters.append(rightDaughter)

            for daughter in calcDaughters:
                # initial pressure gradiant due to viscos effects without gravity
                initialPressureDiff = initialValues[daughter]['Pressure'][1] - initialValues[daughter]['Pressure'][0]
                # update p0 with new p1 from mother including gravity
                p0 = initialValues[leftMother]['Pressure'][1]

                p1 = p0 + initialPressureDiff + self.vessels[daughter].netGravity[0] * self.vessels[daughter].length

                initialValues[daughter]['Pressure'] = [p0, p1]

                if p1 < 0. :
                    raise ValueError("classVascularNetwork.initializeGravityHydrostaticPressure(), \n calculated negative pressure in initialization of vessel {} with inital values {}".format(daughter, [p0, p1]))
        
        return initialValues

    def initializeVenousGravityPressure(self):
        """
        Calculate and initialze the venous pressure depending on gravity for the 2 and 3 element windkessel models
        """
        self.venousSystemCollapse = False
        if self.venousPool is not None:
            venousPoolPressure = self.venousPool.P[0]
        else:
            venousPoolPressure = 5*133.32

        # calculate absolute and relative venous pressure at boundary nodes
        for vesselId in self.boundaryVessels:

            relativeVenousPressure = 0.0 + self.globalFluid['rho'] * self.vessels[vesselId].positionEnd[0][2] * self.gravityConstant - self.vessels[vesselId].externalPressure

            if self.venousPool is not None:
                if self.venousPool.Pmin is not None:
                    if relativeVenousPressure < self.venousPool.Pmin:
                        relativeVenousPressure = self.venousPool.Pmin
                        self.venousSystemCollapse = True
                        self.warning("Venous system showing collapsing dynamics! Calculated Pressure %f".format(relativeVenousPressure), noException= True, quiet=True)

            for bc in self.boundaryConditions[vesselId]:
                # update venous pressure at boundary nodes
                if bc.name in ['_Windkessel-2Elements', 'Windkessel-2Elements', '_Windkessel-3Elements', 'Windkessel-3Elements']:
                    bc.update({'venousPressure':relativeVenousPressure})
        
        # # print(out of method)
        if self.quiet == False and self.venousPool is not None:
            print('\n=============================================================')
            print('_______________Venous Pressures _____________________________')
            print('%s %36.1f' % ('Central venous pressure:', round(self.venousPool.P[0], 2)))

    def print3D(self):

        # # print
        print('==========================================================================================  \n')
        print('__________________________Vessel Id: position, net Gravity________________________________')

        # traverse vascular network
        for vesselId in sorted(self.treeTraverseList):
            # positionStart = self.vessels[vesselId].positionStart
            # positionEnd = self.vessels[vesselId].positionEnd
            # print('Start position  : vessel  {} {:19.3f} {:20.3f} {:21.3f}'.format(vesselId, positionStart[0],   positionStart[1],   positionStart[2]))
            # print('End position    : vessel  {} {:19.3f} {:20.3f} {:21.3f}'.format(vesselId, positionEnd[0],     positionEnd[1],     positionEnd[2]))
            if self.gravitationalField == True:
                print('%s %2i %19.3f' % ('Net gravity     : vessel ', vesselId, self.vessels[vesselId].netGravity[0]))

    def calculate3DpositionsAndGravity(self, nTsteps=None, nSet=None):
        """
        Initializing the position and rotation of each vessel in 3D space
        Initializing netGravity of the vessels.

        Coordinate system is RHS with Z vertical so gravity acts in -Z.

        """
        # TODO: what is this?
        if nSet is not None:
            nTsteps = 0

        for n in range(nTsteps+1):

            if nSet is not None: n = nSet

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
                if rightDaughter is not None:
                    self.vessels[rightDaughter].caculatePositionAndGravity(n, positionEndMother, rotToGlobalSysMother)

                if rightMother is not None:
                    if np.sum(self.vessels[rightMother].positionEnd - self.vessels[leftMother].positionEnd) < 3.e-15:
                        raise NotImplementedError('ERROR: 3d positions of anastomosis {} {} {} is not correct!'.format(leftMother, rightMother, leftDaughter))

    def initializeVenousGravityPressureTime(self, nTsteps):
        """
        Calculate and initialze the venous pressure depending on gravity for the 2 and 3 element windkessel models
        """

        self.venousSystemCollapse = False
        if self.venousPool is not None:
            venousPoolPressure = self.venousPool.P[0]
        else:
            venousPoolPressure = self.centralVenousPressure


        # calculate absolute and relative venous pressure at boundary nodes
        for vesselId in self.boundaryVessels:
            relativeVenousPressure = np.empty(nTsteps+1)
            for n in range(nTsteps+1):

                relativeVP = venousPoolPressure + self.globalFluid['rho'] * self.vessels[vesselId].positionEnd[n][2] * self.gravityConstant - self.vessels[vesselId].externalPressure

                if self.venousPool is not None:
                    if self.venousPool.Pmin is not None:
                        if relativeVP < self.venousPool.Pmin:  
                            relativeVP = self.venousPool.Pmin
                            self.warning("Venous system showing collapsing dynamics!", noException= True, quiet=True)

                relativeVenousPressure[n] = relativeVP

            # update bc
            for bc in self.boundaryConditions[vesselId]:
                # update venous pressure at boundary nodes
                if bc.name in ['_Windkessel-2ElementsDAE', 'Windkessel-2ElementsDAE','_Windkessel-2Elements', 'Windkessel-2Elements', '_Windkessel-3Elements', 'Windkessel-3Elements']:
                    bc.update({'venousPressure':relativeVenousPressure})


# def calculateReflectionCoefficientBifurcations(vascularNetwork, solutionDataSet = None):
#     """
#     This method calculates the reflection coefficients of a vascularNetwork by
#     traversing it as binary tree.
#     Inputparameter:
#         vascularNetwork (instance::classVascularNetwork())
#         solutionDataSet = None (optional)
#         If a solutionDataSet is given, the coefficients are calculated transient taking the 
#         soultion data into account.
#         Else, it the coefficients are calculated with the initial condictions.
#     Retrun:  Rf_transient= {}  == { str(bifurcation) : [Rf_mother, Rf_leftDaughter, Rf_rightDaughter] }
#     """
#     
#     ## firstTraverse the tree and create a list of all Bifurcations
#     bifurcationList = []
#     viz = []
#     #find root
#     root = vascularNetwork.root[0]        
#     # add root to the viz vessels if root has daughters:
#     if vascularNetwork.vessels[root].leftDaugther is not None:
#         viz.append(root)
#     # loop through tree until all daughters are conected
#     while len(viz) != 0:
#         # get the mother vessel (already added) and add its daughters
#         motherVessel = viz.pop(0)
#         # find left daughter
#         leftDaughter  = vascularNetwork.vessels[motherVessel].leftDaugther
#         rightDaughter = vascularNetwork.vessels[motherVessel].rightDaugther
#         #append the mother to the calc list
#         curCalcList = [motherVessel]
#         
#         if leftDaughter is not None:
#             #append the left daughter to the calc list
#             curCalcList.append(leftDaughter)
#            
#             if rightDaughter is not None:
#                 #append the left daughter to the calc list
#                 curCalcList.append(rightDaughter)
#             else:
#                 curCalcList.append(None)
#             # check if leftDaughter has also daughters 
#             if vascularNetwork.vessels[leftDaughter].leftDaugther is not None:
#                 viz.append(leftDaughter)
#                 
#             if rightDaughter is not None:
#                 # check if rightDaughter has also daughters 
#                 if vascularNetwork.vessels[rightDaughter].leftDaugther is not None:
#                     viz.append(rightDaughter)
#         bifurcationList.append(curCalcList)
#           
#     # static // initial condition
#     if solutionDataSet == None:
#         Rf_static = {}
#         rho = vascularNetwork.globalFluid['rho'] # local rho!!
#         print('Rf at bifurcations for initial condition')
#         for bifurcation in bifurcationList:
#             mother = bifurcation[0]
#             leftDaughter = bifurcation[1]
#             rightDaughter = bifurcation[2]
#             # bifurcation    
#             Y_m = C*c# = (vascularNetwork.vessels[mother].r[-1]**2*np.pi) / (rho*vascularNetwork.vessels[mother].c[-1])
#             Y_lD # =  (vascularNetwork.vessels[leftDaughter].r[0]**2*np.pi) / (vascularNetwork.globalFluid[rho]*vascularNetwork.vessels[leftDaughter].c[0])
#             if rightDaughter is not None:
#                 Y_rD # = (vascularNetwork.vessels[rightDaughter].r[0]**2*np.pi) / (vascularNetwork.globalFluid[rho]*vascularNetwork.vessels[rightDaughter].c[0])          
#             else :
#                 Y_rD = 0.0
#             Rf_static[bifurcation] = (Y_m-Y_lD-Y_rD) / (Y_m+Y_lD+Y_rD)
#              
#         return Rf_static
#     else:
#     # over time // after simulation
#         print('Rf at bifurcations transient')
#         Rf_transient = {}
#         rho = vascularNetwork.globalFluid['rho'] # change to local
#         for bifurcation in bifurcationList:
#             mother = bifurcation[0]
#             leftDaughter = bifurcation[1]
#             rightDaughter = bifurcation[2]
#             
#             motherAsol = solutionDataSet['Area'][mother][:,[-1]]
#             motherCsol = solutionDataSet['WaveSpeed'][mother][:,[-1]]
#             leftDaughterAsol = solutionDataSet['Area'][leftDaughter][:,[0]]
#             leftDaughterCsol = solutionDataSet['WaveSpeed'][leftDaughter][:,[0]]
#             
#             Y_m  = (motherAsol) / (rho*motherCsol)
#             Y_lD  =  (leftDaughterAsol) / (rho*leftDaughterCsol)
#             
#             if rightDaughter is not None:
#                 rightDaughterAsol = solutionDataSet['Area'][rightDaughter][:,[0]]
#                 rightDaughterCsol = solutionDataSet['WaveSpeed'][rightDaughter][:,[0]]
#         
#                 Y_rD  =  (rightDaughterAsol) / (rho*rightDaughterCsol)
#             else :
#                 Y_rD = 0.0
#                 
#             Rf_transient[str(bifurcation)] = [(Y_m-Y_lD-Y_rD) / (Y_m+Y_lD+Y_rD),(Y_lD-Y_m-Y_rD) / (Y_m+Y_lD+Y_rD),(Y_rD-Y_lD-Y_m) / (Y_m+Y_lD+Y_rD)]
#             
#         return Rf_transient
