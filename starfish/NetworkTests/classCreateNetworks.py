'''
Created on Jun 28, 2017

@author: fredrik
'''
import starfish.UtilityLib.classStarfishBaseObject as cSBO

from copy import deepcopy


class CreateNetwork(cSBO.StarfishBaseObject):
    """
    Create networks to be used in convergencetest
    """


    def __init__(self, vascularNetwork, quiet=False):
        
        
        
        self.vascularNetwork = vascularNetwork
        self.vascularNetwork.quiet = quiet
        self.vascularNetwork.initialize(initializeForSimulation = True)
        self.lumpedValues = self.vascularNetwork.calcLumpedValues()
        #self.vascularNetwork.initializeNetworkForSimulation()
        self.name = vascularNetwork.name
        self.randomInputManager = vascularNetwork.randomInputManager
        self.globalFluid = vascularNetwork.globalFluid
        self.externalStimuli = vascularNetwork.externalStimuli
        self.baroreceptors = vascularNetwork.baroreceptors
        self.communicators = vascularNetwork.communicators
        self.venousPool = vascularNetwork.venousPool
        self.treeTraverseList = vascularNetwork.treeTraverseList
        self.treeTraverseList_sorted = vascularNetwork.treeTraverseList[:]
        self.treeTraverseList_sorted.sort()
        self.vessels = vascularNetwork.vessels
        
        self.totalTime = vascularNetwork.totalTime
        self.dt = vascularNetwork.dt
        self.description = vascularNetwork.description
        self.CFL = vascularNetwork.CFL
        self.timeSaveBegin = vascularNetwork.timeSaveBegin
        self.minSaveDt = vascularNetwork.minSaveDt
        self.maxMemory = self.vascularNetwork.maxMemory
        self.gravitationalField = self.vascularNetwork.gravitationalField
        self.gravityConstant = self.vascularNetwork.gravityConstant
        self.solvingSchemeField = self.vascularNetwork.solvingSchemeField
        self.rigidAreas = self.vascularNetwork.rigidAreas
        self.simplifyEigenvalues = self.vascularNetwork.simplifyEigenvalues
        self.riemannInvariantUnitBase = self.vascularNetwork.riemannInvariantUnitBase
        self.automaticGridAdaptation = self.vascularNetwork.automaticGridAdaptation
        self.initialsationMethod = self.vascularNetwork.initialsationMethod
        self.initMeanFlow = self.vascularNetwork.initMeanFlow
        self.initMeanPressure = self.vascularNetwork.initMeanPressure
        self.estimateWindkesselCompliance = self.vascularNetwork.estimateWindkesselCompliance
        self.compPercentageWK3 = self.vascularNetwork.compPercentageWK3
        self.compPercentageTree = self.vascularNetwork.compPercentageTree
        self.compTotalSys = self.vascularNetwork.compTotalSys
        self.boundaryConditions = vascularNetwork.boundaryConditions
        self.getVariableValue = vascularNetwork.getVariableValue
        
        self.boundaryVessels = vascularNetwork.boundaryVessels
        self.boundaryVesselsBaseline = vascularNetwork.boundaryVessels
        self.lumpedValues = vascularNetwork.lumpedValues
        
        root = vascularNetwork.root
        for bcTmp in self.boundaryConditions[root]:
            if bcTmp.type == 1:
                
                freq = bcTmp.freq
        self.period = 1./freq
    
    def findMindt(self, Nmin=5, CFL=0.9):
        
        
        minDt = 1000.
        
        for vesselId in self.treeTraverseList_sorted:
            l = self.vessels[vesselId].length
            [U_in, U_out]= self.lumpedValues[vesselId]['Velocity']
            Pm = self.lumpedValues[vesselId]['Pressure'] #['Pressure']
            [c_in, c_out] = self.vessels[vesselId].calcVesselWavespeed_in_out(P=Pm)
            dx = l/Nmin
            
            U_pluss_c = max([abs(U_in) + c_in, abs(U_out) + c_out])
            
            dt = CFL*dx/U_pluss_c
            
            if dt < minDt:
                minDt = dt
        
        return minDt
    
    
    def assignNodes(self, dt, CFL=0.9):
        
        self.vascularNetwork.dt = dt
        for vesselId in self.treeTraverseList_sorted:
            l = self.vessels[vesselId].length
            [U_in, U_out]= self.lumpedValues[vesselId]['Velocity']
            Pm = self.lumpedValues[vesselId]['Pressure'] #['Pressure']
            [c_in, c_out] = self.vessels[vesselId].calcVesselWavespeed_in_out(P=Pm)
            
            U_pluss_c = max([abs(U_in) + c_in, abs(U_out) + c_out])
            
            N = int(round(CFL*l/(dt*U_pluss_c)))
            CFL_actual = dt*U_pluss_c*N/l
            if CFL_actual > CFL + 0.1 or CFL_actual > 1.:
                print "Warning CFL is to big"
            
            self.vessels[vesselId].N = N
        
        
            
            
        