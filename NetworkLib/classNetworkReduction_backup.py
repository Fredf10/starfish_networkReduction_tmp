import sys
import os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur + '/../')

import UtilityLib.classStarfishBaseObject as cSBO

import classVessel as cVes
import classBaroreceptor as cBRX
import classVenousPool as classVenousPool

import UtilityLib.moduleFilePathHandler as mFPH

from VascularPolynomialChaosLib.classRandomInputManager import RandomInputManager
import numpy as np
import math
from scipy import interpolate
import pprint
import h5py
from classBoundaryConditions import *
from scipy.integrate import simps

from UtilityLib import classRuntimeMemoryManager

class NetworkReduction(cSBO.StarfishBaseObject):
    """
    Class representing a vascular Network
    The vascular network consists out of vessels defined in classVessel::Vessel()
    Additional Topology, BoundaryConditions and the SimulationContext are saved.
    """


    solutionMemoryFields    = ["simulationTime", "arterialVolume"]
    solutionMemoryFieldsToSave = ["simulationTime", "arterialVolume"]

    def __init__(self, vascularNetwork):
        
        
        
        self.vascularNetwork = vascularNetwork
        self.vascularNetwork.initialize(initializeForSimulation = True)
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
        
        self.boundaryVessels = vascularNetwork.boundaryVessels
        self.lumpedValues = vascularNetwork.lumpedValues
        
        self.useVesselsImpedance = False
        self.useAverageValues = False
        
        #print self.vessels
        #print self.boundaryConditions
        #self.test()
    
    
    def reduceNetwork(self, truncateFile):
        
        f = open(truncateFile)
        for n, line in enumerate(f):
            if n==1:
                line = line.replace('\n', '')
                strTruncateList = line.split(',')
        print strTruncateList
        #exit()
        if len(strTruncateList) == 1 and strTruncateList[0] == "None":
            print "not reducing anything"
        
        else:     
            for vesselId in strTruncateList:
            
                self.truncate(int(vesselId))
        
        if self.useVesselsImpedance:
            self.setToVesselImpedance()
                    
        


    def checkIfTerminalBif(self, vesselId):
        """ Check if vesselId has left daughter and right daughter who both have terminal Bc's"""
        
        # TODO: change to check if lD and rD are in self.boundaryvessels..
        leftDaughter, rightDaughter = self.vessels[vesselId].leftDaughter, self.vessels[vesselId].rightDaughter
        
        terminalBif = True
        
        if leftDaughter == None:
            terminalBif = False
        else:
            leftDaughterl =  self.vessels[leftDaughter].leftDaughter
            
        if rightDaughter == None:
            terminalBif = False
        else:
            leftDaughterr =  self.vessels[rightDaughter].leftDaughter
        
        if terminalBif:
            if leftDaughterl == None and leftDaughterr == None:
                terminalBif = True
            else:
                terminalBif = False
       
            
        return terminalBif
    
    
    def condenseTerminalBif(self, vesselId):
        """ This method reduces a motherVessel with two daughterVessels who both have a terminal wk3 bc. 
            First the two daughter vessels are reduced to two wk2, then from the values of R1, C1, R2 and C3 
            a single terminal Wk3 is assigned at the distal end of the motherVessel. The method is taken from Epstein S. et al"""
            
        leftDaughter, rightDaughter = self.vessels[vesselId].leftDaughter, self.vessels[vesselId].rightDaughter
        self.reduceTerminal(leftDaughter)
        self.reduceTerminal(rightDaughter)
        self.reduceTerminalBifurcation(vesselId)
    
    
    def truncate(self, vesselId):
        
        toVisit = []
        branchVessels = [] 
        toVisit.append(vesselId)
        branchVessels.append(vesselId)
        keepvisiting = True
        toVisit.append(vesselId)
        while len(toVisit)>0:
            tmpVessel = toVisit[0]
            leftDaughter, rightDaughter = self.vessels[tmpVessel].leftDaughter, self.vessels[tmpVessel].rightDaughter
            if leftDaughter not in branchVessels and leftDaughter !=None:
                branchVessels.append(leftDaughter)
            if rightDaughter not in(branchVessels)and rightDaughter !=None:
                branchVessels.append(rightDaughter)
            if leftDaughter != None:
                toVisit.append(leftDaughter)
            if rightDaughter != None:
                toVisit.append(rightDaughter)
            toVisit.remove(tmpVessel)
        
        print branchVessels
        
        terminalBifurcations = []
        
        for tmpVessel in branchVessels:
            
            terminalbif = self.checkIfTerminalBif(tmpVessel)
            if terminalbif:
                terminalBifurcations.append(tmpVessel)
        
        print terminalBifurcations
        if len(terminalBifurcations) == 1:
            keepTruncating = True
            tmpTerminalBif = terminalBifurcations[0]
            while keepTruncating:
                if tmpTerminalBif == vesselId:
                    self.condenseTerminalBif(tmpTerminalBif)
                    keepTruncating = False
                else:
                    self.condenseTerminalBif(tmpTerminalBif)
                    tmpTerminalBif = self.vessels[tmpTerminalBif].leftMother
        
        
    
    def reduceTerminal(self, vesselId, N=501):
        """ This method reduces a terminal WK3 vessel into a
            lumped Wk2 BC according to Epstein S. et al eq. 6 
            and 7. Rnew and Cnew are stored as vessels[vesselID].Rnew
            and vessels[vesselID].Cnew
        """
        
        
        #print self.boundaryConditions[vesselId]
        Z = self.boundaryConditions[vesselId][0].Z
        R = self.boundaryConditions[vesselId][0].Rc
        
        C = self.boundaryConditions[vesselId][0].C
        
        if self.useAverageValues:
            [p0, p1] = self.lumpedValues[vesselId]['Pressure']
            #print p0, p1
            areaProximal = self.vessels[vesselId].A_nID([p0, p1], 0)
            areaDistal = self.vessels[vesselId].A_nID([p0, p1], -1)
            radiusProximal, radiusDistal = np.sqrt(areaProximal/np.pi), np.sqrt(areaDistal/np.pi)
            #exit()
            
        else:
            radiusProximal = self.vessels[vesselId].radiusProximal
            radiusDistal = self.vessels[vesselId].radiusDistal
        
        h = self.vessels[vesselId].wallThickness
        E = self.vessels[vesselId].youngModulus

        Beta = (4./3)*np.sqrt(np.pi)*E*h
        
        rho = self.vessels[vesselId].rho
        gamma = self.vessels[vesselId].gamma
        my = self.vessels[vesselId].my
        
        r = np.linspace(radiusProximal, radiusDistal, N)
        
        Ad = np.pi*r**2
        
        cd = np.sqrt(Beta/(2*rho*Ad))*Ad**(1./4)
        
        x = np.linspace(0, self.vessels[vesselId].length, N)
        
        f1 = Ad/(cd**2)
        f3 = 1./(Ad**2)
        
        K1 = simps(f1, x)
        K3 = simps(f3, x)
        
        Cv = K1/(rho)
        Rv = 2*(gamma + 2)*np.pi*my*K3
        
        Rnew = Z + R + Rv
        Cnew = (Cv*R + Cv*Z + C*R + Rv*Cv)/Rnew

        
        self.vessels[vesselId].Rnew = Rnew
        self.vessels[vesselId].Cnew = Cnew
        
        
    def reduceTerminalBifurcation(self, vesselId, N=501):
        """ This function reduces a terminal bifurcation into a single 
            vessel with a wk3 BC according to Epstein S et al eq. 8 and 9. 
            Rnew1, Rnew2, Cnew1 and Cnew2 have already been calculated 
            by reducing left and right daughter terminal wk3 vessels.
            A new instance self.boundaaryConditions[vesselID] is created by copying
            the left daughter BC. left and right daughter bc is deleted. left and right daughter vessels
            are deleted, and mother vessels left and right daughter is set to None
        """
        
        leftDaughter, rightDaughter = self.vessels[vesselId].leftDaughter, self.vessels[vesselId].rightDaughter
            
        radiusProximal = self.vessels[vesselId].radiusProximal
        radiusDistal = self.vessels[vesselId].radiusDistal
            
        h = self.vessels[vesselId].wallThickness
        E = self.vessels[vesselId].youngModulus
        Beta = (4./3)*np.sqrt(np.pi)*E*h
        
        rho = self.vessels[vesselId].rho
        
        Ad_distal = np.pi*radiusDistal**2
        c_distal =  np.sqrt(Beta/(2*rho*Ad_distal))*Ad_distal**(1./4)
        
        Rnew1, Cnew1 = self.vessels[leftDaughter].Rnew, self.vessels[leftDaughter].Cnew
        Rnew2, Cnew2 = self.vessels[rightDaughter].Rnew, self.vessels[rightDaughter].Cnew
        
        Rnew = 1./(1./Rnew1 + 1./Rnew2)
        
        Cnew = Cnew1 + Cnew2
        
        Z_0 = rho*c_distal/Ad_distal
        
        R = Rnew - Z_0
#         print "     R2 = {0} ".format(R)
        
        # todo fix so one can reduce all the way down to self.root
        if vesselId in self.boundaryVessels:
            self.boundaryConditions[vesselId].append(self.boundaryConditions[leftDaughter][0])
            self.boundaryConditions[vesselId][1].Z = Z_0
            self.boundaryConditions[vesselId][1].Rc = R
            self.boundaryConditions[vesselId][1].Rtotal = Z_0 + R
            self.boundaryConditions[vesselId][1].C = Cnew
        else:
            self.boundaryConditions[vesselId] = self.boundaryConditions[leftDaughter][:]
            self.boundaryConditions[vesselId][0].Z = Z_0
            self.boundaryConditions[vesselId][0].Rc = R
            self.boundaryConditions[vesselId][0].Rtotal = Z_0 + R
            self.boundaryConditions[vesselId][0].C = Cnew
        self.deleteWK3(leftDaughter)
        self.deleteWK3(rightDaughter)
        self.deleteVessel(leftDaughter)
        self.deleteVessel(rightDaughter)
        
        self.vessels[vesselId].leftDaughter = None
        self.vessels[vesselId].rightDaughter = None
        
        self.treeTraverseList_sorted.remove(leftDaughter)
        self.treeTraverseList_sorted.remove(rightDaughter)
        self.treeTraverseList.remove(leftDaughter)
        self.treeTraverseList.remove(rightDaughter)
        
        self.boundaryVessels.remove(leftDaughter)
        self.boundaryVessels.remove(rightDaughter)
        
        self.boundaryVessels.append(vesselId)

    def setToVesselImpedance(self):
        
        for vesselId in self.treeTraverseList_sorted:
            if vesselId in self.boundaryVessels:
                bc = self.boundaryConditions[vesselId] 
                if len(bc)>1:
                    bc = bc[1]
                else:
                    bc = bc[0]
            
                bc.Z = 'VesselImpedance'
                bc.Rtotal = None
        
    def deleteVessel(self, inputId):
        """
        Remove vessel from network and delete it
        """
        try:
            del self.vessels[inputId]
        except Exception:
            self.warning("networkReduction.deleteVessel(): vessel with Id {} does not exist! Could not remove vessel".format(inputId))
            print "ERROR networkReduction.deleteVessel(): vessel with Id {} does not exist! Could not remove vessel".format(inputId)


    def deleteWK3(self, inputId):
        """
        Remove vessel from network and delete it
        """
        try:
            del self.boundaryConditions[inputId]
        except Exception:
            self.warning("networkReduction.deleteWK3(): boundaryvessel with Id {} does not exist! Could not remove vessel".format(inputId))
            print "ERROR networkReduction.deleteWK3(): boundaryvessel with Id {} does not exist! Could not remove vessel".format(inputId)

    def getVariableValue(self, variableName):
        """
        Returns value of variable with name : variableName
        States Error if not such variable
        """
        try:
            return self.__getattribute__(variableName)
        except Exception:
            self.warning("vascularNetwork.getVariable() : VascularNetwork has no variable {}".format(variableName))
    
    

    
        
            
            
            