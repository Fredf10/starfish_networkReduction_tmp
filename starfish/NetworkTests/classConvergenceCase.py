'''
Created on Jun 28, 2017

@author: fredrik
'''
from scipy import interpolate
import h5py, pickle
import numpy as np
import os
import matplotlib.pylab as plt

from starfish.VascularPolynomialChaosLib.testBaseClass import TestBaseClass
import starfish.UtilityLib.moduleXML as mXML

import starfish.NetworkTests.classCreateNetworks as cCreateNetwork
import starfish.UtilityLib.moduleFilePathHandler as mFPH
import starfish.VascularPolynomialChaosLib.moduleBatchSimulationManager as mBatchSim
import starfish.UtilityLib.moduleLogFile as mLog

class ConvergenceCase(TestBaseClass):
    
    externVariables      = { 'networkName'              : TestBaseClass.ExtValue(str, strCases = ['anything']),
                             'createNetworks'           : TestBaseClass.ExtValue(bool),
                             'solveNetworks'            : TestBaseClass.ExtValue(bool),
                             'processNetworks'          : TestBaseClass.ExtValue(bool),          
                             'multiprocessing'          : TestBaseClass.ExtValue(bool),  
                             'numberOfProcessors'       : TestBaseClass.ExtValue(int),
                             'numberOfRefinements'      : TestBaseClass.ExtValue(int),
                             'vesselOfInterest'         : TestBaseClass.ExtValue(int),
                             'dataType'                 : TestBaseClass.ExtValue(str, strCases = ['P', 'Q'])
                             
                           } 
    
    externXmlAttributes  = []
    
    externXmlElements    = ['networkName',
                            'createNetworks',
                            'solveNetworks',
                            'processNetworks',
                            'createNetworks',
                            'multiprocessing',
                            'numberOfProcessors',
                            'numberOfRefinements',
                            'vesselOfInterest',
                            'dataType'
                            
                            ]
    
    def __init__(self):
        
        
        self.networkName = None
        self.createNetworks = True
        self.solveNetworks = True
        self.processNetworks = True
        self.multiprocessing = False
        self.numberOfProcessors = 1
        self.numberOfRefinements = 4
        self.vesselOfInterest = 1
        self.dataType = 'P'
    
    def createEvaluationCaseFiles(self, convergencePickleCaseFile=None):
        
        if self.createNetworks:
            vascularNetwork = mXML.loadNetworkFromXML(self.networkName, dataNumber = "xxx")
            vascularNetwork.quiet = True

            
            New_network = cCreateNetwork.CreateNetwork(vascularNetwork, quiet=True)
            self.period = New_network.period

            CFL = 0.8
            dt = New_network.findMindt(CFL=CFL)
            New_network.assignNodes(dt, CFL=CFL)
            dataNumber = 001
            batchDataList = []
            for it in range(self.numberOfRefinements):
                batchData = {}
                
                New_network = cCreateNetwork.CreateNetwork(vascularNetwork, quiet=True)
                
                New_network.assignNodes(dt, CFL=CFL)
                networkXMLFile = mFPH.getFilePath('networkXmlFileConvergence', self.networkName, str(dataNumber), 'write')
    
                mXML.writeNetworkToXML(New_network, dataNumber=str(dataNumber), networkXmlFile=networkXMLFile)
    
                solutionFilehd5 = mFPH.getFilePath('solutionFileConvergence', self.networkName, str(dataNumber), 'write')
                solutionFileDirectory = mFPH.getDirectory('solutionFileConvergenceDirectory', self.networkName, str(dataNumber), 'read')
                batchData['pathSolutionDataDirectory'] = solutionFileDirectory
                batchData['networkName'] = self.networkName
                batchData['networkXmlFileLoad'] = networkXMLFile
                batchData['networkXmlFileSave'] = networkXMLFile
                batchData['pathSolutionDataFilename'] = solutionFilehd5
                batchData['dataNumber'] = str(dataNumber)
                batchData['period'] = self.period
                
                batchData['dt'] = dt
                batchData['CFL'] = CFL 
                
                batchDataList.append(batchData)
                
                dataNumber += 1
                dt *= 0.5
                del New_network
            
            pickle.dump(batchDataList, open(convergencePickleCaseFile, 'wb'))
            
        else:
            batchDataList = pickle.load(open(convergencePickleCaseFile, 'r'))
            self.period = batchDataList[0]['period']
        if self.solveNetworks:
            mBatchSim.runBatchAsMultiprocessing(batchDataList)
            pickle.dump(batchDataList, open(convergencePickleCaseFile, 'wb'))
        if self.processNetworks:
            self.solutionFileDirectory = mFPH.getDirectory('solutionFileConvergenceDirectory', self.networkName, '100', 'read')
            epsilonList = self.compareSolutions(batchDataList)
            
            cClog = mLog.ConvergenceLogFile(convergencePickleCaseFile.replace('.p', '.tex'), self.networkName, epsilonList, batchDataList, self.solutionFileDirectory)
            cClog.writeConvergenceLogfile(compileLogFile=True, deleteAuxiliary=True)
            
    def compareSolutions(self, batchDataList, roundTo=3):
        
        solutionFileRef = batchDataList[-1]['pathSolutionDataFilename']
        dataNumberRef = batchDataList[-1]['dataNumber']
        errorList = []
        
        if os.path.isdir(self.solutionFileDirectory + "/fig/"):
            pass
        else:
            os.mkdir(self.solutionFileDirectory + "/fig/")
        for batchData in batchDataList[:-1]:
            
            solutionFile = batchData['pathSolutionDataFilename']
            dataNumber = batchData['dataNumber']
            epsilon = self.compareWithOtherSolution(solutionFile, solutionFileRef, dataNumber, dataNumberRef, batchData)
    
            errorList.append(round(epsilon*100, roundTo))
        
        print errorList
        
        return errorList
        
    def compareWithOtherSolution(self, solutionFile, solutionFileRef, dataNumber, dataNumberRef, batchData):
        
        hdf5File = h5py.File(solutionFile, 'r')
        hdf5File2 = h5py.File(solutionFileRef, 'r')
        
        time = hdf5File['VascularNetwork']['simulationTime'][:]
        time2 = hdf5File2['VascularNetwork']['simulationTime'][:]
        dt = time[1] - time[0]
        #freq = 1./0.8
        period = self.period
        N = int(round(period/dt))
        
        nCycles = int((time[-1]/period))
        nCycles2 = int((time2[-1]/period))
        
        nCycles = min([nCycles, nCycles2])
        
        t_start = period*(nCycles - 1)
        t_end = period*nCycles
        
        epsilonMax = 0
        time_compare = np.linspace(t_start, t_end, N + 1)

        for vesselName in hdf5File['vessels'].keys():

            vesselId = vesselName.split(' - ')[-1]
            vesselId = int(vesselId)

            startP_all  = hdf5File['vessels'][vesselName][self.dataType + 'sol'][:, 0]
            startP2_all  = hdf5File2['vessels'][vesselName][self.dataType + 'sol'][:, 0]

            tck = interpolate.splrep(time, startP_all)
            tck2 = interpolate.splrep(time2, startP2_all)
        
            
            startP = interpolate.splev(time_compare, tck)
            startP2 = interpolate.splev(time_compare, tck2)

            
            epsilon = self.calcEpsilonAvg(startP2, startP, data_type=self.dataType)
            if vesselId == self.vesselOfInterest:
                epsilonMax = epsilon
                
                plt.figure()
                plt.plot(time_compare, startP, 'r')
                plt.plot(time_compare, startP2, 'k--')
                plt.xlabel('t')
                plt.ylabel(self.dataType)
                plt.legend(['dataN' + dataNumber, 'dataN' + dataNumberRef])
                plt.title("vesselID: {0}".format(vesselId))
                figFile = self.solutionFileDirectory + "/fig/" + dataNumber + "_" + dataNumberRef + ".pdf"
                plt.savefig(figFile)
                batchData['figFile'] = figFile
        
        return epsilonMax
    

    def calcEpsilonAvg(self, refData, numData, data_type="P"):
        
        if data_type == "P":
            RMS = np.sum(np.abs((numData - refData)/refData))/len(refData)
        elif data_type == "Q":
            RMS = np.sum(np.abs((numData - refData)/np.amax(refData)))/len(refData)
            
        return RMS