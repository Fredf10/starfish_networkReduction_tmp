import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))
sys.path.append(cur+'/../')

from VascularPolynomialChaosLib.testBaseClass import TestBaseClass 
import VascularPolynomialChaosLib.moduleBatchSimulationManager as mBatchSim


import moduleFilePathHandlerVNR as mFPH_VNR
import VascularNetworkReductionLib.classNetworkReduction as cNred

import UtilityLib.moduleXML as mXML
import UtilityLib.progressBar as cPB
import pickle


import time
import numpy as np
import multiprocessing

class ReductionCase(TestBaseClass):
    
    externVariables      = { 'caseName'                 : TestBaseClass.ExtValue(str, strCases = ['anything']),
                             'dataNumber'               : TestBaseClass.ExtValue(str, strCases = ['anything']),
                             'simulationDescription'    : TestBaseClass.ExtValue(str, strCases = ['anything']), 
                             'createNetworks'           : TestBaseClass.ExtValue(bool),
                             'solveNetworks'            : TestBaseClass.ExtValue(bool),
                             'batchDataFile'            : TestBaseClass.ExtValue(str, strCases = ['anything']),
                             'batchDataStart'           : TestBaseClass.ExtValue(int),
                             'batchDataEnd'             : TestBaseClass.ExtValue(int),      
                             'postProcessing'           : TestBaseClass.ExtValue(bool),               
                             'multiprocessing'          : TestBaseClass.ExtValue(bool),  
                             'numberOfProcessors'       : TestBaseClass.ExtValue(int),
                             'useGeneralMethod'         : TestBaseClass.ExtValue(bool),
                             'useAverageValues'         : TestBaseClass.ExtValue(bool),
                             'useVesselsImpedance'      : TestBaseClass.ExtValue(bool),
                             'useLumpedValues'          : TestBaseClass.ExtValue(bool),
                             'optimizeParams'           : TestBaseClass.ExtValue(bool),
                             'optimizeParamsDataFile'   : TestBaseClass.ExtValue(str, strCases = ['anything']), # path to file
                             'Wkoptimize'               : TestBaseClass.ExtValue(str, strCases = ['anything']), # ["Wk2", "Wk3" ,"Wk4p", "Wk4s"]
                             'params'                   : TestBaseClass.ExtValue(str, strCases = ['anything'])  # ['R1LCR2', 'R1LC', 'R1L', 'R1', 'LCR2', 'LC', 'L', 'CR2', 'C', 'R2']
                             
                           } 
    
    externXmlAttributes  = []
    
    externXmlElements    = ['caseName',
                            'dataNumber',
                            'simulationDescription',
                            'createNetworks',
                            'solveNetworks',
                            'batchDataFile',
                            'batchDataStart',
                            'batchDataEnd',
                            'postProcessing',                 
                            'multiprocessing',
                            'numberOfProcessors',
                            'useGeneralMethod',
                            'useAverageValues',
                            'useVesselsImpedance',
                            'useLumpedValues',
                            'optimizeParams',
                            'optimizeParamsDataFile',
                            'Wkoptimize',
                            'params'
                            
                            ]
    
    
    def __init__(self, CPUTimeFile=None):
        
        self.networkName = None
        self.dataNumber  = None
        
        self.CPUTimeFile = CPUTimeFile
        
        self.caseName = None
        self.simulationDescription = None
        self.createNetworks     = True
        self.solveNetworks     = False
        # ceate simulation case files       
        
        self.batchDataFile = None # path to pickle file containing name of network and info about where to truncate
        self.batchDataStart = 0
        self.batchDataEnd = -1
        #  run simulations 
        self.multiprocessing        = True
        self.numberOfProcessors     = 8
        self.postProcessing  = True
        
        self.useGeneralMethod = False # if true truncatelist only contain the actual vessels to truncate, otherwise the systematic way to truncate
        self.useAverageValues = False
        self.useVesselsImpedance = False
        self.useLumpedValues = False
        
        self.optimizeParams = False
        self.optimizeParamsDataFile = None
        self.Wkoptimize = 'Wk3'
        self.params = 'R1LCR2'
        
    
    
    def loadBatchDatalist(self):
        
        filepath = self.batchDataFile
        batchDataList = pickle.load(open(filepath, 'rb'))
        self.batchDataList = batchDataList[self.batchDataStart: self.batchDataEnd]

        
        
    def loadoptimizeParamsDataFile(self):
        
        if self.optimizeParams:
            filepath = self.optimizeParamsDataFile
            print filepath
            print os.path.isfile(filepath)
            optParamsDict = pickle.load(open(filepath, 'rb'))
        
            self.optParamsDict = optParamsDict
        else:
            self.optParamsDict = None
    
    def createEvaluationCaseFiles(self): 
        '''
        
        batchDataList <list> := with data for each batch job [batchData1, batchData .. ]
            batchData <dict> := dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: 
                                            'reductionNetworkName': , truncateList}
        
        '''
        if self.useLumpedValues and self.optimizeParams:
            print "Error: both useLumpedValues and optimizeParams is set to True"
        
        #TODO: replace create evaluation files or save them to disc!!!
        if self.createNetworks == True:
            
            sampleSize = len(self.batchDataList)
            
            progressBar = cPB.ProgressBar(35, sampleSize)
            
            for simulationIndex, batchData in enumerate(self.batchDataList):
                
                networkName = batchData['networkName']
                reductionNetworkName = batchData['reductionNetworkName']
                truncateList = batchData['truncateList']
                dataNumber = self.dataNumber
                
                vascularNetwork = mXML.loadNetworkFromXML(reductionNetworkName, dataNumber = "xxx")
                vascularNetwork.quiet = True
                
                vascularNetwork.update({'description':self.simulationDescription,
                                        'dataNumber' :self.dataNumber})
                
                New_network = cNred.NetworkReduction(vascularNetwork, quiet=True)
                
                New_network.initialize(useAverageValues=self.useAverageValues, 
                                       useVesselsImpedance=self.useVesselsImpedance, 
                                       useLumpedValues=self.useLumpedValues,
                                       optimizeParams=self.optimizeParams,
                                       optParamsDict=self.optParamsDict,
                                       Wkoptimize=self.Wkoptimize,
                                       params=self.params)
                
                if self.useGeneralMethod:
                    New_network.reduceNetworkFromListGen(truncateList)
                else:
                    New_network.reduceNetworkFromList(truncateList)
                
                New_network.name = networkName
                newNetworkXmlFile =  mFPH_VNR.getFilePath('reductionNetworkXmlFileXXX', networkName, "xxx", 'write', reductionNetworkName=reductionNetworkName, reductionNetworkCase=self.caseName)
                
                mXML.writeNetworkToXML(New_network, dataNumber = dataNumber, networkXmlFile=newNetworkXmlFile)
            

                
                solutionFileXML = mFPH_VNR.getFilePath('reductionNetworkXmlFileSim', networkName, dataNumber, 'write', reductionNetworkName=reductionNetworkName, reductionNetworkCase=self.caseName)
                solutionFilehd5 = mFPH_VNR.getFilePath('reductionSolutionFile', networkName, dataNumber, 'write', reductionNetworkName=reductionNetworkName, reductionNetworkCase=self.caseName)
                
                batchData['networkXmlFileLoad'] = newNetworkXmlFile
                batchData['networkXmlFileSave'] = solutionFileXML
                batchData['pathSolutionDataFilename'] = solutionFilehd5
                batchData['dataNumber'] = self.dataNumber
                
                del vascularNetwork
                del New_network
            
                progressBar.progress(simulationIndex) 
                
        if self.solveNetworks:
            mBatchSim.runBatchAsMultiprocessing(self.batchDataList, CPUTimeFile=self.CPUTimeFile)
            
        
    