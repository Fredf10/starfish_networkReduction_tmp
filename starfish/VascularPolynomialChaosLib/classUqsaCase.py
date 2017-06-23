from __future__ import print_function, absolute_import
from builtins import range
from future.utils import iteritems, iterkeys, viewkeys, viewitems, itervalues, viewvalues
from builtins import input as input3
import sys,os
from starfish.VascularPolynomialChaosLib.testBaseClass import TestBaseClass 
from starfish.VascularPolynomialChaosLib import moduleFilePathHandlerVPC as mFPH_VPC
import progressbarsimple as cPB
from starfish.VascularPolynomialChaosLib import classLocationOfInterestManager
from starfish.VascularPolynomialChaosLib import classUqsaMethods
from starfish.VascularPolynomialChaosLib import classSampleManager
import shutil
import time
import numpy as np
import multiprocessing

class UqsaCase(TestBaseClass):
    
    externVariables      = { 'createSample'             : TestBaseClass.ExtValue(bool),
                             'createEvaluationXmlFiles' : TestBaseClass.ExtValue(bool),  
                             'simulateEvaluations'      : TestBaseClass.ExtValue(bool),  
                             'preProcessData'           : TestBaseClass.ExtValue(bool),  
                             'postProcessing'           : TestBaseClass.ExtValue(bool),  
                             'localEvaluation'          : TestBaseClass.ExtValue(bool),             
                             'multiprocessing'          : TestBaseClass.ExtValue(bool),  
                             'numberOfProcessors'       : TestBaseClass.ExtValue(int),
                             'simulateEvaluationNumbers': TestBaseClass.ExtValue(int, multiVar=True),
                             'sampleManager'            : TestBaseClass.ExtObject({'sampleManager':classSampleManager.SampleManager}),
                             'uqsaMethods'              : TestBaseClass.ExtDict('uqsaMethod', TestBaseClass.ExtObject({'uqsaMethodPolynomialChaos':classUqsaMethods.UqsaMethodPolynomialChaos,
                                                                                                                       'uqsaMethodMonteCarlo'     :classUqsaMethods.UqsaMethodMonteCarlo,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirLR': classUqsaMethods.UqsaMethodPolynomialChaosDepDirLR,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirQR': classUqsaMethods.UqsaMethodPolynomialChaosDepDirQR,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirQL': classUqsaMethods.UqsaMethodPolynomialChaosDepDirQL,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirLRorder': classUqsaMethods.UqsaMethodPolynomialChaosDepDirLRorder,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirFLR' : classUqsaMethods.UqsaMethodPolynomialChaosDepDirFLR,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirFR' : classUqsaMethods.UqsaMethodPolynomialChaosDepDirFR,
                                                                                                                       'uqsaMethodPolynomialChaosDepDirFL' : classUqsaMethods.UqsaMethodPolynomialChaosDepDirFL,
                                                                                                                       'uqsaMethodMonteCarloParametrizedBootstrapping' : classUqsaMethods.UqsaMethodMonteCarloParametrizedBootstrapping},
                                                                                                                        )),
                             'locationOfInterestManager' : TestBaseClass.ExtObject({'LocationOfInterestManager':classLocationOfInterestManager.LocationOfInterestManager}),
                           } 
    
    externXmlAttributes  = []
    
    externXmlElements    = ['createSample',
                            'createEvaluationXmlFiles',
                            'simulateEvaluations',
                            'preProcessData' ,
                            'postProcessing',
                            'localEvaluation',                 
                            'multiprocessing',     
                            'numberOfProcessors',   
                            'simulateEvaluationNumbers',
                            'sampleManager',
                            'uqsaMethods',
                            'locationOfInterestManager']
    
    def __init__(self):
        self.networkName = None
        self.dataNumber  = None
        ### data read in from file
        ##control variables
        # create samples ( TRUE == create and save, FALSE == load existing)
        self.createSample     = True
        # ceate simulation case files       
        self.createEvaluationXmlFiles = True
        #  run simulations 
        self.simulateEvaluations    = True
        self.localEvaluation        = True #TODO: add functions for server
        self.multiprocessing        = True
        self.numberOfProcessors     = 12
        self.simulateEvaluationNumbers = []
        # pre process data for all quantitiy of interest
        self.preProcessData   = True
        self.preProcessData   = True
        # post processing - uncertainty quantification and sensitivity analysis
        self.postProcessing  = True
        # sample manager for the case
        self.sampleManager = None
        ## uqsa method instance of this case to use
        self.uqsaMethods = None # {}
        ## location of interests i.e. datastructure including all evaluated data y = f(z)
        self.locationOfInterestManager = None
        ### data assoziated during run time
        ## samples of Z
        
    def initialize(self,networkName, dataNumber):
        '''
        Initialize case class
        '''
        self.networkName = networkName
        self.dataNumber  = dataNumber
        
        self.locationOfInterestManager.initialize()
        
        for uqsaMethodName in iterkeys(self.uqsaMethods):
            print("Info classUQSACase 96: running ", uqsaMethodName)
        
    def aquireSamples(self, distributionManager, randomInputsExtDist):
        '''
        Function that envokes either sample creation of loading depending on the defined control variable        
        '''
        randomVariableNames = []
        for randomInput in randomInputsExtDist:
            randomVariableNames.append(randomInput.name) 
        
        maxSampleSize = 0
        abcSample = False
        #find out maximum numbers of samples needed and if ABC sample is needed
        for uqsaMethod in itervalues(self.uqsaMethods):
            sampleSizeCurrent, abcSampleCurrent = uqsaMethod.evaluateSamplesSize(distributionManager.distributionDimension)   
            
            if sampleSizeCurrent > maxSampleSize:
                maxSampleSize = sampleSizeCurrent
            if abcSampleCurrent == True:
                abcSample = True
            print("maxSampleSize", maxSampleSize) 
        if self.createSample == True:
            sampleSaveFile = mFPH_VPC.getFilePath('uqsaSampleFile', self.networkName, self.dataNumber, 
                                                  mode = "write", caseName = self.sampleManager.samplingMethod)
            self.sampleManager.createSamples(distributionManager,randomVariableNames, 
                                             maxSampleSize, abcSample, sampleSaveFile)
            
        else:
            sampleLoadFile = mFPH_VPC.getFilePath('uqsaSampleFile', self.networkName, self.dataNumber, 
                                                  mode = "read", caseName = self.sampleManager.samplingMethod)
            self.sampleManager.evaluateLoadedSamples(sampleLoadFile, maxSampleSize, randomVariableNames)
                            
        self.locationOfInterestManager.sampleSize = self.sampleManager.currentSampleSize
    
    def createEvaluationCaseFiles(self): 
        '''
        
        batchDataList <list> := with data for each batch job [batchData1, batchData .. ]
            batchData <dict> := dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
        
        '''
        self.evaluationCaseFiles = [] # list of dict:  [ caseFileDict1,caseFileDict2 ..]  for each evaluation
        
        #TODO: replace create evaluation files or save them to disc!!!
        if self.simulateEvaluations == True or self.preProcessData == True or self.createEvaluationXmlFiles == True:
        
            sampleSize = self.sampleManager.currentSampleSize
            
            print("Create evaluation case file list")
            #progressBar = cPB.ProgressBar(35, sampleSize)
            
            for simulationIndex in range(sampleSize):
                                        
                networkXmlFileLoad = mFPH_VPC.getFilePath('uqsaEvaluationNetworkXmlFile', self.networkName, 
                                                           self.dataNumber, 'write',
                                                           caseName = self.sampleManager.samplingMethod, 
                                                           evaluationNumber=simulationIndex)
                networkXmlFileSave = networkXmlFileLoad
                pathSolutionDataFilename = mFPH_VPC.getFilePath('uqsaEvaluationSolutionDataFile', 
                                                                self.networkName, self.dataNumber, 'write',
                                                                 caseName = self.sampleManager.samplingMethod, 
                                                                 evaluationNumber=simulationIndex)
                
                caseFileDict1= {'simulationIndex': simulationIndex,
                                'networkName': self.networkName,
                                'dataNumber': self.dataNumber,
                                'networkXmlFileLoad': networkXmlFileLoad,
                                'networkXmlFileSave': networkXmlFileSave,
                                'pathSolutionDataFilename': pathSolutionDataFilename}
                
                self.evaluationCaseFiles.append(caseFileDict1)
            
                #progressBar.progress() 
            
    def getSimulationBatchFileList(self):
        '''
        Returns simulation batch file list, as defined in configs
        '''
        
        startIndex = 0
        endIndex   = int(self.sampleManager.currentSampleSize)
        
        newRange = self.simulateEvaluationNumbers
        if len(newRange) == 2:
            # check if the indices are avaliable
            if all([i in range(self.sampleManager.currentSampleSize) for i in newRange]):
                if newRange[0] < newRange[1]:
                    startIndex = newRange[0]
                    endIndex   = newRange[1]
                    
        batchFileList = self.evaluationCaseFiles[startIndex:endIndex+1]
    
        return batchFileList
        
    def preprocessSolutionData(self):
        '''
        envoke preprocessing of solution data
        '''
        if self.preProcessData == True:
            #TODO rename caseName!!
            caseName = '_'.join([self.sampleManager.samplingMethod])
            preprocessedSolutionData = mFPH_VPC.getFilePath('preprocessedDataFile', self.networkName, self.dataNumber, 
                                                     mode = "write", caseName = caseName )
            
            simulationTimeFileSave = mFPH_VPC.getFilePath('simulationTime', self.networkName, self.dataNumber, 
                                                     mode = "write", caseName =  caseName)
            simulationTimeFileLoad = mFPH_VPC.getFilePath('simulationTime', self.networkName, self.dataNumber, 
                                                     mode = "read", caseName = caseName, exception = 'No')
            self.locationOfInterestManager.preprocessSolutionData(self.evaluationCaseFiles,
                                                                  preprocessedSolutionData,
                                                                  simulationTimeFileSave,
                                                                  simulationTimeFileLoad)
                    
    def quantifyUncertaintyAndAnalyseSensitivtiy(self, distributionManager):
        '''
        invoke uq sa process
        '''
        if self.postProcessing == True:
            
            # copy of preprocessed data file 
            caseName = '_'.join([self.sampleManager.samplingMethod])
            preprocessedSolutionData = mFPH_VPC.getFilePath('preprocessedDataFile', self.networkName, self.dataNumber, 
                                                     mode = "read", caseName = caseName)
            uqsaSolutionDataFile = mFPH_VPC.getFilePath('uqsaSolutionDataFile', self.networkName, self.dataNumber, 
                                                     mode = "write", caseName = caseName)
            shutil.copy(preprocessedSolutionData,uqsaSolutionDataFile)
            # open solution file
            self.locationOfInterestManager.openQuantityOfInterestFile(uqsaSolutionDataFile, mode = 'r+')
            # loop through data objects
            for qoi in self.locationOfInterestManager.getQoiIterator():
                multiprocessingUQSA = False
                if multiprocessingUQSA == True:
                    self.multiprocessingUQSA(qoi, distributionManager)
                else:
                    timeStartTotal = time.time()
                    
                    for uqsaMethodName,uqsaMethod in iteritems(self.uqsaMethods):
                        
                        timeStartBatch = time.time()
                        print("calculate uqsa measure for {}".format(uqsaMethodName))
                        
                        uqsaMeasures = uqsaMethod.calculateStatistics(distributionManager, self.sampleManager, qoi)
                        qoi.addUqsaMeasures(uqsaMethodName, uqsaMeasures)
                        
                        timeBatchJob= time.time()-timeStartBatch
                        minutesBatch = int(timeBatchJob/60.)
                        secsBatch = timeBatchJob-minutesBatch*60.
                        print('=====================================')
                        print('runtime:  {} min {} sec'.format(minutesBatch,secsBatch))
                        print('=====================================')
                        print()
            
                    timeTotal= time.time()-timeStartTotal
                    minutesTotal = int(timeTotal/60.)
                    secsTotal = timeTotal-minutesTotal*60.
                    print('=====================================')
                    print('total runtime:  {} min {} sec'.format(minutesTotal,secsTotal))
                    print('=====================================')
                    print()
            
            self.locationOfInterestManager.closeAndSaveQuantityOfInterestFile()
    
    def multiprocessingUQSA(self,qoi, distributionManager):
        '''
        Run all uqsa methods in as a local multiprocess
        '''
        print("running Multiprocessing UQSA")
        timeStartBatch = time.time()
        
        # create batch list for all jobs    
        batchList = []
        for uqsaMethodName,uqsaMethod in iteritems(self.uqsaMethods):
            batchList.append([uqsaMethodName,uqsaMethod,distributionManager,self.sampleManager,qoi])         
        # run jobs
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.imap(self.batchJobUQSA,batchList)
        pool.close() 
        pool.join()
        
        timeBatchJob= time.time()-timeStartBatch
        minutesBatch = int(timeBatchJob/60.)
        secsBatch = timeBatchJob-minutesBatch*60.
        print('=====================================')
        print('total runtime:  {} min {} sec'.format(minutesBatch,secsBatch))
        print('=====================================')
        print()
    
    def batchJobUQSA(self, args):
        '''
        batch job for local uqsa multiprocessing
        '''
        uqsaMethodName, uqsaMethod, distributionManager, sampleManager, qoi = args
        uqsaMeasures = uqsaMethod.calculateStatistics(distributionManager, sampleManager, qoi)
        qoi.addUqsaMeasures(uqsaMethodName, uqsaMeasures)
