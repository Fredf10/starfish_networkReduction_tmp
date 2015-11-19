import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))
sys.path.append(cur+'/../')

from testBaseClass import TestBaseClass 

import moduleFilePathHandlerVPC as mFPH_VPC

import classLocationOfInterestManager
import classUqsaMethods
import classSampleManager

import shutil

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
                                                                                                                       'uqsaMethodPolynomialChaosDepDir': classUqsaMethods.UqsaMethodPolynomialChaosDepDir})),
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
        
        for uqsaMethodName in self.uqsaMethods.iterkeys():
            print "Info classUQSACase 96: running ", uqsaMethodName
        
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
        for uqsaMethod in self.uqsaMethods.itervalues():
            sampleSizeCurrent,abcSampleCurrent = uqsaMethod.evaluateSamplesSize(distributionManager.distributionDimension)   
            
            if sampleSizeCurrent > maxSampleSize:
                maxSampleSize = sampleSizeCurrent
            if abcSampleCurrent == True:
                abcSample = True
                
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
        
        for simulationIndex in xrange(self.sampleManager.currentSampleSize):
                                    
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
            
        
    def getSimulatioNBatchFileList(self):
        '''
        Returns simulation batch file list, as defined in configs
        '''
    
        startIndex = 0
        endIndex   = int(self.sampleManager.currentSampleSize)
        
        newRange = self.simulateEvaluationNumbers
        if len(newRange) == 2:
            # check if the indices are avaliable
            if all([i in xrange(self.sampleManager.currentSampleSize) for i in newRange]):
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
            
            simulationTimeFileSave  = mFPH_VPC.getFilePath('simulationTime', self.networkName, self.dataNumber, 
                                                     mode = "write", caseName =  caseName)
            simulationTimeFileLoad = mFPH_VPC.getFilePath('simulationTime', self.networkName, self.dataNumber, 
                                                     mode = "read", caseName = caseName, exception = 'No')
            
            self.locationOfInterestManager.preprocessSolutionData(self.evaluationCaseFiles,
                                                                  preprocessedSolutionData,
                                                                  simulationTimeFileSave,
                                                                  simulationTimeFileLoad)
                    
    def quantifyUncertaintyAndAnalyseSensitivtiy(self, distributionManager):
        '''
        evnoke uq sa process
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
                
                import time
                timeStartBatch = time.time()
                import numpy as np
                basis = np.linspace(0.27,0.29,100)
                print "hashDataForGivenBases {}".format(basis)
                qoi.hashDataForGivenBases(basis, self.sampleManager.currentSampleSize)
                
                timeBatchJob= time.time()-timeStartBatch
                minutesBatch = int(timeBatchJob/60.)
                secsBatch = timeBatchJob-minutesBatch*60.
                print '====================================='
                print 'total runtime:  {} min {} sec'.format(minutesBatch,secsBatch)
                print '====================================='
                print
                
                for uqsaMethodName,uqsaMethod in self.uqsaMethods.iteritems():
                    
                    timeStartBatch = time.time()
                    print "calculate uqsa measure for {}".format(uqsaMethodName)
                    
                    uqsaMeasures = uqsaMethod.calculateStatistics(distributionManager, self.sampleManager, qoi)
                    qoi.addUqsaMeasures(uqsaMethodName, uqsaMeasures)
                    
                    timeBatchJob= time.time()-timeStartBatch
                    minutesBatch = int(timeBatchJob/60.)
                    secsBatch = timeBatchJob-minutesBatch*60.
                    print '====================================='
                    print 'total runtime:  {} min {} sec'.format(minutesBatch,secsBatch)
                    print '====================================='
                    print
            
            self.locationOfInterestManager.closeAndSaveQuantityOfInterestFile()
    