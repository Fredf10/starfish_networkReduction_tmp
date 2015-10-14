import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))
sys.path.append(cur+'/../')

from testBaseClass import TestBaseClass 

import moduleFilePathHandlerVPC as mFPH_VPC

import classLocationOfInterestManager
import classUqsaMethods

import h5py

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
                             'uqsaMethod'               : TestBaseClass.ExtObject({'uqsaMethodPolynomialChaos':classUqsaMethods.uqsaMethodPolynomialChaos,
                                                                                   'uqsaMethodMonteCarlo'     :classUqsaMethods.uqsaMethodMonteCarlo}),
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
                            'uqsaMethod',
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
        self.simulateEvaluationNumbers      = []
        # pre process data for all quantitiy of interest
        self.preProcessData   = True
        # post processing - uncertainty quantification and sensitivity analysis
        self.postProcessing  = True
        
        ## uqsa method instance of this case to use
        self.uqsaMethod                = None
        
        ## location of interests i.e. datastructure including all evaluated data y = f(z)
        self.locationOfInterestManager = None
        
        ### data assoziated during run time
        ## samples of Z
        self.randomVariableNames = []
        self.samples = None
        self.samplesSize = None
        self.samplesDependent = None
        
    def initialize(self,networkName, dataNumber):
        '''
        Initialize case class
        '''
        self.networkName = networkName
        self.dataNumber  = dataNumber
        
        self.locationOfInterestManager.initialize()
        
        print "Info 88: uqsaCase running ", self.uqsaMethod.name()
        
    def aquireSamples(self, distributionManager, randomInputsExtDist):
        '''
        Function that envokes either sample creation of loading depending on the defined control variable        
        '''
        self.randomVariableNames = []
        for randomInput in randomInputsExtDist:
            self.randomVariableNames.append(randomInput.name) 
        
        if self.createSample == True:
            self.samples,self.samplesDependent = self.uqsaMethod.createSamples(distributionManager)
            self.samplesSize = len(self.samples)
            if len(self.randomVariableNames) != len(self.samples[0]):
                raise ValueError("Created sample matrix does not match with defined random variable vector {}".format(self.randomVariableNames))
            sampleFile = mFPH_VPC.getFilePath('uqsaSampleFile', self.networkName, self.dataNumber, mode = "write", caseName=self.uqsaMethod.name())
            self.saveSamples(sampleFile)
        else:
            sampleFile = mFPH_VPC.getFilePath('uqsaSampleFile', self.networkName, self.dataNumber, mode = "read", caseName=self.uqsaMethod.name())
            self.loadSamples(sampleFile)

        self.locationOfInterestManager.sampleSize = self.samplesSize
    
    def loadSamples(self, sampleFile):
        '''
        load the current sample to disc so it is available for postprocessing or
        sequencielle working process
        for generation gPCE the sample nodes corresponding to the data are needed.
        '''
        f = h5py.File(sampleFile,'r')
        dset = f['sampleSpace']
        self.samples = dset[:]
        self.samplesSize    = dset.attrs.get('samplesSize')
        randomVariableNames = dset.attrs.get('randomVariableNames')
        
        if randomVariableNames != self.randomVariableNames:
            raise ValueError("Loaded randomVariableNames {} for samples does not match with defined random variable vector {}".format(randomVariableNames, self.randomVariableNames))
                    
        if len(self.randomVariableNames) != len(self.samples[0]):
            raise ValueError("Created sample matrix does not match with defined random variable vector {}".format(self.randomVariableNames))
                
        if 'sampleSpaceDependent' in f.keys():
            dset = f['sampleSpaceDependent']
            self.samplesDependent = dset[:]
            
        f.close()
      
    def saveSamples(self, sampleFile):
        '''
        save the current sample to disc so it is available for postprocessing or
        sequencielle working process
        for generation gPCE the sample nodes corresponding to the data are needed.
        '''        
        f = h5py.File(sampleFile,'w')
        dset = f.create_dataset("sampleSpace", data=self.samples)
        dset.attrs.create('samplesSize', data=self.samplesSize)
        dset.attrs.create('randomVariableNames', data = self.randomVariableNames)
        if self.samplesDependent != None:
            dset = f.create_dataset("sampleSpaceDependent", data=self.samplesDependent)
        
        f.flush()
        f.close()  
        
    def createEvaluationCaseFiles(self): 
        '''
        
        batchDataList <list> := with data for each batch job [batchData1, batchData .. ]
            batchData <dict> := dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
        
        '''
        self.evaluationCaseFiles = [] # list of dict:  [ caseFileDict1,caseFileDict2 ..]  for each evaluation
        
        for simulationIndex in xrange(self.samplesSize):
                                    
            networkXmlFileLoad = mFPH_VPC.getFilePath('uqsaEvaluationNetworkXmlFile', self.networkName, self.dataNumber, 'write',
                                                       caseName=self.uqsaMethod.name(), evaluationNumber=simulationIndex)
            networkXmlFileSave = networkXmlFileLoad
            pathSolutionDataFilename = mFPH_VPC.getFilePath('uqsaEvaluationSolutionDataFile', self.networkName, self.dataNumber, 'write',
                                                             caseName=self.uqsaMethod.name(), evaluationNumber=simulationIndex)
            
            caseFileDict1= {'simulationIndex': simulationIndex,
                            'networkName': self.networkName,
                            'dataNumber': self.dataNumber,
                            'networkXmlFileLoad': networkXmlFileLoad,
                            'networkXmlFileSave': networkXmlFileSave,
                            'pathSolutionDataFilename': pathSolutionDataFilename}
            
            self.evaluationCaseFiles.append(caseFileDict1)
            
    def getSample(self, sampleIndex):
        '''
        Returns sample for a certain sampleIndex
        
        dependenSample or independenSample is chooses dependent on the case definitions
        '''
        if sampleIndex in xrange(self.samplesSize):
            if self.uqsaMethod.dependentCase == False:
                sample = self.samples[sampleIndex]
            else:
                sample = self.samplesDependent[sampleIndex]
        else:
            raise ValueError('sampleIndex {} out of range 0:{}'.format(sampleIndex,self.samplesSize))
                             
        return sample
    
    def getSimulatioNBatchFileList(self):
        '''
        Returns simulation batch file list, as defined in configs
        '''
    
        startIndex = 0
        endIndex   = int(self.samplesSize)
        
        newRange = self.simulateEvaluationNumbers
        if len(newRange) == 2:
            # check if the indices are avaliable
            if all([i in xrange(self.samplesSize) for i in newRange]):
                if newRange[0] < newRange[1]:
                    startIndex = newRange[0]
                    endIndex   = newRange[1]
                    
        batchFileList = self.evaluationCaseFiles[startIndex:endIndex+1]
    
        return batchFileList
    
    def saveCaseSolutionDataFile(self):
        '''
        save hdf solution data file
        '''
        uqsaSolutionDataFile = mFPH_VPC.getFilePath('uqsaSolutionDataFile', self.networkName, self.dataNumber, mode = "write", caseName = self.uqsaMethod.name() )
        self.locationOfInterestManager.saveQuantitiyOfInterestData(uqsaSolutionDataFile)
    
    def preprocessSolutionData(self):
        '''
        envoke preprocessing of solution data
        '''
        if self.preProcessData == True:
            self.locationOfInterestManager.preprocessSolutionData(self.evaluationCaseFiles)
            self.saveCaseSolutionDataFile()
            
            
    def quantifyUncertaintyAndAnalyseSensitivtiy(self, distributionManager):
        '''
        evnoke uq sa process
        '''
        if self.postProcessing == True:
        
            # loop through data objects
            for qoi in self.locationOfInterestManager.getQoiIterator():
                stats = self.uqsaMethod.calculateStatistics(distributionManager, self.samples, self.samplesDependent, qoi.getData(),qoi.confidenceAlpha)
                qoi.update(stats)
            self.saveCaseSolutionDataFile()
    