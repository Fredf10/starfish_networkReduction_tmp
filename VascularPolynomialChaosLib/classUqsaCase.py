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
        
    def aquireSamples(self, distributionManager):
        '''
        Function that envokes either sample creation of loading depending on the defined control variable        
        '''
        if self.createSample == True:
            self.samples,self.samplesDependent = self.uqsaMethod.createSamples(distributionManager)
            self.samplesSize = len(self.samples)
            
            sampleFile = mFPH_VPC.getFilePath('uqsaSampleFile', self.networkName, self.dataNumber, mode = "write", caseName=self.uqsaMethod.name())
            self.saveSamples(sampleFile)
        else:
            sampleFile = mFPH_VPC.getFilePath('uqsaSampleFile', self.networkName, self.dataNumber, mode = "read", caseName=self.uqsaMethod.name())
            self.loadSamples(sampleFile)

    
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
        if self.samplesDependent != None:
            dset = f.create_dataset("sampleSpaceDependent", data=self.samplesDependent)
        
        f.flush()
        f.close()  
        