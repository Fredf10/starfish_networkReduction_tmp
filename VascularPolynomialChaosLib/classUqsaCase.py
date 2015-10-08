import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))
sys.path.append(cur+'/../')

from testBaseClass import TestBaseClass 

import classLocationOfInterestManager
import classUqsaMethods

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
        ## dependent samples
        self.sampleDependent = None
        
        
        
        
        