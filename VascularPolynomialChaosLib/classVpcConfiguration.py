#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys,os
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append('/'.join([cur,'..','UtilityLib']))
import moduleXML
import moduleFilePathHandlerVPC as mFPH_VPC

class VpcConfiguration(object):
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    def __init__(self, networkName, dataNumber):
        '''
        define all variables with comments here
        
        try to:
        update variables from vpcConfig-file (networkName,dataNumber)
        '''
        ### network name and datanumber
        self.networkName = networkName
        self.dataNumber  = dataNumber
        
        ### 0. set up polychaos variables
        self.createDistributions = False  #(TRUE == create and save, FALSE == load existing)
        ##  0.1 orthogonal polynoms ( TRUE == create and save, FALSE == load existing)
        self.createOrthoPoly  = False
        ##  0.2 collocation method ( TRUE == create and save, FALSE == load existing)
        self.createSample     = False
        
        ### 1.step genrealized polynomial chaos evaluations + data storing
        self.runSimulations   = False
        
        ### 2.step Construct generalized polynomial chaos expansion and  pre process data
        self.calculateGPCE    = False
        #control variables
        # 2.1 is the data already preprocessed?
        self.preProcessData   = True
        # 2.2 create plots for min max points
        self.plotMinMaxPoints = True
            
        ### 3.step post processing - sensitivity analysis
        self.postProcessing  = True
        self.plotMeanSTD     = True
        self.plotPeaks       = True
        #plotSensitiviy  = False
        
        #----POLYNOMIAL CHAOS DEFINITIONS -------------------------------------------------------------#
        ####    
        
        #polynomialOrders of the polynomial chaos expansion || if more then one given they are processed consecutevely
        self.polynomialOrders = [3]
        # method of the spares grid collocation 
        self.sampleMethod = 'M'
        #Parameters
        #----------
        #sample : str
        #         Normal sampling schemes
        #         Key     Name                Nested
        #         ----    ----------------    ------
        #         "K"     Korobov             no
        #         "R"     (Pseudo-)Random     no
        #         "L"     Latin hypercube     no
        #         "S"     Sobol               yes
        #         "H"     Halton              yes
        #         "M"     Hammersley          yes
        #     
        #         Grided sampling schemes
        #         Key     Name                Nested
        #         ----    ----------------    ------
        #         "C"     Chebyshev nodes     no
        #         "NC"    Nested Chebyshev    yes
        #         "G"     Gaussian quadrTrueature no
        #         "RG"    Regular grid        no
        # 
        #  
        
        
        #----CHOICE OF INVESTIGATION POINTS FOR WHICH THE POLYNOMIAL CHAOS IS CALCULATED-----------#

        
        # [[vesselId, node],[vesselId, node],[vesselId, node]]
        self.locationsToEvaluate = [[0,0] ]#   , [0,99]  , [0,199],[1,99],[1,199]]
        self.locationNames = ['A' ]#,'B','C','D','E']
        
        #### WAVE SPLITTING
        # linear or non linear wave splitting ( non linear =  default)
        self.linearWaveSplit = True
        # velocity profile coefficient (1 = linear, 2= peuiseulle ..)
        self.velocityProfileCoefficient = 2.0
        
        #### PEACK EVALUATION
        # fine tune min max function delta to find all the peaks
        self.delta ={'A': 
                {'Pressure' :1.0,'Flow':1.0e-8,
                'Pressure_f':1.0,'Flow_f':1.0e-8,
                'Pressure_b':1.0,'Flow_b':1.0e-8} ,
                'B':
                {'Pressure' :2.7,'Flow':0.5e-7,
                'Pressure_f':2.7,'Flow_f':0.5e-7,
                'Pressure_b':2.7,'Flow_b':0.5e-7} , 
                'C':
                {'Pressure' :2.7,'Flow':0.5e-7,
                'Pressure_f':2.7,'Flow_f':0.5e-7,
                'Pressure_b':2.7,'Flow_b':0.5e-7} ,
                'D':
                {'Pressure' :2.7,'Flow':0.5e-7,
                'Pressure_f':2.7,'Flow_f':0.5e-7,
                'Pressure_b':2.7,'Flow_b':0.5e-7} ,
                'E':
                {'Pressure' :2.7,'Flow':0.5e-7,
                'Pressure_f':2.7,'Flow_f':0.5e-7,
                'Pressure_b':2.7,'Flow_b':0.5e-7} }
        
        #check extrema in minMaxPointsand write number of peaks you want compare starting with 0 in the brackets
        self.peaksToEvaluate =    {'A': 
                            {'extremaPressure':[0],'extremaFlow':[0],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[0],'extremaFlow_b':[1]} ,
                            'B':
                            {'extremaPressure':[0],'extremaFlow':[0,1],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[0],'extremaFlow_b':[0]},  
                            'C':
                            {'extremaPressure':[0],'extremaFlow':[0],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[0],'extremaFlow_b':[0]} , 
                            'D':
                            {'extremaPressure':[0],'extremaFlow':[0],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[],'extremaFlow_b':[]} , 
                            'E':
                            {'extremaPressure':[0],'extremaFlow':[0],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[],'extremaFlow_b':[]}  }
        
        #--POSTPROCESSING: PLOTTING-------------------------------------------------------------------------------------#
        ### 
        
        self.plotDirectory = 'PolynomialChaos'
        ## load additional data for postprocessing
        # load additional simulation results from deterministic simulations
        # (stored in the same directory as the xml file
        self.deterministicDataSetNumbers = [] # [200,199,201] # not in use
        # load additional polynomials with a different order then defined above 
        # if 0 then the current order is respected only
        self.polynomsToPlotOrder = [0]
        
        ##########################################################################################
        ## MEAN and STD plots
        # bool for confidence interval area for 100-plotMeanConfidenceAlpha % confidence
        self.plotMeanConfidenceInterval = False 
        self.plotMeanConfidenceAlpha    = 1.0 
        # bool for sigma interval area
        self.plotMeanSigmaInterval = False
            
        ##########################################################################################
        ## Plot of Peak Analysis
        
        # do peak analysis and save it if True // else: load peak file 
        self.peakAnalysis = False 
        # Confidence interval with 100-plotPeaksConfidenceAlpha % confidence
        self.plotPeaksConfidenceAlpha = 1.0
        
        ### analytic bar plots
        # compare occ.time and amplitude to analytic values (e.g. Sherwin Bifurcation) in bar plots
        self.plotPeaksAnalyticSensitivity = False
        
        ### Mean and STD Box plots 
        # save plots for each parameters in a seperate file if True // else: save one countious
        self.plotPeaksMeanSTDBoxPlotsSingle = False
        
        ##########################################################################################  
        #---Additional hidden fine tuning parameters---
        # set start and end value for the time-line, if plotxTimeEnd = 'totalTime' simulation-endtime is taken
        self.totalTime = 'vascularNetwork.totalTime'
        
        self.plotxTimeStart = 0.0 
        self.plotxTimeEnd   = 'totalTime'
        self.plotXLablesTime = ['0','0.2','0.4','0.6','0.8']
        self.plotYLabelNumber = 5
            
        self.startEvaluation = 0 ## start the polychaos evaluations from number xxx
        
        self.plotFileType = '.pdf'#'.png'
        ## MEAN and STD plots
        self.latexUnits = False
        if self.latexUnits: self.plotMeanStdQuantities = {'Pressure':r'$[mmHg]$', 'Pressure_f':r'$[mmHg]$', 'Pressure_b':r'$[mmHg]$',
                                                          'Flow':r'$[\frac{ml}{s}]$','Flow_f':r'$[\frac{ml}{s}]$','Flow_b':r'$[\frac{ml}{s}]$'}
        else: self.plotMeanStdQuantities = {'Pressure':'[mmHg]'} #, 'Pressure_f':'[mmHg]', 'Pressure_b':r'[mmHg]', 
                                            #'Flow':'[ml/s]','Flow_f':r'[ml/s]','Flow_b':r'[ml/s]'}
        
        #fine tune y-axis for Mean  and STD plots if [0,0] use automatic-limit-caluclation
        self.limits= {'Pressure':  [ 0,0],'Flow':  [ 0,0],
                      'Pressure_b':[ 0,0],'Flow_b':[ 0,0],
                      'Pressure_f':[ 0,0],'Flow_f':[ 0,0]}
         
         
        #--Updating: data from file------------------------------------------------------------------------------------#
        
        vpcConfigXmlFile =  mFPH_VPC.getFilePath('vpcConfigXmlFile', networkName, dataNumber, 'read')
    
        self.update(moduleXML.loadPolyChaosXML(vpcConfigXmlFile))
        
        
        
    def update(self, Dict):
        '''
        updates the class data using a dictionary in from of 
        dataDict = {'variableName': value}
        '''
        for key,value in Dict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except:
                print "ERROR VpcConfiguration.updateData (vesselId {}) Wrong key: {}, could not update varibale".format(self.Id, key)
    