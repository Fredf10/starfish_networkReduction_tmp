#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys,os
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur+'/../')

from testBaseClass import TestBaseClass 


class VpcConfiguration(TestBaseClass):
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    #----External Variables -------------------------------------------------------------#
    externVariables = {  'createSample'             : TestBaseClass.ExtValue(bool),
                         'createEvaluationXmlFiles' : TestBaseClass.ExtValue(bool),  
                         'simulateEvaluations'      : TestBaseClass.ExtValue(bool),  
                         'preProcessData'           : TestBaseClass.ExtValue(bool),  
                         'postProcessing'           : TestBaseClass.ExtValue(bool),  
                         'localEvaluation'          : TestBaseClass.ExtValue(bool),             
                         'multiprocessing'          : TestBaseClass.ExtValue(bool),  
                         'numberOfProcessors'       : TestBaseClass.ExtValue(int),
                         'evaluationNumbers'        : TestBaseClass.ExtValue(int, multiVar=True),
                         'runPolynomialChaos'       : TestBaseClass.ExtValue(bool),
                         'polynomialOrders'         : TestBaseClass.ExtValue(int, multiVar=True), 
                         'sampleMethod'             : TestBaseClass.ExtValue(str, strCases = ['K','R','L','S','H','M','C','NC','G','RG']),
                         'runMonteCarlo'            : TestBaseClass.ExtValue(bool)
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
                            'evaluationNumbers',
                            'runPolynomialChaos',
                            'polynomialOrders', 
                            'sampleMethod',
                            'runMonteCarlo']
    
    def __init__(self):
        '''
        
        '''
        
        ### network name and datanumber
        #self.networkName = networkName
        #self.dataNumber  = dataNumber
        
        #control variables
        ##  0.2 collocation method ( TRUE == create and save, FALSE == load existing)
        self.createSample     = True
        
        ### 1.step genrealized polynomial chaos evaluations + data storing
        self.createEvaluationXmlFiles = True
        
        self.simulateEvaluations    = True
        self.localEvaluation        = True #TODO: add functions for server
        self.multiprocessing        = True
        self.numberOfProcessors     = 12
        self.evaluationNumbers      = []
        
        # 2.1 pre process data for all locations of interest
        self.preProcessData   = True
        
        ### 3.step post processing - orthogonal polynomials, gpce expansion, uncertainty quantification and sensitivity analysis
        self.postProcessing  = True
        
        #----POLYNOMIAL CHAOS DEFINITIONS -------------------------------------------------------------#
        self.runPolynomialChaos = True
        #polynomialOrders of the polynomial chaos expansion || if more then one given they are processed consecutevely
        self.polynomialOrders = [2,3]
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
        
        #----MONTE CARLO DEFINITIONS -------------------------------------------------------------#
        self.runMonteCarlo = False
        