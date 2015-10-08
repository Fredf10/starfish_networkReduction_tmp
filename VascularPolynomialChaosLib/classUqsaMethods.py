#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys,os
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur+'/../')

from testBaseClass import TestBaseClass 


class uqsaMethodPolynomialChaos(TestBaseClass):
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    #----External Variables -------------------------------------------------------------#
    externVariables = {  'polynomialOrder'         : TestBaseClass.ExtValue(int), 
                         'samplingMethod'          : TestBaseClass.ExtValue(str, strCases = ['K','R','L','S','H','M','C','NC','G','RG']),
                         'sampleFactor'            : TestBaseClass.ExtValue(int)
                         }
                
    externXmlAttributes  = []
    
    externXmlElements    = ['polynomialOrder', 
                            'samplingMethod',
                            'sampleFactor']
    
    def __init__(self):
        '''
        
        '''
        #----POLYNOMIAL CHAOS DEFINITIONS -------------------------------------------------------------#
        self.sampleFactor = 2
        #polynomialOrders of the polynomial chaos expansion || if more then one given they are processed consecutevely
        self.polynomialOrder = 3
        # method of the spares grid collocation 
        self.samplingMethod = 'M'
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
  
  
class uqsaMethodMonteCarlo(TestBaseClass):
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    #----External Variables -------------------------------------------------------------#
    externVariables = {  'sensitivityAnalysis'   : TestBaseClass.ExtValue(bool), 
                         'samplingMethod'        : TestBaseClass.ExtValue(str, strCases = ['K','R','L','S','H','M','C','NC','G','RG']),
                         'sampleSize'            : TestBaseClass.ExtValue(int)
                         }
                
    externXmlAttributes  = []
    
    externXmlElements    = ['sampleSize', 
                            'samplingMethod',
                            'sensitivityAnalysis']
    
    def __init__(self):
        '''
        
        '''
        #
        self.sampleSize = 10
        #
        self.sensitivityAnalysis = True
        # method of the spares grid collocation 
        self.samplingMethod = 'M'
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