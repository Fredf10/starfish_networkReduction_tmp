#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys,os
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur+'/../')

from testBaseClass import TestBaseClass 

import chaospy as cp
import numpy as np

class uqsaMethodPolynomialChaos(TestBaseClass):
    
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    #----External Variables -------------------------------------------------------------#
    externVariables = {  'polynomialOrder' : TestBaseClass.ExtValue(int), 
                         'samplingMethod'  : TestBaseClass.ExtValue(str, strCases = ['K','R','L','S','H','M','C','NC','G','RG']),
                         'sampleFactor'    : TestBaseClass.ExtValue(int),
                         'dependentCase'   : TestBaseClass.ExtValue(bool), 
                         }
                
    externXmlAttributes  = []
    
    externXmlElements    = ['polynomialOrder', 
                            'samplingMethod',
                            'sampleFactor',
                            'dependentCase']
    
    def __init__(self):
        '''
        
        '''
        
        
        self.methodType = 'polynomialChaos'
        
        #----POLYNOMIAL CHAOS DEFINITIONS -------------------------------------------------------------#
        self.dependentCase = False
        
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
        
    def name(self):
        '''
        returns name of the case
        '''
        return '_'.join(['polynomialChaos',self.samplingMethod,'Ord',str(self.polynomialOrder),'sf',str(self.sampleFactor)])
  
    def createSamples(self, distributionManager):
        '''
        create samples for the defined distribution for given samplesSize and sampleMethod
        using the chaospy toolbox
        Input
            samplesSize : int,array_like
        
            sampleMethod : str
                (from chaospy)
                Alternative sampling techniques
            
                Normal sampling schemes
                Key     Name                Nested
                ----    ----------------    ------
                "K"     Korobov             no
                "R"     (Pseudo-)Random     no
                "L"     Latin hypercube     no
                "S"     Sobol               yes
                "H"     Halton              yes
                "M"     Hammersley          yes
            
                Grided sampling schemes
                Key     Name                Nested
                ----    ----------------    ------
                "C"     Chebyshev nodes     maybe
                "G"     Gaussian quadrature no
                "E"     Gauss-Legendre      no
            
            expansionOrder: float
                calculate optimal samplesSize for gPC with the following rule
                    samplesSize =  2* number gPC-expansion terms
            
        '''
        # calculate samplesSize from expansion order 
        samplesSize = int(self.sampleFactor*cp.terms(self.polynomialOrder, distributionManager.distributionDimension))
        # creating samples
        samples = distributionManager.jointDistribution.sample(samplesSize,self.samplingMethod).transpose()   
        # reshape samples
        if distributionManager.distributionDimension == 1:
            samples = samples.reshape(samplesSize,1)
        # create dependent samples if correlation exists
        samplesDependent = None
        if self.dependentCase == True:
            samplesDependent = self.createDependentSamples(distributionManager, samples)
        return samples, samplesDependent
  
    def createDependentSamples(self,distributionManager, samples):
        '''
        Create dependen samples if dependen case = True and dependent random variables true
        '''
        if distributionManager.jointDistributionDependent != None:
            samplesDependent = distributionManager.jointDistributionDependent.inv(distributionManager.jointDistribution.fwd(samples))
        else:
            raise ValueError("uqsaMethodPolynomialChaos.createDependentSamples(), cannot create dependent samples as distributionManager.jointDistributionDependent is not defined!")
        return samplesDependent
                
    def calculateOrthogonalPolynomials(self,distributionManager):
        '''
        Method to calculate orthogonal polynomials
        '''
        
        self.orthogonalPolynomials = cp.orth_ttr(self.expansionOrder,distributionManager.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_gs(self.expansionOrder,self.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_chol(self.expansionOrder,self.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_svd(self.expansionOrder,self.jointDistribution)
  
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
        self.methodType = 'MonteCarlo'        
        #
        self.dependentCase = False # not implemented yet!
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
    def name(self):
        '''
        returns name of the case
        '''
        return '_'.join(['MonteCarlo',self.samplingMethod,str(self.sampleSize)])
    
    def createSamples(self, distributionManager):
        '''
        create samples
        
        samplesA, samplesB and samplesC
        '''
        
        distDim = distributionManager.distributionDimension
        
        ## if sensitivityAnalysis is false, only UQ is done with one matrix
        if self.sensitivityAnalysis == False:
            samples = distributionManager.jointDistribution.sample(self.sampleSize,self.samplingMethod).transpose()   
            # reshape samples
            if distDim == 1:
                samples = samples.reshape(self.sampleSize,1)
                
            return samples,None
        else:
            #TODO: creation of proper A,B,C matrices
            #TODO create Hash for A,B and Ci
            
            # creating samples
            print self.sampleSize*2
            samples = distributionManager.jointDistribution.sample(self.sampleSize*2,self.samplingMethod).transpose()
            
            samplesA = samples[0:self.sampleSize]
            samplesB = samples[self.sampleSize::]
            if distDim == 1:
                samplesA = samplesA.reshape(self.sampleSize,1)
                samplesB = samplesB.reshape(self.sampleSize,1)
                
            samplesTest = np.sum((samplesA-samplesB).ravel())
            if samplesTest == 0:
                print "WARNING: samplesA and samplesB are the same!"
        
            samplesC = np.empty((distDim,self.sampleSize,distDim))
            for i in xrange(distDim):
                samplesC[i,:,:] = samplesB.copy()
                samplesC[i,:,i] = samplesA[:,i].copy()
            
            return samples,None
        