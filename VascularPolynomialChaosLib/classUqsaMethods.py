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
            samplesDependent = distributionManager.jointDistributionDependent.inv(distributionManager.jointDistribution.fwd(samples.T)).transpose()
        else:
            raise ValueError("uqsaMethodPolynomialChaos.createDependentSamples(), cannot create dependent samples as distributionManager.jointDistributionDependent is not defined!")
        return samplesDependent
                
    def calculateOrthogonalPolynomials(self,distributionManager):
        '''
        Method to calculate orthogonal polynomials
        '''
        
        return cp.orth_ttr(self.polynomialOrder,distributionManager.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_gs(self.expansionOrder,self.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_chol(self.expansionOrder,self.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_svd(self.expansionOrder,self.jointDistribution)
  
    def calculateStatistics(self, distributionManager, samples, samplesDependent, data, confidenceAlpha):
        '''
        Function which calculates the gPCExpansion for the given data
        '''
                
        orthogonalPolynomials = self.calculateOrthogonalPolynomials(distributionManager)
                
        # polynomial chaos expansion
        gPCExpansion = cp.fit_regression(orthogonalPolynomials, samples.T, data)
                 
        # statistics
        
        statsDict = {}
                     
        statsDict['expectedValue']  = cp.E(gPCExpansion, distributionManager.jointDistribution)
        statsDict['variance']       = cp.Var(gPCExpansion, distributionManager.jointDistribution)
        statsDict['standardDeviation']   = np.sqrt(statsDict['variance'])
        conficenceInterval  = cp.Perc(gPCExpansion, [confidenceAlpha/2., 100-confidenceAlpha/2.], distributionManager.jointDistribution)
        statsDict['conficenceInterval']  =  conficenceInterval.reshape(2,len(np.atleast_1d(statsDict['expectedValue'])))
        
        # conditional expected values  and sensitivity coefficients
        distributionDimension = len(distributionManager.jointDistribution)
        if distributionDimension > 1:
            # test dependecy or not
            if self.dependentCase == False:
                # independent case: use analytic expression from polynomial chaos expansion
                conditionalExpectedValue = []
                conditionalVariance      = []
                # conditional mean and variance
                for rvIndex in xrange(distributionDimension):
                    currDistMean = cp.E(distributionManager.jointDistribution)
                    currDistMean[rvIndex] = np.nan
                    # reduce polynomials
                    currPolynomTime = gPCExpansion(*currDistMean)
                    conditionalExpectedValue.append(cp.E(currPolynomTime,distributionManager.jointDistribution))
                    conditionalVariance.append(cp.Var(currPolynomTime,distributionManager.jointDistribution))    
                
                statsDict['conditionalExpectedValue'] = conditionalExpectedValue
                statsDict['conditionalVariance']      = conditionalVariance
                
                # sensitivity indices
                statsDict['firstOrderSensitivities'] = cp.Sens_m(gPCExpansion,distributionManager.jointDistribution)
                statsDict['totalSensitivities']      = cp.Sens_t(gPCExpansion,distributionManager.jointDistribution)
            else:
                # dependent rancom variables
                
                # this method is broken
                #sensindices = cp.Sens_nataf(distributionManager.expansionOrder, distributionManager.jointDistributionDependent, distributionManager.samplesDependent.T, self.data)
                #
                
                statsDict['firstOrderSensitivities'] = cp.Sens_m_nataf(self.polynomialOrder,
                                                               distributionManager.jointDistributionDependent,
                                                               samplesDependent.T,
                                                               data)
                
                statsDict['totalSensitivities']      = cp.Sens_t_nataf(self.polynomialOrder,
                                                               distributionManager.jointDistributionDependent,
                                                               samplesDependent.T,
                                                               data)
        
        return statsDict
            
  
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
        
        if distDim == 1:
            print "WARNING: random distribution dimension == 1, no sensitivity analysis possible!"
            self.sensitivityAnalysis = False
        
        ## if sensitivityAnalysis is false, only UQ is done with one matrix
        if self.sensitivityAnalysis == False:
            samples = distributionManager.jointDistribution.sample(self.sampleSize,self.samplingMethod).transpose()   
            # reshape samples
            if distDim == 1:
                samples = samples.reshape(self.sampleSize,1)
                
            return samples,None
        else:
            # createHash table
            self.createSampleMatixHashTable(distDim)
            
            samples = distributionManager.jointDistribution.sample(self.sampleSize*2,self.samplingMethod).transpose()
            
            samplesA = samples[self.matrixHash['A'][0]:self.matrixHash['A'][1]]
            samplesB = samples[self.matrixHash['B'][0]:self.matrixHash['B'][1]]
                
            samplesTest = np.sum((samplesA-samplesB).ravel())
            if samplesTest == 0:
                print "WARNING: samplesA and samplesB are the same!"
        
            for i in xrange(distDim):
                samplesC = samplesB.copy()
                samplesC[:,i] = samplesA[:,i].copy()
                
                samples = np.vstack([samples,samplesC])
            
            return samples,None
    
    def createSampleMatixHashTable(self, distDim):
        '''
        Creates hash table for the sample matrix
        '''
        self.matrixHash = {}
        self.matrixHash['A'] = [0,self.sampleSize]
        self.matrixHash['B'] = [self.sampleSize,2*self.sampleSize]
        for d in xrange(distDim):
            self.matrixHash[''.join(['C',str(d)])] = [self.sampleSize*(d),self.sampleSize*(d+1),self.sampleSize*(d+2),self.sampleSize*(d+3)]
          
    def calculateStatistics(self, distributionManager, samples, samplesDependent, data, confidenceAlpha):
        '''
        Function which calculates the gPCExpansion for the given data
        '''
                
        # statistics
        
        statsDict = {}
        
        if self.dependentCase == False:
            
            
            # check if samples, and data have the right format
            distDim = distributionManager.distributionDimension
            if len(samples) != (distDim+2)*self.sampleSize:
                print 'WARNING: Not simulated as sensitivity analysis case as len(samples) != (distDim+2)*self.sampleSize'
                self.sensitivityAnalysis = False
                
            if self.sensitivityAnalysis == False:
            
                statsDict['expectedValue']       = np.sum(data,axis=0)/self.sampleSize
                # variance = (data**2-mean)/len(samples)
                statsDict['variance']            = np.sum(data**2-statsDict['expectedValue']**2,axis=0)/len(samples)
                statsDict['standardDeviation']   = np.sqrt(statsDict['variance'])
                
                quantiles = [confidenceAlpha/2, 100.-confidenceAlpha/2]
                statsDict['conficenceInterval']  = np.percentile(data,quantiles, axis=0)
            
            else:
                      
                self.createSampleMatixHashTable(distDim)
                
                dataA = data[self.matrixHash['A'][0]:self.matrixHash['A'][1]]
                dataB = data[self.matrixHash['B'][0]:self.matrixHash['B'][1]]
                dataC = np.empty((distDim,self.sampleSize,len(data[0])))
                for d in xrange(distDim):
                    key = ''.join(['C',str(d)])
                    s = self.matrixHash[key][2]
                    e = self.matrixHash[key][3]
                    dataC[d] = data[s:e]
                
                meanA = np.sum(dataA,axis=0)/self.sampleSize
                meanB = np.sum(dataB,axis=0)/self.sampleSize
                #print np.shape(dataC), dataC
                #meanC = np.mean(dataC,axis= 0)
                #print np.shape(meanC), meanC
                mean = (meanA+meanB)/2.0
                dataAmm = dataA -mean
                dataBmm = dataB -mean
                dataCmm = dataC -mean
                # sensitivity
                f0sq = np.mean(dataAmm*dataBmm, axis=0)
                
                varianceA = np.sum(dataAmm**2,axis=0)/self.sampleSize - f0sq
                varianceB = np.sum(dataAmm**2,axis=0)/self.sampleSize - f0sq
            
                # conditional variance
                conditionalVarianceGivenQ = np.empty((distDim,len(data[0])))
                conditionalVarianceNotQ   = np.empty((distDim,len(data[0])))
            
                Si  =  []
                STi =  []
            
                for i in xrange(distDim):
                    conditionalVarianceGivenQ[i] = np.sum(dataAmm*dataCmm[i],axis=0)/self.sampleSize - f0sq
                    conditionalVarianceNotQ[i]   = np.sum(dataBmm*dataCmm[i],axis=0)/self.sampleSize - f0sq
                    Si.append(conditionalVarianceGivenQ[i]/(varianceA+(varianceA==0))*(varianceA!=0))
                    STi.append(1.- conditionalVarianceNotQ[i]/(varianceA+(varianceA==0))*(varianceA!=0))
                   
                statsDict['expectedValue']      = mean
                statsDict['variance']           = varianceA
                
                if varianceA.all() > 0:
                    statsDict['standardDeviation']  = np.sqrt(varianceA)
                
                quantiles = [confidenceAlpha/2, 100.-confidenceAlpha/2]
                statsDict['conficenceInterval'] = np.percentile(data,quantiles, axis=0) 
                
                statsDict['conditionalExpectedValue'] = None
                statsDict['conditionalVariance']      = np.array(conditionalVarianceGivenQ)
                
                # sensitivity indices
                statsDict['firstOrderSensitivities'] = np.array(Si)
                statsDict['totalSensitivities']      = np.array(STi)
        
        else:
            # TODO: implement dependent dist methods for MC
            raise NotImplementedError("NO MC methods for dependet distributions implemented")
        
        return statsDict  
            