#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys,os
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur+'/../')

from testBaseClass import TestBaseClass 
from classUqsaMeasures import UqsaMeasures

import chaospy as cp
import numpy as np

class UqsaMethodPolynomialChaos(TestBaseClass):
    
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    #----External Variables -------------------------------------------------------------#
    externVariables = {  'polynomialOrder' : TestBaseClass.ExtValue(int), 
                         'sampleFactor'    : TestBaseClass.ExtValue(int),
                         }
                
    externXmlAttributes  = []
    
    externXmlElements    = ['polynomialOrder', 
                            'sampleFactor']
    
    def __init__(self):
        '''
        
        '''                
        self.sampleFactor = 2
        #polynomialOrders of the polynomial chaos expansion || if more then one given they are processed consecutevely
        self.polynomialOrder = 3
                   
    def evaluateSamplesSize(self, distributionDimension):
        '''
        function to evaluate the sample size
        '''
        # calculate samplesSize from expansion order 
        samplesSize = int(self.sampleFactor*cp.terms(self.polynomialOrder, distributionDimension))
        abcSample = False
        return samplesSize,abcSample
                    
    def calculateOrthogonalPolynomials(self,distributionManager):
        '''
        Method to calculate orthogonal polynomials
        '''
        
        return cp.orth_ttr(self.polynomialOrder,distributionManager.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_gs(self.expansionOrder,self.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_chol(self.expansionOrder,self.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_svd(self.expansionOrder,self.jointDistribution)
  
    def calculateStatistics(self, distributionManager, sampleManager, qoi):
        
        sampleSize,abcSample = self.evaluateSamplesSize(distributionManager.distributionDimension)
        samples,samplesDependent = sampleManager.getSampleMatrices(sampleSize,abcSample)
        data                     = qoi.getData(sampleSize,abcSample)
        
        dependentCase = sampleManager.dependentCase
        confidenceAlpha = qoi.confidenceAlpha
        '''
        Function which calculates the gPCExpansion for the given data
        '''
                
        orthogonalPolynomials = self.calculateOrthogonalPolynomials(distributionManager)
                
        # polynomial chaos expansion
        gPCExpansion = cp.fit_regression(orthogonalPolynomials, samples.T, data)
        
        
        # add this polynomial chaos method psudo spectral ...
        # q,w = cp.generate_quadrature(order, dist, rule="g") # 'g' just for independent; dependent use 'c'
        # gPCExpansion = cp.fit_quadrature(orth, nodes, weights, solves, retall, norms)
        
        
        # statistics
        
        statsDict = {}
                     
        statsDict['expectedValue']  = cp.E(gPCExpansion, distributionManager.jointDistribution)
        statsDict['variance']       = cp.Var(gPCExpansion, distributionManager.jointDistribution)
        statsDict['standardDeviation']   = np.sqrt(statsDict['variance'])
        conficenceInterval  = cp.Perc(gPCExpansion, [confidenceAlpha/2., 100-confidenceAlpha/2.], distributionManager.jointDistribution)
        statsDict['conficenceInterval']  =  conficenceInterval.reshape(2,len(np.atleast_1d(statsDict['expectedValue'])))
        statsDict['confidenceAlpha'] = confidenceAlpha
        
        # conditional expected values  and sensitivity coefficients
        distributionDimension = len(distributionManager.jointDistribution)
        if distributionDimension > 1:
            # test dependecy or not
            if dependentCase == False:
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
        
        # statistics
        uqsaMeasures = UqsaMeasures()
        uqsaMeasures.setVariablesDict(statsDict)
        return uqsaMeasures
         
         
class UqsaMethodPolynomialChaosDepDir(TestBaseClass):
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    #----External Variables -------------------------------------------------------------#
    externVariables = {  'polynomialOrder' : TestBaseClass.ExtValue(int), 
                         'sampleFactor'    : TestBaseClass.ExtValue(int),
                         }
                
    externXmlAttributes  = []
    
    externXmlElements    = ['polynomialOrder', 
                            'sampleFactor']
    
    def __init__(self):  
        self.sampleFactor = 2
        #polynomialOrders of the polynomial chaos expansion || if more then one given they are processed consecutevely
        self.polynomialOrder = 3
        
    def evaluateSamplesSize(self, distributionDimension):
        '''
        function to evaluate the sample size
        '''
        # calculate samplesSize from expansion order 
        samplesSize = int(self.sampleFactor*cp.terms(self.polynomialOrder, distributionDimension))
        abcSample = False
        return samplesSize,abcSample
              
    def calculateOrthogonalPolynomials(self,distributionManager):
        '''
        Method to calculate orthogonal polynomials
        '''
        
        return cp.orth_ttr(self.polynomialOrder,distributionManager.jointDistribution)
    
    def calculateStatistics(self, distributionManager, sampleManager, qoi):
        
        
        sampleSize,abcSample = self.evaluateSamplesSize(distributionManager.distributionDimension)
        
        samples,samplesDependent = sampleManager.getSampleMatrices(sampleSize,abcSample)
        data                     = qoi.getData(sampleSize,abcSample)
        
        dependentCase = sampleManager.dependentCase
        confidenceAlpha = qoi.confidenceAlpha
        
        
        q = samples.T
        
        data = data.T
        
        space = True
        
        if space == True:
            nSpacePoints = len(data)
            
            X = np.linspace(0, 0.8, nSpacePoints)
            
            E,V = np.empty((2,nSpacePoints))
            
            # loop over all places x in the trajectory 
            for j,x in enumerate(X[1:-1]):
                print "calculate number {} from {} at position {}".format(j,nSpacePoints-2,x)
                i = j+1
                #z == Xvalue we are right now ... 
                trans = lambda q: \
                    [q[1], q[2],
                     (x-q[0])*(x>q[0]),
                     (q[0]-x)*(x<=q[0])]
                    
                dist = cp.Dist(_length=4)
                dist._mom = cp.momgen(100, distributionManager.jointDistribution, trans=trans, rule="C",
                        composit=[x,.5,.5])
        
                orth = cp.orth_chol(self.polynomialOrder, dist, normed=0)
                
                #y = np.array(map(solver, q.T))
                y = data[i]
                
                approx = cp.fit_regression(orth, trans(q), y,
                        rule="T", order=1)
        
                # save exp and var for this x value
                E[i]   = cp.E(approx, dist)
                V[i]   = cp.Var(approx, dist)
            
        else:
            # TODO!!
            nSpacePoints = len(data)
            
            X = np.linspace(0, 0.8, nSpacePoints)
            
            E,V = np.empty((2,nSpacePoints))
            
            # loop over all places x in the trajectory 
            for j,x in enumerate(X[1:-1]):
                print "calculate number {} from {} at position {}".format(j,nSpacePoints-2,x)
                i = j+1
                #z == Xvalue we are right now ... 
                trans = lambda q: \
                    [q[1], q[2],
                     (x-q[0])*(x>q[0]),
                     (q[0]-x)*(x<=q[0])]
                    
                dist = cp.Dist(_length=4)
                dist._mom = cp.momgen(100, distributionManager.jointDistribution, trans=trans, rule="C",
                        composit=[x,.5,.5])
        
                orth = cp.orth_chol(self.polynomialOrder, dist, normed=0)
                
                #y = np.array(map(solver, q.T))
                y = data[i]
                
                approx = cp.fit_regression(orth, trans(q), y,
                        rule="T", order=1)
        
                # save exp and var for this x value
                E[i]   = cp.E(approx, dist)
                V[i]   = cp.Var(approx, dist)
            
        # collect data for return
        statsDict = {}
        statsDict['expectedValue']      = E
        statsDict['variance']           = V
        
        # statistics
        uqsaMeasures = UqsaMeasures()
        uqsaMeasures.setVariablesDict(statsDict)
        return uqsaMeasures
  
class UqsaMethodMonteCarlo(TestBaseClass):
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    #----External Variables -------------------------------------------------------------#
    externVariables = {  'sensitivityAnalysis'   : TestBaseClass.ExtValue(bool), 
                         'sampleSize'            : TestBaseClass.ExtValue(int)
                         }
                
    externXmlAttributes  = []
    
    externXmlElements    = ['sampleSize', 
                            'sensitivityAnalysis']
    
    def __init__(self):
        '''
        
        '''
        self.sampleSize = 10
        #
        self.sensitivityAnalysis = False
        
         
    def evaluateSamplesSize(self, distributionDimension):
        '''
        function to evaluate the sample size
        '''
        # calculate samplesSize from expansion order 
        abcSample = self.sensitivityAnalysis
        return self.sampleSize,abcSample
         
    def calculateStatistics(self, distributionManager, sampleManager, qoi):
        '''
        Function which calculates the gPCExpansion for the given data
        '''
                
        sampleSize,abcSample = self.evaluateSamplesSize(distributionManager.distributionDimension)
        
        samples,samplesDependent = sampleManager.getSampleMatrices(sampleSize,abcSample)
        data                     = qoi.getData(sampleSize,abcSample)
                
        dependentCase = sampleManager.dependentCase
        confidenceAlpha = qoi.confidenceAlpha
                
        statsDict = {}
        if dependentCase == False:
            
            
            # check if samples, and data have the right format
            distDim = distributionManager.distributionDimension
            if len(samples) != (distDim+2)*self.sampleSize:
                print 'WARNING: Not simulated as sensitivity analysis case as len(samples) != (distDim+2)*self.sampleSize'
                self.sensitivityAnalysis = False
                
            if self.sensitivityAnalysis == False:
            
                statsDict['expectedValue']       = np.sum(data,axis=0)/len(samples)
                # variance = (data**2-mean)/len(samples)
                statsDict['variance']            = np.sum(data**2-statsDict['expectedValue']**2,axis=0)/len(samples)
                
                quantiles = [confidenceAlpha/2, 100.-confidenceAlpha/2]
                statsDict['conficenceInterval']  = np.percentile(data,quantiles, axis=0)
                statsDict['confidenceAlpha'] = confidenceAlpha
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
                statsDict['confidenceAlpha'] = confidenceAlpha
                
                statsDict['conditionalExpectedValue'] = None
                statsDict['conditionalVariance']      = np.array(conditionalVarianceGivenQ)
                
                # sensitivity indices
                statsDict['firstOrderSensitivities'] = np.array(Si)
                statsDict['totalSensitivities']      = np.array(STi)
        
        else:
            # TODO: implement dependent dist methods for MC
            raise NotImplementedError("NO MC methods for dependet distributions implemented")
        
        
        # statistics
        uqsaMeasures = UqsaMeasures()
        uqsaMeasures.setVariablesDict(statsDict)
        return uqsaMeasures
            