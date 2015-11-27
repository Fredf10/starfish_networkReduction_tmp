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
        
        sampleSize,abcSample     = self.evaluateSamplesSize(distributionManager.distributionDimension)
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
        
        statsDict['numberOfSamples'] = sampleSize
        # statistics
        uqsaMeasures = UqsaMeasures()
        uqsaMeasures.setVariablesDict(statsDict)
        return uqsaMeasures
         
         
class UqsaMethodPolynomialChaosDepDirLR(TestBaseClass):
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
        
        trajectoryBasis = qoi.hdf5Group['trajectoryBasis']
        
        peakTime = qoi.hdf5Group['dataBasisMinMaxPeak'][:]
        peakTime = peakTime.T[2][:sampleSize]
                
        dependentCase = sampleManager.dependentCase
        confidenceAlpha = qoi.confidenceAlpha
        
        samples = samples.T        
        data = data.T
        
        # TODO!!
        nTimePoints = len(data)
        E,V = np.empty((2,nTimePoints))
            
        v = cp.orth_ttr(self.polynomialOrder, distributionManager.jointDistribution)
        
        alphas = np.empty(nTimePoints)
        
        # move to jonathans namespace
        t = trajectoryBasis
        Q = distributionManager.jointDistribution
        order = self.polynomialOrder
        t_max = peakTime
        
        # high resolution for tracking
        top = cp.fit_regression(v, samples, peakTime)
        
        # loading bar
        import sys
        write = sys.stdout.write
        loadingBarElementCount = 1
        nElements = 50
        if nTimePoints < nElements:
            nElements = nTimePoints
        spaces = ''.join([' ' for i in xrange(nElements)])
        loadingBar = ''.join(['[',spaces,']'])
        write(loadingBar)
        sys.stdout.flush()
        backspacing = ''.join(['\b' for i in xrange(nElements+1)])
        write(backspacing)
                    
        ## SPLIT WITH LEFT AND RIGHT POLYNOMIALS
        for j in xrange(nTimePoints):
            t_ = t[j]
    
            samples_left = samples[:,t_>t_max]
            samples_right = samples[:,t_<=t_max]
            f_vals_left = data[j][t_>t_max]
            f_vals_right = data[j][t_<=t_max]
    
            N_left = len(f_vals_left)
            N_right = len(f_vals_right)
    
            nu = N_left*1./(N_left+N_right)
            n_left = int(nu*len(v))
            n_right = int((1-nu)*len(v))
    
            if n_left:
                v_left = cp.orth_ttr(int(order*nu)+1, Q)
                poly_left = cp.fit_regression(v_left, samples_left, f_vals_left)
            if n_right:
                v_right = cp.orth_ttr(int(order*(1-nu))+1, Q)
                poly_right = cp.fit_regression(v_right, samples_right, f_vals_right)
    
            if not n_right:
                poly = poly_left
                dist = Q
    
            elif not n_left:
                poly = poly_right
                dist = Q
    
            else:
                trans = lambda q: \
                        np.array([poly_left(*q)*(top(*q)<=t_), poly_right(*q)*(top(*q)>t_)])
                dist = cp.Dist(_length=2)
                dist._mom = cp.momgen(100, Q, trans=trans, rule="C",
                        composit=[t_, t_])
    
                orth = cp.orth_chol(1, dist, normed=0)
                poly = cp.fit_regression(orth, trans(samples), data[j],
                        rule="T", order=1, alpha=None, retall=2)[0]
            
            
            E[j] = cp.E(poly, dist)
            V[j]   = cp.Var(poly, dist)
            
            
            if divmod(j+1,nTimePoints/nElements)[0] == loadingBarElementCount:
                loadingBarElementCount = loadingBarElementCount+1
                write("#")
                sys.stdout.flush()
        
        write("\n")
            
        # collect data for return
        statsDict = {}
        statsDict['expectedValue']      = E
        statsDict['variance']           = V
        statsDict['numberOfSamples'] = sampleSize
        
        # statistics
        uqsaMeasures = UqsaMeasures()
        uqsaMeasures.setVariablesDict(statsDict)
        return uqsaMeasures
    
class UqsaMethodPolynomialChaosDepDirQR(TestBaseClass):
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
        
        
        sampleSize,abcSample = self.evaluateSamplesSize(distributionManager.distributionDimension+1)
        
        samples,samplesDependent = sampleManager.getSampleMatrices(sampleSize,abcSample)
        data                     = qoi.getData(sampleSize,abcSample)
        
        trajectoryBasis = qoi.hdf5Group['trajectoryBasis']
        
        peakTime = qoi.hdf5Group['dataBasisMinMaxPeak'][:]
        peakTime = peakTime.T[2][:sampleSize]
                
        dependentCase = sampleManager.dependentCase
        confidenceAlpha = qoi.confidenceAlpha
        
        samples = samples.T        
        data = data.T
        
        # TODO!!
        nTimePoints = len(data)
        E,V = np.empty((2,nTimePoints))
            
        v = cp.orth_ttr(self.polynomialOrder, distributionManager.jointDistribution)
        
        alphas = np.empty(nTimePoints)
        
        # move to jonathans namespace
        t = trajectoryBasis
        Q = distributionManager.jointDistribution
        order = self.polynomialOrder
        t_max = peakTime
        
        # high resolution for tracking
        top = cp.fit_regression(v, samples, peakTime)
        
        # loading bar
        import sys
        
        print "estimate alphas"
        
        write = sys.stdout.write
        loadingBarElementCount = 1
        nElements = 50
        if nTimePoints < nElements:
            nElements = nTimePoints
        spaces = ''.join([' ' for i in xrange(nElements)])
        loadingBar = ''.join(['[',spaces,']'])
        write(loadingBar)
        sys.stdout.flush()
        backspacing = ''.join(['\b' for i in xrange(nElements+1)])
        write(backspacing)
        
        
        ## FIRST TRY
        for j in xrange(nTimePoints):
             
            t_ = trajectoryBasis[j]
            #trans = lambda q: \
            #        np.array([(t_-top(*q))*(top(*q)>t_), (t_-top(*q))*(top(*q)<=t_)])
            trans = lambda q: \
                    np.array([q[0], q[1],q[2], (t_-top(*q))*(top(*q)<=t_)])
            dist = cp.Dist(_length=4)
            dist._mom = cp.momgen(100, distributionManager.jointDistribution, trans=trans, rule="C",
                    composit=[0.5,0.5,0.5,t_]) #composit=[.5, .5, t_])
     
            orth = cp.orth_chol(self.polynomialOrder, dist, normed=0)
            poly = cp.fit_regression(orth, trans(samples), data[j],
                    rule="T", order=1, alpha=None, retall=2)
            alphas[j] = poly[3]
            #mean[j] = cp.E(poly[0], dist)
     
            if divmod(j+1,nTimePoints/nElements)[0] == loadingBarElementCount:
                loadingBarElementCount = loadingBarElementCount+1
                write("#")
                sys.stdout.flush()
        
        write("\n")
     
        alpha = np.median(alphas)
        print "alpha", alphas
        
        
        print "calculate gpce expansion"
        
        write = sys.stdout.write
        loadingBarElementCount = 1
        nElements = 50
        if nTimePoints < nElements:
            nElements = nTimePoints
        spaces = ''.join([' ' for i in xrange(nElements)])
        loadingBar = ''.join(['[',spaces,']'])
        write(loadingBar)
        sys.stdout.flush()
        backspacing = ''.join(['\b' for i in xrange(nElements+1)])
        write(backspacing)
        
        for j in xrange(nTimePoints):
     
            t_ = trajectoryBasis[j]
            
            trans = lambda q: \
                    np.array([q[0], q[1],q[2], (t_-top(*q))*(top(*q)<=t_)])
            dist = cp.Dist(_length=4)
            dist._mom = cp.momgen(100, distributionManager.jointDistribution, trans=trans, rule="C",
                    composit=[0.5,0.5,0.5,t_])
            
            orth = cp.orth_chol(self.polynomialOrder, dist, normed=0)
            poly = cp.fit_regression(orth, trans(samples), data[j],
                    rule="T", order=1, alpha=alpha)
            E[j] = cp.E(poly, dist)
            V[j]   = cp.Var(poly, dist)
            
            if divmod(j+1,nTimePoints/nElements)[0] == loadingBarElementCount:
                loadingBarElementCount = loadingBarElementCount+1
                write("#")
                sys.stdout.flush()
        
        write("\n")
            
        # collect data for return
        statsDict = {}
        statsDict['expectedValue']      = E
        statsDict['variance']           = V
        statsDict['numberOfSamples'] = sampleSize
        
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
                #print 'WARNING: Not simulated as sensitivity analysis case as len(samples) != (distDim+2)*self.sampleSize'
                self.sensitivityAnalysis = False
                
            if self.sensitivityAnalysis == False:
            
                statsDict['expectedValue']       = np.sum(data,axis=0)/len(samples)
                # variance = (data**2-mean)/len(samples)
                statsDict['variance']            = np.sum(data**2-statsDict['expectedValue']**2,axis=0)/len(samples)
                
                quantiles = [confidenceAlpha/2, 100.-confidenceAlpha/2]
                statsDict['conficenceInterval']  = np.percentile(data,quantiles, axis=0)
                statsDict['confidenceAlpha'] = confidenceAlpha
            else:
                print "Monte Carlo Sensitivity Analysis"  
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
        statsDict['numberOfSamples'] = sampleSize
        uqsaMeasures.setVariablesDict(statsDict)
        return uqsaMeasures
            