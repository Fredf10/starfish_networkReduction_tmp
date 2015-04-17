
import chaospy as cp
import numpy as np

class QuantityOfInterest(object):
    '''
    
    '''
    def __init__(self,quantityName, locationName, confidenceAlpha):
        
        self.locationName = locationName
        self.quantityName = quantityName
        self.data = None
        
        # statistic properties
        self.gPCExpansion = None
        
        self.expectedValue       = None
        self.variance            = None
        
        self.conditionalExpectedValue = None
        self.conditionalVariance      = None
        
        self.conficenceInterval  = None
        self.confidenceAlpha     = confidenceAlpha
        
        self.firstOrderSensitivities = None
        self.totalSensitivities      = None
        
    def calculateStatisticsPolynomialChaos(self, orthogonalPolynomials, samples, distributions):
        '''
        Function which calculates the gPCExpansion for the given data
        '''
        if self.data != None:
            print "    starting the polychaos polynomial calculation from polychaos simulation result!!"
            print self.locationName,'--',self.quantityName
            
            # polynomial chaos expansion
            self.gPCExpansion = cp.fit_regression(orthogonalPolynomials, samples.T, self.data) 
                     
            # statistics
            self.expectedValue       = cp.E(self.gPCExpansion, distributions)
            self.variance            = cp.Var(self.gPCExpansion, distributions)
            self.conficenceInterval  = cp.Perc(self.gPCExpansion, [self.confidenceAlpha/2., 100-self.confidenceAlpha/2.], distributions)
            
            distributionDimension = len(distributions)
            if distributionDimension > 1:
                self.conditionalExpectedValue = []
                self.conditionalVariance      = []
                # conditional mean and variance
                for rvIndex in xrange():
                    currDistMean = cp.E(distributions)
                    currDistMean[rvIndex] = np.nan
                    # reduce polynomials
                    currPolynomTime = self.gPCExpansion(*currDistMean)
                    self.conditionalExpectedValue.append(cp.E(currPolynomTime,distributions))
                    self.conditionalVariance.append(cp.Var(currPolynomTime,distributions))    
            
                # sensitivity indices
                self.firstOrderSensitivities = cp.Sens_m(self.gPCExpansion,distributions)
                self.totalSensitivities      = cp.Sens_t(self.gPCExpansion,distributions)
                         
        else:
            raise ValueError('QuantityOfInterest {}-{} cannot calculate gPCExpansion and statistics as there is no data defined'.format(self.locationName,self.quantityName))
        
    def calculateStatisticsMonteCarlo(self):
        '''
        Function to calculate statistics
        '''
        pass
        
        