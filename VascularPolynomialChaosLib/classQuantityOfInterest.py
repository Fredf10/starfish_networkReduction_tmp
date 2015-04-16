
import chaospy as cp
import numpy as np

class QuantityOfInterest(object):
    '''
    
    '''
    def __init__(self,quantityName, locationName):
        
        self.locationName = locationName
        self.quantityName = quantityName
        self.data = None
        
        # statistic properties
        self.gPCExpansion = None
        
        self.expectedValue       = None
        self.variance            = None
        
        self.conficenceIntervals = None
        self.confidenceAlpha     = 0.8
        
        self.partialSensitivities = None
        self.totalSensitivities   = None
        
    def calculatePolynomialChaosExpansions(self, orthogonalPolynomials, samples):
        '''
        Function which calculates the gPCExpansion for the given data
        '''
        if self.data != None:
            print "    starting the polychaos polynomial calculation from polychaos simulation result!!"
            print self.locationName,'--',self.quantityName
            
            self.gPCExpansion = cp.fit_regression(orthogonalPolynomials, samples.T, self.data) 
                     
#                     for tag,data in sol.iteritems():
#                         # polynoms for the total pressure signal
#                         print "        polynoms for ",str(tag)
#                         if 'extrema' not in tag:
#                              
#                             polynomial = pc.fitter_lr(orthoPoly, sample.T, data)     
#                             polyDict[tag]= polynomial
#                              
#                         else:
#                             print sample.shape
#                              
#                             #print orthoPoly.dim
#                             #print orthoPoly.shape
#                              
#                             #print sample.ravel().shape
#                             #print len(data[0])
#                             polynomialTime = pc.fitter_lr(orthoPoly, sample.T, data[0])    
#                             polynomialAmp  = pc.fitter_lr(orthoPoly, sample.T, data[1])    
#                              
#                             extremaDict = {'Time':polynomialTime,'Amp':polynomialAmp}
#                             polyDict[tag]= extremaDict
#                              
#                     endTime = time.clock()
#                     polynomsT.append(polyDict)
                     
        
        else:
            raise ValueError('QuantityOfInterest {}-{} cannot calculate gPCExpansion as there is no data defined'.format(self.locationName,self.quantityName))
        
    def calculateStatistics(self,distributions):
        '''
        Function to calculate statistics
        '''
        pass
        
        