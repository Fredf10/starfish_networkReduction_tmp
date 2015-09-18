
import chaospy as cp
import numpy as np

class QuantityOfInterest(object):
    '''
    
    '''
    def __init__(self,quantityName, locationName, confidenceAlpha, data = None):
        
        self.queryLocation = locationName
        self.quantityName = quantityName
        self.data = data
        
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
        
        # use saltellis MC method evaluating the gpc-expansion
        self.firstOrderSensitivitiesMC = None
        self.totalSensitivitiesMC = None
        
    def calculateStatisticsPolynomialChaos(self,distributionManager):
        '''
        Function which calculates the gPCExpansion for the given data
        '''
        assert self.data != None, 'QuantityOfInterest {} : {} cannot calculate gPCExpansion and statistics as there is no data defined'.format(self.queryLocation,self.quantityName)
                   
        print "    starting the polychaos polynomial calculation from polychaos simulation result!!"
        print self.queryLocation,'--',self.quantityName
                
        # polynomial chaos expansion
        self.gPCExpansion = cp.fit_regression(distributionManager.orthogonalPolynomials, distributionManager.samples.T, self.data)
                 
        # statistics
        self.expectedValue       = cp.E(self.gPCExpansion, distributionManager.jointDistribution)
        self.variance            = cp.Var(self.gPCExpansion, distributionManager.jointDistribution)
        self.conficenceInterval  = cp.Perc(self.gPCExpansion, [self.confidenceAlpha/2., 100-self.confidenceAlpha/2.], distributionManager.jointDistribution)
        self.conficenceInterval =  self.conficenceInterval.reshape(2,len(np.atleast_1d(self.expectedValue)))
        
        # conditional expected values  and sensitivity coefficients
        distributionDimension = len(distributionManager.jointDistribution)
        if distributionDimension > 1:
            # test dependecy or not
            if distributionManager.dependentCase == False:
                # independent case: use analytic expression from polynomial chaos expansion
                self.conditionalExpectedValue = []
                self.conditionalVariance      = []
                # conditional mean and variance
                for rvIndex in xrange(distributionDimension):
                    currDistMean = cp.E(distributionManager.jointDistribution)
                    currDistMean[rvIndex] = np.nan
                    # reduce polynomials
                    currPolynomTime = self.gPCExpansion(*currDistMean)
                    self.conditionalExpectedValue.append(cp.E(currPolynomTime,distributionManager.jointDistribution))
                    self.conditionalVariance.append(cp.Var(currPolynomTime,distributionManager.jointDistribution))    
            
                # sensitivity indices
                self.firstOrderSensitivities = cp.Sens_m(self.gPCExpansion,distributionManager.jointDistribution)
                self.totalSensitivities      = cp.Sens_t(self.gPCExpansion,distributionManager.jointDistribution)
            else:
                # dependent rancom variables
                
                # this method is broken
                #sensindices = cp.Sens_nataf(distributionManager.expansionOrder, distributionManager.jointDistributionDependent, distributionManager.samplesDependent.T, self.data)
                #
                
                self.firstOrderSensitivities = cp.Sens_m_nataf(distributionManager.expansionOrder,
                                                               distributionManager.jointDistributionDependent,
                                                               distributionManager.samplesDependent.T,
                                                               self.data)
                
                self.totalSensitivities      = cp.Sens_t_nataf(distributionManager.expansionOrder,
                                                               distributionManager.jointDistributionDependent,
                                                               distributionManager.samplesDependent.T,
                                                               self.data)
            
    def calculateStatisticsMonteCarlo(self):
        '''
        Function to calculate statistics
        '''
        pass
                
    def saveQuantitiyOfInterestData(self, quantitiyGroup):
        '''
        Method to save data and statistics to file
        '''
        quantitiyGroup.create_dataset("data", data=self.data)
        #quantitiyGroup.create_dataset("gPCExpansion", data=self.gPCExpansion)
        variablesToSave = ["expectedValue",
                           "variance",
                           "conficenceInterval",
                           "conditionalExpectedValue",
                           "conditionalVariance",
                           "firstOrderSensitivities",
                           "totalSensitivities"]
        
        for variableName in variablesToSave:
            variableValue = self.getVariableValue(variableName)
            if variableValue != None: 
                quantitiyGroup.create_dataset(variableName, data=variableValue)
                    
    def getVariableValue(self, variableName):
        '''
        Returns value of variable with name : variableName
        States Error if not such variable
        '''
        try:
            return self.__getattribute__(variableName)
        except: 
            print "ERROR QuantityOfInterest.getVariable() : QuantityOfInterest has no variable {}".format(variableName)
        
    def update(self, data):
        '''
        updates the QuantityOfInterest data using a dictionary in form of 
        QuantityOfInterest = {'variableName': value}
        '''
        for key, value in data.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key, value)
            except: 
                print 'WARNING QuantityOfInterest.update(): wrong key: %s, could not update QuantityOfInterest' % key   
        