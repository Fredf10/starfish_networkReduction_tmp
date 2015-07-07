
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
        
    def calculateStatisticsPolynomialChaos(self, orthogonalPolynomials, samples, distributions, dependentCase):
        '''
        Function which calculates the gPCExpansion for the given data
        '''
        assert self.data != None, 'QuantityOfInterest {} : {} cannot calculate gPCExpansion and statistics as there is no data defined'.format(self.locationName,self.quantityName)
                   
        print "    starting the polychaos polynomial calculation from polychaos simulation result!!"
        print self.locationName,'--',self.quantityName
        
        # polynomial chaos expansion
        self.gPCExpansion = cp.fit_regression(orthogonalPolynomials, samples.T, self.data) 
                 
        # statistics
        self.expectedValue       = cp.E(self.gPCExpansion, distributions)
        self.variance            = cp.Var(self.gPCExpansion, distributions)
        self.conficenceInterval  = cp.Perc(self.gPCExpansion, [self.confidenceAlpha/2., 100-self.confidenceAlpha/2.], distributions)
        self.conficenceInterval =  self.conficenceInterval.reshape(2,len(self.expectedValue))
        
        # conditional expected values  and sensitivity coefficients
        distributionDimension = len(distributions)
        if distributionDimension > 1:
            # test dependecy or not
            if dependentCase == False:
                # independent case: use analytic expression from polynomial chaos expansion
                self.conditionalExpectedValue = []
                self.conditionalVariance      = []
                # conditional mean and variance
                for rvIndex in xrange(distributionDimension):
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
                # dependent: use saltellies monte carlo method with polynomial chaos as function for evaluations
                self.firstOrderSensitivities = cp.Sens_m_sample(self.gPCExpansion,distributions,int(1e4))
                self.totalSensitivities      = cp.Sens_t_sample(self.gPCExpansion,distributions,int(1e4))
                #samplesSalelli      
                
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
        