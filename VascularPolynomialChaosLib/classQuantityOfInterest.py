
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
        self.expectedValue       = None
        self.variance            = None
        self.standardDeviation   = None 
        
        self.conditionalExpectedValue = None
        self.conditionalVariance      = None
        
        self.conficenceInterval  = None
        self.confidenceAlpha     = confidenceAlpha
        
        self.firstOrderSensitivities = None
        self.totalSensitivities      = None       
    
    def getData(self):
        '''
        checks if data defined and returns it
        '''
        assert self.data != None, 'QuantityOfInterest {} : {} cannot calculate gPCExpansion and statistics as there is no data defined'.format(self.queryLocation,self.quantityName)
        return self.data
    
    def saveQuantitiyOfInterestData(self, quantitiyGroup):
        '''
        Method to save data and statistics to file
        '''
        quantitiyGroup.create_dataset("data", data=self.data)
        #quantitiyGroup.create_dataset("gPCExpansion", data=self.gPCExpansion)
        variablesToSave = ["expectedValue",
                           "variance",
                           "standardDeviation",
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
        