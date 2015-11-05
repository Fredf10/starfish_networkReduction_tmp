

from testBaseClass import TestBaseClass 
import numpy as np

class UqsaMeasures(TestBaseClass):
    '''
    
    '''
    def __init__(self):
        
        # statistic properties
        self.expectedValue            = None
        self.variance                 = None
        self.standardDeviation        = None 
        
        self.conditionalExpectedValue = None
        self.conditionalVariance      = None
        
        self.conficenceInterval       = None
        self.confidenceAlpha          = None
        
        self.firstOrderSensitivities  = None
        self.totalSensitivities       = None       
        
    
    # --- move to base class
    def saveDataHdf5(self, hdf5SaveGroup):
        '''
        Method to save data and statistics to file
        '''
        
        variablesToSave = ["expectedValue",
                           "variance",
                           "standardDeviation",
                           "conficenceInterval",
                           "conditionalExpectedValue",
                           "conditionalVariance",
                           "firstOrderSensitivities",
                           "totalSensitivities"]
        
        for variableName in variablesToSave:
            variableValue = self.getVariable(variableName)
            if variableValue != None: 
                hdf5SaveGroup.create_dataset(variableName, data=variableValue)
                
    def loadDataHdf5(self, hdf5SaveGroup):
        '''
        method to load data from hdf5
        '''
        loadedData = {}
        ## pure variables
        # in class definition
        variablesToLoad = ["expectedValue",
                           "variance",
                           "standardDeviation",
                           "conficenceInterval",
                           "conditionalExpectedValue",
                           "conditionalVariance",
                           "firstOrderSensitivities",
                           "totalSensitivities"]
        # functionality
        for variableName in variablesToLoad:
            if variableName in hdf5SaveGroup.keys(): 
                variableData = (hdf5SaveGroup[variableName])
                # TODO: use externalVariable definitions dictionary ... #
                # check for shape 
                if np.shape(variableData)== ():
                    loadedData[variableName] = variableData
                else:
                    loadedData[variableName] = variableData[:]
                
        self.setVariablesDict(loadedData)