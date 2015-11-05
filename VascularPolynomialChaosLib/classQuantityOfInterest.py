

from testBaseClass import TestBaseClass 

import numpy as np

class QuantityOfInterest(TestBaseClass):
    '''
    
    '''
    def __init__(self,quantityName, locationName, confidenceAlpha, data = None):
        
        self.queryLocation = locationName
        self.quantityName = quantityName
        self.confidenceAlpha = confidenceAlpha
        self.data = data
        # for large MC cases
        self.data2 = None
        self.data3 = None
        self.data4 = None
        # trajectory object
        self.trajectoryData = None
        
        self.uqsaMeasures = None
    
    def getData(self, sampleSize, abcSample):
        '''
        checks if data defined and returns it
        '''
        assert self.data != None, 'QuantityOfInterest {} : {} cannot calculate gPCExpansion and statistics as there is no data defined'.format(self.queryLocation,self.quantityName)
        
        if abcSample == False:
            return self.data[:sampleSize]
        else:
            #TODO: implement abc sampling hash if monte carlo is used
            raise NotImplementedError("MC sensitivity, abc-sample hash case not implemented yet, exit")
        
    def addUqsaMeasures(self,uqsaMeasureName,uqsaMeasure):
        '''
        Method append a uqsaMeasure class to the dictionary: self.uqsaMeasures
        '''
        
        if self.uqsaMeasures == None:
            self.uqsaMeasures = {}
        if uqsaMeasureName not in self.uqsaMeasures:
            self.uqsaMeasures[uqsaMeasureName] = uqsaMeasure
        else:
            print " did not add {} as it already exist in dict".format(uqsaMeasureName)
    
    
    # --- move to base class
    def saveDataHdf5(self, hdf5SaveGroup):
        '''
        Method to save data and statistics to file
        '''
        # pure variables
        
        # in class definition
        variablesToSave = ["data"]
        # functionality
        for variableName in variablesToSave:
            variableValue = self.getVariable(variableName)
            if variableValue != None: 
                hdf5SaveGroup.create_dataset(variableName, data=variableValue)
        
        # dict of objects
        
        dictsToSave = ["uqsaMeasures"]
        # functionality
        for dictName in dictsToSave:
            dictValue = self.getVariable(dictName)
            dictGroup = hdf5SaveGroup.create_group(dictName)
            if dictValue != None or dictValue == {}:
                # if dict is contains objects
                for dictObjectName,dictObject in dictValue.iteritems():
                    dictEntryGroup = dictGroup.create_group(dictObjectName)
                    dictObject.saveDataHdf5(dictEntryGroup)
        
    def loadDataHdf5(self, hdf5SaveGroup):
        '''
        method to load data from hdf5
        '''
        loadedData = {}
        ## pure variables
        # in class definition
        variablesToLoad = ['data']
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
        
        ## dictionary to load
        dictsToLoad = ["uqsaMeasures"]
        
        # functionality
        for dictName in dictsToLoad:
            # check if dictionary not none and not empty
            dictL = self.getVariable(dictName) 
            if dictL != {} and dictL != None:
                # serach for dict group node for the dictionary
                if dictName in hdf5SaveGroup.keys():
                    dictGroup = hdf5SaveGroup[dictName]
                    for dictObjectName,dictObject in dictL.iteritems():
                        if dictObjectName in dictGroup.keys():
                            dictObjectGroup = dictGroup[dictObjectName]
                            dictObject.loadDataHdf5(dictObjectGroup)
                            
        self.setVariablesDict(loadedData)       
        