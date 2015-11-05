


import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))


from classLocationOfInterest import LocationOfInterest

#sys.path.append(''.join([cur,'/../UtilityLib']))
#import moduleXML as mXML

from UtilityLib import moduleXML


import numpy as np
import h5py

from testBaseClass import TestBaseClass 

class LocationOfInterestManager(TestBaseClass):
    '''
    
    '''
    externVariables      = {'locationsOfInterest' : TestBaseClass.ExtDict('locationOfInterest',TestBaseClass.ExtObject({'LocationOfInterest':LocationOfInterest})),
                            'evaluateSimulationTime': TestBaseClass.ExtValue(bool)
                           } 
    externXmlAttributes  = []
    externXmlElements    = ['evaluateSimulationTime',
                            'locationsOfInterest']
    
    def __init__(self):
        
        self.locationsOfInterest = {}
        self.sampleSize = None
            
        self.evaluateSimulationTime = False
        self.simulationTime = None
            
    def initialize(self):
        '''
        Function to initialize all locations of interest which creates quantity of intertes objects
        '''
        for locationOfInterest in self.locationsOfInterest.itervalues():
            locationOfInterest.initialize()
                
    def addLocationOfInterest(self,locationId, locationName, quantitiesOfInterestToProcess, xVal, confidenceAlpha):
        '''
        
        '''
        self.locationsOfInterest[locationId] = LocationOfInterest(locationName,quantitiesOfInterestToProcess, xVal, confidenceAlpha)
                        
    def loadQuantitiyOfInterestData(self,vpcQuantityOfInterestFile, simulationTimeOnly = False):
        '''
        
        '''
        print "loadQuantitiyOfInterestData() not implemented yet"
        
        hdf5SaveGroup = h5py.File(vpcQuantityOfInterestFile,'r')
        
        #TODO: add subgroup
        
        loadedData = {}
        
        ## pure variables
        # in class definition
        variablesToLoad = ["simulationTime"]
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
        
        if simulationTimeOnly == False:
            
            ## dictionary to load
            dictsToLoad = ["locationsOfInterest"]
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
        hdf5SaveGroup.close()
        
    def saveQuantitiyOfInterestData(self,vpcQuantityOfInterestFile, simulationTimeOnly = False):
        '''
        
        '''
        hdf5SaveGroup = h5py.File(vpcQuantityOfInterestFile,'w')
        
        # pure variables
        
        # in class definition
        variablesToSave = ["simulationTime"]
        # functionality
        for variableName in variablesToSave:
            variableValue = self.getVariable(variableName)
            if variableValue != None: 
                hdf5SaveGroup.create_dataset(variableName, data=variableValue)
        
        if simulationTimeOnly == False:
            # dict of objects
            
            dictsToSave = ["locationsOfInterest"]
            # functionality
            for dictName in dictsToSave:
                dictValue = self.getVariable(dictName)
                dictGroup = hdf5SaveGroup.create_group(dictName)
                if dictValue != None or dictValue == {}:
                    # if dict is contains objects
                    for dictObjectName,dictObject in dictValue.iteritems():
                        dictEntryGroup = dictGroup.create_group(dictObjectName)
                        dictObject.saveDataHdf5(dictEntryGroup)
                    
        
        hdf5SaveGroup.flush()
        hdf5SaveGroup.close()
            
    def preprocessSolutionData(self,
                               evaluationCaseFiles,
                               uqsaSolutionDataFileSave,
                               simulationTimeFileSave,
                               simulationTimeaFileLoad):
        '''
        load all samples and pass data to locations of interest
        
        find simulation time array with largest dt
        
        invoke time cropping and post processing for the data        
        '''
        # loop through all solution data files 
        # find time array
        #solutionTime min
        
        numberOfTimePoints = []
        timeStart = []
        timeEnd = []
        
        if self.evaluateSimulationTime == True:
            print "estimate simulation time of all simulations"
            for batchData in evaluationCaseFiles:
                networkName              = batchData['networkName']
                dataNumber               = batchData['dataNumber']
                networkXmlFileLoad       = batchData['networkXmlFileLoad']
                pathSolutionDataFilename = batchData['pathSolutionDataFilename']
                
                vascularNetwork = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = networkXmlFileLoad, pathSolutionDataFilename = pathSolutionDataFilename)
                vascularNetwork.linkSolutionData()
                
                numberOfTimePoints.append(len(vascularNetwork.simulationTime))
                timeStart.append(min(vascularNetwork.simulationTime))
                timeEnd.append(max(vascularNetwork.simulationTime))
                            
                vascularNetwork.solutionDataFile.close()
                del vascularNetwork
            
            self.simulationTime = np.linspace(max(timeStart), min(timeEnd), min(numberOfTimePoints))
            # save the simulationTime
            self.saveQuantitiyOfInterestData(simulationTimeFileSave, True) 
            
        else:
            print "load estimated simulation time of all simulations"
            if simulationTimeaFileLoad != None:
                self.loadQuantitiyOfInterestData(simulationTimeaFileLoad)
            else:
                raise ValueError('simulationTime hdf5 file for case does not exist! {}'.format(simulationTimeaFileLoad))
        
        print "load data for quantities of interest"
        # pass the data to the locationsOfInterests which will load the information needed
        for batchData in evaluationCaseFiles:
            simulationIndex          = batchData['simulationIndex']
            networkName              = batchData['networkName']
            dataNumber               = batchData['dataNumber']
            networkXmlFileLoad       = batchData['networkXmlFileLoad']
            pathSolutionDataFilename = batchData['pathSolutionDataFilename']
            
            vascularNetwork = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = networkXmlFileLoad, pathSolutionDataFilename = pathSolutionDataFilename)
            vascularNetwork.linkSolutionData()
            for locationOfInterest in self.locationsOfInterest.values():
                locationOfInterest.preprocessSolutionData(vascularNetwork,self.simulationTime, self.sampleSize, simulationIndex)
        
        # second postprocessing find extrema if needed also for variables defined over space
        for locationOfInterest in self.locationsOfInterest.values():
            locationOfInterest.preprocessSolutionDataExtremaAndInflectionPoints(self.simulationTime, self.sampleSize)
    
        for locationOfInterest in self.locationsOfInterest.values():
            locationOfInterest.preprocessSolutionDataTrajectory(self.simulationTime, self.sampleSize)
        
        self.saveQuantitiyOfInterestData(uqsaSolutionDataFileSave)
        
    def getQoiIterator(self):
        '''
        Function that creates and iterator element which iterates through all stored qoi data elements
        '''
        qoiList = []
        for locationOfInterest in self.locationsOfInterest.values():
            for qoi in locationOfInterest.quantitiesOfInterest.itervalues():
                qoiList.append(qoi)
        return qoiList
    
    