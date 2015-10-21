


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
                           } 
    externXmlAttributes  = []
    externXmlElements    = ['locationsOfInterest']
    
    def __init__(self):
        
        self.locationsOfInterest = {}
        self.sampleSize = None
                                
        self.testDict = {}
            
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
                        
    def loadQuantitiyOfInterestData(self):
        '''
        
        '''
        print "loadQuantitiyOfInterestData() not implemented yet"
        
    def saveQuantitiyOfInterestData(self,vpcQuantityOfInterestFile):
        '''
        
        '''
        saveFile = h5py.File(vpcQuantityOfInterestFile,'w')
        # add simulation time
        saveFile.create_dataset('simulationTime', data=self.simulationTime)
        for locationId,locationOfInterest in self.locationsOfInterest.iteritems():        
            # add information of each quantity in each location
            locationGroup = saveFile.create_group(locationId)
            locationOfInterest.saveQuantitiyOfInterestData(locationGroup)
        saveFile.flush()
        saveFile.close()
            
    def preprocessSolutionData(self,evaluationCaseFiles):
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
        
    def getQoiIterator(self):
        '''
        Function that creates and iterator element which iterates through all stored qoi data elements
        '''
        qoiList = []
        for locationOfInterest in self.locationsOfInterest.values():
            for qoi in locationOfInterest.quantitiesOfInterest.itervalues():
                qoiList.append(qoi)
        return qoiList
    
    