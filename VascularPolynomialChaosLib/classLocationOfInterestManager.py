


import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))

from classLocationOfInterest import LocationOfInterest

sys.path.append(''.join([cur,'/../UtilityLib']))
import moduleXML
import numpy as np
class LocationOfInterestManager(object):
    '''
    
    '''
    def __init__(self, sampleSize):
        
        self.locationOfInterests = []
        self.sampleSize = sampleSize
        
    def addLocationOfInterest(self,locationName, quantitiesOfInterestToProcess, xVals):
        
        self.locationOfInterests.append(LocationOfInterest(locationName,quantitiesOfInterestToProcess, xVals))    
            
    def loadQuantitiyOfInterestData(self):
        '''
        
        '''
        print "loadQuantitiyOfInterestData() not implemented yet"
        
    def saveQuantitiyOfInterestData(self):
        '''
        
        '''
        print "saveQuantitiyOfInterestData() not implemented yet" 
            
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
        
        for networkName,dataNumber,vpcNetworkXmlEvaluationFile,vpcEvaluationSolutionDataFile in evaluationCaseFiles:
            vascularNetwork = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = vpcNetworkXmlEvaluationFile, pathSolutionDataFilename = vpcEvaluationSolutionDataFile)
            vascularNetwork.linkSolutionData()
            
            numberOfTimePoints.append(len(vascularNetwork.simulationTime))
            timeStart.append(min(vascularNetwork.simulationTime))
            timeEnd.append(max(vascularNetwork.simulationTime))
            
            vascularNetwork.solutionDataFile.close()
            del vascularNetwork
        
        self.simulationTime = np.linspace(max(timeStart), min(timeEnd), min(numberOfTimePoints))
        
        # pass the data to the locationsOfInterests which will load the information needed
        for sampleIndex,[networkName,dataNumber,vpcNetworkXmlEvaluationFile,vpcEvaluationSolutionDataFile] in enumerate(evaluationCaseFiles):
            vascularNetwork = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = vpcNetworkXmlEvaluationFile, pathSolutionDataFilename = vpcEvaluationSolutionDataFile)
            vascularNetwork.linkSolutionData()
            for locationOfInterest in self.locationOfInterests:
                locationOfInterest.loadSolutionData(vascularNetwork,self.simulationTime, self.sampleSize, sampleIndex)
                        
    def calculatePolynomialChaosExpansions(self,orthogonalPolynomials, samples):
        '''
        
        '''
        for locationOfInterest in self.locationOfInterests:
            for quantity in locationOfInterest.quantitiesOfInterestToProcess:
                locationOfInterest.quantitiesOfInterest[quantity].calculatePolynomialChaosExpansions(orthogonalPolynomials, samples)
                
    
    def calculateStatistics(self):
        '''
        
        '''
        pass
      
                        