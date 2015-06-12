


import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))
sys.path.append(cur+'/../')


from classLocationOfInterest import LocationOfInterest
import moduleFilePathHandlerVPC as mFPH_VPC

#sys.path.append(''.join([cur,'/../UtilityLib']))
from UtilityLib import moduleXML

import numpy as np
import h5py

class LocationOfInterestManager(object):
    '''
    
    '''
    def __init__(self, sampleSize):
        
        self.locationOfInterests = []
        self.sampleSize = sampleSize
        
    def addLocationOfInterest(self,locationName, quantitiesOfInterestToProcess, xVals, confidenceAlpha):
        
        self.locationOfInterests.append(LocationOfInterest(locationName,quantitiesOfInterestToProcess, xVals, confidenceAlpha))    
            
    def loadQuantitiyOfInterestData(self):
        '''
        
        '''
        print "loadQuantitiyOfInterestData() not implemented yet"
        
    def saveQuantitiyOfInterestData(self,networkName,dataNumber,gPCEmethod,gPCEorder):
        '''
        
        '''
        vpcQuantityOfInterestFile = mFPH_VPC.getFilePath('vpcSolutionDataFile', networkName, dataNumber, mode = "write", gPCEmethod=gPCEmethod, gPCEorder=gPCEorder)
        saveFile = h5py.File(vpcQuantityOfInterestFile,'w')
        for locationOfInterest in self.locationOfInterests:        
            # add information of each quantity in each location
            locationOfInterest.saveQuantitiyOfInterestData(saveFile)
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
                locationOfInterest.preprocessSolutionData(vascularNetwork,self.simulationTime, self.sampleSize, sampleIndex)
                        
    def calculateStatisticsPolynomialChaos(self,orthogonalPolynomials, samples, distributions):
        '''
        Calculate statisitcs for each quantity of interest at each location of interest based
        on the generalized polynomial chaos expansion.
        '''
        for locationOfInterest in self.locationOfInterests:
            for quantity in locationOfInterest.quantitiesOfInterestToProcess:
                locationOfInterest.quantitiesOfInterest[quantity].calculateStatisticsPolynomialChaos(orthogonalPolynomials, samples, distributions)
                
    
    def calculateStatisticsMonteCarlo(self):
        '''
        Calculate statisitcs for each quantity of interest at each location of interest based
        on a Monte Carlo simulation
        '''
        pass
        for locationOfInterest in self.locationOfInterests:
            for quantity in locationOfInterest.quantitiesOfInterestToProcess:
                locationOfInterest.quantitiesOfInterest[quantity].calculateStatisticsMonteCarlo()
                        
