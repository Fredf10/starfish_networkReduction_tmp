

from classQuantityOfInterest import QuantityOfInterest

from copy import copy as copy

import numpy as np

import UtilityLib.processing as mProc

import matplotlib.pyplot as plt

class LocationOfInterest(object):
    '''
    
    '''
    def __init__(self,locationName, quantitiesOfInterestToProcess, xval, confidenceAlpha):
        
        self.locationName = locationName
        self.xval = xval #position in x 
        self.quantitiesOfInterestToProcess = quantitiesOfInterestToProcess
        self.quantitiesOfInterest = {}
        self.confidenceAlpha = confidenceAlpha
                
        # check if vessel and add missing qoi to the list if needed
        
        quantitiesOfInterestToCreate = copy(quantitiesOfInterestToProcess)
        
        for quantity in quantitiesOfInterestToProcess:
            if "ward" in quantity:
                quantitiesOfInterestToCreate.extend(['Flow','Pressure','Area','WaveSpeed'])
                quantitiesOfInterestToCreate = list(set(quantitiesOfInterestToCreate))
                
            if 'Extrema' in quantity:
                quantityNamePure = quantity.split('Extrema')[-1]
                quantitiesOfInterestToCreate.extend([quantityNamePure])
                quantitiesOfInterestToCreate = list(set(quantitiesOfInterestToCreate))
                
            if 'InflectionPoint' in quantity:
                quantityNamePure = quantity.split('InflectionPoint')[-1]
                quantitiesOfInterestToCreate.extend([quantityNamePure])
                quantitiesOfInterestToCreate = list(set(quantitiesOfInterestToCreate))
        
        for quantity in quantitiesOfInterestToCreate:
            self.quantitiesOfInterest[quantity] = QuantityOfInterest(quantity,locationName, confidenceAlpha)
        
    def preprocessSolutionData(self, vascularNetwork, simulationTime, sampleSize, sampleIndex):
        '''
        
        '''
    
        if "vessel" in self.locationName:
            vesselId = int(self.locationName.split('_')[-1])
            
            dataDict = vascularNetwork.getSolutionData(vesselId, self.quantitiesOfInterest.keys(), simulationTime, [self.xval])
            for quantitiyName,quantityObject in self.quantitiesOfInterest.iteritems():
                # normal quantity of interest
                if 'Extrema' not in quantitiyName and 'InflectionPoint' not in quantitiyName:    
                    if quantityObject.data == None: quantityObject.data = np.empty((sampleSize,len(simulationTime)))
                    quantityObject.data[sampleIndex] = dataDict[quantitiyName][:,0]
                
        elif "baroreceptor" in self.locationName:
            baroId = int(self.locationName.split('_')[-1])
            ## TODO: Jacob add here your stuff
                         
        else:
            ## TODO: add more locations as necessary baroreceptor etc.
            print "class LoacationOfInterest: location {} is not supported yet".format(self.locationName)
     
    def preprocessSolutionDataExtremaAndInflectionPoints(self, simulationTime, sampleSize):
        '''
        
        '''
        if "vessel" in self.locationName: 
            vesselId = int(self.locationName.split('_')[-1])
            
            for quantitiyName,quantityObject in self.quantitiesOfInterest.iteritems():
                
                if 'Extrema' in quantitiyName:
                    
                    quantityNamePure = quantitiyName.split('Extrema')[-1]
                    print quantityNamePure
                    print quantitiyName,quantityObject
                    
                    dataPure = self.quantitiesOfInterest[quantityNamePure].data
                                        
                    delta = 0.005
                    myData = []
                    numberOfPointsFirst = None
                    for sampleIndex in xrange(sampleSize):
                        dataPureCurr = dataPure[sampleIndex]
                        minMaxPoints = mProc.minMaxFunction(dataPureCurr,timeValues=simulationTime,delta=delta, seperateMinMax = False ) 
                        myData.append(minMaxPoints)
                        
                        numberOfPointsCurrent = len(minMaxPoints[0])
                        if numberOfPointsFirst == None: numberOfPointsFirst = numberOfPointsCurrent
                        
                        colors = ['r','g','b','m','k','c']
                        
                        for c,t,x in zip(colors,minMaxPoints[1],minMaxPoints[0]):
                            plt.plot(t,x,'o',color = c, markersize = 5)
                        
                        if numberOfPointsCurrent == numberOfPointsFirst:
                            plt.plot(simulationTime,dataPureCurr,'g', alpha = 0.2)
                        elif numberOfPointsCurrent < numberOfPointsFirst:
                            plt.plot(simulationTime,dataPureCurr,'m', alpha = 1.0)
                        elif numberOfPointsCurrent > numberOfPointsFirst:
                            plt.plot(simulationTime,dataPureCurr,'k', alpha = 1.0)
                            
                    #plt.show()
                    #print myData
                    #quantityObject.data = np.empty((sampleSize,len(simulationTime)))
                    
                    
                elif 'InflectionPoint' in quantitiyName:
                    print quantitiyName,quantityObject
            
        
            
    def saveQuantitiyOfInterestData(self, saveFile):
        '''
        method to save the data of the location and all quantities of interest to file
        '''
        locationGroup = saveFile.create_group(self.locationName)
        locationGroup.attrs.create('quantitiesOfInterestToProcess', data=self.quantitiesOfInterestToProcess)
        locationGroup.attrs.create('xval', data=self.xval)
        locationGroup.attrs.create('confidenceAlpha', data=self.confidenceAlpha)
        
        for quantitiyName,quantityObject in self.quantitiesOfInterest.iteritems():
                quantitiyGroup = locationGroup.create_group(quantitiyName)
                quantityObject.saveQuantitiyOfInterestData(quantitiyGroup)
                
    def loadQuantitiyOfInterestData(self, saveFile):
        '''
        
        '''
        