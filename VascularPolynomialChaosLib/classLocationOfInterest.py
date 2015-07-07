

from classQuantityOfInterest import QuantityOfInterest

from copy import copy as copy

import numpy as np

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
        
        for quantitiy in quantitiesOfInterestToProcess:
            if "ward" in quantitiy:
                quantitiesOfInterestToCreate.extend(['Flow','Pressure','Area','WaveSpeed'])
                quantitiesOfInterestToCreate = list(set(quantitiesOfInterestToCreate))
                break
        
        for quantity in quantitiesOfInterestToCreate:
            self.quantitiesOfInterest[quantity] = QuantityOfInterest(quantity,locationName, confidenceAlpha)
        
    def preprocessSolutionData(self, vascularNetwork, simulationTime, sampleSize, sampleIndex):
        '''
        
        '''
    
        if "vessel" in self.locationName:
            vesselId = int(self.locationName.split('_')[-1])
            dataDict = vascularNetwork.getSolutionData(vesselId, self.quantitiesOfInterest.keys(), simulationTime, [self.xval])
            
            for quantitiyName,quantityObject in self.quantitiesOfInterest.iteritems():
                if quantityObject.data == None: 
                    quantityObject.data = np.empty((sampleSize,len(simulationTime)))
                quantityObject.data[sampleIndex] = dataDict[quantitiyName][:,0]
            ##
            # TODO: peak detection and saving of amplitude and timing if wanted
            
            
        elif "baroreceptor" in self.locationName:
            baroId = int(self.locationName.split('_')[-1])
            dataDict = {}
            for quantitiyName in self.quantitiesOfInterest.keys():
                dataDict[quantitiyName] = vascularNetwork.baroreceptors[baroId].getVariableValue(quantitiyName)
            
            for quantitiyName,quantityObject in self.quantitiesOfInterest.iteritems():
                if quantityObject.data == None: 
                    quantityObject.data = np.empty((sampleSize,len(simulationTime)))
                quantityObject.data[sampleIndex] = dataDict[quantitiyName][:]
            ## TODO: Jacob add here your stuff
                         
        else:
            ## TODO: add more locations as necessary baroreceptor etc.
            print "class LoacationOfInterest: location {} is not supported yet".format(self.locationName)
            
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
        