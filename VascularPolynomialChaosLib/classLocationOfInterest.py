

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
     
     
    def decideOnPoints(self):
        pass
    
    def preprocessSolutionDataExtremaAndInflectionPoints(self, simulationTime, sampleSize):
        '''
        
        '''
        if "vessel" in self.locationName: 
            vesselId = int(self.locationName.split('_')[-1])
            
            for quantitiyName,quantityObject in self.quantitiesOfInterest.iteritems():
                
                if 'Extrema' in quantitiyName:
                    print "cKOI80, calculating Extrema now"
                    quantityNamePure = quantitiyName.split('Extrema')[-1]
                    print quantityNamePure
                    print quantitiyName,quantityObject
                    
                    dataPure = self.quantitiesOfInterest[quantityNamePure].data
                          
                    allDesieredPointsDetected = False
                    delta = 0.005
                    
                    colordef  = ['r','g','b','m','k','c']
                    markerdef = ['o','d','x','s','*','v']
                    
                    colors       = colordef*len(markerdef)
                    markerstyles = [item for sublist in [ [i]*len(colordef) for i in markerdef] for item in sublist]
                            
                    
                    fig = plt.figure(1)
                    
                    while allDesieredPointsDetected == False:              
                        myData = []
                        numberOfPointsFirst = None
                        for sampleIndex in xrange(sampleSize):
                            
                            dataPureCurr = dataPure[sampleIndex]
                            minMaxPoints = mProc.minMaxFunction(dataPureCurr,timeValues=simulationTime,delta=delta, seperateMinMax = False ) 
                            
                            numberOfPointsCurrent = len(minMaxPoints[0])
                            if numberOfPointsFirst == None: numberOfPointsFirst = numberOfPointsCurrent
                            
                            deltaAltered = delta
                            while numberOfPointsCurrent != numberOfPointsFirst:
                                if numberOfPointsCurrent < numberOfPointsFirst:
                                    deltaAltered -= delta*0.1
                                elif numberOfPointsCurrent > numberOfPointsFirst:
                                    deltaAltered += delta*0.1
                                dataPureCurr = dataPure[sampleIndex]
                                minMaxPoints = mProc.minMaxFunction(dataPureCurr,timeValues=simulationTime,delta=deltaAltered, seperateMinMax = False ) 
                                numberOfPointsCurrent = len(minMaxPoints[0])
                                                    
                            myData.append(minMaxPoints)
                            
                            for c,m,t,x in zip(colors,markerstyles,minMaxPoints[1],minMaxPoints[0]):
                                plt.plot(t,x,'o',color = c,marker=m, markersize = 5)
                            
                            if numberOfPointsCurrent == numberOfPointsFirst:
                                plt.plot(simulationTime,dataPureCurr,'g', alpha = 0.2)
                            elif numberOfPointsCurrent < numberOfPointsFirst:
                                plt.plot(simulationTime,dataPureCurr,'m', alpha = 1.0)
                            elif numberOfPointsCurrent > numberOfPointsFirst:
                                plt.plot(simulationTime,dataPureCurr,'k', alpha = 1.0)
                                
                        plt.title(quantityNamePure)
                        plt.show(block = False)
                        
                        print '\n    Extrema detection all simulations are evaluated:'
                        print '      [0] All desiered peaks are detected and not more or less'
                        print '      [1] Not all desiered peaks are detected (redo and adjust search-delta)'
                        suggestions = [str(i) for i in xrange(2)]
                        answer = "nothing"
                        while answer not in suggestions:
                            answer = raw_input("What do you want to do: ")
                        if answer == '0':
                            allDesieredPointsDetected = True
                            print "\n    please enter the numbers of the points you want to analyse:"
                            for i in xrange(numberOfPointsFirst):
                                print "      [{}] - color: {}, marker style {}".format(i,colors[i],markerstyles[i])
                            
                            
                            
                        elif answer == '1':
                            answer2 = "haha"
                            convertableToFloat = False
                            while convertableToFloat == False:
                                answer2 = raw_input("please redefine delta, current delta value == {} ".format(delta))
                                try:
                                    delta = float(answer2)
                                    convertableToFloat = True
                                except ValueError as v:
                                    print "{} not convertable to float!!".format(v)
                        plt.close()
                    #print myData
                    #quantityObject.data = np.empty((sampleSize,len(simulationTime)))
                    
                elif 'InflectionPoint' in quantitiyName:
                    print quantitiyName,quantityObject
            
            print "LOI199 \n"
                        
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
        