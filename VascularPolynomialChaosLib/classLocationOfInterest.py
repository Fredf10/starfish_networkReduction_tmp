

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
        
        self.pointEvaluationQuantities = []
                        
        # check if vessel and add missing qoi to the list if needed
        quantitiesOfInterestToCreate = copy(quantitiesOfInterestToProcess)
        
        for quantity in quantitiesOfInterestToProcess:
            if 'Extrema' in quantity:
                quantityNamePure = quantity.split('Extrema')[-1]
                quantitiesOfInterestToCreate.extend([quantityNamePure])
                quantitiesOfInterestToCreate = list(set(quantitiesOfInterestToCreate))
                self.pointEvaluationQuantities.append(quantity)
                
            if 'InflectionPoint' in quantity:
                quantityNamePure = quantity.split('InflectionPoint')[-1]
                quantitiesOfInterestToCreate.extend([quantityNamePure])
                quantitiesOfInterestToCreate = list(set(quantitiesOfInterestToCreate))
                self.pointEvaluationQuantities.append(quantity)
        
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
        for quantityName in self.pointEvaluationQuantities:
            quantityObject = self.quantitiesOfInterest[quantityName]
            
            # TODO: outsource peak selection method
            
            if 'InflectionPoint' in quantityName:
                quantityNamePure = quantityName.split('InflectionPoint')[-1]
                searchPointOfInflection = True                
            
            elif 'Extrema' in quantityName:
                quantityNamePure = quantityName.split('Extrema')[-1]
                searchPointOfInflection = False
                
                
            dataPure = self.quantitiesOfInterest[quantityNamePure].data
                  
            allDesieredPointsDetected = False
            delta = 0.005
            
            colordef  = ['r','g','b','m','k','c']
            markerdef = ['o','d','x','s','*','v']
            colors       = colordef*len(markerdef)
            markerstyles = [item for sublist in [ [i]*len(colordef) for i in markerdef] for item in sublist]
                
                                            
            while allDesieredPointsDetected == False:   
                amplitudesList = []
                timingList = []
                numberOfPointsFirst = None
                for sampleIndex in xrange(sampleSize):
                    
                    dataPureCurr = dataPure[sampleIndex]
                    
                    if searchPointOfInflection == True:
                        
                        pointInfAmp,pointInfTime = mProc.calculatePointOfInflection(dataPureCurr, simulationTime, simulationTime[0])
                        minMaxPoints = [[pointInfAmp],[pointInfTime]]
                        numberOfPointsCurrent = len(minMaxPoints[0])
                        if numberOfPointsFirst == None: numberOfPointsFirst = numberOfPointsCurrent
                        
                    else:
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
                    
                                                                        
                    amplitudesList.append(minMaxPoints[0])
                    timingList.append(minMaxPoints[1])
                    
                    for c,m,t,x in zip(colors,markerstyles,minMaxPoints[1],minMaxPoints[0]):
                        plt.plot(t,x,'o',color = c,marker=m, markersize = 5)
                    
                    if numberOfPointsCurrent == numberOfPointsFirst:
                        plt.plot(simulationTime,dataPureCurr,'g', alpha = 0.2)
                    elif numberOfPointsCurrent < numberOfPointsFirst:
                        plt.plot(simulationTime,dataPureCurr,'m', alpha = 1.0)
                    elif numberOfPointsCurrent > numberOfPointsFirst:
                        plt.plot(simulationTime,dataPureCurr,'k', alpha = 1.0)
                
                # break for inflection point
                if searchPointOfInflection == True:
                    allDesieredPointsDetected = True
                    break
                        
                plt.title(quantityNamePure)
                plt.show(block = False)
                
                print '\n    Extrema detection all simulations are evaluated:'
                print '      [0] - Analyse all points'
                print '      [1] - Choose points for analysis'
                print '      [2] - Not all desiered peaks are detected (redo and adjust search-delta)'
                suggestions = [str(i) for i in xrange(3)]
                answer = "nothing"
                while answer not in suggestions:
                    answer = raw_input("What do you want to do: ")
                    
                if answer == '0':
                    answerSelectionSplit = [str(i) for i in xrange(numberOfPointsFirst)]
                    allDesieredPointsDetected = True
                                    
                if answer == '1':
                    allDesieredPointsDetected = True
                    print "\n    points found for analysis:"
                    for i in xrange(numberOfPointsFirst):
                        print "      [{}] - color: {}, marker style {}".format(i,colors[i],markerstyles[i])
                    answerSelectionSplit = ['']
                    while sum([item in [str(i) for i in xrange(numberOfPointsFirst)] for item in answerSelectionSplit]) != len(answerSelectionSplit):
                        answerSelection = raw_input("    please enter the numbers of the points (separated with comma): ")
                        if ',' in answerSelection:
                            answerSelectionSplit = answerSelection.split(',')
                        elif ' ' in answerSelection:
                            answerSelectionSplit = answerSelection.split(' ')
                        else:
                            answerSelectionSplit = [answerSelection]
                                       
                elif answer == '2':
                    answer2 = "haha"
                    convertableToFloat = False
                    while convertableToFloat == False:
                        answer2 = raw_input("please redefine delta, current delta value == {} ".format(delta))
                        try:
                            delta = float(answer2)
                            convertableToFloat = True
                        except ValueError as v:
                            print "{} not convertable to float!!".format(v)
                            
                #TODO: save the figures and close ..
                plt.close()
                
            # get the data of the selected points
            selectedPoints = np.array([int(i) for i in answerSelectionSplit])    
            amplitudeData =  np.array(amplitudesList)
            timingData    =  np.array(timingList)
            
            for selectedPoint in selectedPoints:
                # create Names
                qOINameAmplitude = ''.join([quantityNamePure,'Extremum',str(selectedPoint),'Amplitude'])
                qOINameTiming    = ''.join([quantityNamePure,'Extremum',str(selectedPoint),'Timing'])
                # rename if point of inflection
                if searchPointOfInflection == True:              
                    qOINameAmplitude = ''.join([quantityNamePure,'InflectionPointAmplitude'])
                    qOINameTiming    = ''.join([quantityNamePure,'InflectionPointTiming'])
                # create 2 new quantities of interest
                quantityOfInterestAmplitude = QuantityOfInterest(qOINameAmplitude, quantityObject.locationName, quantityObject.confidenceAlpha, data = amplitudeData[:,selectedPoint])
                quantityOfInterestTiming = QuantityOfInterest(qOINameTiming, quantityObject.locationName, quantityObject.confidenceAlpha, data = timingData[:,selectedPoints])
                # append them in the quqntity of interest vector
                self.quantitiesOfInterest[qOINameAmplitude] = quantityOfInterestAmplitude
                self.quantitiesOfInterest[qOINameTiming] = quantityOfInterestTiming
                # append them to the quantitiy to process list
                self.quantitiesOfInterestToProcess.append(qOINameAmplitude)
                self.quantitiesOfInterestToProcess.append(qOINameTiming)
            
            # remove the old one        
            del self.quantitiesOfInterest[quantityName]
            self.quantitiesOfInterestToProcess.remove(quantityName)
            
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
        