

from copy import copy as copy
import numpy as np

from classQuantityOfInterest import QuantityOfInterest
import UtilityLib.processing as mProc
import matplotlib.pyplot as plt

from testBaseClass import TestBaseClass 

class LocationOfInterest(TestBaseClass):
    '''
    
    '''
    
    # defined external data
    externVariables      = {'quantitiesOfInterestToProcess' : TestBaseClass.ExtValue(str, strCases = ['anything'], multiVar = True),
                             'queryLocation'                : TestBaseClass.ExtValue(str, strCases = ['anything']),
                             'xVal'                         : TestBaseClass.ExtValue(float,  unit = 'm'),
                             'confidenceAlpha'              : TestBaseClass.ExtValue(float)}
    
    externXmlAttributes  = []
    
    externXmlElements    = ['queryLocation',
                             'xVal',
                             'quantitiesOfInterestToProcess',
                             'confidenceAlpha']
            
        
#         self.externHdf5 = ['queryLocation',
#                            'xVal',
#                            'quantitiesOfInterestToProcess']
#         
#         self.externCSV = ['queryLocation',
#                           'xVal',
#                           'quantitiesOfInterestToProcess']
        
    def __init__(self):
        
        self.queryLocation                 = "queryLocation"
        self.xVal                          = 0 #xVal #position in x 
        self.quantitiesOfInterestToProcess = [] #quantitiesOfInterestToProcess
        self.quantitiesOfInterest          = {}
        self.confidenceAlpha               = 1 #confidence Alpha
        
        self.pointEvaluationQuantities = []
        self.trajectoryEvaluationQuantities = []
        
    def initialize(self):
                        
        # check if vessel and add missing qoi to the list if needed
        quantitiesOfInterestToCreate = copy(self.quantitiesOfInterestToProcess)
        
        for quantity in self.quantitiesOfInterestToProcess:
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
                
            if 'Trajectory' in quantity:
                self.trajectoryEvaluationQuantities.append(quantity)
        
        for quantity in quantitiesOfInterestToCreate:
            self.quantitiesOfInterest[quantity] = QuantityOfInterest(quantity, self.queryLocation, self.confidenceAlpha)
        
    def preprocessSolutionData(self, vascularNetwork, simulationTime, sampleSize, sampleIndex):
        '''
        
        '''
    
        if "vessel" in self.queryLocation:
            vesselId = int(self.queryLocation.split('_')[-1])
            
            dataDict = vascularNetwork.getSolutionData(vesselId, self.quantitiesOfInterest.keys(), simulationTime, [self.xVal])
            for quantitiyName,quantityObject in self.quantitiesOfInterest.iteritems():
                # normal quantity of interest
                if 'Extrema' not in quantitiyName and 'InflectionPoint' not in quantitiyName and 'Trajectory' not in quantitiyName:    
                    if quantityObject.data == None: quantityObject.data = np.empty((sampleSize,len(simulationTime)))
                    quantityObject.data[sampleIndex] = dataDict[quantitiyName][:,0]
                
                if 'Trajectory' in quantitiyName:
                    # allocates data for a vessel over space and saves the raw data in quantityObject.dataSpace
                    #TODO: create own class tha xVal is not missued as number of space points per vessel
                    vesselIds = [int(i) for i in self.queryLocation.split('vessel')[-1].split('_') if i != '']
                    data = None
                    
                    quantitiyPure = quantitiyName.split('Trajectory')[-1] 
                    
                    if "FB_" in quantitiyPure:
                        
                        quantitiyPure = quantitiyPure.split('FB_')[-1]
                        
                        maxLength = self.xVal
                        
                        maxNumberPoints = 20
                        nPointsUsed = 0
                        position = 0
                        
                        lastVessel = False
                        trajectory = None                   
                        
                        if len(vesselIds) == 1:
                            vesselId = vesselIds[0]
                            length = vascularNetwork.vessels[vesselId].length
                            
                            maxLength = length*2
                            # do it twice once forward and once backward
                            for direction in ['Forward','Backward']:
                            
                                
                                remainingLength = maxLength-position
                                if length < remainingLength:
                                    # get number of points the solution in x
                                    nPoints = int(maxNumberPoints*length/maxLength)
                                    xValEnd = length
                                    
                                if length >= remainingLength:
                                    xValEnd = remainingLength
                                    lastVessel = True
                                    nPoints = maxNumberPoints - nPointsUsed
                                
                                xvals = np.linspace(0,xValEnd,nPoints)
                                
                                dataTemp = vascularNetwork.getSolutionData(vesselId, [''.join([direction,quantitiyPure])], simulationTime, xvals = xvals)[''.join([direction,quantitiyPure])]
                                
                                if data == None: data = dataTemp.T
                                else:
                                    if direction == 'Forward':  
                                        data = np.vstack([data,dataTemp.T])
                                    else:
                                        data = np.vstack([data,dataTemp.T[::-1]])
                                
                                if trajectory == None: trajectory = xvals
                                else: 
                                    #trajectory = np.append(trajectory,xvals+position)
                                    
                                    if direction == 'Forward':                                    
                                        trajectory = np.append(trajectory,xvals+position)
                                    else:
                                        trajectory = np.append(trajectory,xvals[::-1])
                                
                                position = position+xValEnd
                                nPointsUsed = nPointsUsed + nPoints
                                
                                if lastVessel: 
                                    #print "reached end of trajectory breaking"
                                    break
                    
                    
                    else:                    
                    
                        maxLength = self.xVal
                        
                        maxNumberPoints = 20
                        nPointsUsed = 0
                        position = 0
                        
                        lastVessel = False
                        trajectory = None                   
                        
                        
                        
                        for vesselId in vesselIds:
                                
                            length = vascularNetwork.vessels[vesselId].length
                            remainingLength = maxLength-position
                           
                            if length < remainingLength:
                                # get number of points the solution in x
                                nPoints = int(maxNumberPoints*length/maxLength)
                                
                            if length >= remainingLength:
                                length = remainingLength
                                lastVessel = True
                                nPoints = maxNumberPoints - nPointsUsed
                                
                            xvals = np.linspace(0,length,nPoints)
                            dataTemp = vascularNetwork.getSolutionData(vesselId, [quantitiyPure], simulationTime, xvals = xvals)[quantitiyPure]
                            
                            if data == None: data = dataTemp.T
                            else: data = np.vstack([data,dataTemp.T])
                            
                            if trajectory == None: trajectory = xvals
                            else: trajectory = np.append(trajectory,xvals+position)
                            
                            position = position+length
                            nPointsUsed = nPointsUsed + nPoints
                            
                            if lastVessel: 
                                print "reached end of trajectory breaking"
                                break
                       
                    
                       
                    if quantityObject.trajectoryData == None:
                        quantityObject.trajectoryData =  np.empty((sampleSize,maxNumberPoints))   
                    quantityObject.trajectoryData[sampleIndex] = trajectory
                          
                    print sampleSize,maxNumberPoints,len(simulationTime)
                    if quantityObject.data == None: quantityObject.data = np.empty((sampleSize,maxNumberPoints,len(simulationTime)))
                    quantityObject.data[sampleIndex] = data
                    
                
        elif "baroreceptor" in self.queryLocation:
            baroId = int(self.queryLocation.split('_')[-1])
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
            print "class LoacationOfInterest: location {} is not supported yet".format(self.queryLocation)
     
     
    
    def preprocessSolutionDataExtremaAndInflectionPoints(self, simulationTime, sampleSize):
        '''
        
        '''
        for quantityName in self.pointEvaluationQuantities:
            quantityObject = self.quantitiesOfInterest[quantityName]
            
            if 'InflectionPoint' in quantityName:
                quantityNamePure = quantityName.split('InflectionPoint')[-1]
                searchPointOfInflection = True                
            
            elif 'Extrema' in quantityName:
                quantityNamePure = quantityName.split('Extrema')[-1]
                searchPointOfInflection = False
                
            dataPure = self.quantitiesOfInterest[quantityNamePure].data
                  
            selectedPoints, amplitudeData, timingData = self.extremaFinderFunction(dataPure, quantityNamePure, sampleSize, simulationTime, searchPointOfInflection)
                        
            for selectedPoint in selectedPoints:
                # create Names
                qOINameAmplitude = ''.join([quantityNamePure,'Extremum',str(selectedPoint),'Amplitude'])
                qOINameTiming    = ''.join([quantityNamePure,'Extremum',str(selectedPoint),'Timing'])
                # rename if point of inflection
                if searchPointOfInflection == True:              
                    qOINameAmplitude = ''.join([quantityNamePure,'InflectionPointAmplitude'])
                    qOINameTiming    = ''.join([quantityNamePure,'InflectionPointTiming'])
                # create 2 new quantities of interest
                quantityOfInterestAmplitude = QuantityOfInterest(qOINameAmplitude, quantityObject.queryLocation, quantityObject.confidenceAlpha, data = amplitudeData[:,selectedPoint])
                quantityOfInterestTiming = QuantityOfInterest(qOINameTiming, quantityObject.queryLocation, quantityObject.confidenceAlpha, data = timingData[:,selectedPoints])
                # append them in the quqntity of interest vector
                self.quantitiesOfInterest[qOINameAmplitude] = quantityOfInterestAmplitude
                self.quantitiesOfInterest[qOINameTiming] = quantityOfInterestTiming
                # append them to the quantitiy to process list
                self.quantitiesOfInterestToProcess.append(qOINameAmplitude)
                self.quantitiesOfInterestToProcess.append(qOINameTiming)
            
            # remove the old one        
            del self.quantitiesOfInterest[quantityName]
            self.quantitiesOfInterestToProcess.remove(quantityName)
            
    def preprocessSolutionDataTrajectory(self, simulationTime, sampleSize):
        
        
        import matplotlib.pyplot as plt
                
        for quantityName in self.trajectoryEvaluationQuantities:
            quantityObject = self.quantitiesOfInterest[quantityName]
            
            trajectoryData = quantityObject.trajectoryData
            
            data = np.swapaxes(quantityObject.data, 0, 1)
            
            useMax = True
            
            
            if useMax == True:
                
                
                # find maxima of the data
                maxValue = np.amax(data,axis=2)
                ## currently not used as qoI for gpc
                
                # find indices where the data has maximums
                maxIndices =  np.argmax(data, axis=2)
                
                # get time points of the maximums
                maxTimes = simulationTime[maxIndices].T
                
                fig = plt.figure()
                fig.canvas.set_window_title(''.join([quantityName,'PC']))
                for x,t in zip(trajectoryData,maxTimes):
                    plt.plot(t,x)       
                #plt.show()
                
                fixSpace = False
                
                if fixSpace == True:
                
                    xN = 200
                    xInt = np.linspace(0, self.xVal, xN)
                    maxTimeMatched = np.empty((sampleSize,xN))
                    
                    
                    for i in xrange(sampleSize): 
                        maxTimeMatched[i] = np.interp(xInt, trajectoryData[i],maxTimes[i])
                    quantityObject.trajectoryData = xInt
                    quantityObject.data = maxTimeMatched
                    print "Approximated error in forward wave jump: as a multiple of timesteps"
                    print np.abs(np.min(maxTimes[:,1::]-maxTimes[:,0:-1], axis=1))/(simulationTime[1]-simulationTime[0])
                    
                                         
                    fig = plt.figure()
                    fig.canvas.set_window_title(''.join([quantityName,'PC_interp']))
                    for t in maxTimeMatched:
                        plt.plot(xInt,t)
                
                else:
                    #fix time
                    
                    #find the minimum of all return times
                    startTimes = np.min(maxTimes,axis=1)
                    maxStartTime = np.max(startTimes)
                    
                    endTimes = np.max(maxTimes,axis=1)
                    minEndTime = np.min(endTimes)
                    
                    tN = 200
                    tInt = np.linspace(maxStartTime,minEndTime,tN)
                    
                    spacePointMatched = np.empty((sampleSize,tN))
                    
                    for i in xrange(sampleSize): 
                        spacePointMatched[i] = np.interp(tInt, maxTimes[i], trajectoryData[i])
                        
                    quantityObject.trajectoryData = tInt
                    quantityObject.data = spacePointMatched
                    
                    fig = plt.figure()
                    fig.canvas.set_window_title(''.join([quantityName,'PC_interp']))
                    for x in spacePointMatched:
                        plt.plot(tInt,x)
                    
                    
                plt.show()
                
            else:
                # find peaks                
                for dataOfSpacePoint in xrange(len(data)):
                    self.extremaFinderFunction(dataOfSpacePoint, quantityName, sampleSize, simulationTime, searchPointOfInflection = False)
                      
              
    def extremaFinderFunction(self, dataPure, quantityNamePure, sampleSize, simulationTime, searchPointOfInflection):
        
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
            
        return selectedPoints, amplitudeData, timingData
            
    def saveDataHdf5(self, hdf5SaveGroup):
        '''
        method to save the data of the location and all quantities of interest to file
        '''
        # pure variables
        
        # in class definition
        variablesToSave = ["quantitiesOfInterestToProcess",
                           'xVal',
                           'confidenceAlpha']
        # functionality
        for variableName in variablesToSave:
            variableValue = self.getVariable(variableName)
            if variableValue != None: 
                hdf5SaveGroup.create_dataset(variableName, data=variableValue)
        
        # dict of objects
        
        dictsToSave = ["quantitiesOfInterest"]
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
        variablesToLoad = ["quantitiesOfInterestToProcess",
                           'xVal',
                           'confidenceAlpha']
        # functionality
        for variableName in variablesToLoad:
            if variableName in hdf5SaveGroup.keys(): 
                variableData = (hdf5SaveGroup[variableName])
                # TODO: use externalVariable definitions dictionary ... #
                # check for shape 
                if np.shape(variableData)== ():
                    loadedData[variableName] = variableData[()] # for scalar data sets!!
                else:
                    loadedData[variableName] = variableData[:]
        
        ## dictionary to load
        dictsToLoad = ["quantitiesOfInterest"]
        
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
        