import numpy as np 

import pprint
import subprocess
import math
import numpy


import sys,os,time
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../')
#sys.path.append(cur+'/../Visualisation')

from VisualisationLib.classRealTimeVisualisation import realTimeVisualisation

class CommunicatorBaseClass(object):
    '''
    Base-class for all boundary conditions    
    
    including variables:
        communicatorId
        filename
        blockingWait
    
    including functions:
        update
        readCommunicatorFile
        writeCommunicatorFile
    '''
    communicatorId = 0
    filename = None
    blockingWait = False
        
    def update(self,comDict):
        '''
        updates the updateCommunicatorDict data using a dictionary in form of 
        comDict = {'variableName': value}
        '''
        for key,value in comDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except: 
                print 'ERROR communicator.update(): wrong key: %s, could not set up communicator' %key
                
    def readCommunicatorFile(self, samplingTime = 0.001):
        '''
        This function reads a data string out of a file with given filename
        
        input:    filename <string>
                  blockingWaitRead <bool>
                  timeWait = 0.5 sec <float> time function waits between reading
                  maxTime  = 1.0 ??  <float> max time until connection is stopped 
                                
        return: dataString, None
        '''
        dataString = None
        
        while dataString == None:
            try:
                with open(''.join([cur,'/',self.filenameRead]),'r') as dataFile:
                    for dataLine in dataFile:
                        dataString = dataLine
                        break
                    os.remove(''.join([cur,'/',self.filenameRead]))
            except: 
                if self.blockingWaitRead == False:
                    break
            time.sleep(samplingTime)
                        
        return dataString
    
    def writeCommunicatorFile(self, dataString, samplingTime = 0.0001, maxTime = 0.1):
        '''
        This function reads a data string out of a file with given filename
        
        input:    dataString <string>
        '''
        
        fileExist = False
        if self.blockingWaitWrite == True:
            while fileExist == False:
                try:
                    open(''.join([cur,'/',self.filenameWrite])) 
                except:
                    fileExist = True
                    with open(''.join([cur,'/',self.filenameWrite]),'w') as dataFile:
                        dataFile.write(dataString)         
        else:
            with open(''.join([cur,'/',self.filenameWrite]),'w') as dataFile:
                dataFile.write(dataString)
                
class CommunicatorRealTimeViz(CommunicatorBaseClass):
    def __init__(self, comDict): 
        '''
        Communicator class envokes and communicates with realTimeVisualisation
        '''
        # def variables to set with comDict          
        self.currentMemoryIndex = [0] # n
        self.currentTimeIndex   = [0]
        self.currentTimeStep    = None
        self.dn                 = 1   #dn
        self.dt                 = 0.1 #dt
        self.node               = 0   # node
        self.quantitiesToPlot   = []  # quantitiesToPlot  #['Pressure','Flow']
        self.data               = {}  # data # dictionary with data to plot! data should be acsessable with [quantity][time n][node]
        self.vesselId           = 0
        self.comType            = ''
        self.comId              = 0
        
        ### to delete
        self.boundaryCondition  = 0
        
        # update variables with values in comDict      
        self.update(comDict)
        
        # fixed variables
        self.blockingWaitRead   = False
        self.blockingWaitWrite  = True
        self.count              = 0
        self.filenameWrite      = str(self.communicatorId).join(['../.realtimeWrite','.viz'])
        self.filenameRead       = str(self.communicatorId).join(['../.realtimeRead','.viz'])
        self.unitDict           = {'Pressure': 1./133.32, 'Flow': 1.e6, 'Area': 1000., 'InOutVolum':1.e6}
        
        #if 'bloodVolume' in self.quantitiesToPlot:
        
        
        # variables depending on init variables
        self.inititalValues     = []
        for quantity in self.quantitiesToPlot:
            self.inititalValues.append(str(self.data[quantity][0][self.node]*self.unitDict[quantity]))
                
        try: os.remove(''.join([cur,'/',self.filename]))
        except: pass
    
    def startRealTimeVisualisation(self):
        '''
        This function starts the realtime visualisation in a subprocess
        '''
        # start the visualisation
        visualisationProcess = ' '.join(['python',cur+'/../VisualisationLib/classRealTimeVisualisation.py',
                                         '-f',str(self.filenameWrite),
                                         '-t',str(self.dt*self.dn),
                                         '-q','-'.join(self.quantitiesToPlot),
                                         '-i','='.join(self.inititalValues)]) 
        self.realTimeVisualisationProcess = subprocess.Popen(visualisationProcess, shell = True )
                 
    def __call__(self):
        '''
        call function of communicator realtimevisualisation
        saves data in file
        '''
        n = self.currentMemoryIndex[0]
        self.count += 1
        
        if self.count == self.dn:        
            dataDict = {}
            for quantity in self.quantitiesToPlot:
                dataDict[quantity] = self.data[quantity][n][self.node]*self.unitDict[quantity]
            
            self.writeCommunicatorFile(str(dataDict))            
            self.count = 0

    def terminateRealtimeViz(self):
        '''
        Terminates subprocess realTimeVisualisation
        '''
        self.realTimeVisualisationProcess.terminate()

    def stopRealtimeViz(self):
        '''
        Stops realTimeVisualisation.update method
        called after simulation, if this is called the graphs can not be updated again
        '''
        self.writeCommunicatorFile('STOP')


class CommunicatorBaroreceptor(CommunicatorBaseClass):
    def __init__(self, comDict):
        '''
        
        '''
        # def variables to set with comDict         
         
        self.currentTimeStep    = [0] # n
        self.currentMemoryIndex = [0]
        self.dn                 = 1  #dn
        self.dt                 = 0.1 #dt
        self.node               = 0   # node
        self.data               = {}  # data # dictionary with data to plot! data should be acsessable with [quantity][time n][node]
        self.vesselId           = 0
        self.boundaryCondition  = 0
        #self.comType            = ''
        self.comType            = 'CommunicatorBaroreceptor'
        self.comId              = 0
        
        ### to delete
        self.quantitiesToPlot   = []
        
        # update variables with values in comDict
        self.update(comDict)
        
        
        # fixed variables
        self.blockingWaitRead   = False
        self.blockingWaitWrite  = True
        self.count              = 0
        self.filenameRead       = str(self.communicatorId).join(['../baroreceptorRead','.baro'])
        self.filenameWrite      = str(self.communicatorId).join(['../baroreceptorWrite','.baro'])
        self.unitDict           = {'Pressure': 1./133.32, 'Flow': 1.e6, 'Area': 1000.}
        
        # variables depending on init variables
        self.inititalValues     = []
        self.Ao                 = self.data['Area'][0]
        
        # Area, radius and strain
        self.Ro = numpy.power(self.Ao/math.pi,0.5)
        self.A = self.Ao
        self.R = self.Ro
        
        sizeEpsilon = numpy.shape(self.A)
        
        self.epsilon = numpy.zeros(sizeEpsilon)
        self.epsMean = 0
        
                
        try:    os.remove(''.join([cur,'/',self.filenameRead]))
        except: pass
        try:    os.remove(''.join([cur,'/',self.filenameWrite]))
        except: pass
        
    def __call__(self):
        '''
        
        '''
        n = self.currentTimeStep[0]
        self.count += 1
        
        ## calculate strain and save to data dictionary
        
        self.A = self.data['Area'][n]
        self.R = numpy.power(self.A/math.pi,0.5)
        
        self.epsilon = (self.R - self.Ro)/self.Ro
        self.epsMean = numpy.mean(self.epsilon)
        
        self.data['Strain'][n] = self.epsilon
        self.data['MStrain'][n] = self.epsMean
        print 'epsilon'
        print self.epsMean
        
                       
        time1 = 5.
        time2 = 10
        a = True
        if a == True:
            if n*self.dt > time1-self.dt/2. and n*self.dt < time1+self.dt/2.: 
                
                ## read heart rate and apply it to boundaryConditions
                Tnew = 0.8 ## read from file
                updateTime = n*self.dt
                        
                self.boundaryCondition.updatePeriodRuntime(Tnew,updateTime)
                
    
            if n*self.dt > time2-self.dt/2. and n*self.dt < time2+self.dt/2.: 
                
                Tnew = 3.0 ## read from file
                updateTime = n*self.dt
                        
                self.boundaryCondition.updatePeriodRuntime(Tnew,updateTime)

if __name__ == '__main__':
    
    com = CommunicatorRealTimeViz()
    while True:
        com()
