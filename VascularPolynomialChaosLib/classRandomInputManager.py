

from classRandomInput import RandomInput
import numpy as np
import moduleFilePathHandlerVPC as mFPH_VPC

import itertools

class RandomInputManager(object):
    '''
    
    '''
    def __init__(self):        
        self.map = {} # hash map with {randomInput.location : randomInput.randomInputId}
        self.randomInputs = [] # randomInput as they stand in xml
        
        self.randomInputVector = [] # randomInput which have a external distribution assosiated
        
        
    def addRandomInput(self,randomInputDict):
        '''
        Function which adds a random variable to randomInputs
        '''
        randomInputLocation = randomInputDict['location']
        
        if randomInputLocation not in self.map.keys():
            randomInput = RandomInput()
            randomInputId = len(self.randomInputs)
            
            randomInputDict['randomInputId'] = randomInputId
            randomInput.update(randomInputDict)
            
            self.randomInputs.append(randomInput)
            self.map[randomInputLocation] = randomInputId
        else:
            print "Random input with location {} is already defined!".format(randomInputLocation)
                    
    def __call__(self, index = None):
        '''
        call function returns the random vector
        '''
        if index == None:
            return self.randomInputs
        else:
            try:
                return self.randomInputs[index]
            except:
                return []
            
    def linkRandomInputUpdateFunctions(self, vascularNetwork):
        '''
        link update functions and create randomInputvector
        
        input vascularNetwork 
        '''
        
        randomInputMap = {}
        for randomInput in self.randomInputs:
            if randomInput.type == 'parametricRandomInput':
                # check distribution
                dist = randomInput.distributionType 
                loc = randomInput.location.split('_')
                objType = loc[0]
                
                if objType == "boundaryCondition":
                    for bc in vascularNetwork.boundaryConditions[int(loc[2])]:
                        if bc.getVariableValue('name') == loc[1]:
                            randomInput.updateMethods = {randomInput.variableName[0]:
                                                         bc.update}
            
                elif objType == "vessel":
                    randomInput.updateMethods = {randomInput.variableName[0]:
                                                 vascularNetwork.vessels[int(loc[1])].update}
                else: break
                if randomInput.updateMethods == {}: break
                
                if dist in ['Normal','Uniform']:
                    # append random input to random input vector
                    self.randomInputVector.append(randomInput)
                else:
                    if dist not in randomInputMap:
                        randomInputMap[dist] = {randomInput.randomInputId : randomInput.passRealisationToAssosiatedObj}
                    else:
                        randomInputMap[dist][randomInput.randomInputId] = randomInput.passRealisationToAssosiatedObj
        
        for dist,updateMethods in randomInputMap.iteritems():
            # find random input with dist:
            for randomInputLocation in self.map.keys():
                if dist in randomInputLocation:
                    generalRandomInput = self(self.map[randomInputLocation])
                    generalRandomInput.updateMethods = updateMethods
                    generalRandomInput.variableName = updateMethods.keys()
                    self.randomInputVector.append(generalRandomInput)
         
    def update(self, dataDict):
        '''
        updates the data using a dictionary in from of 
        dataDict = {'variableName': value}
        '''
        for key,value in dataDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except:
                print "ERROR RandomInputManager.updateData Wrong key: {}, could not update varibale".format( key)
           
    def saveRealisationLog(self, networkName, dataNumber, gPCEmethod, gPCEorder):
        '''
        method to save a log file with the realisations passed for each sample iteration
        '''
        
        logData = np.array([randomInput.updateLog for randomInput in self.randomInputs]).transpose()
        evaluationLogFile = mFPH_VPC.getFilePath('evaluationLogFile', networkName, dataNumber, mode = "write", gPCEmethod=gPCEmethod, gPCEorder=gPCEorder)
        
        logfile = open(evaluationLogFile, "wb")
        logfile.write(''.join(['Stochastic simulation ',str(networkName),' DataNumber ',dataNumber,'\n','\n']))
        logfile.write(''.join(['gPCEorder :',str(gPCEorder),'  gPCEmethod: ',gPCEmethod,'  number of evaluations: ',str(len(logData)),'\n','\n']))
        
        for info in self.generateInfo():
            logfile.write(''.join([info,'\n']))
        randomInputIds = [randomInput.randomInputId for randomInput in self.randomInputs]
        logfile.write(''.join(['\n{:<10}'.format(''), '{:20}'.format('RandomVariableId'),'\n','\n']))
        logfile.write(''.join(['{:<10}'.format('EvalNr'), ''.join(['{:20}'.format(j) for j in ['{:<5}'.format(k) for k in randomInputIds]]),'\n',]))     
        for i,logLine in enumerate(logData):
            logfile.write(''.join(['{:<10}'.format(i), ''.join(['{:20}'.format(j) for j in ['{:<.5}'.format(k) for k in logLine]]),'\n']))
        
        
    def printOutInfo(self):
        '''
        Function to print out random variable informations
        '''
        for info in self.generateInfo():
            print info
            
    def generateInfo(self):
        '''
        Function which generates info of the random inputs
        and returns it in a list
        '''
        randomInputManagerInfo = []
        randomInputManagerInfo.append("\n Defined Random Inputs\n")
        randomInputManagerInfo.append( '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","distribution"))
        randomInputManagerInfo.append( "-------------------------------------------------------------------- \n")
        randomInputInfos = list(itertools.chain.from_iterable([randomInput.generateInfo() for randomInput in self.randomInputs]))
        for info in randomInputInfos:
            randomInputManagerInfo.append( info)
        randomInputManagerInfo.append( "\n Defined Random Variables\n")
        randomInputManagerInfo.append( '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","distribution"))
        randomInputManagerInfo.append( "-------------------------------------------------------------------- \n")
        randomInputInfos = list(itertools.chain.from_iterable([randomInput.generateInfo() for randomInput in self.randomInputs]))
        for info in randomInputInfos:
            randomInputManagerInfo.append(info)
        return randomInputManagerInfo