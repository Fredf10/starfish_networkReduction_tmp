 
 

from classRandomInput import RandomInput
import numpy as np

import itertools

class RandomInputManager(object):
    '''
    
    '''
    def __init__(self):        
        self.map = {} # hash map with {randomInput.location : randomInput.randomInputId}
        self.randomInputs = [] # randomInput as they stand in xml
        
        self.randomInputVector = [] # randomInput which have a external distribution assosiated
        self.randomInputDimension = 0
        
    def addRandomInput(self,randomInputDict):
        '''
        Function which adds a random variable to randomInputs
        '''
        randomInputLocation = randomInputDict['location']
        
        if randomInputLocation == '' or randomInputLocation == None:
            if randomInputDict['randomInputType'] == 'generalRandomInput':
                randomInputLocation = ''.join([randomInputDict['randomInputType'],'_',randomInputDict['name']])
                randomInputDict['location'] = randomInputLocation
        
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
            except IndexError:
                print "cRIM49: error index does not exist!"
                return []
                    
    def linkRandomInputUpdateFunctions(self, vascularNetwork):
        '''
        link update functions and create randomInputvector
        
        input vascularNetwork 
        
        only 1 level of nested randomInputs is currently implemented
        
        parametricRandomInputs can be in the form of:
            a+b*Uniform(0,1) / a+b*Normal(0,1)
            or
            a+b*Z1 where Z1 is a general random input defined as 
        
        '''
        
        randomInputMap = {}
        for randomInput in self.randomInputs:
            dist = randomInput.distributionType 
            if randomInput.randomInputType == 'parametricRandomInput':
                # check distribution
                loc = randomInput.location.split('_')
                objType = loc[0]
                
                if randomInput.variableName == []:
                    randomInput.variableName.append(loc[-1])
                
                if objType == "boundaryCondition":
                    for bc in vascularNetwork.boundaryConditions[int(loc[2])]:
                        if bc.getVariableValue('name') == loc[1]:
                            randomInput.updateMethods = {randomInput.variableName[0]:
                                                         bc.update}
            
                elif objType == "vessel":
                    randomInput.updateMethods = {randomInput.variableName[0]:
                                                 vascularNetwork.vessels[int(loc[1])].update}
                
                elif objType == "baroreceptor":
                    randomInput.updateMethods = {randomInput.variableName[0]:
                                                 vascularNetwork.baroreceptors[int(loc[1])].update}
                    
                    print "cRIM81, linked baroreceptor",int(loc[1]), 'to ',vascularNetwork.baroreceptors[int(loc[1])].update
                
                else: break
                if randomInput.updateMethods == {}: break
                
            # find links between random inputs
            if dist in ['Normal','Uniform']:
                # append random input to random input vector
                self.randomInputVector.append(randomInput)
            else:
                if dist not in randomInputMap:
                    randomInputMap[dist] = {randomInput.randomInputId : randomInput.passRealisationToAssosiatedObj}
                else:
                    randomInputMap[dist][randomInput.randomInputId] = randomInput.passRealisationToAssosiatedObj
        
        # link general random inputs 
        for dist,updateMethods in randomInputMap.iteritems():
            # find random input with dist:
            for randomInputLocation in self.map.keys():
                if dist in randomInputLocation:
                    generalRandomInput = self(self.map[randomInputLocation])
                    generalRandomInput.updateMethods = updateMethods
                    generalRandomInput.variableName = updateMethods.keys()
                    #self.randomInputVector.append(generalRandomInput)
        
        
        
        self.randomInputDimension = len(self.randomInputVector)
         
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
           
    def saveRealisationLog(self, evaluationLogFile, networkName, dataNumber, caseName): 
        '''
        method to save a log file with the realisations passed for each sample iteration
        '''
        
        logData = np.array([randomInput.updateLog for randomInput in self.randomInputs]).transpose()
        
        logfile = open(evaluationLogFile, "wb")
        logfile.write(''.join(['Stochastic simulation ',str(networkName),' DataNumber ',dataNumber,'\n','\n']))
        logfile.write(''.join(['uqsaCase :', caseName,'  number of evaluations: ',str(len(logData)),'\n','\n']))
        
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
        randomInputInfos = list(itertools.chain.from_iterable([randomInput.generateInfo() for randomInput in self.randomInputVector]))
        for info in randomInputInfos:
            randomInputManagerInfo.append(info)
        return randomInputManagerInfo
    
    def deleteAllRandomInputs(self):
        '''
        Function removes all existing random inputs
        '''
        for rI in self.randomInputs:
            del rI
        self.randomInputs = []
        self.randomInputVector = []
        self.map = {}
        self.randomInputDimension = 0
                