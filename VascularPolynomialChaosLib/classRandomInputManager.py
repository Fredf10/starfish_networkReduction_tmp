 
 


import numpy as np

import itertools

from testBaseClass import TestBaseClass 

import classRandomInput 

class RandomInputManager(TestBaseClass):
    '''
    
    '''
    
    externVariables      = {'randomInputs' : TestBaseClass.ExtDict('randomInput',TestBaseClass.ExtObject({'ParametricRandomInput':classRandomInput.ParametricRandomInput,
                                                                                                          'GeneralRandomInput':classRandomInput.GeneralRandomInput})),
                            'correlationMatrix': TestBaseClass.ExtValue([float], multiVar=True, strCases = ['anything'])
                           } 
    externXmlAttributes  = []
    externXmlElements    = ['randomInputs',
                            'correlationMatrix']
    
    def __init__(self):        
        
        
        self.randomInputs = {} # randomInput as they stand in xml
        
        self.correlationMatrix = None
        
        self.randomInputsExtDist = [] # randomInput which have a external distribution assosiated
        self.randomInputDimension = 0
                                                
    def initialize(self, vascularNetwork):
        '''
        Method to initialize the random inputs and
        defined them into the different types
        '''
        
        randomInputParameters = []
        
        for name,randomInput in self.randomInputs.iteritems():
            # set name            
            randomInput.name = name
            
            if randomInput.parameter not in randomInputParameters:
                randomInputParameters.append(randomInput.parameter)
            else:
                raise ValueError("assoziated parameter <{}> of randomInput {} is doubled defined!".format(randomInput.location,randomInput.name))
            
        self.linkRandomInputUpdateFunctions(vascularNetwork)
            
    def linkRandomInputUpdateFunctions(self, vascularNetwork):
        '''
        link update functions and create randomInputvector
        
        input vascularNetwork 
        
        only 1 level of nested self.randomInputs.values() is currently implemented
        
        parametricRandomInputs can be in the form of:
            a+b*Uniform(0,1) / a+b*Normal(0,1)
            or
            a+b*Z1 where Z1 is a general random input defined as 
        
        '''
        
        randomInputMap = {}
        for randomInputName,randomInput in self.randomInputs.iteritems():
            dist = randomInput.distributionType 
            
            if randomInput.randomInputType == 'parametricRandomInput':
                # check distribution
                loc = randomInput.parameter.split('_')
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
                
                else: raise ValueError("RandomInputManager: Parameter {} of randomInput {} is linkable.".format(randomInput.location,randomInputName))
                if randomInput.updateMethods == {}: break
                
            # find links between random inputs
            if dist in ['Normal','Uniform']:
                # append random input to random input vector
                self.randomInputsExtDist.append(randomInput)
            else:
                if dist not in randomInputMap:
                    randomInputMap[dist] = {randomInputName : randomInput.passRealisationToAssosiatedObj}
                else:
                    randomInputMap[dist][randomInputName] = randomInput.passRealisationToAssosiatedObj
        
        # link general random inputs 
        for dist,updateMethods in randomInputMap.iteritems():
            # find random input with dist:
            for randomInputName in self.randomInputs.iterkeys():
                if dist == randomInputName:
                    generalRandomInput = self.randomInputs[randomInputName]
                    generalRandomInput.updateMethods = updateMethods
                    generalRandomInput.variableName = updateMethods.keys()
                    print updateMethods.keys(), 
                    #self.randomInputsExtDist.append(generalRandomInput)
        
        self.randomInputDimension = len(self.randomInputsExtDist)
                    
    def saveRealisationLog(self, evaluationLogFile, networkName, dataNumber, caseName): 
        '''
        method to save a log file with the realisations passed for each sample iteration
        '''
        
        logData = np.array([randomInput.updateLog for randomInput in self.randomInputs.values()]).transpose()
        
        logfile = open(evaluationLogFile, "wb")
        logfile.write(''.join(['Stochastic simulation ',str(networkName),' DataNumber ',dataNumber,'\n','\n']))
        logfile.write(''.join(['uqsaCase :', caseName,'  number of evaluations: ',str(len(logData)),'\n','\n']))
        
        for info in self.generateInfo():
            logfile.write(''.join([info,'\n']))
        randomInputName = [randomInput.name for randomInput in self.randomInputs.values()]
        logfile.write(''.join(['\n{:<10}'.format(''), '{:20}'.format('RandomVariableId'),'\n','\n']))
        logfile.write(''.join(['{:<10}'.format('EvalNr'), ''.join(['{:20}'.format(j) for j in ['{:<5}'.format(k) for k in randomInputName]]),'\n',]))     
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
        randomInputInfos = list(itertools.chain.from_iterable([randomInput.generateInfo() for randomInput in self.randomInputs.values()]))
        for info in randomInputInfos:
            randomInputManagerInfo.append( info)
        randomInputManagerInfo.append( "\n Defined Random Variables\n")
        randomInputManagerInfo.append( '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","distribution"))
        randomInputManagerInfo.append( "-------------------------------------------------------------------- \n")
        randomInputInfos = list(itertools.chain.from_iterable([randomInput.generateInfo() for randomInput in self.randomInputsExtDist]))
        for info in randomInputInfos:
            randomInputManagerInfo.append(info)
        return randomInputManagerInfo
    
    def deleteAllRandomInputs(self):
        '''
        Function removes all existing random inputs
        '''
        for rI in self.randomInputsList:
            del rI
        self.randomInputsList = []
        self.randomInputsExtDist = []
        self.map = {}
        self.randomInputDimension = 0
                