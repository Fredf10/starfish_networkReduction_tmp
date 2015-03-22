

from classRandomInput import RandomInput


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
           
        
    def printOutInfo(self):
        '''
        Function to print out random variable informations
        '''
        print "\n Defined Random Inputs\n"
        
        print '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","distribution")
        print "---------------------------------------------------------- \n"
        for randomInput in self.randomInputs:
            randomInput.printOutInfo()
            
        print "\n Defined Random Variables\n"
        print '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","distribution")
        print "---------------------------------------------------------- \n"
        for randomInput in self.randomInputVector:
            randomInput.printOutInfo()