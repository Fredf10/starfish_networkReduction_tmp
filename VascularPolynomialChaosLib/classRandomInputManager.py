

from classRandomInput import RandomInput


class RandomInputManager(object):
    '''
    
    '''
    def __init__(self):        
        self.map = {} # hash map with {randomInput.location : randomInput.randomInputId}
        self.randomInputs = [] # randomInput as they stand in xml
        
        self.randomInputVector = [] # randomInput which have a external distribution assosiated
        
        
    def addRandomInput(self,randomInputLocation, randomInputDict):
        '''
        Function which adds a random variable to randomInputs
        '''
        
        if randomInputLocation not in self.map.keys():
            randomInput = RandomInput()
            randomInputId = len(self.randomInputs)
            
            randomInputDict['randomInputId'] = randomInputId
            randomInputDict['location'] = randomInputLocation
            randomInput.update(randomInputDict)
            
            self.randomInputs.append(randomInput)
            self.map[randomInputLocation] = randomInputId
        
                    
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
        
    def printOutInfo(self):
        '''
        Function to print out random variable informations
        '''
        
        print "\n Defined Random Inputs\n"
        
        print '{:3} | {:20} | {:20} | {}'.format("Id","variableName","location","distribution")
        print "---------------------------------------------------------- \n"
        for randomInput in self.randomInputs:
            randomInput.printOutInfo()