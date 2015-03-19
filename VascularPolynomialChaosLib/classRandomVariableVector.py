

from classRandomVariable import RandomVariable


class RandomVariableVector(object):
    '''
    
    '''
    def __init__(self):        
        self.mapOfRandomVectors = {} # hash map with {randomVariable.location : randomVariable.randomVariableId}
        self.randomVector = [] # random variables as they stand in xml
        
        self.randomVectorCP = [] # random variables as they are used in vascular polynomial chaos
        
        self.jointDistributions = []
        
        
    def addRandomVariable(self,randomVariableLocation, randomVariableDict):
        '''
        Function which adds a random variable to randomVector
        '''
        
        randomVariable = RandomVariable()
        randomVariableId = len(self.randomVector)
        
        randomVariableDict['randomVariableId'] = randomVariableId
        randomVariableDict['location'] = randomVariableLocation
        randomVariable.update(randomVariableDict)
        
        self.randomVector.append(randomVariable)
        self.mapOfRandomVectors[randomVariableLocation] = randomVariableId
        
        
    def __call__(self, index = None):
        '''
        call function returns the random vector
        '''
        if index == None:
            return self.randomVector
        else:
            try:
                return self.randomVector[index]
            except:
                return []
        
if __name__ == '__main__':
    
    data1 = {'distributionType' : 'Uniform',
            'a'                : 2,  
            'b'                : 12,
            'variableName'     : 'radius',
            'updateMethod'     : None}
    
    name1 = "vessel_1"
    
    data2 = {'distributionType' : 'Uniform',
            'a'                : 2,  
            'b'                : 12,
            'variableName'     : 'beta',
            'updateMethod'     : None}
    
    name2 = "vessel_2"
    
    Z = RandomVariableVector()
    Z.addRandomVariable(name1, data1)
    Z.addRandomVariable(name2, data2)