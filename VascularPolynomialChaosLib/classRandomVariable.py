


class RandomVariable(object):
    '''
    class decription of a random variable of STARFiSh
    '''
    def __init__(self, dataDict = None):
        ## meta information # location of the deterministic variable
        self.location         = None # location e.g. vessel_1 or Sinus2_1
        ##        
        self.variableName     = None # name of the deterministic variable
        self.randomVariableId = None # position in the random vector
        ##
        self.distributionType = None # distribution type: Normal(a,b), Uniform(a,b), another random variable Z. a+bZ
        self.a                = None # 
        self.b                = None #
        self.cpDist           = None # chaospy distribution
        ##
        self.updateMethod     = None 
        self.sample           = None # realisation of the random sample
        
        if dataDict != None:
            self.update(dataDict)
        
    def setSample(self,sample):
        '''
        methods update connected deterministic variable value from given
        'sample' drawn from random vector 
        
        input:
            sample:  <float> 'sample' drawn from random vector 
        '''
        self.sample = sample
        if sample != None:
            if self.updateMethod != None:
                self.updateMethod({self.variableName:sample})
    
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
                print "ERROR RandomVariable.updateData (randomVariableId {}) Wrong key: {}, could not update varibale".format(self.randomVariableId, key)

    def getVariableValue(self, variableName):
            '''
            Returns value of variable with name : variableName
            States Error if not such variable
            '''
            try:
                return self.__getattribute__(variableName)
            except: 
                print "ERROR RandomVariable.getVariable() : RandomVariable has no variable {}".format(variableName)
        

if __name__ == '__main__':
    
    data = {'distributionType' : 'Uniform',
            'a'                : 2,  
            'b'                : 12,
            'variableName'     : 'beta',
            'updateMethod'     : None}
    
    rv = RandomVariable(data)
    