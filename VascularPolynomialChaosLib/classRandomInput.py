

class RandomInput(object):
    '''
    class decription of a random input of STARFiSh
    '''
    def __init__(self, dataDict = None):
        ##
        self.type             = None # 'generalRandomInput' 'parametricRandomInput'
        ## meta information # location of the deterministic variable
        self.location         = None # location e.g. vessel_1 or Sinus2_1
        ##        
        self.name             = None # generalRandomInputName
        self.variableName     = [] # name of the deterministic variable
        self.randomInputId    = None # position in the random vector
        ##
        self.distributionType = None # distribution type: Normal(a,b), Uniform(a,b), another random input Z. a+bZ
        self.a                = 0 # 
        self.b                = 0 #
        ##
        self.updateMethods    = {} 
        self.updateLog        = [] # list to save the passed realisations
        
        if dataDict != None:
            self.update(dataDict)
        
        self.printOutDistributions = {'Uniform': 'U(0,1)',
                                      'Normal' : 'N(0,1)'}
        
    def evaluateRealisationFromSample(self,sample_i):
        '''
        calculate the realisation of the random variable defined as
        realisation =  a + b sample_i
        
        where the sample_i comes from the governing distribution: eg. Uniform, Normal or another RandomInput 
        '''
        realisation = self.a + self.b * sample_i
        return realisation
        
    def passRealisationToAssosiatedObj(self,input):
        '''
        methods update connected deterministic variable value from given
        'sample' drawn from random vector 
        
        input:
            sample:  <float> 'sample' drawn from random vector 
        '''
        sample_i = None      
        # if input comes from other random Input
        if isinstance(input,dict):
            sample_i = input[self.randomInputId]
        elif isinstance(input,float):        
            # if input comes from distribution handler
            sample_i = input        
        if sample_i != None:
            realisation = self.evaluateRealisationFromSample(sample_i)
            if self.updateMethods != {}:
                for i,loc in enumerate(self.location.split('_')):
                    if i == 0:
                        print '{:3} | {:20} | {:21} | {:.4} '.format(self.randomInputId,self.variableName,loc, realisation)
                    else:
                        print '{:3} | {:20} | {:21} | '.format(' ',' ',loc)
                print  
                for variableIdentifier,updateMethod in self.updateMethods.iteritems():                    
                    updateMethod({variableIdentifier : realisation})
                self.updateLog.append(realisation)
    
    
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
                print "ERROR RandomInput.updateData (randomInputId {}) Wrong key: {}, could not update varibale".format(self.randomInputId, key)

    def getVariableValue(self, variableName):
        '''
        Returns value of variable with name : variableName
        States Error if not such variable
        '''
        try:
            return self.__getattribute__(variableName)
        except: 
            print "ERROR RandomVariable.getVariable() : RandomVariable has no variable {}".format(variableName)
                
    def generateInfo(self):
        '''
        Function to print random input information to console
        '''
        if self.distributionType in self.printOutDistributions:
            dist = self.printOutDistributions[self.distributionType]
        else:
            dist = self.distributionType
        
        randomInputInfo = []
        for i,loc in enumerate(self.location.split('_')):
            if i == 0:
                info = '{:3} | {:20} | {:21} | {:3} + {:3} {} '.format(self.randomInputId,self.variableName,loc, self.a, self.b, dist)
            else:
                info = '{:3} | {:20} | {:21} | '.format(' ',' ',loc)
            randomInputInfo.append(info)
        randomInputInfo.append('\n')
        return randomInputInfo
        
        
        
