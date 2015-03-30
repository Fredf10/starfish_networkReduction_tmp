

import chaospy as cp



class DistributionManager(object):
    
    def __init__(self, randomInputVector = None):
        
        self.randomInputVector     = randomInputVector
        self.marginalDistributions = []
        self.jointDistribution     = None
                
        self.toolboxType = 'chaospy'
                                
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
                print "ERROR DistributionManager.updateData Wrong key: {}, could not update varibale".format(self.randomInputId, key)
    
    def createRandomVariables(self):
        '''
        Alias for the outside world
        '''
        if  self.toolboxType == 'chaospy':
            self.createRandomVariablesChaospy()
    
    def setToolbox(self,toolbox):
        '''
        method to choose a toolbox for creating distributions
        
        currently only chaospy toolbox is supported
        '''
        self.toolboxType = 'chaospy'
                           
    def createRandomVariablesChaospy(self):
        '''
        create a random variable vector from
        for the random input variables in the random input variable vector
        and the joint distribution 
        '''
        if self.randomInputVector == None:
            print "WARNING: DistributionManager.createRandomVariablesChaospy() no randomInputVector are defined"
            return
        
        evalDistDict = {'Uniform': cp.Uniform, 
                        'Normal' : cp.Normal}
        
        # create marignal distributions
        for randomInput in self.randomInputVector:
            distType = randomInput.distributionType
            if distType in evalDistDict.keys():
                marginalDistribution = evalDistDict[distType]() #eval(distType,{"__builtins__":None},evalDistDict)()
                self.marginalDistributions.append(marginalDistribution)
            #self.jointDistribution = cp.J(self.jointDistribution, marginalDistribution)
        
        # create joint distributions
        self.jointDistribution = cp.J(*self.marginalDistributions)
        
    def passRealisation(self, sample):
        '''
        Function to pass sample of the random variables to
        the random inputs
        '''
        if len(sample) == len(self.randomInputVector):
            for randomInput,sample_i in zip(self.randomInputVector,sample):
                randomInput.passRealisationToAssosiatedObj(sample_i)
        