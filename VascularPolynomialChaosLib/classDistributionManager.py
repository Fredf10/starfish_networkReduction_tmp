

import chaospy as cp



class DistributionManager(object):
    
    def __init__(self, randomInputVector = None):
        # distribution
        self.randomInputVector     = randomInputVector
        self.marginalDistributions = []
        self.jointDistribution     = None
        self.distributionDimension = None
        # samples
        self.expansionOrder = 0
        self.samples         = None
        self.samplesSize     = 0
        self.sampleMethod   = None
        
    def passRealisation(self, sampleIndex):
        '''
        Function to pass samples of the random variables to
        the random inputs
        '''
        sample = self.samples[sampleIndex]
        if len(sample) == len(self.randomInputVector):
            for randomInput,sample_i in zip(self.randomInputVector,sample):
                randomInput.passRealisationToAssosiatedObj(sample_i)
    
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
    
    def loadSamples(self):
        '''
        
        '''
        print "loadSampleSpace not jet implemented"
#         print "number of simulations",len(samples)
#         #save file
#         saveFile = open(sampleFile,"wb")       
#         cPickle.dump(samples,saveFile,protocol=2)
#         saveFile.close()
#         print ".. done"
#     
    def saveSamples(self):
        '''
        
        '''
        print "saveSampleSpace not jet implemented"
#         try:
#             print " load samples space "
#             loadFile = open(sampleFile,"rb")
#             # load pickle
#             samples = cPickle.load(loadFile)
#             loadFile.close()
#             print ".. done"
#         except:
#             print 'File does not exits:'
#             print sampleFile
#             exit()
              
    ## methods created by the toolbox-child class implementation
    def createRandomVariables(self):
        '''
        create a random variable vector from
        for the random input variables in the random input variable vector
        and the joint distribution 
        '''
        pass
    
    def createSamples(self):
        '''
        
        '''
        pass
    
    
                     
class DistributionManagerChaospy(DistributionManager):
    
    def __init__(self, randomInputVector = None):
        super(DistributionManagerChaospy, self).__init__(randomInputVector)
        
        
    def createRandomVariables(self):
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
        self.distributionDimension = len(self.jointDistribution)
    
        
    def createSamples(self, sampleMethod, sampleSize = 1, expansionOrder = None):
        '''
        create samples for the defined distribution for given samplesSize and sampleMethod
        using the chaospy toolbox
        Input
            samplesSize : int,array_like
        
            sampleMethod : str
                (from chaospy)
                Alternative sampling techniques
            
                Normal sampling schemes
                Key     Name                Nested
                ----    ----------------    ------
                "K"     Korobov             no
                "R"     (Pseudo-)Random     no
                "L"     Latin hypercube     no
                "S"     Sobol               yes
                "H"     Halton              yes
                "M"     Hammersley          yes
            
                Grided sampling schemes
                Key     Name                Nested
                ----    ----------------    ------
                "C"     Chebyshev nodes     maybe
                "G"     Gaussian quadrature no
                "E"     Gauss-Legendre      no
            
            expansionOrder: float
                calculate optimal samplesSize for gPC with the following rule
                    samplesSize =  2* number gPC-expansion terms
            
        '''
        self.sampleMethod   = sampleMethod
        
        self.samplesSize = sampleSize
        # calculate samplesSize from expansion order if given
        if expansionOrder != None:
            self.expansionOrder = expansionOrder
            self.samplesSize = 2*cp.terms(expansionOrder,self.distributionDimension)
        self.samples = self.jointDistribution.sample(self.samplesSize,sampleMethod).transpose()   
             
        if self.distributionDimension == 1:
            self.samples = self.samples.reshape(self.samplesSize,1)
                