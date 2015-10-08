
import chaospy as cp
import h5py


class DistributionManager(object):
    
    def __init__(self, randomInputVector = None):
        # distribution
        self.randomInputVector     = randomInputVector
        self.marginalDistributions = []
        self.jointDistribution     = None
        self.jointDistributionDependent = None
        self.distributionDimension = None
        self.dependentCase = False
        # samples
        self.expansionOrder  = 0
        self.samples         = None
        self.samplesDependent = None
        self.samplesSize     = 0
        self.sampleMethod    = None
        # orthogonalPolynomials
        self.orthogonalPolynomials = None
        
    def passRealisation(self, sampleIndex):
        '''
        Function to pass samples of the random variables to
        the random inputs
        '''
        if self.dependentCase == False:
            sample = self.samples[sampleIndex]
        else:
            sample = self.samplesDependent[sampleIndex]
            
        print sample         
        if len(sample) == len(self.randomInputVector):
            print "\nSample number {}".format(sampleIndex)
            print '{:3} | {:20} | {:21} | {}'.format("Id","variableName","location","realisation")
            print "--------------------------------------------------------------------"          
            for randomInput,sample_i in zip(self.randomInputVector,sample):
                print "random variable {} with realisation {}".format(randomInput.randomInputId,sample_i)
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
            except StandardError:
                print "WARNING DistributionManager.updateData Wrong key: {}, could not update varibale".format(self.randomInputId, key)
    
 
        
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
    
        
    def createDependentDistribution(self,CorrelationMatrix):
        '''
        Method creates a dependent distribution with Nataf transformation 
        based on a correlation matrix of the dependent random variables
        '''
        # TODO: check size of correlation and distributions
        self.dependentCase = True
        self.jointDistributionDependent = cp.Nataf(self.jointDistribution, CorrelationMatrix)
    
         
        