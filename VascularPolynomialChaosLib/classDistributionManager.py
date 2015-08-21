
import moduleFilePathHandlerVPC as mFPH_VPC

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
        # samples
        self.expansionOrder  = 0
        self.samples         = None
        self.samplesDependent = None
        self.samplesSize     = 0
        self.sampleMethod    = None
        # orthogonalPolynomials
        self.orthogonalPolynomials = None
        
    def passRealisation(self, sampleIndex, dependentCase = False):
        '''
        Function to pass samples of the random variables to
        the random inputs
        '''
        if dependentCase == False:
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
    
    def loadSamples(self, networkName, dataNumber, gPCEmethod, gPCEorder):
        '''
        load the current sample to disc so it is available for postprocessing or
        sequencielle working process
        for generation gPCE the sample nodes corresponding to the data are needed.
        '''
        vpcSampleFile = mFPH_VPC.getFilePath('vpcSampleFile', networkName, dataNumber, mode = "read", gPCEmethod=gPCEmethod, gPCEorder=gPCEorder)
        f = h5py.File(vpcSampleFile,'r')
        dset = f['sampleSpace']
        self.samples = dset[:]
        self.samplesSize    = dset.attrs.get('samplesSize')
        self.sampleMethod   = dset.attrs.get('sampleMethod')
        self.expansionOrder = dset.attrs.get('expansionOrder')
        
        if 'sampleSpaceDependent' in f.keys():
            dset = f['sampleSpaceDependent']
            self.samplesDependent = dset[:]
            
        f.close()
        # check if data is accordingly to the case
        if self.sampleMethod != gPCEmethod:
            raise ValueError("loadSamples for {} - {}: wrong '{}' saved in hdf5 file: {} (file) != {} (vpc)".fromat(networkName, dataNumber,'sampleMethod',self.sampleMethod,gPCEmethod))
        if self.expansionOrder != gPCEorder:
            raise ValueError("loadSamples for {} - {}: wrong '{}' saved in hdf5 file: {} (file) != {} (vpc)".fromat(networkName, dataNumber,'expansionOrder',self.expansionOrder,gPCEorder))
        if self.samplesSize != len(self.samples):
            raise ValueError("loadSamples for {} - {}: wrong '{}' saved in hdf5 file: {} (file) != {} (vpc)".fromat(networkName, dataNumber,'samplesSize',self.samplesSize,len(self.samples)))
      
      
    def saveSamples(self, networkName, dataNumber, gPCEmethod, gPCEorder):
        '''
        save the current sample to disc so it is available for postprocessing or
        sequencielle working process
        for generation gPCE the sample nodes corresponding to the data are needed.
        '''
        vpcSampleFile = mFPH_VPC.getFilePath('vpcSampleFile', networkName, dataNumber, mode = "write", gPCEmethod=gPCEmethod, gPCEorder=gPCEorder)
        f = h5py.File(vpcSampleFile,'w')
        dset = f.create_dataset("sampleSpace", data=self.samples)
        dset.attrs.create('samplesSize', data=self.samplesSize)
        dset.attrs.create('sampleMethod', data=self.sampleMethod)
        dset.attrs.create('expansionOrder', data=self.expansionOrder)
        if self.samplesDependent != None:
            dset = f.create_dataset("sampleSpaceDependent", data=self.samplesDependent)
        
        f.flush()
        f.close()
        
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
    
        
    def createSamples(self,networkName, dataNumber, sampleMethod, sampleSize = 1, expansionOrder = None):
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
        
        # create dependent samples if correlation exists
        if self.jointDistributionDependent != None:
            self.createDependentSamples()
    
    def createDependentDistribution(self,CorrelationMatrix):
        '''
        Method creates a dependent distribution with Nataf transformation 
        based on a correlation matrix of the dependent random variables
        '''
        self.jointDistributionDependent = cp.Nataf(self.jointDistribution, CorrelationMatrix)
        
    def createDependentSamples(self):
        '''
                
        '''
        self.samplesDependent = self.jointDistributionDependent.inv(self.jointDistribution.fwd(self.samples))
                
        
    def calculateOrthogonalPolynomials(self):
        '''
        Method to calculate orthogonal polynomials
        
        TODO: saving/loading og orthogonalPolynomials
        '''
        self.orthogonalPolynomials = cp.orth_ttr(self.expansionOrder,self.jointDistribution)
        
        #self.orthogonalPolynomils = pc.orth_gs(self.expansionOrder,self.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_chol(self.expansionOrder,self.jointDistribution)
        #self.orthogonalPolynomils = pc.orth_svd(self.expansionOrder,self.jointDistribution)
         
        