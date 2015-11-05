

from testBaseClass import TestBaseClass 

class RandomInput(TestBaseClass):
    '''
    class decription of a random input of STARFiSh
    '''
    
    def __init__(self):
        
        self.variableName     = [] # list of variable names (e.g. betaHayashi, Z1) which are assoziated with random input
        ##        
        self.name             = None # name of the random input
        ##
        self.distributionType = None # distribution type: Normal(a,b), Uniform(a,b), another random input Z. a+bZ
        self.a                = 0 # 
        self.b                = 0 #
        ##
        self.updateMethods    = {} 
        self.updateLog        = [] # list to save the passed realisations
                
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
        
    def passRealisationToAssosiatedObj(self, input):
        '''
        methods update connected deterministic variable value from given
        'sample' drawn from random vector 
        
        input:
            sample:  <float> 'sample' drawn from random vector 
        '''
        sample_i = None      
        # if input comes from other random Input
        if isinstance(input,dict):
            sample_i = input[self.name]
        elif isinstance(input,float):        
            # if input comes from distribution handler
            sample_i = input        
        if sample_i != None:
            realisation = self.evaluateRealisationFromSample(sample_i)
            if self.updateMethods != {}:
                #for i,loc in enumerate(self.parameter.split('_')):
                #    if i == 0:
                #        print '{:3} | {:20} | {:21} | {:.4} '.format(self.name,self.variableName,loc, realisation)
                #    else:
                #        print '{:3} | {:20} | {:21} | '.format(' ',' ',loc)
                #print  
                for variableIdentifier,updateMethod in self.updateMethods.iteritems():                    
                    updateMethod({variableIdentifier : realisation})
                self.updateLog.append(realisation)
                    
    def generateInfo(self):
        '''
        Function to print random input information to console
        '''
        if self.distributionType in self.printOutDistributions:
            dist = self.printOutDistributions[self.distributionType]
        else:
            dist = self.distributionType
        
        randomInputInfo = []
        
        for i,loc in enumerate(self.parameter.split('_')):
            if i == 0:
                info = '{:3} | {:20} | {:21} | {:3} + {:3} {} '.format(self.name,self.variableName,loc, self.a, self.b, dist)
            else:
                info = '{:3} | {:20} | {:21} | '.format(' ',' ',loc)
            randomInputInfo.append(info)
        randomInputInfo.append('\n')
        return randomInputInfo
        
        
class ParametricRandomInput(RandomInput):
    
    externVariables      = {'parameter': TestBaseClass.ExtValue(str,strCases = ['anything']),
                            'distributionType':TestBaseClass.ExtValue(str,strCases = ['anything']),
                            'a':TestBaseClass.ExtValue(float),
                            'b':TestBaseClass.ExtValue(float),
                            } 
    
    externXmlAttributes  = []
    externXmlElements    = ['parameter',
                            'distributionType',
                            'a',
                            'b'
                            ]    
    
    def __init__(self):
        super(ParametricRandomInput,self).__init__()
        
        self.randomInputType = 'parametricRandomInput'
        self.parameter  = None
        
        
class GeneralRandomInput(RandomInput):
    
    externVariables   = {'distributionType':TestBaseClass.ExtValue(str,strCases = ['anything']),
                         'a':TestBaseClass.ExtValue(float),
                         'b':TestBaseClass.ExtValue(float),
                        } 
    
    externXmlAttributes  = []
    externXmlElements    = ['distributionType',
                            'a',
                            'b'
                            ]
    
    def __init__(self):
        super(GeneralRandomInput,self).__init__()
        
        self.randomInputType = 'generalRandomInput'
        self.parameter = '-'
    