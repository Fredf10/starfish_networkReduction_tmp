
class QuantityOfInterest(object):
    '''
    
    '''
    def __init__(self,quantityName):
        
        self.quantityName = quantityName
        self.data = None

        # statistic properties
        self.gPCExpansion = None
        
        self.expectedValue       = None
        self.variance            = None
        
        self.conficenceIntervals = None
        self.confidenceAlpha     = 0.8
        
        self.partialSensitivities = None
        self.totalSensitivities   = None
        
    def asd(self):
        '''
        
        '''
        
    