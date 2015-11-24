

from testBaseClass import TestBaseClass 

import numpy as np

class QuantityOfInterest(TestBaseClass):
    '''
    
    '''
    # in class definition
    ## pure variables
    variablesHdf5Memory = []
    
    ## dictionary with objects to load
    objectDictsHdf5Memory = ['uqsaMeasures']
        
    def __init__(self,quantityName, locationName, confidenceAlpha):
        
        self.queryLocation = locationName
        self.quantityName = quantityName
        self.confidenceAlpha = confidenceAlpha
        self.hdf5Group = None
        # trajectory object
        self.simulationTime = None # link to simulation time
        #
        self.uqsaMeasures = None
        
    
    def getData(self, sampleSize, abcSample):
        '''
        checks if data defined and returns it
        '''
        if abcSample == False:
            # check for trajectory stuff TODO: do it nicer
            if 'trajectoryData' in self.hdf5Group.keys():
                # get min max of all basis saved
                return self.hdf5Group['trajectoryData'][:sampleSize]
            
            elif 'data' in self.hdf5Group.keys():
                return self.hdf5Group['data'][:sampleSize]
        else:
            #TODO: implement abc sampling hash if monte carlo is used
            raise NotImplementedError("MC sensitivity, abc-sample hash case not implemented yet, exit")
        
    def hashDataForGivenBases(self, basis, sampleSize):
        '''
        
        '''
        self.hdf5Group.create_dataset('trajectoryData', (sampleSize,len(basis)) , dtype='float64')
        self.hdf5Group.create_dataset('trajectoryBasis', data = basis, dtype='float64')
        
        import sys
        write = sys.stdout.write
        loadingBarElementCount = 1
        nElements = 50
        spaces = ''.join([' ' for i in xrange(nElements)])
        loadingBar = ''.join(['[',spaces,']'])
        write(loadingBar)
        sys.stdout.flush()
        backspacing = ''.join(['\b' for i in xrange(nElements+1)])
        write(backspacing)
        
        for n in xrange(sampleSize):
            
            data      = self.hdf5Group['data'][n]
            dataBasis = self.hdf5Group['dataBasis'][n]
            
            self.hdf5Group['trajectoryData'][n] = np.interp(basis, dataBasis, data)
        
            if divmod(n,sampleSize/nElements)[0] == loadingBarElementCount:
                loadingBarElementCount = loadingBarElementCount+1
                write("#")
                sys.stdout.flush()
        
        write("\n")
        
    def addUqsaMeasures(self,uqsaMeasureName,uqsaMeasure):
        '''
        Method append a uqsaMeasure class to the dictionary: self.uqsaMeasures
        '''
        
        if self.uqsaMeasures == None:
            self.uqsaMeasures = {}
        if uqsaMeasureName not in self.uqsaMeasures:
            self.uqsaMeasures[uqsaMeasureName] = uqsaMeasure
        else:
            print " did not add {} as it already exist in dict".format(uqsaMeasureName)
    
    
    def retainDsetRowData(self,dsetName, dataRow, sampleIndex, totalDataShape):
        '''
        
        '''
        if dsetName not in self.hdf5Group.keys():
            self.hdf5Group.create_dataset(dsetName, totalDataShape, dtype='float64')
        self.hdf5Group[dsetName][sampleIndex] = dataRow
        
    
        