import sys, os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur + '/NetworkLib')
from classVessel import Vessel

class DataHandler(object):
    
    def __init__(self, currentTimeStep, nTsteps, vessels, currentMemoryIndex, memoryArraySize):
        '''
        Initialize DataHandler numerical object for solver
    
        Input:
           currentTimeStep    := the current position of the solver relative to the total solution
           nTsteps            := the total number of time steps required for the solution
           vessels            := a dictionary of all the vessel composing the vascular network
           currentMemoryIndex := the index in vessel memory buffer corresponding to currentTimeStep
           memoryArraySize    := the number of time steps to be stored in the runtime memory
        '''
        
        self.nTsteps = nTsteps
        
        self.chunkCount = 0
        
        self.memoryArraySizeTime = memoryArraySize      
        
        self.currentTimeStep = currentTimeStep
        
        # the temporal position of beginning of the data currently stored in memory
        self.memoryOffset = [0]
        
        # the temporal position of the current time point in memory
        self.currentMemoryIndex = currentMemoryIndex
        
        self.vessels = vessels
        
    def __call__(self):
        '''
        call function for DataHandler to save the data for each vessel in the network
        '''
        currentMemoryIndex = self.currentMemoryIndex[0]
        currentTimeStep    = self.currentTimeStep[0]
        
        if currentMemoryIndex == self.memoryArraySizeTime - 2:
            # ## after the memory is filled
            # # initiate flush of memoryArrays
            for vessel in self.vessels.itervalues():
                vessel.flushMemory(self.chunkCount, self.memoryOffset[0])
                
            self.chunkCount += 1
            self.memoryOffset[0] = (self.memoryArraySizeTime - 1) * self.chunkCount

        elif currentTimeStep == self.nTsteps - 1:
            #### if the simulation is finished but memory not filled
            # # initiate flush of memoryArrays
            for vessel in self.vessels.itervalues():
                vessel.flushMemory(self.chunkCount, self.memoryOffset[0])
