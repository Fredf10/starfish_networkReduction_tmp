import sys, os
import math
# set the path relative to THIS file not the executing file!
cur = os.path.dirname(os.path.realpath(__file__))
#sys.path.append(cur + '/NetworkLib')
class RuntimeMemoryManager(object):
    def __init__(self, nSaveBegin, nSaveEnd, nTsteps, network):
        self.nTSteps = nTsteps
        self.memoryArraySizeTime = None
        self.totalDataPointsPerTimeStep = 0
        self.registeredData = []
        
        self.chunkCount = 0
        
        self.network = network
        self.nDCurrent = network.nDCurrent 
        self.nSaveBegin= network.nSaveBegin
        self.nSaveEnd = network.nSaveEnd
        self.nSkipShift = network.nSkipShift
        self.nSaveSkip = network.nSaveSkip
        self.maxMemory = network.maxMemory
        
    def registerDataSize(self,dataSizes):
        """
        Args:
            dataSizes := a tuple with the number of values stored at a single time point
        """
        self.totalDataPointsPerTimeStep += sum(dataSizes)
        
        estimatedMemorySolutionDataSpace = self.totalDataPointsPerTimeStep*8
        self.memoryArraySizeTime = int(math.floor(self.maxMemory * 1024.*1024. / estimatedMemorySolutionDataSpace))
        # Don't allocate more memory than needed
        if self.memoryArraySizeTime > (self.nTSteps + 1):
            self.memoryArraySizeTime = self.nTSteps + 1
    
    def registerSimulationData(self, solutionMemory, dataBuffers):
        """
        Args:
            solutionMemory := a list of numpy array objects to be registered with the memory manager
            dataBuffers := a list corresponding to the solutinMemoryList with either a dataBuffer or
             the value None if the corresponding vector doesn't need to be saved.
        """
        self.registeredData.extend(zip(solutionMemory,dataBuffers))
        # Each simulation object registers itself here by passing in a tuple of spatial dimensions
        # The
    
    def registerSaveData(self, solutionMemory,  dataBuffers):
        pass
    
    def rollSimulationData(self):
        for solutionMemory, dataBuffer in self.registeredData:
            solutionMemory[0] = solutionMemory[-1]
        
    
    def flushSolutionMemory(self):
       
        """
        saving utility function to determine if solution data needs to be sent to the output file,
        and to calculate the correct indices between solution memory and the data output file.
       
        Explanation of index variables
        nCB,nCE where the beginning and end of the current solution data in memory would lie
         in a full time history of the solution. These are position indices
        nMB,nME what indices of the current memory mark the beginning and end of what should be saved
         nME is a slice index, i.e. position + 1
        nSB,nSE where the beginning and end of the current data to save would lie in the whole of the
         simulation. nSE is a slice index, i.e. position + 1
        nDB,nDE where does the current selection of data belong in the whole of the saved data
        """

        memoryArraySize = self.memoryArraySizeTime;
        offset = (memoryArraySize - 1) * self.chunkCount
        # indices mapping beginning and end of memory to the absolute number time steps in solution
        nCB = offset+1
        nCE = (memoryArraySize - 1) * (self.chunkCount + 1)# nCE == currentTimeStep+1
        nDB = None
        nDE = None
        nSB = None
        nSE = None
        # check if we need to save
        saving = not(nCE < self.nSaveBegin or nCB > self.nSaveEnd) # not(not saving)

        if saving:
            
            nMB = 1   + self.nSkipShift
            nSB = nCB + self.nSkipShift
            if self.nSaveBegin > nCB:
                nMB = 1 + (self.nSaveBegin-nCB)
                nSB = self.nSaveBegin
                
            #determine length to write
            # 1. assume we write out through the end of memory
            # 1.a Accounting for skipping write first value, then as many values as remain divisible by the skip
            nME = memoryArraySize
            nSE = nCE + 1 # Add one to get the last slice index... 
            
            
            
            # correct this if save index is less than the current time step
            if self.nSaveEnd < nCE:
                # set the index to end saving
                nME -= (nCE - self.nSaveEnd)
                nSE = self.nSaveEnd + 1
            
            numWrittenAfterFirst = ((nME - nMB) -1)//self.nSaveSkip
            lengthToWrite = 1 + numWrittenAfterFirst
            
            # How many to skip on the next round?
            self.nSkipShift = (self.nSaveSkip - (nME-nMB)%self.nSaveSkip)%self.nSaveSkip
            
            # Where in the data file to put the data    
            nDB = self.nDCurrent
            nDE = self.nDCurrent + lengthToWrite
            self.nDCurrent += lengthToWrite
            self.chunkCount += 1
        
        
        for solutionMemory, dataBuffer in self.registeredData:
            
            if saving and dataBuffer:
                dataBuffer[nDB:nDE] = solutionMemory[nMB:nME]

            solutionMemory[0] = solutionMemory[-1]
    
    
    
class DataHandler(object):
    
    def __init__(self, currentTimeStep, nTsteps, network, currentMemoryIndex, memoryArraySize):
        """
        Initialize DataHandler numerical object for solver
    
        Args:
           currentTimeStep    := the current position of the solver relative to the total solution
           nTsteps            := the total number of time steps required for the solution
           network            := the vascular network object
           currentMemoryIndex := the index in vessel memory buffer corresponding to currentTimeStep
           memoryArraySize    := the number of time steps to be stored in the runtime memory
        """
                
        self.nTsteps = nTsteps
        
        self.chunkCount = 0
        
        self.memoryArraySizeTime = memoryArraySize      
        
        self.currentTimeStep = currentTimeStep
        
        # the temporal position of beginning of the data currently stored in memory
        self.memoryOffset = [0]
        
        # the temporal position of the current time point in memory
        self.currentMemoryIndex = currentMemoryIndex
        
        self.network = network
    
    def emergencyFlush(self):
        # self.flushSolutionMemory(self.currentTimeStep[0], self.currentMemoryIndex[0],self.chunkCount)
        self.network.flushSolutionMemory( self.currentTimeStep[0], self.currentMemoryIndex[0],self.chunkCount)

    def __call__(self):
        '''
        call function for DataHandler to save the data for each vessel in the network
        '''
        currentMemoryIndex = self.currentMemoryIndex[0]
        currentTimeStep    = self.currentTimeStep[0]
        
        if currentMemoryIndex == self.network.memoryArraySizeTime - 2:
            
            ## after the memory is filled
            # initiate flush of memoryArrays
            self.network.flushSolutionMemory(currentTimeStep, currentMemoryIndex,self.chunkCount)
            #self.flushSolutionMemory(self.currentTimeStep[0], self.currentMemoryIndex[0],self.chunkCount)
            self.chunkCount += 1
            self.memoryOffset[0] = (self.memoryArraySizeTime - 1) * self.chunkCount

        elif currentTimeStep == self.nTsteps - 1:
            ## if the simulation is finished but memory not filled
            # initiate flush of memoryArrays
            self.network.flushSolutionMemory(currentTimeStep, currentMemoryIndex, self.chunkCount)
            # self.flushSolutionMemory(self.currentTimeStep[0], self.currentMemoryIndex[0],self.chunkCount)
            
        
        
        