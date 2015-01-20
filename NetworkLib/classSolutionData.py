import h5py
import numpy as np

class SolutionDataVessel(object):
    '''
    class for vessel solution data handling:
    
    Pressure
    Flow
    Area
    '''
    def __init__(self):
    
        self.P = None
        self.Q = None
        self.A = None
        
        self.nDCurrent = None
        
        self.save = True
        self.memoryArraySizeTime = None
        
        self.dsetGroup = None
        
        self.dsetP = None
        self.dsetQ = None
        self.dsetA = None
        
        self.nSaveBegin = None
        self.nSaveEnd = None
                
    def allocateMemory(self, memoryArraySizeTime, dsetGroup, nSaveBegin, nSaveEnd, nTsteps, numberOfGridPoints):
        '''
        allocates memory for the solution data to store data of a simulation
        Input:
            initialSolution     := intial Solution 
            memoryArraySize     := size of one array in memory
            dsetGroup           := data set group in the hdf5 file
            nSaveBegin          := time index when to start saving
            nSaveEnd            := time index when to stop saving
            numberOfGridPoints  := number of spacial grid points of the vessel
        '''
        
        self.P = np.ones((memoryArraySizeTime,numberOfGridPoints))
        self.Q = np.zeros((memoryArraySizeTime,numberOfGridPoints))
        self.A = np.zeros((memoryArraySizeTime,numberOfGridPoints))
        
        self.memoryArraySizeTime = memoryArraySizeTime
                
        if self.save == True:
            
            self.nSaveBegin = nSaveBegin
            self.nSaveEnd   = nSaveEnd 
        
            savedArraySize = nSaveEnd-nSaveBegin+1
            
            self.dsetP = dsetGroup.create_dataset("Pressure", (savedArraySize,numberOfGridPoints), dtype='float64')
            self.dsetQ = dsetGroup.create_dataset("Flow", (savedArraySize,numberOfGridPoints), dtype='float64')
            self.dsetA = dsetGroup.create_dataset("Area", (savedArraySize,numberOfGridPoints), dtype='float64')
            
            self.nDCurrent = 0
            
    def saveInitialState(self):
        '''
        saves the initial values (first time step)
        '''
        # check if everything from start should be saved 
        if self.nSaveBegin == 0 and self.save == True:
            
            self.dsetP[0] = self.P[0] 
            self.dsetQ[0] = self.Q[0]
            self.dsetA[0] = self.A[0]
        
            self.nDCurrent = 1
        
    def linkDataSets(self, dsetGroup):
        '''
        Links the dataset files of a loaded file to the appropriate pointer        
        '''
        self.dsetP = dsetGroup['Pressure']
        self.dsetQ = dsetGroup['Flow']
        self.dsetA = dsetGroup['Area']
        
    def loadDataInMemory(self):
        '''
        Loads everything into the memory which is saved in the dataset-files
        '''
        self.P = self.dsetP[:]
        self.Q = self.dsetQ[:]
        self.A = self.dsetA[:]        
        
    def flushMemory(self, chunkCount, offset):
        '''
        save a data from the memory chunk in the dset's of the hdf5 file
        
        Input variables:
        
        chunkCount := current chunk number to save
        offset := offset in the time point counter
        '''
        if self.save == True:
            # stored data saving indices
            nCB = offset+1
            nCE = offset+self.memoryArraySizeTime-1 # nCE == currentTimeStep+1
            
            # check if we need to save
            if not(nCE < self.nSaveBegin or nCB > self.nSaveEnd): # if not(not saving)
                
                ## memory indices
                nMB = 1
                if (self.nSaveBegin-nCB)>0:
                    nMB = 1+(self.nSaveBegin-nCB)
                
                #length to write
                lengthToWrite = self.memoryArraySizeTime - nMB 
                    
                nME = None
                if (self.nSaveEnd-nCE)<0:
                    nME = (self.nSaveEnd-nCE)
                    lengthToWrite += nME # -(-nME) as nME is negative
                                   
                nDB = self.nDCurrent
                nDE = self.nDCurrent+lengthToWrite
                                            
                self.dsetP[nDB:nDE] = self.P[nMB:nME]
                self.dsetQ[nDB:nDE] = self.Q[nMB:nME]
                self.dsetA[nDB:nDE] = self.A[nMB:nME]
                
                self.nDCurrent += lengthToWrite
                
            # re initialize
            self.P[0] = self.P[-1]
            self.Q[0] = self.Q[-1]
            self.A[0] = self.A[-1]