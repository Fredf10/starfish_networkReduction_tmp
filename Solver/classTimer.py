import numpy as np


class TimerBaseClass(object):
    """
    Base class for a timer object.
    intended to implement "timed" events in the simulation, e.g. Valsalva maneuver, Hemorrhage, motion,...
    """
    def __init__(self):

        self.currentTimeStep = 0
        self.currentMemoryIndex = 0
        self.nTsteps = 0
        self.dt = 0
        self.Tstart = 0 # start time of the event
        self.Tend = 0 # end time of the event
        self.type = '' # type of the event, e.g. Valsalva
        
    def update(self,TimeDict):
        '''
        updates the data using a dictionary in form of 
        TimeDict = {'variableName': value}
        '''
        for key,value in TimeDict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except: 
                print 'ERROR communicator.update(): wrong key: %s, could not set up communicator' %key
        
    

class Valsalva(TimerBaseClass):
    
    """
    definition of a Valsalva maneuver, a clinical test during which the intrathoracic pressure is changed
    """
    

    def __init__(self,TimerDict):
        
        TimerBaseClass.__init__(self)
        
        self.VesselsToModify = 0 # Vessels on which the external pressure is changed
        self.deltaP = -30*133.32 # pressure difference
        self.InitialexternalPressure = 0 # the initial external pressure
        
        self.update(TimerDict) # updated method
        
        self.newexternalPressure = 0 # array for the new external pressure values
        self.calculateNewExternalPressure() # calculate the new external pressure values
        
        
        import matplotlib.pyplot as plt
        
        plt.plot(self.newexternalPressure)
        plt.show()
        
                                                          
    def calculateNewExternalPressure(self):
        """
        calculate the new external pressures given the time where the Valsalva maneuver starts and where it ends and the pressure change
        the function will be linear (like a pressure ramp)
        """
        start = int(self.Tstart/self.dt) # time step where pressure ramp starts 
        end = int(self.Tend/self.dt) # time step where pressure ramp ends
        
        startPressure = self.InitialexternalPressure * np.ones(start) # pressure at the beginning
        endPressure = (self.InitialexternalPressure+self.deltaP)*np.ones(self.nTsteps+1-end) # pressure after Valsalva
        
        changingPressure = np.linspace(0,self.deltaP,num = end-start) # pressure values during ramp
        self.newexternalPressure = np.append(startPressure,np.append(changingPressure,endPressure)) # concatenate arrays to have the whole pressue history in one array
        
    
    def __call__(self):
        
        """
        updates the external pressure for each timestep in the Compliance instance of the concerned vessels
        does work for Compliance models, which include the external pressure, e.g. Hayashi
        """
        
        n = self.currentTimeStep[0]
        
        for key,vessel in self.VesselsToModify.iteritems():
            
            vessel.compliance.update({"externalPressure": self.newexternalPressure[n+1]})


class Hemorrhage(TimerBaseClass):
    
    """
    class to implement a hemorrhagic event (blood loss)
    """
    
    def __init__(self):
        
        """
        constructor will depend on the implementation
        e.g. outflow terminal with resistance that is set to a small value or a sink term in the 
        system equations
        """
    