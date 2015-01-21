import numpy as np 
#from numpy.linalg import solve
#from scipy.optimize import fsolve
#from numpy.linalg import inv

import copy

import pprint
import math
import sys,os


# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')

class Baroreceptor(object):
    """
    Mother class for all baroreceptor models
    """
    
    def __init__(self, BaroDict):
        '''
        Baroreceptor model initialisation
        
        '''
        #System and Vessel Variables
        self.dt = 0
        self.n = 0 # current time step
        self.Tsteps = 0
        self.data = {}
        self.modelName = ''
        
        # the boundary conditions of type 1 and type 2
        self.boundaryCondition = 0
        self.boundaryConditionII = 0
        
        # Model from CellML or hardcoded
        self.cellMLBaroreceptorModel = False
        
          
        # update Time for proximal boundary
        self.oldUpdateTime = 0
        self.newUpdateTime = self.boundaryCondition.Tperiod/self.dt
        
        #initialize the CellML Baroreceptor model
        if self.cellMLBaroreceptorModel == True:
            
            # import the model from the file which defines it
            cellMLimport = "import" + " " + self.modelName + " " + " as " + " baroreceptorCellML "
            exec(cellMLimport)
            
            # input to te CellML model and output from it
            self.CellMLinputID = baroreceptorCellML.inputID
            self.CellMLoutputArray = baroreceptorCellML.outputArray
            self.CellMLoutputID = baroreceptorCellML.outputID
            
             
            (iniStates, self.constants) = baroreceptorCellML.initConsts()
            timeArray = np.linspace(0,self.dt,2)
            self.voi, self.states, self.algebraic = baroreceptorCellML.solver2(timeArray,iniStates,self.constants)
             
               
        else:
            print "Error: No CellML Baroreceptor Model provided!"

            
    ### solving of imported CellML model
    def solveCellML(self, input):
        """
        solve CellML model
        
        """
        nbElements = np.shape(input)[0]            
        #timeArray = np.linspace(0,(nbElements-2)*self.dt,(nbElements-1))
        
        #for it in xrange(0,(nbElements-1)):
        for it in xrange(0,(nbElements)):
                
            self.constants[self.cellMLinputID] = input[it]
            self.voi, self.states, self.algebraic = baroreceptorCellML.solver2([0,self.dt],self.states[-1][:],self.constants)
        
        
        return self.voi, self.states, self.algebraic
    
    
    def calcAndupdatePeriodTypeIcellML(self,n,input):
         """
         Function to calculate new period
         Calls self.solveCellML()
         is used by self.__call__()
         """
         
         if n == (round(self.newUpdateTime-2)):
            
            # solve the cellML system using the function defined above
            self.voi, self.states, self.algebraic = self.solveCellML(input)
            
            Tperiod  = self.algebraic[-1][self.CellMLoutputID]
            #print "BR128 - TPeriod"
            #print self.algebraic[-1][self.CellMLoutputID]
            
            self.boundaryCondition.updatePeriodRuntime(Tperiod,(self.newUpdateTime-2)*self.dt)
            
            self.oldUpdateTime = self.newUpdateTime
            self.newUpdateTime = self.oldUpdateTime + self.boundaryCondition.Tperiod/self.dt
        
            
    def update(self,baroDict):
            '''
            updates the updateBaroreceptorDict data using a dictionary in form of 
            baroDict = {'variableName': value}
            '''
            for key,value in baroDict.iteritems():
                try:
                    self.__getattribute__(key)
                    self.__setattr__(key,value)
                except: 
                    print 'ERROR baroreceptor.update(): wrong key: %s, could not set up baroreceptor' %key
        

           
class AorticBaroreceptor(Baroreceptor):
    
    '''
    for models of the AorticBaroreceptors
    Aortic Baroreceptor models with strain input and period of the heart cycle as output
    '''
    
    def __init__(self,BaroDict):
        
        # intialize with mother class constructor
        super().__init__(self)
        self.VesselID = 0
        
        # boundary condition that will be adjusted (type 1 condition)
        self.boundaryCondition  = 0
        
        # initial Compliance information for estimation of unstretched radius
        self.initialCompliance1 = 0
        self.initialCompliance2 = 0
        
        self.Ao_1 = self.data['Area1'][0] 
        self.Ao_2 = self.data['Area2'][0]
        
        self.Po_1 = self.data['Pressure1'][0]
        self.Po_2 = self.data['Pressure2'][0]
        
        
        # update with data from dictionary
        self.update(BaroDict)
        
        #initial area and unstretched radius
        self.Ro_1 = np.zeros(np.size['Area1'][0])
        self.Ro_2 = np.zeros(np.size['Area2'][0])
        
        # estimate the unstretched Radii R0 for the calculation of the strain
        self.estimateUnstretchedRadius(self)
        
        # update Time for the period of the inflow function
        self.oldUpdateTime = 0
        self.newUpdateTime = self.boundaryCondition.Tperiod/self.dt
        
            
    def estimateUnstretchedRadius(self):
        """
        Function to estimate the unstretched radius of the Vessels of the Aortic Arch for the calculation of strain.
        Takes information about compliance and initial vessel geometry.
        Takes the initial compliance and assumes it as constant, uses the initial Area to calculate 
        the Area at zero pressure, where the vessel is in its unstressed state, so that the true
        strain can be calculated.
        """
        
        compliance1 = self.initialCompliance1.C(self.Po_1)
        compliance2 = self.initialCompliance2.C(self.Po_2)
        
        deltaA1 = compliance1 * self.Po_1
        deltaA2 = compliance2 * self.Po_2
        
        self.Ro_1 = np.power((self.Ao_1 - deltaA1)/math.pi,0.5)
        self.Ro_2 = np.power((self.Ao_2 - deltaA2)/math.pi,0.5)
        
    
       
    def __call__(self):
        
        n = self.n[0]
        
        ## read area, calculate radius and strain and save to data dictionary
        A1 = self.data['Area1'][n]
        R1 = np.power(A1/math.pi,0.5)
        
        A2 = self.data['Area2'][n]
        R2 = np.power(A2/math.pi,0.5)
        
        epsilon1 = (R1 - self.Ro_1)/self.Ro_1
        epsilon2 = (R2 - self.Ro_2)/self.Ro_2
        
        # concatenate the strain arrays of the two vessels of the Aortic Arch
        epsilon = np.concatenate(epsilon1,epsilon2)
        epsMean = np.mean(epsilon)
        
        # update the data dictionary with the strain and the mean strain
        self.data['Strain'][n+1]  = epsilon
        self.data['MStrain'][n+1] = epsMean
        
        
        # use CellML model to calculate new period, the method calcAndupdatePeriodcellML
        # does this on a beat-to-beat basis
        if self.cellMLBaroreceptorModel == True:
            
            self.calcAndupdatePeriodTypeIcellML(self,n,self.data['MStrain'][self.oldUpdateTime:(self.newUpdateTime-2)])
            
            # update the data dictionary with current period
            if self.CellMLoutputArray == 'algebraic':
                
                self.data['T'][n+1] = self.algebraic[-1][self.CellMLoutputID]
            
            elif self.CellMLoutputArray == 'states':
                
                self.data['T'][n+1] = self.states[-1][self.CellMLoutputID]
            
        else:
             print 'Error: currently no hardcoded models for the Aortic BR available'
            
        
        

   
class CarotidBaroreceptor(object):

    '''
    for models of the Carotid Baroreceptors
    '''
    
    def __init__(self,BaroDict):
        
        """
        constructor for the CarotidBaroreceptor
        """
        #System and Vessel Variables
        super().__init(self) # from mother class
        self.vesselIdLeft = 0
        self.VesselIDright = 0
        self.terminalBoundaries = 0
        
        # get data from dictionary
        self.update(BaroDict)
        
        
        #######################################
        """
        Parameters for the Ursino model of the carotid baroreflex
        """
        #######################################
        
        # model parameters afferent part of the Ursino model - Ursino 1999
        
        self.pn = 12266.
        self.ka = 1567.
        self.fmin = 2.51
        self.tau_z= 6.37
        self.fmax = 47.78
        self.tau_p= 2.076
        
        # model parameters efferent part of the Ursino model  - Ursino 1999
        
        self.fe_inf = 2.10
        self.fe_0 = 16.11
        self.ke = 0.0675
        self.fe_min = 2.66
        
        # model parameters for the TPR effector part of the Ursino model  - Ursino 1999
        
        self.cR = 42263190800.
        self.tauR = 6.
        self.DR = 2.0
        #self.R0 = 81326644700.
    
        # model parameters for the EmaxLV effector part of the Ursino model  - Ursino 1999
        
        self.cE = 0.475 * 133.322368
        self.tauE = 8.
        self.DE = 2.0
        #self.E0 = 2.392*133.322368
        
        # model parameters for the Vusv effector part of the Ursino model  - Ursino 1999
        
        self.cVusv = -0.625*10e-3
        self.tauVusv = 20.0
        self.DVusv = 3.0
        #self.Vusv0 = 3.213e-3
        
        
        
        # states of the Ursino model
        
        self.PtildLeft = np.zeros(self.Tsteps)
        self.PtildLeft[0] = 12266. # initial value of the dynamic block (has the dimension of a pressure)
        
        self.PtildRight = np.zeros(self.Tsteps)
        self.PtildRight[0] = 12266.
        
        self.F_cs_left = np.zeros(self.Tsteps)
        self.F_cs_right = np.zeros(self.Tsteps)
        
        self.F_cs = np.zeros(self.Tsteps)
        self.F_efferent = np.zeros(self.Tsteps)
        
        self.delta_TPR = np.zeros(self.Tsteps)
        self.delta_Emax = np.zeros(self.Tsteps)
        self.delta_Vusv = np.zeros(self.Tsteps)
        
    #############################################################################
    
    """
    methods defining the different parts of the Ursino Baroreceptor model
    - Afferent
    - Efferent
    - Effector
    
    """
    
    def UrsinoAfferent(self,n):
        """
        Afferent part of the Baroreceptor model as defined by Ursino 1999
        The firing rate of the left and right carotid sinus are averaged to give a global firing rate
        the firing rate of the carotid sinus will be processed by the efferent part
        """
        
        pLeft = np.mean(self.data['pressureLeft'][n])
        pRight = np.mean(self.data['pressureRight'][n])
        
        if n == 0:
            
            dpLeft =  0.
            dpRight = 0.
            
        else:
            
            dpLeft = (pLeft - np.mean(self.data['pressureLeft'][n-1]))/self.dt
            dpRight = (pRight - np.mean(self.data['pressureRight'][n-1]))/self.dt
                
        self.PtildLeft[n+1] = (pLeft + self.tau_z*dpLeft-self.PtildLeft[n])*self.dt/self.tau_p+self.PtildLeft[n]
        self.PtildRight[n+1] = (pRight + self.tau_z*dpRight-self.PtildRight[n])*self.dt/self.tau_p+self.PtildRight[n]
        
        self.F_cs_left[n+1] = (self.fmin + self.fmax*math.exp((self.PtildLeft[n+1]-self.pn)/self.ka))/(1+math.exp((self.PtildLeft[n+1]-self.pn)/self.ka))
        self.F_cs_right[n+1] = (self.fmin + self.fmax*math.exp((self.PtildRight[n+1]-self.pn)/self.ka))/(1+math.exp((self.PtildRight[n+1]-self.pn)/self.ka))
        
        self.F_cs[n+1] = 0.5*(self.F_cs_left[n+1]+self.F_cs_right[n+1])
    
    
    
              
    def UrsinoEfferent(self,n):
        """
        Efferent part of the Baroreceptor model as defined by Ursino 1999
        Calculation of Efferent firing rate resp. sympathetic activity
        """
        self.F_efferent[n+1] = self.fe_inf + (self.fe_0-self.fe_inf)*math.exp(-self.ke*self.F_cs[n+1])
        
        
    
            
    def UrsinoResistanceEffector(self,n):
        """
        Effector part for TPR of the Baroreceptor model as defined by Ursino 1999
        Calculation of effect efferent firing rate (sympathetic activity) on system quantity 
        """
        delay = round(self.DR/self.dt)
        vR = 0.0
        
        if n <= delay:
            
            vR = 0.
            
        elif n > delay:
            deltaF = self.F_efferent[n+1-delay] > self.fe_min
            
            if deltaF >= 0.0:
                vR = self.cR * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                vR = 0.0
        
        self.delta_TPR[n+1] = self.dt/self.tauR *(-self.delta_TPR[n] + vR) + self.delta_TPR[n]
        
        
    def UrsinoEmaxLVEffector(self,n):
        """
        Effector of maximal contractility of left ventricle
        """    
        delay = round(self.DE/self.dt)
        vE = 0.0
        
        if n <= delay:
            vE = 0.
            
        elif n > delay:
            
            deltaF = self.F_efferent[n+1-delay] > self.fe_min
            
            if deltaF >= 0.0:
                vE = self.cE * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                vE = 0.0
        
        self.delta_Emax[n+1] = self.dt/self.tauE *(-self.delta_Emax[n] + vE) + self.delta_Emax[n]
    
    
    def UrsinoVusvEffector(self,n):
        """
        Effector for venous unstressed volume
        """
        delay = round(self.DVusv/self.dt)
        vVusv = 0.0
        
        if n <= delay:
            
            vVusv = 0.
            
        elif n > delay:
            
            deltaF = self.F_efferent[n+1-delay] > self.fe_min
            
            if deltaF >= 0.0:
                
                vVusv = self.cVusv * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                
                vVusv = 0.0
        
        self.delta_Vusv[n+1] = self.dt/self.tauVusv *(-self.delta_Vusv[n] + vVusv) + self.delta_Vusv[n]
        
    
    
    def updateBC2(self,n):
        """
        update Boundary condition of type 2
        TPR and Emax for Varying Elastance heart model
        """
        # loop through all boundary conditions of type 2 and update the resistance at every timestep
        # update the heart properties at completion of the beat
        
        deltaWK_resistance = deltaTPR[n+1]*self.terminalBoundaries
        
        for bcId,bcs in self.boundaryConditionII:
            
                
            if bcs.name == 'Windkessel-3Elements':    
                bcs.Rtotal = bcs.Rtotal + deltaWK_resistance
                    
            elif bcs.name == 'Resistance':
                bcs.Rc = bcs.Rc + deltaWK_resistance
                    
            elif bcs.name == 'Windkessel-2Elements':    
                bcs.Rc = bcs.Rc + deltaWK_resistance
                    
            elif bcs.name == 'VaryingElastanceHeart':
                if bcs.newCycle == True:
                    bcs.Emax = bcs.Emax + self.delta_Emax[n+1]
                
            

       
#     #def updateVenousSide(self):
#       #  """
#       #  to update the Venous unstretched Volume
#       #  """
#               
     
    def UrsinoBRmodel(self,n):
        """
        Calculate Baroreflex as described by Ursino_1999
        Effected quantities are TPR, Emax of the left heart, Vusv
        """
          
        self.UrsinoAfferent(n)
        self.UrsinoEfferent(n)
        self.UrsinoResistanceEffector(n)
        self.UrsinoEmaxLVEffector(n) 
        self.UrsinoVusvEffector(n)
        #self.updateBC2()
        #self.updateVenousSide()

     
    def __call__(self):
        
        n = self.n[0]
        
        if self.modelName == 'Ursino':
            
            self.UrsinoBRmodel(n)
            self.updateBCII(n)
            print self.F_efferent[n]
            print self.delta_TPR[n]
            

