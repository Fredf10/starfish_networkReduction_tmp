import numpy as np 

import copy

import pprint
import math
import sys,os


# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')
sys.path.append(cur+'/cellMLBaroreflexModels')

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
        self.currentTimeStep = 0
        self.currentMemoryIndex = 0
        self.nTsteps = 0
        self.receptorType = ''
        self.modelName = ''
        
        # Model from CellML or hardcoded
        self.cellMLBaroreceptorModel = False
        self.cellMLimport = '' # for the import of the CellML file
        
        self.boundaryCondition = 0
        self.boundaryConditionII = 0
        self.boundaryConditionIIout = {}
        
        
        self.update(BaroDict)
                
        #initialize the CellML Baroreceptor model
        if self.cellMLBaroreceptorModel == True:
            
            # import the model from the file which defines it
            self.cellMLimport = "import" + " " + self.modelName + " " + " as "+ "baroreceptorCellML "            
            exec(self.cellMLimport)
            
            # input to te CellML model and output from it
            self.cellMLinputID = baroreceptorCellML.inputID
            self.cellMLoutputArray = baroreceptorCellML.outputArray
            self.cellMLoutputID = baroreceptorCellML.outputID
            
             
            (iniStates, self.constants) = baroreceptorCellML.initConsts()
            timeArray = np.linspace(0,self.dt,2)
            self.voi, self.states, self.algebraic = baroreceptorCellML.solver2(timeArray,iniStates,self.constants)
             
               
        else:
            print "No CellML Baroreceptor Model provided!"

            
    ### solving of imported CellML model
    def solveCellML(self, input):
        """
        solve CellML model
        """
        exec(self.cellMLimport)
        
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
            
            Tperiod  = self.algebraic[-1][self.cellMLoutputID]
            
            self.boundaryCondition.updatePeriodRuntime(Tperiod,(self.newUpdateTime-2)*self.dt)
            
            self.oldUpdateTime = self.newUpdateTime
            self.newUpdateTime = self.oldUpdateTime + round(self.boundaryCondition.Tperiod/self.dt)
        
        
    def calcAndupdatePeriodTypeIIcellML(self,input):
        """
        Function to calculate new period
        Calls self.solveCellML()
        is used by self.__call__()
        """
        if self.boundaryConditionII.name == 'VaryingElastanceSimple':
            self.voi, self.states, self.algebraic = self.solveCellML(input)
                
            if self.boundaryConditionII.newCycle == True:
                self.boundaryConditionII.T = self.algebraic[-1][self.cellMLoutputID]
                self.boundaryConditionII.Tpeak = 0.43*self.boundaryConditionII.T
                self.newCycles = np.append(self.newCycles,self.currentTimeStep)
                
                print "BR122"
                print self.boundaryConditionII.T
        else: 
            print "not a varying elastance heart"
                 
         
         
            
    def update(self,baroDict):
            '''
            updates the Baroreceptor using a dictionary in form of 
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
        Baroreceptor.__init__(self,BaroDict)
        self.vesselId = 0
        
        # boundary conditions
        #self.boundaryCondition  = 0
        #self.boundaryConditionII = 0
        
        self.Area1 = 0
        self.Area2 = 0
        self.Pressure1 = 0
        self.Pressure2 = 0
        self.Strain = 0
        self.MStrain = 0
        self.T = 0
        
        if self.modelName == 'bugenhagenAorticBR':
            self.n = np.ones(self.nTsteps+1)*self.algebraic[0][3]
            self.Tsym = np.ones(self.nTsteps+1)*self.algebraic[0][7]
            self.Tparasym = np.ones(self.nTsteps+1)*self.algebraic[0][8]
            self.c_nor = np.ones(self.nTsteps+1)*self.states[0][3]
            self.c_ach = np.ones(self.nTsteps+1)*self.states[0][4]
        
        
        elif self.modelName == 'pettersenAorticBR':
            self.n = np.ones(self.nTsteps+1)*self.algebraic[0][5]
            self.Tsym = np.ones(self.nTsteps+1)*self.algebraic[0][9]
            self.Tparasym = np.ones(self.nTsteps+1)*self.algebraic[0][10]
            self.c_nor = np.ones(self.nTsteps+1)*self.states[0][2]
            self.c_ach = np.ones(self.nTsteps+1)*self.states[0][3]
        
        
        # initial Compliance information for estimation of unstretched radius
        #self.initialCompliance1 = 0
        #self.initialCompliance2 = 0
        
        # update with dictionary
        self.update(BaroDict)
        
        self.newCycles = np.zeros(1)
        
        self.Ao_1 = self.Area1[0] 
        self.Ao_2 = self.Area2[0]
        
        self.Po_1 = self.Pressure1[0]
        self.Po_2 = self.Pressure2[0]
            
        #initial area and unstretched radius
        self.Ro_1 = np.zeros(np.shape(self.Area1)[0])
        self.Ro_2 = np.zeros(np.shape(self.Area2)[0])
        
        # estimate the unstretched Radii R0 for the calculation of the strain
        self.estimateUnstretchedRadius()
        
        # update Time for the period of the inflow function
        self.oldUpdateTime = 0
        if self.boundaryCondition != 0:
            self.newUpdateTime = round(self.boundaryCondition.Tperiod/self.dt)
            print "BR160"
            print self.newUpdateTime
        
            
    def estimateUnstretchedRadius(self):
        """
        Function to estimate the unstretched radius of the Vessels of the Aortic Arch for the calculation of strain.
        Takes information about compliance and initial vessel geometry.
        Takes the initial compliance and assumes it as constant, uses the initial Area to calculate 
        the Area at zero pressure, where the vessel is in its unstressed state, so that the true
        strain can be calculated.
        """
        
        #compliance1 = self.initialCompliance1.C(self.Po_1)
        #compliance2 = self.initialCompliance2.C(self.Po_2)
        
        #deltaA1 = compliance1 * self.Po_1
        #deltaA2 = compliance2 * self.Po_2

        #self.Ro_1 = np.power((self.Ao_1 - deltaA1)/math.pi,0.5)
        #self.Ro_2 = np.power((self.Ao_2 - deltaA2)/math.pi,0.5)
        
        
        f = 1.3 # estimated value, probably somewhere between 1.35 and 1.65
        
        self.Ro_1 = np.power((self.Ao_1)/math.pi,0.5)/f
        self.Ro_2 = np.power((self.Ao_2)/math.pi,0.5)/f
        
    
       
    def __call__(self):
        
        n = self.currentTimeStep[0]
        n_mem = self.currentMemoryIndex[0]
        
        ## read area, calculate radius and strain
        A1 = self.Area1[n_mem]
        R1 = np.power(A1/math.pi,0.5)
       
        A2 = self.Area2[n_mem]
        R2 = np.power(A2/math.pi,0.5)
        
        epsilon1 = (R1 - self.Ro_1)/self.Ro_1
        epsilon2 = (R2 - self.Ro_2)/self.Ro_2
        
        # concatenate the strain arrays of the two vessels of the Aortic Arch
        epsilon = np.concatenate((epsilon1,epsilon2),axis=0)
        epsMean = np.mean(epsilon)
        
        
        # update the strain and the mean strain
        if n < self.nTsteps-1:
            self.Strain[n+1]  = epsilon
            self.MStrain[n+1] = epsMean
            
            
        # use CellML model to calculate new period, the method calcAndupdatePeriodcellML
        # does this on a beat-to-beat basis
        if self.cellMLBaroreceptorModel == True:
            
            if self.boundaryCondition != 0:
                self.calcAndupdatePeriodTypeIcellML(n,self.MStrain[self.oldUpdateTime:(self.newUpdateTime-2)])
                print "BR255"
                print self.algebraic[-1][self.cellMLoutputID]
                print self.boundaryCondition.Tperiod
                
            else:
                self.calcAndupdatePeriodTypeIIcellML(self.MStrain[n:(n+1)])
            
            # update the current period
            if self.cellMLoutputArray == 'algebraic':
                
                self.T[n+1] = self.algebraic[-1][self.cellMLoutputID]
            
            
            elif self.cellMLoutputArray == 'states':
                
                self.T[n+1] = self.states[-1][self.cellMLoutputID]
            
            
            ### for Bugenhagen###
            self.n[n+1] = self.algebraic[-1][3]
            self.Tsym[n+1] = self.algebraic[-1][7]
            self.Tparasym[n+1] = self.algebraic[-1][8]
            self.c_nor[n+1] = self.states[-1][3]
            self.c_ach[n+1] = self.states[-1][4]
            
            
            ### for Pettersen###
            #self.n[n+1] = self.algebraic[-1][5]
            #self.Tsym[n+1] = self.algebraic[-1][9]
            #self.Tparasym[n+1] = self.algebraic[-1][10]
            #self.c_nor[n+1] = self.states[-1][2]
            #self.c_ach[n+1] = self.states[-1][3]
            
            
        else: pass
            #print 'Error: currently no hardcoded models for the Aortic BR available'
            

   
class CarotidBaroreceptor(Baroreceptor):
    '''
    for models of the Carotid Baroreceptors
    '''
    def __init__(self,BaroDict):
        
        """
        constructor for the CarotidBaroreceptor
        """
        #System and Vessel Variables
        Baroreceptor.__init__(self,BaroDict) # from mother class
        self.vesselIdLeft = 0
        self.vesselIdRight = 0
        self.terminalBoundaries = 0
        self.VenousPool = 0
        
        self.pressureLeft = 0
        self.pressureRight = 0
        
        
        # update with dictionary
        self.update(BaroDict)
        
        
        ### total resistance
        self.Res0 = {} # Windkssel resistances of all boundaries --> for updating the resistances
        self.ResTot0 = 0.0 # total resistance of the windkessels
        
        for key in self.boundaryConditionIIout:
            self.Res0[key] = self.boundaryConditionIIout[key].Rtotal
            self.ResTot0 = self.ResTot0 + 1/self.boundaryConditionIIout[key].Rtotal
                
            
        # total resisitance of the Windkessels    
        self.ResTot0 = 1.0/self.ResTot0
        
        
        ### for update of period in type 1 boundary condition
        self.oldUpdateTime = 0
        if self.boundaryCondition != 0:
            self.newUpdateTime = round(self.boundaryCondition.Tperiod/self.dt)
        
        
        #######################################
        """
        Parameters for the Ursino model of the carotid baroreflex
        """
        #######################################
        
        ### model parameters afferent part of the Ursino model - Ursino 1999
        
        self.pn = 135*133.32#92*133.32 #12266. adjusted from 92mmHg, 85 was found in literature
        self.ka = 1567.
        self.fmin = 2.51
        self.tau_z= 6.37
        self.fmax = 47.78
        self.tau_p= 2.076
        
        # model parameters efferent part of the Ursino model  - Ursino 1999
        
        self.fe_inf = 2.10 #
        self.fe_0 = 16.11 #4.12
        self.ke = 0.0675 #
        self.fe_min = 2.66 #
        
        # model parameters for the TPR effector part of the Ursino model  - Ursino 1999
        
        self.cR = 0.317*133.32e6 # 
        self.tauR = 6 #
        self.DR = 2.0 # 
        self.R0 = 0.61*133.32e6 #
        
        
        
        # model parameters for the EmaxLV effector part of the Ursino model  - Ursino 1999
        
        self.cE = 0.475 * 133.322368e6
        self.tauE = 8.
        self.DE = 2.0
        self.E0 = 2.392*133.322368e6
        
        # model parameters for the T effector part of the Ursino model  - Ursino 1999
        
        self.cT = -0.08
        self.tauT = 2.0
        self.DT = 2.0
        self.T0 = 0.58 
        
        
        # model parameters for the Vusv effector part of the Ursino model  - Ursino 1999
        
        self.cVusv = -625.0e-6#-125e-6 
        self.tauVusv = 20.0
        self.DVusv = 3.0
        self.Vusv0 = 3213e-6
        
        
        ### initial values for affected quantities
        
        self.dTPRin = 5.85e7 #5.939e7 #5.85e7
        self.dTin = -0.1114 #-0.11259 #-0.1114
        self.dEmaxin = 87.65e6 #88.96815e6 #87.65e6
        self.dVusvin = -877.54e-6 #-875.0e-6
        self.Ptildin = 10600. #11680.0 #10600.
        
        self.ratio = (1.53792e8)/(self.R0 + self.dTPRin) #1.099887574 
        
        ### states of the Ursino model###
        # afferent
        self.PtildLeft = np.zeros(self.nTsteps+1)
        self.PtildLeft[0] = self.Ptildin
        
        self.PtildRight = np.zeros(self.nTsteps+1)
        self.PtildRight[0] = self.Ptildin
        
        self.F_cs_left = np.zeros(self.nTsteps+1)
        self.F_cs_right = np.zeros(self.nTsteps+1)
        
        self.F_cs = np.zeros(self.nTsteps+1) # mean value of afferent firing (mean of left and right carotid sinus)
        
        
        # efferent
        self.F_efferent = np.zeros(self.nTsteps+1)
        
        
        # effectors
        self.delta_TPR = np.ones(self.nTsteps+1)*self.dTPRin
        self.delta_Emax = np.ones(self.nTsteps+1)*self.dEmaxin
        self.delta_Vusv = np.ones(self.nTsteps+1)*self.dVusvin
        self.delta_T = np.ones(self.nTsteps+1)*self.dTin
        
        
        ### quantities in the left ventricle --> to extract for postprocessing
        self.newCycles = np.zeros(1)
        self.Pheart = np.zeros(self.nTsteps+1)
        self.Vheart = np.zeros(self.nTsteps+1)
        
        
    #############################################################################
    
    """
    methods defining the different parts of the Ursino Baroreceptor model
    - Afferent
    - Efferent
    - Effector
    """
    
    def UrsinoAfferent(self,n,n_mem):
        """
        Afferent part of the Baroreceptor model as defined by Ursino 1999
        The firing rate of the left and right carotid sinus are averaged to give a global firing rate
        the firing rate of the carotid sinus will be processed by the efferent part
        """
        
        pLeft = np.mean(self.pressureLeft[n_mem])
        pRight = np.mean(self.pressureRight[n_mem])
        
        if n == 0:
            
            dpLeft =  0.
            dpRight = 0.
            
        else:
            
            dpLeft = (pLeft - np.mean(self.pressureLeft[n_mem-1]))/self.dt
            dpRight = (pRight - np.mean(self.pressureRight[n_mem-1]))/self.dt
            
        
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
            vR = self.dTPRin
            
        elif n > delay:
            deltaF = self.F_efferent[n+1-delay] - self.fe_min
            
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
        vE = 0
        
        if n <= delay:
            vE = self.dEmaxin
            
        elif n > delay:
            
            deltaF = self.F_efferent[n+1-delay] - self.fe_min
            
            if deltaF >= 0.0:
                vE = self.cE * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                vE = 0.0
        
        self.delta_Emax[n+1] = self.dt/self.tauE *(-self.delta_Emax[n] + vE) + self.delta_Emax[n]
    
    
    
    def UrsinoTEffector(self,n):
        """
        Effector of maximal contractility of left ventricle
        """    
        delay = round(self.DT/self.dt)
        vT = 0
        
        if n <= delay:
            vT = self.dTin
            
        elif n > delay:
            
            deltaF = self.F_efferent[n+1-delay] - self.fe_min
            
            if deltaF >= 0.0:
                vT = self.cT * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                vT = 0.0
        
        self.delta_T[n+1] = self.dt/self.tauT *(-self.delta_T[n] + vT) + self.delta_T[n]
    
    
    
    def UrsinoVusvEffector(self,n):
        """
        Effector for venous unstressed volume
        """
        delay = round(self.DVusv/self.dt)
        vVusv = 0.0
        
        if n <= delay:
            vVusv = self.dVusvin
            
        elif n > delay:
            deltaF = self.F_efferent[n+1-delay] - self.fe_min
            
            if deltaF >= 0.0:
                vVusv = self.cVusv * math.log(deltaF + 1)
                
            elif deltaF < 0.0:
                vVusv = 0.0
        
        self.delta_Vusv[n+1] = self.dt/self.tauVusv *(-self.delta_Vusv[n] + vVusv) + self.delta_Vusv[n]
        
    
    
    def updateBC1(self,n):
        """
        update period in BC of type 1 at inflow boundaries
        uses the method updatedPeriodRuntime of the BoundaryCondition class
        uses newUpdateTime
        """
        
        if n == (round(self.newUpdateTime-2)):
             
            newPeriod = self.T0 + self.delta_T[n+1]
            
            self.boundaryCondition.updatePeriodRuntime(newPeriod,(self.newUpdateTime-2)*self.dt)
            
            self.oldUpdateTime = self.newUpdateTime
            self.newUpdateTime = self.oldUpdateTime + round(self.boundaryCondition.Tperiod/self.dt)
    
    
    def updateBC2(self,n):
        """
        update Boundary condition of type 2
        TPR and Emax for Varying Elastance heart model
        loop through all boundary conditions of type 2 and update the resistance at every timestep
        --> the single resistances are updated with respect to their initial value dRtot/Rtot0 = dRi/Ri0
        update the heart properties at completion of the beat
        """
        
        newTotalResistance = self.ratio*(self.R0 + self.delta_TPR[n+1])
         
        for key in self.boundaryConditionIIout:    
            if self.boundaryConditionIIout[key].name == 'Windkessel-3Elements':
                self.boundaryConditionIIout[key].Rtotal = (newTotalResistance)/self.ResTot0*self.Res0[key]  
                
            ## for Resistance outlets      
            elif self.boundaryConditionIIout[key].name == 'Resistance':
                self.boundaryConditionIIout[key].Rc = (newTotalResistance)/self.ResTot0*self.Res0[key]
                
            ## for WK 2 outlets       
            elif self.boundaryConditionIIout[key].name == 'Windkessel-2Elements':    
                self.boundaryConditionIIout[key].Rc = (newTotalResistance)/self.ResTot0*self.Res0[key]
            
        ## update inlet BC's of type 2 --> Varying Elastance heart           
        if self.boundaryConditionII != 0:
            if self.boundaryConditionII.newCycle == True:
                self.boundaryConditionII.Emax = self.E0 + self.delta_Emax[n+1]
                self.boundaryConditionII.T = self.T0 + self.delta_T[n+1]
                self.newCycles = np.append(self.newCycles,n)
                self.Pheart[n] = self.boundaryConditionII.pressure[n]
                self.Vheart[n] = self.boundaryConditionII.volume[n]
                
                #print "BR538"
                #print self.boundaryConditionII.Emax
                #print self.boundaryConditionII.T
            
            else:
                self.Pheart[n] = self.boundaryConditionII.pressure[n]
                self.Vheart[n] = self.boundaryConditionII.volume[n]
                

    def updateVenousSide(self,n):
         """
         to update the Venous unstretched Volume
         """
         ### first update when the delay in the effector of Vusv is complete
         if self.currentTimeStep*self.dt > self.DVusv:
             self.VenousPool.Vusv = self.Vusv0 + self.delta_Vusv[n]
     
    def UrsinoBRmodel(self,n,n_mem):
        """
        Calculate Baroreflex as described by Ursino_1999
        Effected quantities are TPR, Emax of the left heart, Vusv
        Update the affected quantities in their classes
        """
          
        self.UrsinoAfferent(n,n_mem) # calculate the afferent signals
        self.UrsinoEfferent(n) # calculate the efferent signal
        
        self.UrsinoResistanceEffector(n) # the different effector parts
        self.UrsinoEmaxLVEffector(n) 
        self.UrsinoVusvEffector(n)
        self.UrsinoTEffector(n)
        
        #self.updateBC1(n)
        self.updateBC2(n) # update BC type 2
        self.updateVenousSide(n) # update the venous side with Vusv

     
    def __call__(self):
        
        """
        call function
        takes the Baroreceptor model (Ursino_1999) and evaluates it
        """
        
        n = self.currentTimeStep[0]
        n_mem = self.currentMemoryIndex[0]
        
        if self.modelName == 'Ursino':
            
            self.UrsinoBRmodel(n,n_mem)
            #print "BR646"
            #print self.F_efferent[n]
            #print self.delta_TPR[n]
            #print self.delta_T[n]
            #print self.delta_Emax[n]
            #print self.delta_Vusv[n]

    
