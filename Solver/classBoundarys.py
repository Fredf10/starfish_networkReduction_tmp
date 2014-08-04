import numpy as np 
#from numpy.linalg import solve
#from scipy.optimize import fsolve
#from numpy.linalg import inv

import copy

import pprint

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')
from classBoundaryConditions import *


class Boundary():
    def __init__(self, vessel, boundaryConditions, rigidArea, dt, n, Tsteps, systemEquation):
        '''
        Constructor of Boundary
        
        Initializes one Boundary of a vessel
        input: - boundaryConditions:
                  a list of all boundaryConditions of the vessel at one position i.g 0 or -1 

        Note: if network consits of a single vessel make sure the function is called twice 
              once for the start and once for the end boundary at 0 and -1!
              
        variables (most important):
        self.position        :    position of the boundary
        self.ID              :    Vessel ID to which the boundary belongs
        self.duFunction      :    function which gives du-vector back
        self.bcType1         :    list of all type1 functions found in the given boundaryConditions
        self.omegaInput      :    function which gives the _omega of omega_ back depending on self.position
        self.omegaFunctio    :    function with boundaryCondtion of type2 calculating the omega-vector
        '''
        self.position = None
        
        self.name = ' '.join(['Boundary',str(vessel.Id)])
        self.type = ''
        
        #Function which gives du-vector back
        self.duFunction = None
        self.duVector = np.zeros((Tsteps,2))
        self._du_ = np.empty(2)
        # list of all type1 functions found in the given boundaryConditions, which determine 
        # the du-Function
        self.bcType1 = []
        #Function which gives the _omega of omega_ back depending on self.position
        self.omegaInput = None
        #Function with boundaryCondtion of type2 calculating the omega-vector
        self.omegaFunction = None
        #list of all type2 functions found in the given boundaryConditions
        self.bcType2 = []
        
        #System and Vessel Variables
        self.z  = vessel.z
        self.dt = dt
        self.n  = n
        self.systemEquation = systemEquation
        
        #SolutionVariables
        self.P = vessel.Psol
        self.Q = vessel.Qsol
        self.A = vessel.Asol
        
        self.BloodFlowSep  = np.zeros((Tsteps,2)) ## [ dbloodVolume in, dbloodVolume out] blood volume from one timestep to the next
        self.BloodVolumen  = np.zeros(2)
                
        posTemp = []
        for bC in boundaryConditions:            
            if bC.type == 1:
                self.bcType1.append(bC)
            elif bC.type == 2:
                self.bcType2.append(bC)
            posTemp.append(bC.position)
            
        #find out position of the Boundary + check if all conditions are defined at the same side
        if sum(posTemp) == 0: self.position = 0
        elif sum(np.array(posTemp)+1) == 0: self.position = -1
        else: raise ValueError("One position of one boundaryCondition is not correct")
        
        #initialize solution variable
        self.BloodFlowSep[-1][self.position] = vessel.Qsol[0][self.position]
        
        # 1. Set duFunction 
        if len(self.bcType1) == 0:    
            # in for callMacCormack normal
            self.duFunction = self.duFunctionZero
            self.duRuntimeEvaluation = False
        elif len(self.bcType1) == 1:
            self.duFunction = self.duFunctionSingle
            self.bcType1ConditionQuantity = self.bcType1[0].conditionQuantity
            self.duRuntimeEvaluation      = self.bcType1[0].runtimeEvaluation
            precribeTotalValues           = self.bcType1[0].prescribeTotalValues
        elif len(self.bcType1) > 1:
            self.duFunction = self.duFunctionMulti
            self.bcType1ConditionQuantity = self.bcType1[0].conditionQuantity
            self.duRuntimeEvaluation      = self.bcType1[0].runtimeEvaluation
            precribeTotalValues           = self.bcType1[0].prescribeTotalValues
            print "WARNING: Multiple Type1 Boundary Conditions are prescribed as influx condition!"
        # prepare duVector from boundaryConditions:
        self.duEvaluateVector()
        # set bcType1 to type instead of list
        try: self.bcType1 = self.bcType1[0]
        except:pass
        
        # 3. Set Condition there should only one! if None apply the one condition depending on the side
        if len(self.bcType2) == 1:
            self.omegaFunction = self.bcType2[0]
            self.type = self.bcType2[0].name
            
        elif len(self.bcType2) == 0:    
            # set calculation type for Standard Boundary
            self.type = self.bcType1.name
            try:
                if precribeTotalValues == False:
                    self.omegaFunction = PrescribedInflux()
                elif precribeTotalValues == True and self.bcType1ConditionQuantity == 'Flow':
                    self.omegaFunction = PrescribedTotalFlow()
                elif precribeTotalValues == True and self.bcType1ConditionQuantity == 'Pressure':
                    self.omegaFunction = PrescribedTotalPressure()
            except:
                self.omegaFunction = PrescribedInflux()
            self.omegaFunction.setPosition(self.position)
        else:
            print "ERROR classBoundary: Too many type2-boundary Conditions defined!"
        
        # 4. Define the output of A, dependend if rigidArea
        if rigidArea == True:
            self.AFunction = self.AFunctionSys0
        else:
            self.A_nID = vessel.A_nID
            self.AFunction = self.AFunctionSys1
        
        ## 5. Define the call function depending on the solving Scheme
        #if solvingScheme == "MacCormack_Field": 
        #    self.__call__ = self.callMacCormackField

    ### Function which calculated du
    def duFunctionZero(self,n,dt):
        '''
        Determine the du-vector = [0,0] if no type1 boundaryConditions are given
        '''
        return np.zeros(2)
    
    def duFunctionSingle(self,n,dt):
        '''
        Determine the du-vector with the values of the given
        type1 boundaryCondition 
        '''     
        return self.bcType1.calculateDu(n,dt)
    
    def duFunctionMulti(self,n,dt):
        '''
        Determine the summized du-vector with the values of all given
        type1 boundaryConditions 
        '''
        du = np.zeros(2)
        for bc in self.bcType1:
            du = du + bc.calculateDu(n,dt)
        return du
    
    def duEvaluateVector(self):
        '''
        Pre-Calculate the BoundaryConditions duVector of Type1
        '''
        for bc in self.bcType1:
            self.duVector = self.duVector+bc.calculateDuVector(len(self.duVector),self.dt)
        
    def duFunctionRuntime(self,n,dt,Z1,Z2,position):   
        
        if self.duRuntimeEvaluation == True:
            duPrescribed = self.duFunction(n,dt) 
        else:
            duPrescribed = self.duVector[n]
        # check if influx Values should be prescribed
        if self.bcType1 != []:
            # check if Flow is prescribed
            if self.bcType1ConditionQuantity == 'Flow':
                #position 0: pf = Zf*Qf
                if position == 0:
                    duPrescribed[0]  = Z1*duPrescribed[1]
                #position -1: pb = -Zb*Qb
                elif position == -1:
                    duPrescribed[0] = -Z2*duPrescribed[1]
            # check if Pressure is prescribed
            elif self.bcType1ConditionQuantity == 'Pressure':
                #position 0: Qf = pf/Zf
                if position == 0:
                    duPrescribed[1]  = duPrescribed[0]/Z1
                #position -1: Qb = - pf/Zb
                elif position == -1:
                    duPrescribed[1] = -duPrescribed[0]/Z2

        return duPrescribed
    
    ## Function to define the output of A, dependend on the characteristic system 0.1
    
    def AFunctionSys0(self,A,p):
        return A
    
    def AFunctionSys1(self,A,p):
        return self.A_nID(p,self.position)
    
    def __call__(self): #callMacCormackField(self):    
        '''
        new Boundary method calculates the values at the boundary
        and applies it to the new values
        '''
        
        # create local variables for this timestep
        dt = self.dt
        n = self.n[0]
        
        P = self.P[n]
        Q = self.Q[n]
        A = self.A[n]
        
        z = self.z
        position = self.position
        
        # calculate need values
        self.systemEquation.updateLARL(P,Q,A,idArray=[position])
        L    = self.systemEquation.L[position]
        R    = self.systemEquation.R[position]
        LMBD = self.systemEquation.LAMBDA[position]
        Z1   = self.systemEquation.Z[position][0]
        Z2   = self.systemEquation.Z[position][1]
        _omega_field = self.systemEquation.domega[position]
        
        #calculate the du_vector using given boundaryConditions of type 1
        duPrescribed = self.duFunctionRuntime(n, dt, Z1, Z2, position)
        
        #calculate the dp and dq using given boundaryConditions of type 2
        dPQ_calc,dBloodVolumen = self.omegaFunction(_omega_field, duPrescribed, R, L, n, dt, P[position], Q[position], A[position],Z1, Z2)      

        # calculate new values for p and q
        P_calc = dPQ_calc[0]+ P[position]
        Q_calc = dPQ_calc[1]+ Q[position]
                
        # check new p value
        if P_calc < 0:
            print "ERROR: {} calculated negativ pressure at time {} (n {},dt {}), exit system".format(self.name,n*dt,n,dt)
            print P_calc
            exit()
        
        # calculate new a value
        A_calc =  self.AFunction(A[position],[P_calc])
                
        # apply values to solution array
        self.P[n+1][position] = P_calc
        self.Q[n+1][position] = Q_calc
        self.A[n+1][position] = A_calc
        
        try: self.BloodFlowSep[n] = self.BloodFlowSep[n-1]+dBloodVolumen
        except: print "passed bloodflow integration"
        try:    self.BloodVolumen = self.BloodVolumen + 0.5*(self.BloodFlowSep[n-1]+self.BloodFlowSep[n])*dt
        except: self.BloodVolumen = self.BloodVolumen + self.BloodFlowSep[n]*dt
        
        #if n > 3680: print "cB285 ", abs(position), self.BloodVolumen[0]*1.e6, self.BloodVolumen[1]*1.e6
        