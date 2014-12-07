import numpy as np 
#from numpy.linalg import solve
#from scipy.optimize import fsolve
#from numpy.linalg import inv

import copy

import pprint

import sys,os

from CellMLBaroReceptorModel import *

# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')

class BaroReceptor(object):
    def __init__(self, BaroDict):
        '''
        Baroreceptor model
        
        '''
        #System and Vessel Variables
        
        self.dt = 0
        self.n = 0
        self.Tsteps = 0
        #self.vesselNi = vesselNi
        
        self.cellMLBaroreceptorModel = False
        
        self.updateBaro(BaroDict)
        
        self.epsilon = np.zeros(self.Tsteps)
        self.HeartRate = np.zeros(self.Tsteps)
        
        if self.cellMLBaroreceptorModel == True:
            
            #initialize the CellML Baroreceptor model
            
            (self.states, self.constants) = initConsts()
            timeArray = linspace(0,self.dt,2)
            self.voi, self.states, self.algebraic = solver2(timeArray,self.states,self.constants)
             
        else:
            print "No CellML Baroreceptor Model provided"
        
        
    def solveCellML(self,eps):
                           
        timeArray = linspace(0,self.dt,2)
        self.constants[34] = eps
        self.voi, self.states, self.algebraic = solver2(timeArray,self.states[-1],self.constants)
        #print "heart rate"    
        #print self.algebraic[-1][11]
        #print self.algebraic[-1][10]
        #print self.algebraic[-1][9]
        return self.voi, self.states, self.algebraic, self.algebraic[-1][11], self.constants[34]
        
        
        
    #def solve(self,eps):
           
      #  self.epsilon[self.n] = eps
        
      #  if self.cellMLBaroreceptorModel == True:
       #         
         #   self.solveCellML()

     
    def  __call__(self,eps):
          
        (t, s, a, hr, c) = self.solveCellML(eps)
        self.HeartRate[self.n] = hr
    
    
    def updateBaro(self,baroDict):
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
        
        #self.A = A
        #self.A0 = copy.deepcopy(A[vesselNi])
        
        # constants of baroreceptor
        #self.F0    = 0 #called L0 in article
        #self.g     = 0
        #self.alpha = 0
        #self.tau1  = 0
        #self.tau2  = 0
        
        #self.F1 = F0
        #self.F2 = F0
        
        #self.bc = bc
        
           
    #def loadCellMLBaroModel(self):
      #   '''
      #   create the python class for baroreceptor model from cellML file
      #   '''
         
      #   self.filename     
         
      #   self.cellMlBaroreceporModel = baroReceptorBuganhage
        
   # def __call__(self):
        
        # create local variables for this timestep
      #  dt = self.dt
      #  n = self.n[0]
        
      #  if n%50 == 0:
        
            # 1. calculate strain
        #    An = self.A[n][self.vesselNi]
        #    A0 = self.A0
          #  epsilon = (An-A0)/A0
            
            ## call cellML baor model
<<<<<<< .mine
          #  self.cellMlBaroreceporModel.solve(epsilon, dt)
=======
            self.cellMlBaroreceporModel.solve(P, dt)
>>>>>>> .r890
            
          #  heartRate = self.cellMlBaroreceporModel.heartRate
            
            ## update hartRate in VaryingElastance model
                   
          #  self.bc.update({'T':1.0/heartRate})
        
        
        
    #def baroReceptorBuganhage(self,epsilon,dt):        
         
        ## 2. solve F1,F2  (== L1 L2) functions
        # explicit euler
       # self.F1 = self.F1 + dt/self.tau1 (epsilon-self.F1)
       # self.F2 = self.F2 + dt/self.tau2 (epsilon-self.F2)
         
        # 3. calculate L(L1,L2)
       # F = (self.alpha * self.F1 - self.F2)*(alpha-1)
         
        # 4. calculate firing rate phi
      #  phi = self.phi0 + self.g*(F-self.F0)* 0.5 * (np.sign((F-self.F0)) + 1)
         
        # 5. calcualte phi_sn phi_pn
         
        ## should go into the classBoundaryConditions.py
        # 6. (adrelanin ODE) phi_sn -> calculate heart rate // apply to boundaryConditions
         
      #  return phi
 
            
    
    