import numpy as np 
#from numpy.linalg import solve
#from scipy.optimize import fsolve
#from numpy.linalg import inv

import copy

import pprint
import math

import sys,os

import CellMLBaroReceptorModel as baroreceptorCellML

# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')

class BaroReceptor(object):
    def __init__(self, BaroDict):
        '''
        Baroreceptor model initialisation
        
        '''
        #System and Vessel Variables
        self.dt = 0
        self.n = 0
        self.Tsteps = 0
        self.vesselId = 0
        self.data = {}
        self.boundaryCondition  = 0
        self.CellML = False
        
        # Model from CellML
        self.cellMLBaroreceptorModel = False
        
        
        # update with the class with information from dictionary
        self.update(BaroDict)
        
        
        # initial area, remains unchanged
        self.Ao = self.data['Area'][0]
        
        # Area, radius and strain
        self.Ro = np.power(self.Ao/math.pi,0.5)
        
        sizeEpsilon = np.shape(self.Ao)
        
        self.epsilon = np.zeros(sizeEpsilon)
        self.epsMean = 0
                
        
        if self.cellMLBaroreceptorModel == True:
            
            #initialize the CellML Baroreceptor model
            
            (self.states, self.constants) = baroreceptorCellML.initConsts()
            timeArray = np.linspace(0,self.dt,2)
            self.voi, self.states, self.algebraic = baroreceptorCellML.solver2(timeArray,self.states,self.constants)
             
        else:
            print "No CellML Baroreceptor Model provided"
        
        
    #def solveCellML(self,eps):
    def solveCellML(self, epsMean):
                           
        timeArray = np.linspace(0,self.dt,2)
        self.constants[34] = epsMean
        self.voi, self.states, self.algebraic = baroreceptorCellML.solver2(timeArray,self.states[-1],self.constants)
        
        return self.voi, self.states, self.algebraic, self.algebraic[-1][11], self.constants[34]
        
        
        
        
    def __call__(self):
        
        n = self.n[0]
        
        ## read area, calculate radius and strain and save to data dictionary
        
        A = self.data['Area'][n]
        R = np.power(A/math.pi,0.5)
        
        epsilon = (R - self.Ro)/self.Ro
        epsMean = np.mean(epsilon)
        
        # update the data dictionary with the strain and the mean strain
        self.data['Strain'][n+1]  = epsilon
        self.data['MStrain'][n+1] = epsMean
        
        # print the calculated mean strain
        print 'cBR95: epsilon'
        print epsMean
         
        if n%20 == 0:
            # solve the cellML system using the function defined above
            self.solveCellML(epsMean)
        
        
        print "heart rate"    
        print self.algebraic[-1][11]
        print self.algebraic[-1][10]
        print self.algebraic[-1][9]
        
        # update the heart rate and the data dictionary with the current heart rate
        self.data['HR'][n+1] = self.algebraic[-1][11]
        
        # update boundary conditions with new heart rate / period
    
    
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
        
            