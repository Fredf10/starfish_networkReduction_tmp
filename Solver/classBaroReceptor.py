import numpy as np 
#from numpy.linalg import solve
#from scipy.optimize import fsolve
#from numpy.linalg import inv

import copy

import pprint
import math

import sys,os

from CellMLBaroReceptorModel import *

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
        
        
        # Model from CellML
        self.cellMLBaroreceptorModel = False
        
        
        # update with the class with information from dictionary
        self.updateBaro(BaroDict)
        
        
        # initial area, remains unchanged
        self.Ao = self.data['Area'][0]
        
        # Area, radius and strain
        self.Ro = np.power(self.Ao/math.pi,0.5)
        self.A = self.Ao
        self.R = self.Ro
        
        sizeEpsilon = np.shape(self.A)
        
        self.epsilon = np.zeros(sizeEpsilon)
        self.epsMean = 0
        
        # heart rate
        self.HeartRate = np.zeros(self.Tsteps)    
        
        
        if self.cellMLBaroreceptorModel == True:
            
            #initialize the CellML Baroreceptor model
            
            (self.states, self.constants) = initConsts()
            timeArray = linspace(0,self.dt,2)
            self.voi, self.states, self.algebraic = solver2(timeArray,self.states,self.constants)
             
        else:
            print "No CellML Baroreceptor Model provided"
        
        
    #def solveCellML(self,eps):
    def solveCellML(self):
                           
        timeArray = linspace(0,self.dt,2)
        self.constants[34] = self.epsMean
        self.voi, self.states, self.algebraic = solver2(timeArray,self.states[-1],self.constants)
        
        return self.voi, self.states, self.algebraic, self.algebraic[-1][11], self.constants[34]
        
        
        
        
    def __call__(self):
        
        
        ## read area, calculate radius and strain and save to data dictionary
        
        self.A = self.data['Area'][self.n]
        self.R = np.power(self.A/math.pi,0.5)
        
        self.epsilon = (self.R - self.Ro)/self.Ro
        self.epsMean = np.mean(self.epsilon)
        
        # update the data dictionary with the strain and the mean strain
        self.data['Strain'][self.n] = self.epsilon
        self.data['MStrain'][self.n] = self.epsMean
        
        # print the calculated mean strain
        print 'epsilon'
        print self.epsMean
         
        # solve the cellML system using the function defined above
        self.solveCellML()
        
        print "heart rate"    
        print self.algebraic[-1][11]
        print self.algebraic[-1][10]
        print self.algebraic[-1][9]
        
        # update the heart rate and the data dictionary with the current heart rate
        self.HeartRate[self.n] = self.algebraic[-1][11]
        self.data['HR'][self.n] = self.algebraic[-1][11]
        
        # update boundary conditions with new heart rate / period
    
    
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
        
            