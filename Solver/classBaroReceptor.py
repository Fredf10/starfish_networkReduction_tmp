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
                
        # update Time for the heart rate
        self.oldUpdateTime = 0
        self.newUpdateTime = self.boundaryCondition.Tperiod/self.dt
        
        
        if self.cellMLBaroreceptorModel == True:
            
            #initialize the CellML Baroreceptor model
            
            (self.states, self.constants) = baroreceptorCellML.initConsts()
            timeArray = np.linspace(0,self.dt,2)
            self.voi, self.states, self.algebraic = baroreceptorCellML.solver2(timeArray,self.states,self.constants)
            
            print "SHAPE"
            print np.shape(self.states)
            
        else:
            print "No CellML Baroreceptor Model provided"
        
        
    #def solveCellML(self,eps):
    #def solveCellML(self, epsMean):
                           
        #timeArray = np.linspace(0,self.dt,2)
        #self.constants[34] = epsMean
        #self.voi, self.states, self.algebraic = baroreceptorCellML.solver2(timeArray,self.states[-1],self.constants)
        
        #return self.voi, self.states, self.algebraic, self.algebraic[-1][11], self.constants[34]

    # function to calculate new heart rate on a beat to beat basis
    def solveCellML(self, epsMean):
        
        nbElements = np.shape(epsMean)[0]            
        timeArray = np.linspace(0,(nbElements-2)*self.dt,(nbElements-1))
        
        for it in xrange(0,(nbElements-1)):
            
            self.constants[34] = epsMean[it]
            self.voi, self.states, self.algebraic = baroreceptorCellML.solver2(timeArray[it:(it+1)],self.states[-1][:],self.constants)
        
        print "SHAPE"
        print np.shape(self.states)
        return self.voi, self.states, self.algebraic
        
        
        
        
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
        #print 'cBR95: epsilon'
        #print epsMean
         
        #if n%20 == 0:
        if n == (round(self.newUpdateTime-2)):
            
            # solve the cellML system using the function defined above
            self.voi, self.states, self.algebraic = self.solveCellML(self.data['MStrain'][self.oldUpdateTime:(self.newUpdateTime-2)])
            
            Tperiod  = self.algebraic[-1][12]
            #print BR122
            #print 'T -BC'
            #print self.boundaryCondition.Tperiod
            #print 'T'
            #print Tperiod
            #print 'HR'
            print "BR128 - TPeriod"
            print self.algebraic[-1][12]
            
            self.boundaryCondition.updatePeriodRuntime(Tperiod,(self.newUpdateTime-2)*self.dt)
            print self.boundaryCondition.Tperiod
            
            print self.constants[34]
            
            
            self.oldUpdateTime = self.newUpdateTime
            self.newUpdateTime = self.oldUpdateTime + self.boundaryCondition.Tperiod/self.dt
        
        # update the heart rate and the data dictionary with the current heart rate
        self.data['HR'][n+1] = self.algebraic[-1][11]
        
        
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
        
            