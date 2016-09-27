import os,sys
import h5py
import numpy as np
from scipy.integrate import simps


# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))
# sys.path.append(cur+'/../'+'/NetworkLib')

import NetworkLib.classVascularNetwork as cVascNw
from NetworkLib.classBoundaryConditions import *

import moduleFilePathHandler as mFPH

class NetworkLogFile:
    
    def __init__(self, vascularNetwork, dataNumber = "xxx", networkLogFile = None):
        
        self.networkName = vascularNetwork.getVariableValue('name')
        self.dataNumber = dataNumber
        if networkLogFile == None:
            self.networkLogFile =  mFPH.getFilePath('networkLogFile', self.networkName, dataNumber, 'write')
            #print networkLogFile
        else:
            self.networkLogFile = networkLogFile
        TemplateLogFileDirectory = mFPH.getDirectory('networkXmlFileTemplateDirectory', 'singleVessel_template', dataNumber, 'read')
        #print TemplateLogFileDirectory
        self.templateFile = ''.join([TemplateLogFileDirectory,'/logFile_template.tex'])
        
        self.solutionFile = mFPH.getFilePath('solutionFile', self.networkName, dataNumber, 'read')
        
        self.vessels = vascularNetwork.vessels
        self.boundaryConditions = vascularNetwork.boundaryConditions
        
        treeTraverseList_sorted = vascularNetwork.treeTraverseList[:]
        treeTraverseList_sorted.sort()
        self.treeTraverseList_sorted = treeTraverseList_sorted
        self.boundaryVessels = vascularNetwork.boundaryVessels
        self.vascularNetwork = vascularNetwork
        print vascularNetwork.tsol
    
    def writeNetworkLogfile(self):
        
        fTemplate = open(self.templateFile, 'r')
        fLogFile = open(self.networkLogFile, 'w')
        for line in fTemplate:
            
            if "% Fill in title here" in line:
                self.writeTitle(fLogFile)
                            
            elif "% Fill model parameters here" in line:
                self.writeModelParameters(fLogFile)

            elif "% Fill lumped values here" in line:
                if self.vascularNetwork.initialsationMethod not in ['ConstantPressure']:
                    #self.writeLumpedValues(fLogFile)
                    self.findAverageValuesFromSolution()
                    self.writeLumpedAndSolutionValues(fLogFile)
            else:
                fLogFile.write(line)
        
        fTemplate.close()
        fLogFile.close()
        

    def writeTitle(self, fLogFile):
        titleLine = ''.join([r'\section{Network : ',self.networkName, ', Datanumber: ', self.dataNumber, '}'])
        fLogFile.write(titleLine)
        fLogFile.write('\n')        
    

    def writeModelParameters(self, fLogFile):
        
        for n, vesselId in enumerate(self.treeTraverseList_sorted):
            
            if vesselId in self.boundaryVessels:
                bc = self.boundaryConditions[vesselId]
                if len(bc)>1:
                    bc = bc[1]
                else:
                    bc = bc[0]
                try:
                    R = round(bc.Rc/1e9, 3)
                    Z = round(bc.Z/1e9, 3)
                    C = round(bc.C*1e10, 3)
                except:
                    print "Warning ( moduleLogfile line 50): found bc but could not extract WK3 parameters"
            
            else:
                R = ""
                Z = ""
                C = ""
                
            name = self.vessels[vesselId].name
            l = round(self.vessels[vesselId].length*100, 3)
            rProx = round(self.vessels[vesselId].radiusProximal*100, 3)
            rDist = round(self.vessels[vesselId].radiusDistal*100, 3)
            
            tableLine = " {0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} ".format(vesselId, name, l, rProx, rDist, Z, R, C)
            
            fLogFile.write(tableLine)
            fLogFile.write(r'\\')
            fLogFile.write('\n')
            
            
            if n == 54 and len(self.treeTraverseList_sorted)>55:
                startNewtableLine1 = r"\begin{tabular}{ll*6{c}}"
                startNewtableLine2 = r"\hline"
                startNewtableLine3 = r"Id & Name & Length & p. Radius & d. Radius & R1/Z & R2 & C \\"
                startNewtableLine4 = r" & & [$cm$] & [$cm$] & [$cm$] & [$Pa\, s\, m^{-3}$] & [$Pa\, s\, m^{-3}$] &  [$m^{3}\, {Pa}^{-1}$] \\"
                startNewtableLine5 = r" & &  &  &  & $\cdot 10^{-9}$& $\cdot 10^{-9}$ & $\cdot 10^{10}$  \\"
                startNewtableLine6 = r"\hline"
                
                spesificLines = [startNewtableLine1, startNewtableLine2, startNewtableLine3, startNewtableLine4, startNewtableLine5, startNewtableLine6]
                
                self.endAndStartTable(fLogFile, spesificLines)
                

    def writeLumpedValues(self, fLogFile):
        initialValues = self.vascularNetwork.initialValues
        Rcum = self.vascularNetwork.Rcum
        root = self.vascularNetwork.root
        Q_root = initialValues[root]['Flow'][0]
        for n, vesselId in enumerate(self.treeTraverseList_sorted):
            [Pin, Pout] = initialValues[vesselId]['Pressure'] 
            [Qin, Qout] = initialValues[vesselId]['Flow']
            R = Rcum[vesselId]
            
            name = self.vessels[vesselId].name
            Qin_ml_s = round(Qin*1e6, 3)
            Qin_perc = round(Qin/Q_root, 3)
            Pin, Pout = round(Pin/133.32, 3), round(Pout/133.32, 3)
            R = round(R*1e-9, 3)
            
            tableLine1 = " {0} & {1} & {2} & {3} & {4} & {5} & {6} ".format(vesselId, name, Qin_ml_s, Qin_perc, Pin, Pout, R)
            tableLine2 = "& & & & &"
            fLogFile.write(tableLine1)
            fLogFile.write(tableLine2)
            fLogFile.write(r'\\')
            fLogFile.write('\n')
            
            
            if n == 54 and len(self.treeTraverseList_sorted)>55:
                 
                startNewtableLine1 = r"\begin{tabular}{l|*5{c}|*5{c}}"
                startNewtableLine2 = r"\hline"
                startNewtableLine3 = r"Id & $Q$ & $Q$ & $P_{in}$ & $P_{out}$ & $R_{cum}$ & $Q$ & $Q$ & $P_{in}$ & $P_{out}$ & $R_{cum}$ \\"
                startNewtableLine4 = r" & [$ml \, s^{-1}$] & [$\%$] & [$mmHg$] & [$mmHg$] & [$Pa\, s\, m^{-3}$] & [$ml \, s^{-1}$] & [$\%$] & [$mmHg$] & [$mmHg$] & [$Pa\, s\, m^{-3}$]  \\"
                startNewtableLine5 = r" &  &  &  & & $\cdot 10^{-9}$ &  &  &  & & $\cdot 10^{-9}$\\"
                startNewtableLine6 = r"\hline"

                spesificLines = [startNewtableLine1, startNewtableLine2, startNewtableLine3, startNewtableLine4, startNewtableLine5, startNewtableLine6]
                
                self.endAndStartTable(fLogFile, spesificLines)
    
             
    def writeLumpedAndSolutionValues(self, fLogFile):
        lumpedValues = self.vascularNetwork.lumpedValues
        averageValues = self.averageValues
        Rcum = self.vascularNetwork.Rcum
        root = self.vascularNetwork.root
        Q_root = lumpedValues[root]['Flow'][0]
        Q_root_avg = averageValues[root]['Qin']
        for n, vesselId in enumerate(self.treeTraverseList_sorted):
            [Pin, Pout] = lumpedValues[vesselId]['Pressure'] 
            [Qin, Qout] = lumpedValues[vesselId]['Flow']
            R = Rcum[vesselId]
            
            Qin_ml_s = round(Qin*1e6, 3)
            Qin_perc = round(100*Qin/Q_root, 3)
            Pin, Pout = round(Pin/133.32, 3), round(Pout/133.32, 3)
            R = round(R*1e-9, 3)
            
            tableLine = " {0} & {1} & {2} & {3} & {4} & {5} ".format(vesselId, Qin_ml_s, Qin_perc, Pin, Pout, R)
            
            Pin, Pout = averageValues[vesselId]['Pin'], averageValues[vesselId]['Pout']
            Qin = averageValues[vesselId]['Qin']
            R = averageValues[vesselId]['R']
            
            Qin_ml_s = round(Qin*1e6, 3)
            Qin_perc = round(100*Qin/Q_root_avg, 3)
            Pin, Pout = round(Pin/133.32, 3), round(Pout/133.32, 3)
            R = round(R*1e-9, 3)
            
            tableLine2 = "& {0} & {1} & {2} & {3} & {4} ".format(Qin_ml_s, Qin_perc, Pin, Pout, R)
            
            fLogFile.write(tableLine)
            fLogFile.write(tableLine2)
            fLogFile.write(r'\\')
            fLogFile.write('\n')
            
            
            if n == 54 and len(self.treeTraverseList_sorted)>55:
                 
                startNewtableLine1 = r"\begin{tabular}{l|*5{c}|*5{c}}"
                startNewtableLine2 = r"\hline"
                startNewtableLine3 = r"Id & $Q$ & $Q$ & $P_{in}$ & $P_{out}$ & $R_{cum}$ & $Q$ & $Q$ & $P_{in}$ & $P_{out}$ & $R_{cum}$ \\"
                startNewtableLine4 = r" & [$ml \, s^{-1}$] & [$\%$] & [$mmHg$] & [$mmHg$] & [$Pa\, s\, m^{-3}$] & [$ml \, s^{-1}$] & [$\%$] & [$mmHg$] & [$mmHg$] & [$Pa\, s\, m^{-3}$]  \\"
                startNewtableLine5 = r" &  &  &  & & $\cdot 10^{-9}$ &  &  &  & & $\cdot 10^{-9}$\\"
                startNewtableLine6 = r"\hline"

                spesificLines = [startNewtableLine1, startNewtableLine2, startNewtableLine3, startNewtableLine4, startNewtableLine5, startNewtableLine6]
                
                self.endAndStartTable(fLogFile, spesificLines)
    
                
    def endAndStartTable(self, fLogFile, spesificLines):
        lineList = [r"\end{tabular}", r"\end{table}", r"\begin{table}"]
        
        for line in spesificLines:
            lineList.append(line)
        
        for line in lineList:
            fLogFile.write(line)
            fLogFile.write("\n")
        
    
    def findAverageValuesFromSolution(self):
        
        hdf5File = h5py.File(self.solutionFile, 'r')
        time = hdf5File['VascularNetwork']['simulationTime'][:] - hdf5File['VascularNetwork']['simulationTime'][0]
        dt = time[1] - time[0]
        root = self.vascularNetwork.root
        freq = self.boundaryConditions[root][0].freq
        period = 1./freq
        N = int(period/dt)
        N = len(time) - N - 1
        time = time[N:]
        
        print "starting to load solutiondata"
        print "averaging between t_start = {0} and t_end = {1}".format(time[0], time[-1])
        averageValues = {}
        for vesselName in hdf5File['vessels'].keys():
            tmpValues = {}
            vesselId = vesselName.split(' - ')[-1]
            vesselId = int(vesselId)
            
            startP  = hdf5File['vessels'][vesselName]['Psol'][N:, 0]
            endP =  hdf5File['vessels'][vesselName]['Psol'][N:, -1]
            
            
            startQ  = hdf5File['vessels'][vesselName]['Qsol'][N:, 0]
            #endQ =  hdf5File['vessels'][vesselName]['Qsol'][-N, -1]
            
            Pin, Pout = self.findMean(time, startP), self.findMean(time, endP)
            Qin = self.findMean(time, startQ)
            R = (Pin - self.vascularNetwork.venousPool.P[0])/Qin
            tmpValues['Pin'], tmpValues['Pout'] = Pin, Pout
            tmpValues['Qin'] = Qin
            tmpValues['R'] = R
            averageValues[vesselId] = tmpValues
        
        self.averageValues = averageValues
            
            
    def findMean(self, x, f):
        
        F = simps(f, x)
        
        f_mean = F/(x[-1] - x[0])
        
        return f_mean
        
        




    



