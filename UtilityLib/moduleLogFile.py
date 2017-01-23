import os,sys
import h5py
import numpy as np
from scipy.integrate import simps



# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))


import moduleFilePathHandler as mFPH

class NetworkLogFile:
    
    def __init__(self, vascularNetwork, dataNumber = "xxx", networkLogFile = None, dt=None, CpuTimeInit=[None, None], CpuTimeSolve=[None, None]):
        
        self.networkName = vascularNetwork.getVariableValue('name')
        self.dataNumber = dataNumber
        if networkLogFile == None:
            self.networkLogFile =  mFPH.getFilePath('networkLogFile', self.networkName, dataNumber, 'write')
            #print networkLogFile
        else:
            self.networkLogFile = networkLogFile
        
        self.solutionFileDirectory = mFPH.getDirectory('solutionFileDirectory', self.networkName, dataNumber, 'r')
        
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
        
        self.CpuTimeInit = CpuTimeInit
        self.CpuTimeSolve = CpuTimeSolve

        self.totalTime = vascularNetwork.totalTime
        self.timeSaveBegin = vascularNetwork.timeSaveBegin
        self.globalFluid = vascularNetwork.globalFluid
        self.dt = dt
        print self.dt
        self.CpuTimeInit = CpuTimeInit
        self.CpuTimeSolve = CpuTimeSolve
        self.venousPressure = vascularNetwork.venousPool.P[0]

    
    def writeNetworkLogfile(self, compileLogFile=False, deleteAuxiliary=False):
        
        fTemplate = open(self.templateFile, 'r')
        fLogFile = open(self.networkLogFile, 'w')
        for line in fTemplate:
            
            if "% Fill in title here" in line:
                self.writeTitle(fLogFile)

            elif "% Fill solving data here" in line:
                self.writeSolvingData(fLogFile)

            elif "% Fill numerical data here" in line:
                self.writeNumericalData(fLogFile)
                
            elif "% Fill total vlues here" in line:
                self.writeTotalValues(fLogFile)
                
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
        
        if compileLogFile:
            self.compilePdf(deleteAuxiliary=deleteAuxiliary)
    
    def compilePdf(self, deleteAuxiliary=False):
        
        filesInDirBefore = os.listdir(self.solutionFileDirectory)
        print filesInDirBefore
        compileString = ''.join(['pdflatex -output-directory=',self.solutionFileDirectory, ' ', self.networkLogFile])
        os.system(compileString)
        #subprocess.Popen(compileString, shell=True)
        filesInDirAfter = os.listdir(self.solutionFileDirectory)
        
        if deleteAuxiliary:
            for filex in filesInDirAfter:
                if filex not in filesInDirBefore and 'pdf' not in filex:
                    completeFilePath = ''.join([self.solutionFileDirectory, '/', filex])
                    os.remove(completeFilePath)
                
        
    def writeTitle(self, fLogFile):
        #titleLine = ''.join([r'\section{Network : ',self.networkName, ', Datanumber: ', self.dataNumber, '}'])
        titleLine =''.join([r'\Large Network : ',self.networkName, ', Datanumber: ', self.dataNumber, r'\\'])
        fLogFile.write(titleLine)
        fLogFile.write('\n')


    def writeSolvingData(self, fLogFile):
        
        fLogFile.write("dt (ms): & {0}".format(self.dt*1000))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("Number of domains: & {0}".format(len(self.treeTraverseList_sorted)))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("Simulation time: & {0}".format(self.totalTime))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("Time save begin: & {0}".format(self.timeSaveBegin))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("CPU time solving: & {0} min, {1} sec".format(self.CpuTimeSolve[0], self.CpuTimeSolve[1]))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("CPU time initialization: & {0} min, {1} sec".format(self.CpuTimeInit[0], self.CpuTimeInit[1]))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("Velocity profile parameter: & {0}".format(self.globalFluid['gamma']))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("Density ($Kg/m^3$): & {0}".format(self.globalFluid['rho']))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("Blood viscosity (mPa s): & {0}".format(self.globalFluid['my']*1000))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("Outflow Windkessel pressure (Pa): & {0}".format(self.venousPressure))
        fLogFile.write('\n')
        

    def writeNumericalData(self, fLogFile):
        
        self.findNumericalValuesFromSolution()
        
        maxdz = 0
        mindz = 100
        maxNode = 0
        minNode = 10000
        maxCFL = 0
        minCFL = 1000
        minCFLmin = 1000
        maxc = 0
        minc = 1000
        
        maxdzID = None
        mindzID = None
        maxNodeID = None
        minNodeID = None
        maxCFLID = None
        minCFLID = None
        minCFLminID = None
        maxcID = None
        mincID = None
        
        for vesselID in self.treeTraverseList_sorted:
            
            cmin_in, cmax_in = self.numericalValues[vesselID]['cin']
            cmin_out, cmax_out = self.numericalValues[vesselID]['cout']
            dz = self.numericalValues[vesselID]['dz']
            Nodes = self.numericalValues[vesselID]['N']
            startCFLmax = self.numericalValues[vesselID]['CFLin']
            endCFLmax = self.numericalValues[vesselID]['CFLout']
            
            cValues = [cmin_in, cmax_in, cmin_out, cmax_out]
            cminVessel = min(cValues)
            cmaxVessel = max(cValues)
            
            CFLminVessel = min([startCFLmax, endCFLmax])
            CFLmaxVessel = max([startCFLmax, endCFLmax])
            
            if dz > maxdz:
                maxdz = dz
                maxdzID = vesselID
            if dz < mindz:
                mindz = dz
                mindzID = vesselID
                
            if Nodes > maxNode:
                maxNode = Nodes
                maxNodeID = vesselID
            if Nodes < minNode:
                minNode = Nodes
                minNodeID = vesselID

            if cmaxVessel > maxc:
                maxc = cmaxVessel
                maxcID = vesselID
            if cminVessel < minc:
                minc = cminVessel
                mincID = vesselID
                
            if CFLmaxVessel > maxCFL:
                maxCFL = CFLmaxVessel
                maxCFLID = vesselID
            if CFLmaxVessel < minCFL:
                minCFL = CFLmaxVessel
                minCFLID = vesselID
            if CFLminVessel < minCFLmin:
                minCFLmin = CFLminVessel
                minCFLminID = vesselID
            
        fLogFile.write("max CFL Values: & {0}, {1} ({2}, {3})".format(round(minCFL, 3), round(maxCFL, 3), minCFLID, maxCFLID))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("min CFL: & {0} ({1})".format(round(minCFLmin, 3), minCFLminID))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("max c (m/s): & {0} ({1})".format(round(maxc, 3), maxcID))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("min c (m/s): & {0} ({1})".format(round(minc, 3), mincID))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("max dz (mm): & {0} ({1})".format(round(maxdz*1000, 3), maxdzID))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("min dz (mm): & {0} ({1})".format(round(mindz*1000, 3), mindzID))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("max Nodes: & {0} ({1})".format(maxNode, maxNodeID))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write("min Nodes: & {0} ({1})".format(minNode, minNodeID))
        fLogFile.write('\n')

    
    def writeTotalValues(self, fLogFile):
        root = self.vascularNetwork.root
        R_total = self.vascularNetwork.Rcum[root]
        _R_total_bc = 0
        
        C_total = 0
        C_bc = 0
        C_vessels = 0
        
        
        for n, vesselId in enumerate(self.treeTraverseList_sorted):
            
            if vesselId in self.boundaryVessels:
                bc = self.boundaryConditions[vesselId]
                if len(bc)>1:
                    bc = bc[1]
                else:
                    bc = bc[0]
                R = bc.Rc
                Z = bc.Z
                C = bc.C
                
                Rt = R + Z
                C_total += C
                C_bc += C
                
                _R_total_bc += 1./Rt
            
            Cv = self.vessels[vesselId].Cv
            C_vessels += Cv
            C_total += Cv
        
        R_total_bc = 1./_R_total_bc
        
        fLogFile.write(r"Total arterial resistance ($Pa \, s \, m^{-3}$): & " + str(R_total))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write(r"Total peripheral resistance ($Pa \, s \, m^{-3}$): & " + str(R_total_bc))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write(r"Total arterial Compliance ($ m^{3} \, s^{-1}$): &  " + str(C_total))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write(r"Total peripheral Compliance ($ m^{3} \, s^{-1}$): &  " + str(C_bc))
        fLogFile.write(r' \\')
        fLogFile.write('\n')
        fLogFile.write(r"Total vessel Compliance ($ m^{3} \, s^{-1}$): &  " + str(C_vessels))
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
            
            
            Rv = round(self.vessels[vesselId].resistance/1e9, 5)
            Cv = round(self.vessels[vesselId].Cv*1e10, 3)
            cd_in = round(self.vessels[vesselId].cd_in, 2)
            cd_out =round(self.vessels[vesselId].cd_out, 2)

                
            name = self.vessels[vesselId].name
            l = round(self.vessels[vesselId].length*100, 2)
            rProx = round(self.vessels[vesselId].radiusProximal*1000, 2)
            rDist = round(self.vessels[vesselId].radiusDistal*1000, 2)
            
            
            tableLine = "{0} & {1} & {2} & {3} $\\rightarrow$ {4} & {5} $\\rightarrow$ {6} & {7} & {8} & {9} & {10} & {11}"
            tableLine = tableLine.format(vesselId, name, l, rProx, rDist, cd_in, cd_out, Cv, Rv, Z, R, C)
            fLogFile.write(tableLine)
            fLogFile.write(r'\\')
            fLogFile.write('\n')
            
            
            if n == 54 and len(self.treeTraverseList_sorted)>55:
                startNewtableLine1 = r"\begin{tabular}{ll*8{c}}"
                startNewtableLine2 = r"\hline"
                startNewtableLine3 = r"Id & Name & Length & $R_{in} \rightarrow R_{out} $  & $c_{in} \rightarrow c_{out} $ & $C_v$ & $R_v$ & R1/Z & R2 & C \\"
                startNewtableLine4 = r" & & [$cm$] & [$mm$] & [$m/s$] &  [$m^{3}\, {Pa}^{-1}$] & [$Pa\, s\, m^{-3}$] & [$Pa\, s\, m^{-3}$] & [$Pa\, s\, m^{-3}$] &  [$m^{3}\, {Pa}^{-1}$] \\"
                startNewtableLine5 = r" & &  &  &  & $\cdot 10^{10}$ & $\cdot 10^{-9}$ & $\cdot 10^{-9}$ & $\cdot 10^{-9}$ & $\cdot 10^{10}$  \\"
                startNewtableLine6 = r"\hline"
                
                spesificLines = [startNewtableLine1, startNewtableLine2, startNewtableLine3, startNewtableLine4, startNewtableLine5, startNewtableLine6]
                
                self.endAndStartTable(fLogFile, spesificLines, resizeTable=True)
                

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
                
                self.endAndStartTable(fLogFile, spesificLines, resizeTable=True)
    
             
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
                
                self.endAndStartTable(fLogFile, spesificLines, resizeTable=False)
    
                
    def endAndStartTable(self, fLogFile, spesificLines, resizeTable=False):
        if resizeTable==True:
            lineList = [r"\end{tabular}}", r"\end{table}", r"\begin{table}", r"\resizebox{\textwidth}{!}{"]
        else:
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
        

    def findNumericalValuesFromSolution(self):
        
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
        numericalValues = {}
        for vesselName in hdf5File['vessels'].keys():
            tmpValues = {}
            vesselId = vesselName.split(' - ')[-1]
            vesselId = int(vesselId)
            
            startP  = hdf5File['vessels'][vesselName]['Psol'][N:, 0]
            endP =  hdf5File['vessels'][vesselName]['Psol'][N:, -1]
            
            startA  = hdf5File['vessels'][vesselName]['Asol'][N:, 0]
            endA =  hdf5File['vessels'][vesselName]['Asol'][N:, -1]
            
            startc = self.vessels[vesselId].waveSpeed(startA, self.vessels[vesselId].C_nID([startP], 0))
            endc = self.vessels[vesselId].waveSpeed(endA, self.vessels[vesselId].C_nID([endP], -1))
            
            startQ  = hdf5File['vessels'][vesselName]['Qsol'][N:, 0]
            endQ =  hdf5File['vessels'][vesselName]['Qsol'][-N, -1]
            
            startU = startQ/startA
            endU = endQ/endA
            
            dz = self.vessels[vesselId].dz[0]
            Nodes = self.vessels[vesselId].N
            
            startCFL = self.dt*(startU + startc)/dz
            endCFL = self.dt*(endU + endc)/dz
            startCFLmax = np.amax(startCFL)
            endCFLmax = np.amax(endCFL)
            
            tmpValues['cin']= [np.amin(startc), np.amax(startc)]
            tmpValues['cout'] = [np.amin(endc), np.amax(endc)]
            tmpValues['dz'] = dz
            tmpValues['N'] = Nodes
            tmpValues['CFLin']= startCFLmax
            tmpValues['CFLout']= endCFLmax
            
            numericalValues[vesselId] = tmpValues
        
        self.numericalValues = numericalValues
            
            
    def findMean(self, x, f):
        
        F = simps(f, x)
        
        f_mean = F/(x[-1] - x[0])
        
        return f_mean
        
        




    



