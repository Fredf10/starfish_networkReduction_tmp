import os,sys
import h5py
import numpy as np
from scipy.integrate import simps
from scipy import interpolate



# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))


import moduleFilePathHandler as mFPH

class NetworkLogFile:
    
    def __init__(self, vascularNetwork, dataNumber = "xxx", dataNumberCompare=None,
                 networkLogFile=None, solutionFileDirectory=None, solutionFile = None,
                 dt=None, CpuTimeInit=[None, None], CpuTimeSolve=[None, None]):
        
        self.networkName = vascularNetwork.getVariableValue('name')
        self.root = vascularNetwork.root
        self.dataNumber = dataNumber
        if networkLogFile == None:
            self.networkLogFile =  mFPH.getFilePath('networkLogFile', self.networkName, dataNumber, 'write')
            #print networkLogFile
        else:
            self.networkLogFile = networkLogFile
        if solutionFileDirectory == None:
            self.solutionFileDirectory = mFPH.getDirectory('solutionFileDirectory', self.networkName, dataNumber, 'r')
        else:
            self.solutionFileDirectory = solutionFileDirectory
            
        TemplateLogFileDirectory = mFPH.getDirectory('networkXmlFileTemplateDirectory', 'singleVessel_template', dataNumber, 'read')
        #print TemplateLogFileDirectory
        self.templateFile = ''.join([TemplateLogFileDirectory,'/logFile_template.tex'])
        
        if solutionFile == None:
            self.solutionFile = mFPH.getFilePath('solutionFile', self.networkName, dataNumber, 'read')
        
        else:
            self.solutionFile = solutionFile
        if dataNumberCompare != None:
            self.solutionFileCompare = mFPH.getFilePath('solutionFile', self.networkName, dataNumberCompare, 'read')
            self.compare = True
        else:
            self.compare = False
            
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
        #print self.dt
        self.CpuTimeInit = CpuTimeInit
        self.CpuTimeSolve = CpuTimeSolve
        self.venousPressure = vascularNetwork.centralVenousPressure

    
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
                try:
                    self.writeTotalValues(fLogFile)
                except:
                    pass
            
            elif "% Fill periodic values here" in line:
                self.findPeriodicConvergenceValuesFromSolution()
                self.writePeriodicConvergenceValues(fLogFile)
                fLogFile.write(r'\end{tabular}')
                fLogFile.write('\n')
                fLogFile.write(r'\end{table}')
                fLogFile.write('\n')
                if self.compare:
                    self.compareWithOtherSolution()
                    fLogFile.write(r'\begin{table}[h!]')
                    self.writeCompareWithOtherSolutionValues(fLogFile)
                    
                    fLogFile.write(r'\end{tabular}')
                    fLogFile.write('\n')
                    fLogFile.write(r'\end{table}')
                    fLogFile.write('\n')
                    
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
                    for bcTmp in self.boundaryConditions[vesselId]:
                        if bcTmp.type != 1:
                            
                            bc = bcTmp
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
    
    def writePeriodicConvergenceValues(self, fLogFile, roundTo=3):
        
        convergenceValues = self.convergenceValues
        N = len(convergenceValues[self.vascularNetwork.root]['deltaP'])
        
        fLogFile.write(r'\begin{tabular}{')
        for n in range(N + 1):
            fLogFile.write('l')
        fLogFile.write(r'}')
        fLogFile.write('\n')
        
        fLogFile.write('$\Delta{P}_{sys}^1$')
        fLogFile.write(r'    &    ')
        for deltaP in convergenceValues[self.vascularNetwork.root]['deltaP']:
            
            deltaPround = round(deltaP/133.32, roundTo)
            
            if deltaPround > 1.0:
                fLogFile.write(r"\textcolor{red}{" + str(deltaPround) + "}")
            elif deltaPround > 0.1:
                fLogFile.write(r"\textcolor{orange}{" + str(deltaPround) + "}")
            else:
                fLogFile.write(r"\textcolor{green}{" + str(deltaPround) + "}")
            if deltaP != convergenceValues[self.root]['deltaP'][-1]:
                fLogFile.write(r'    &    ')
        fLogFile.write(r'\\')
        fLogFile.write('\n')

        fLogFile.write('$\Delta{P}_{sys}^{max}$')
        fLogFile.write(r'    &    ')
        for n in range(N):
            deltaP_max = 0
            for vesselId in self.treeTraverseList_sorted:
                deltaP = convergenceValues[vesselId]['deltaP'][n]
                
                if deltaP > deltaP_max:
                    deltaP_max = deltaP
                    vesselId_max = vesselId
                    
            deltaPround = round(deltaP_max/133.32, roundTo)
            if deltaPround > 1.0:
                fLogFile.write(r"\textcolor{red}{" + str(deltaPround) + "}")
                fLogFile.write(" ({0})".format(vesselId))
            elif deltaPround > 0.1:
                fLogFile.write(r"\textcolor{orange}{" + str(deltaPround) + "}")
                fLogFile.write(" ({0})".format(vesselId))
            else:
                fLogFile.write(r"\textcolor{green}{" + str(deltaPround) + "}")
                fLogFile.write(" ({0})".format(vesselId))
            if n != N - 1:
                fLogFile.write(r'    &    ')
        fLogFile.write(r'\\')
        fLogFile.write('\n')
        fLogFile.write(r'\hline')
        fLogFile.write('\n')


        fLogFile.write('$\epsilon{P}_{avg}^1$')
        fLogFile.write(r'    &    ')
        for epsilon in convergenceValues[self.vascularNetwork.root]['epsilonAvgP']:
            
            epsilonRound = round(100*epsilon, roundTo)
            
            if epsilonRound > 1.0:
                fLogFile.write(r"\textcolor{red}{" + str(epsilonRound) + "}")
            elif epsilonRound > 0.1:
                fLogFile.write(r"\textcolor{orange}{" + str(epsilonRound) + "}")
            else:
                fLogFile.write(r"\textcolor{green}{" + str(epsilonRound) + "}")
            if epsilon != convergenceValues[self.vascularNetwork.root]['epsilonAvgP'][-1]:
                fLogFile.write(r'    &    ')
        fLogFile.write(r'\\')
        fLogFile.write('\n')

        fLogFile.write('$\epsilon{P}_{avg}^{max}$')
        fLogFile.write(r'    &    ')
        for n in range(N):
            epsilon_max = 0
            for vesselId in self.treeTraverseList_sorted:
                epsilon = convergenceValues[vesselId]['epsilonAvgP'][n]
                
                if epsilon > epsilon_max:
                    epsilon_max =epsilon
                    
            epsilonRound = round(100*epsilon_max, roundTo)
            if epsilonRound > 1.0:
                fLogFile.write(r"\textcolor{red}{" + str(epsilonRound) + "}")
            elif epsilonRound > 0.1:
                fLogFile.write(r"\textcolor{orange}{" + str(epsilonRound) + "}")
            else:
                fLogFile.write(r"\textcolor{green}{" + str(epsilonRound) + "}")
            if n != N - 1:
                fLogFile.write(r'    &    ')
        fLogFile.write(r'\\')
        fLogFile.write('\n')


    def writeCompareWithOtherSolutionValues(self, fLogFile, roundTo=3):
        
        compareOtherSolutionValues = self.compareOtherSolutionValues
        
        fLogFile.write(r'\begin{tabular}{')
        
        fLogFile.write('ll')
        fLogFile.write(r'}')
        fLogFile.write('\n')
        
        vesselToCompare = 1
        fLogFile.write('$\epsilon{P}_{avg}^{' + str(vesselToCompare) + '}$')
        fLogFile.write(r'    &    ')
        
        fLogFile.write(str(round(self.compareOtherSolutionValues[vesselToCompare]['epsilonAvgP']*100, roundTo)))
        fLogFile.write(r'\\')
        fLogFile.write('\n')


        epsilonP_max = 0
        vesselMax = -1
        for vesselId in self.treeTraverseList_sorted:
            epsilonP = compareOtherSolutionValues[vesselId]['epsilonAvgP']
            
            if epsilonP > epsilonP_max:
                epsilonP_max = epsilonP
                vesselMax = vesselId

        fLogFile.write('$\epsilon{P}_{avg}^{max=' + str(vesselMax) + '}$')
        fLogFile.write(r'    &    ')
        epsilonPround = round(epsilonP_max*100, roundTo)
        fLogFile.write(str(epsilonPround))
            
        fLogFile.write(r'\\')
        fLogFile.write('\n')
        fLogFile.write(r'\hline')
        fLogFile.write('\n')
        
        
    def writeModelParameters(self, fLogFile):
        
        for n, vesselId in enumerate(self.treeTraverseList_sorted):
            
            if vesselId in self.boundaryVessels:
                bc = self.boundaryConditions[vesselId]
                if len(bc)>1:
                    for bcTmp in self.boundaryConditions[vesselId]:
                        if bcTmp.type != 1:
                            
                            bc = bcTmp
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

                
            name = self.vessels[vesselId].name.replace('_', ' ')
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
        meanValueSum = 0
        varianceSum = 0
        
        epsilonMax = 0
        epsilonSum = 0

        epsilonMax_Q = 0
        epsilonSum_Q = 0
        
        round_to = 5
        for n, vesselId in enumerate(self.treeTraverseList_sorted):
            [Pin, Pout] = lumpedValues[vesselId]['Pressure'] #['MeanValues']
            [Qin, Qout] = lumpedValues[vesselId]['Flow'] # ['MeanValues']
            R = Rcum[vesselId]
            
            Qin_ml_s_lump = round(Qin*1e6, round_to)
            Qin_perc_lump = round(100*Qin/Q_root, round_to)
            Pin_lump, Pout_lump = round(Pin/133.32, round_to), round(Pout/133.32, round_to)
            R_lump = round(R*1e-9, round_to)
            
           
            
            Pin, Pout = averageValues[vesselId]['Pin'], averageValues[vesselId]['Pout']
            Qin = averageValues[vesselId]['Qin']
            Qin_max = averageValues[vesselId]['Qin_max']
#             R_in, R_out = averageValues[vesselId]['R']
#             lumpedValues[vesselId]['R'] = [R_in, R_out]
#             lumpedValues[vesselId]['R_new'] = [R_in, R_out]
            Qin_max = round(Qin_max*1e6, round_to)
            Qin_ml_s = round(Qin*1e6, round_to)
            Qin_perc = round(100*Qin/Q_root_avg, round_to)
            Pin, Pout = round(Pin/133.32, round_to), round(Pout/133.32, round_to)
            #R_in = round(R_in*1e-9, 3)
            
            P_rel = 100*round(Pin_lump/Pin, round_to)
            
            epsilon = abs(Pin_lump - Pin)/Pin
            epsilonSum += ((Pin_lump - Pin)/Pin)**2

            epsilon_Q = abs((Qin_ml_s_lump - Qin_ml_s)/Qin_max)
            epsilonSum_Q += ((Qin_ml_s_lump - Qin_ml_s)/Qin_max)**2
            
            if epsilon > epsilonMax:
                epsilonMax = epsilon

            if epsilon_Q > epsilonMax_Q:
                epsilonMax_Q = epsilon_Q
            
            meanValueSum += P_rel
            varianceSum += (P_rel - 100)**2
            
            tableLine = " {0} & {1} & {2} & {3} ({4}) & {5} & {6} ".format(vesselId, Qin_ml_s_lump, Qin_perc_lump, Pin_lump, P_rel, Pout_lump, R_lump)
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

#         print 'epsilonMax: ', epsilonMax
#         print 'ST: ', np.sqrt(epsilonSum/(n + 1))
# 
#         print 'epsilonMax_Q: ', epsilonMax_Q
#         print 'ST_Q: ', np.sqrt(epsilonSum_Q/(n + 1))
#         
#         raw_input('press enter to continue')

    
                
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
        for bc in self.vascularNetwork.boundaryConditions[root]:
            if bc.type == 1:
                freq = bc.freq

        period = 1./freq
        
        N = int((period/dt))
        
        nCycles = int((time[-1]/period))
        t_start = period*(nCycles - 1)
        t_end = period*nCycles
        time_compare = np.linspace(t_start, t_end, N + 1)
        
        print "starting to load solutiondata for calculating averageValues"
        print "averaging between t_start = {0} and t_end = {1}".format(t_start, t_end)
        averageValues = {}
        for vesselName in hdf5File['vessels'].keys():
            tmpValues = {}
            vesselId = vesselName.split(' - ')[-1]
            vesselId = int(vesselId)
            
            startP  = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Psol'][:, 0], time_compare)
            endP =  self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Psol'][:, -1], time_compare)
            
            
            startQ  = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Qsol'][:, 0], time_compare)
            endQ =  self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Qsol'][:, -1], time_compare)
            
            Pin, Pout = self.findMean(time_compare, startP), self.findMean(time_compare, endP)
            Qin = self.findMean(time_compare, startQ)
            Qout = self.findMean(time_compare, endQ)
            R_in = (Pin - self.venousPressure)/Qin
            R_out = (Pout - self.venousPressure)/Qout
            tmpValues['Pin'], tmpValues['Pout'] = Pin, Pout
            tmpValues['Qin'] = Qin
            tmpValues['Qin_max'] = np.max(startQ)
            tmpValues['R'] = [R_in, R_out]
            tmpValues['Pressure'] = [Pin, Pout]
            tmpValues['Flow'] = [Qin, Qout]
            
            averageValues[vesselId] = tmpValues
        
        self.averageValues = averageValues
    

    def findNumericalValuesFromSolution(self):
        
        hdf5File = h5py.File(self.solutionFile, 'r')
        time = hdf5File['VascularNetwork']['simulationTime'][:] - hdf5File['VascularNetwork']['simulationTime'][0]
        dt = time[1] - time[0]
        root = self.vascularNetwork.root
        for bc in self.vascularNetwork.boundaryConditions[root]:
            if bc.type == 1:
                freq = bc.freq

        period = 1./freq
        N = int(period/dt)
        N = len(time) - N - 1
        time = time[N:]
        
        print "starting to load solutiondata"
        print "averaging between t_start = {0} and t_end = {1}".format(time[0], time[-1])
        numericalValues = {}
        SystolicDiastolicDict = {'Systolic':{}, 'Diastolic':{}}
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
            SystolicDiastolicDict['Systolic'][vesselId] = {}
            SystolicDiastolicDict['Diastolic'][vesselId] = {}
            SystolicDiastolicDict['Systolic'][vesselId]['Pressure'] = [np.amax(startP), np.amax(endP)]
            SystolicDiastolicDict['Diastolic'][vesselId]['Pressure'] = [np.amin(startP), np.amin(endP)]
            
            createFig = False
            if createFig:
                import matplotlib.pylab as plt
                initialValues = self.vascularNetwork.initialValues
                Wkname = initialValues[vesselId]['WK_name']
                t_Wk = initialValues[vesselId]['tWK']
                P_Wk = initialValues[vesselId]['PWK']
                Q_Wk = initialValues[vesselId]['QWK']
                t_num = hdf5File['VascularNetwork']['simulationTime'][:] - hdf5File['VascularNetwork']['simulationTime'][0]
                P_num =hdf5File['vessels'][vesselName]['Psol'][:, -1]
                Q_num =hdf5File['vessels'][vesselName]['Qsol'][:, -1]
                fig, ax1 = plt.subplots()
                ax1.plot(t_Wk, P_Wk/133.32, 'b')
                ax1.plot(t_num, P_num/133.32, 'k--')
                ax1.set_xlabel('time (s)')
                # Make the y-axis label, ticks and tick labels match the line color.
                ax1.set_ylabel('P', color='b')
                
                ax2 = ax1.twinx()
                ax2.plot(t_Wk, 1e6*Q_Wk, 'b')
                ax2.plot(t_num, 1e6*Q_num, 'k--')
                ax2.set_ylabel('Q', color='r')
                plt.tight_layout()
                plt.savefig('/home/fredrik/starfish_working_directory/Full55ModelDev/fig/{0}_{1}.png'.format(Wkname, vesselId))
                plt.close()
        
            #import pickle
            #pickle.dump(SystolicDiastolicDict, open("/home/fredrik/Documents/Backup_old_desktop/Documents/git/NTNU_KCL/apps/dirAppDev/data/Full55ModelDev/systolicDiastolicPressures.p", "wb"))
        self.numericalValues = numericalValues

    def findPeriodicConvergenceValuesFromSolution(self):
        
        hdf5File = h5py.File(self.solutionFile, 'r')
        time = hdf5File['VascularNetwork']['simulationTime'][:] - hdf5File['VascularNetwork']['simulationTime'][0]
        dt = time[1] - time[0]


        root = self.vascularNetwork.root
        for bcTmp in self.boundaryConditions[root]:
            if bcTmp.type == 1:
                
                freq = bcTmp.freq
        period = 1./freq
        N = int((period/dt))
        
        nCycles = int(((time[-1] - time[0])/period))

        convergenceValues = {}
        for vesselName in hdf5File['vessels'].keys():
            tmpValues = {}
            vesselId = vesselName.split(' - ')[-1]
            vesselId = int(vesselId)
            maxValues = []
            errorValues = []
            startP_all  = hdf5File['vessels'][vesselName]['Psol'][:, 0]
            startP_list = []
            t_start = time[0]
            t_end = t_start + period
            tck = interpolate.splrep(time, startP_all)
            for n in range(nCycles):
            
                time_tmp = np.linspace(t_start, t_end, N + 1)
                startP = interpolate.splev(time_tmp, tck)

                maxValues.append(np.max(startP))
                
                startP_list.append(startP)
                
                t_start += period
                t_end += period
            
            for n in range(nCycles - 1):
                e = self.calcEpsilonAvg(startP_list[-1], startP_list[n], data_type='P')
                errorValues.append(e)
            deltaP = np.abs(np.array(maxValues[1:]) - np.array(maxValues[:-1]))
            
            tmpValues['maxValues'] = maxValues
            tmpValues['deltaP'] = deltaP
            tmpValues['epsilonAvgP'] = errorValues
            convergenceValues[vesselId] = tmpValues
        
        self.convergenceValues = convergenceValues
    
    def compareWithOtherSolution(self):
        
        hdf5File = h5py.File(self.solutionFile, 'r')
        hdf5File2 = h5py.File(self.solutionFileCompare, 'r')
        
        time = hdf5File['VascularNetwork']['simulationTime'][:]
        time2 = hdf5File2['VascularNetwork']['simulationTime'][:]
        dt = time[1] - time[0]
        root = self.vascularNetwork.root
        freq = self.boundaryConditions[root][0].freq
        period = 1./freq
        N = int((period/dt))
        
        nCycles = int((time[-1]/period))
        nCycles2 = int((time2[-1]/period))
        
        #nCycles = min([nCycles, nCycles2])
        
        t_start = period*(nCycles - 1)
        t_end = period*nCycles

        t_start_2 = period*(nCycles2 - 1)
        t_end_2 = period*nCycles2
        
        time_compare = np.linspace(t_start, t_end, N + 1)
        time_compare2 = np.linspace(t_start_2, t_end_2, N + 1)
        compareOtherSolutionValues = {}
        for vesselName in hdf5File['vessels'].keys():
            tmpValues = {}
            vesselId = vesselName.split(' - ')[-1]
            vesselId = int(vesselId)

            startP_all  = hdf5File['vessels'][vesselName]['Psol'][:, 0]
            startP2_all  = hdf5File2['vessels'][vesselName]['Psol'][:, 0]

            tck = interpolate.splrep(time, startP_all)
            tck2 = interpolate.splrep(time2, startP2_all)
        
            
            startP = interpolate.splev(time_compare, tck)
            startP2 = interpolate.splev(time_compare2, tck2)

            
            e = self.calcEpsilonAvg(startP2, startP, data_type='P')
                

            tmpValues['epsilonAvgP'] = e
            compareOtherSolutionValues[vesselId] = tmpValues
        
        self.compareOtherSolutionValues = compareOtherSolutionValues
        

    def mySplineInterpolater(self, x, y, x_new):
        tck = interpolate.splrep(x, y, s=0)

        y_new = interpolate.splev(x_new, tck, der=0)
        
        return y_new
            
    def findMean(self, x, f, splineInterpolate=False):
        
        if splineInterpolate:
            from scipy import interpolate
            N = len(x)*10
            
            tck = interpolate.splrep(x, f, s=0)
            x = np.linspace(x[0], x[-1], N)
            f = interpolate.splev(x, tck, der=0)
        F = simps(f, x)
        
        f_mean = F/(x[-1] - x[0])
        
        return f_mean
    

    def calcEpsilonAvg(self, refData, numData, data_type="P"):
        
        if data_type == "P":
            RMS = np.sum(np.abs((numData - refData)/refData))/len(refData)
        elif data_type == "Q":
            RMS = np.sum(np.abs((numData - refData)/np.amax(refData)))/len(refData)
            
        return RMS
        

class ConvergenceLogFile:
    
    
    def __init__(self, networkLogFile, networkName, epsilonList, batchDataList, solutionFileDirectory):
        

        TemplateLogFileDirectory = mFPH.getDirectory('networkXmlFileTemplateDirectory', 'singleVessel_template', '200', 'read')
        #print TemplateLogFileDirectory
        self.templateFile = ''.join([TemplateLogFileDirectory,'/convergence_template.tex'])
        
        self.networkLogFile = networkLogFile
        
        self.networkName = networkName
        self.epsilonList = epsilonList
        self.batchDataList = batchDataList
        self.solutionFileDirectory = solutionFileDirectory
    
    def writeConvergenceLogfile(self, compileLogFile=False, deleteAuxiliary=False):
        
        fTemplate = open(self.templateFile, 'r')
        fLogFile = open(self.networkLogFile, 'w')
        for line in fTemplate:
            
            if "% Fill in title here" in line:
                self.writeTitle(fLogFile)
            elif "% Fill convergence data here" in line:
                self.writeConvergenceData(fLogFile)
            else:
                fLogFile.write(line)

        
        fTemplate.close()
        fLogFile.close()
        
        if compileLogFile:
            self.compilePdf(deleteAuxiliary=deleteAuxiliary)

    def writeTitle(self, fLogFile):
        #titleLine = ''.join([r'\section{Network : ',self.networkName, ', Datanumber: ', self.dataNumber, '}'])
        titleLine =''.join([r'\Large Network : ',self.networkName, r'\\'])
        fLogFile.write(titleLine)
        fLogFile.write('\n')

    def writeConvergenceData(self, fLogFile):
        
        N = len(self.batchDataList)
        fLogFile.write(r'\begin{table}[h!]')
        fLogFile.write('\n')
        fLogFile.write(r'\begin{tabular}{l')
        for n in range(N):
            fLogFile.write('l')
        
        fLogFile.write(r'}')
        fLogFile.write('\n')
        fLogFile.write(r'\hline')
        fLogFile.write('\n')
        fLogFile.write('dataN')
        fLogFile.write("    &    ")
        for n in range(N):
            fLogFile.write(self.batchDataList[n]['dataNumber'])
            if n != N -1:
                fLogFile.write("    &    ")
        
        fLogFile.write(r'    \\')
        fLogFile.write('\n')
        fLogFile.write(r'\hline')
        fLogFile.write('\n')
        fLogFile.write('dt [ms]')
        fLogFile.write("    &    ")
        for n in range(N):
            fLogFile.write(str(round(self.batchDataList[n]['dt']*1000, 3)))
            if n != N -1:
                fLogFile.write("    &    ")
        
        fLogFile.write(r'    \\')

        fLogFile.write('\n')
        fLogFile.write(r'\hline')
        fLogFile.write('\n')

        fLogFile.write('CFL')
        fLogFile.write("    &    ")
        for n in range(N):
            fLogFile.write(str(self.batchDataList[n]['CFL']))
            if n != N -1:
                fLogFile.write("    &    ")
        
        fLogFile.write(r'    \\')

        fLogFile.write('\n')
        fLogFile.write(r'\hline')
        fLogFile.write('\n')

        fLogFile.write(r'$\epsilon_P [\%]$')
        fLogFile.write("    &    ")
        for n in range(N - 1):
            fLogFile.write(str(self.epsilonList[n]))

            fLogFile.write("    &    ")
        fLogFile.write(r'\\')
        fLogFile.write('\n')
        fLogFile.write(r'\hline')
        fLogFile.write('\n')

        fLogFile.write(r'\end{tabular}')
        fLogFile.write('\n')
        fLogFile.write(r'\end{table}')
        
        imagecount = 0
        fLogFile.write('\n')
        self.writeStartFigure(fLogFile)
        for n in range(N - 1):
            if imagecount == 2:
                self.writeEndFigure(fLogFile)
                self.writeStartFigure(fLogFile)
                imagecount = 0
            self.writeFigure(fLogFile, self.batchDataList[n]['figFile'])
            imagecount += 1
        self.writeEndFigure(fLogFile)
        
        
    
    def writeStartFigure(self, fLogFile):
        fLogFile.write(r'\begin{figure}[h]')
        fLogFile.write('\n')
        fLogFile.write(r'\centering')
        fLogFile.write('\n')
    
    def writeEndFigure(self, fLogFile):
        fLogFile.write(r'\end{figure}')
        fLogFile.write('\n')
    
    def writeFigure(self, fLogFile, figurePath):
        fLogFile.write(r'\includegraphics[height=0.25\textheight,width=0.48\textwidth]{')
        fLogFile.write(figurePath)
        fLogFile.write(r'}')
        fLogFile.write('\n')
        
        
        
        
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

