from __future__ import print_function, absolute_import

from scipy import interpolate
from scipy.integrate import simps
import sys, os, h5py
import numpy as np
import math
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import starfish.UtilityLib.moduleXML as mXML
import starfish.SolverLib.class1DflowSolver as c1dFS
import starfish.VascularNetworkReductionLib.classNetworkReduction as cNred

def calcRMS(refData, numData, data_type="P"):
    
    if data_type == "P":
        RMS = np.sum(np.abs((numData - refData)/refData))/len(refData)
    elif data_type == "Q":
        RMS = np.sum(np.abs((numData - refData)/np.amax(refData)))/len(refData)
        
    return RMS*100

def findMean(x, f):
    
    F = simps(f, x)
    
    f_mean = F/(x[-1] - x[0])
    
    return f_mean

def loadData(networkFile, vesselId, period=0.8, node=0, medical=True, allTime=False):
    
    hdf5File = h5py.File(networkFile,'r')
    time = hdf5File['VascularNetwork']['simulationTime'][:] - hdf5File['VascularNetwork']['simulationTime'][0]
    
    #dt = time[1] - time[0]

    dt = 0.005 #time[1] - time[0]
    
    vesselDict = {}

    N = int(round(period/dt))
    
    nCycles = int(round(time[-1]/period))
    
    
    t_start = period*(nCycles - 1)
    t_end = period*nCycles
    

    
    time_compare = np.linspace(t_start, t_end, N + 1)
    
    Nvessel = 0
    P = None
    Q = None
    for vesselName in hdf5File['vessels'].keys():
        
        tmpvesselId = vesselName.split(' - ')[-1]
        tmpvesselId = int(tmpvesselId)
        
        Nvessel += 1
        if tmpvesselId == vesselId:
            #print "yolo", vesselId, vesselName
            P  = mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Psol'][:, node], time_compare)
            Q  = mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Qsol'][:, node], time_compare)

            vName = vesselName.split(' - ')[0]
            
    if medical and P is not None:
        
        P = P/133.32

        Q = Q*1e6
    
    vesselDict['N'] = Nvessel
    #vesselDict['name'] = Name
    vesselDict['time'] = time_compare
    vesselDict['P'] = P
    vesselDict['Q'] = Q

    return vesselDict


def mySplineInterpolater(x, y, x_new):
    tck = interpolate.splrep(x, y, s=0)

    y_new = interpolate.splev(x_new, tck, der=0)
    
    return y_new


def test_networkReduction():
 
    simulationDescription = None

    useAverageValues = False
    useVesselsImpedance = False
    useLumpedValues = False
        
    optimizeParams = False
    Wkoptimize = 'Wk3'
    params = 'R1LCR2'

    networkName = "AoBif"
    reductionNetworkName = "Full96Model"
    truncateList = [2, 3]


    dataNumber = "100"

    networkXmlFileLoad = cur + "/Full96Model/Full96Model.xml"
    newNetworkXmlFileLoad = cur + "/AoBif.xml"
    

    # Temporary files for saving data
    networkXmlFileSave = cur + "/Full96Model/Full96Model_sol.xml"
    pathSolutionDataFilename = cur + "/Full96Model/Full96Model_sol.hdf5"
    newPathSolutionDataFilenameRef = cur + "/AoBif_sol.hdf5"
    newPathSolutionDataFilename = cur + "/AoBif_tmpSol.hdf5"
    
    vascularNetwork = mXML.loadNetworkFromXML(reductionNetworkName, dataNumber,
                                              networkXmlFile=networkXmlFileLoad,
                                              pathSolutionDataFilename=pathSolutionDataFilename)
    vascularNetwork.quiet = True
    
    vascularNetwork.update({'description':simulationDescription,
                            'dataNumber' :dataNumber})
    
    
    New_network = cNred.NetworkReduction(vascularNetwork, quiet=True)
    
    New_network.initialize(useAverageValues=useAverageValues, 
                           useVesselsImpedance=useVesselsImpedance, 
                           useLumpedValues=useLumpedValues,
                           optimizeParams=optimizeParams,
                           optParamsDict=None,
                           Wkoptimize=Wkoptimize,
                           params=params)
    

    
    New_network.reduceNetworkFromListGen(truncateList)
    New_network.name = networkName
    mXML.writeNetworkToXML(New_network, dataNumber = dataNumber, networkXmlFile=newNetworkXmlFileLoad)
    
    vascularNetworkNew = mXML.loadNetworkFromXML(networkName, dataNumber,
                                          networkXmlFile=newNetworkXmlFileLoad,
                                          pathSolutionDataFilename=newPathSolutionDataFilename)
    flowSolver = c1dFS.FlowSolver(vascularNetworkNew, quiet=True)
    flowSolver.solve()
    vascularNetworkNew.saveSolutionData()
    
    vesselList = [1, 2, 3]
    threshold = 0.01
    TooHighError = False
    for vesselId in vesselList:
        vesselDictBaseline = loadData(pathSolutionDataFilename, vesselId, period=0.8, node=0, medical=True, allTime=False)    
        vesselDictReduced = loadData(newPathSolutionDataFilename, vesselId, period=0.8, node=0, medical=True, allTime=False)
        vesselDictReducedRef = loadData(newPathSolutionDataFilenameRef, vesselId, period=0.8, node=0, medical=True, allTime=False)
        
        P_baseline = vesselDictBaseline['P']
        P_Reduced = vesselDictReduced['P']
        P_Reduced_Ref = vesselDictReducedRef['P']
        
        Q_baseline = vesselDictBaseline['Q']
        Q_Reduced = vesselDictReduced['Q']
        Q_Reduced_Ref = vesselDictReducedRef['Q']
        
        
        RMS1 = calcRMS(P_baseline, P_Reduced_Ref, data_type="P")
        RMS2 = calcRMS(P_baseline, P_Reduced, data_type="P")
        
            
        if abs(RMS1 - RMS2)> threshold :
            TooHighError = True
            print("Error was found to be too high for Vessel {}, key {}, with value {} deviating from {} with more than threshold {}".format(vesselId, 'P', RMS1, RMS2, threshold))

    assert(not TooHighError)
    if not TooHighError:
        print("\nAll values below error threshold")
        print("Test Successful!")
       


def test_SteadyStateSolver():

    networkName = "Full96Model"
    truncateList = [2, 3]


    dataNumber = "100"
    networkXmlFileLoad = cur + "/Full96Model/Full96Model.xml"
    

    # Temporary files for saving data
    networkXmlFileSave = cur + "/Full96Model/Full96Model_sol.xml"
    pathSolutionDataFilename = cur + "/Full96Model/Full96Model_sol.hdf5"
    pathSolutionDataFilenameTmp = cur + "/Full96Model/Full96Model_Tmpsol.hdf5"
    
    vascularNetwork = mXML.loadNetworkFromXML(networkName, dataNumber,
                                              networkXmlFile=networkXmlFileLoad,
                                              pathSolutionDataFilename=pathSolutionDataFilenameTmp)
    vascularNetwork.quiet = True
    
    vascularNetwork.update({'description':'None',
                            'dataNumber' :dataNumber})
    

    flowSolver = c1dFS.FlowSolver(vascularNetwork, quiet=True)

    epsilonMaxP = 0
    epsilonMaxQ = 0
    N = 0
    threshold = 0.1
    TooHighError = False
    for vesselId in vascularNetwork.lumpedValues:
        P_SS_inlet = vascularNetwork.lumpedValues[vesselId]['Pressure'][0]

        Q_SS_inlet = vascularNetwork.lumpedValues[vesselId]['Flow'][0]
        
        vesselDict = loadData(pathSolutionDataFilename, vesselId, period=0.8, node=0, medical=False, allTime=False)
        t = vesselDict['time']
        P = vesselDict['P']
        Q = vesselDict['Q']

        P_avg_inlet = np.mean(P)
        Q_avg_inlet = np.mean(Q)

        epsilonP = abs((P_SS_inlet - P_avg_inlet)/P_avg_inlet)
        epsilonQ = abs((Q_SS_inlet - Q_avg_inlet)/Q_avg_inlet)

        
        if epsilonP > epsilonMaxP:
            epsilonMaxP = epsilonP
            vesselMax = vesselId
        if epsilonQ > epsilonMaxQ:
            epsilonMaxQ = epsilonQ
    
    if abs(epsilonMaxP*100 - 2.39814666085) > threshold :
        TooHighError = True
        print("Error was found to be too high for Vessel {}, key {}, with value {} threshold {}".format(vesselMax, 'P', epsilonMaxP, threshold))

    assert(not TooHighError)
    if not TooHighError:
        print("\nAll values below error threshold")
        print("Test Successful!")
        
if __name__ == "__main__":
    test_networkReduction()
    test_SteadyStateSolver()
