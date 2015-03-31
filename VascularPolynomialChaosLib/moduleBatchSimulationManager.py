

import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))

sys.path.append(''.join([cur,'/../NetworkLib']))
from classVascularNetwork import VascularNetwork

sys.path.append(''.join([cur,'/../Solver']))
from class1DflowSolver import FlowSolver

sys.path.append(''.join([cur,'/../UtilityLib']))
import moduleXML

import gc,time

## module to run simulations as a batch jobs

def runBatchSingleProcess(batchDataList, quiet = False):
    '''
    Run a set of simulations on one core without multiprocessing
    
    Input:
        batchDataList <list> := with data for each batch job [batchData1, batchData .. ]
            batchData <list> := list with [networkName, dataNumber, networkXmlFile, pathSolutionDataFilename]
        
    '''
    for batchData in batchDataList:
        timeStart = time.clock()
        networkName,dataNumber,networkXmlFile,pathSolutionDataFilename = batchData
        #print networkXmlFile
        #print pathSolutionDataFilename
        vascularNetworkTemp = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = networkXmlFile, pathSolutionDataFilename = pathSolutionDataFilename)
        vascularNetworkTemp.quiet = True
        flowSolver = FlowSolver(vascularNetworkTemp, quiet=True)
        flowSolver.solve()
        vascularNetworkTemp.saveSolutionData()
        del flowSolver
        gc.collect()
        
        if quiet == False:
            timeSolverSolve = time.clock()-timeStart
            minutesSolve = int(timeSolverSolve/60.)
            secsSolve = timeSolverSolve-minutesSolve*60.
            print '====================================='
            print '__________Batch {:5} time ___________'.format(batchDataList.index(batchData))
            print 'Runtime:        {} min {} sec'.format(minutesSolve,secsSolve)
            #print '====================================='