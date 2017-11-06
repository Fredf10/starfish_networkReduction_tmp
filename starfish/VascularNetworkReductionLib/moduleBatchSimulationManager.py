from __future__ import print_function, absolute_import
from future.utils import iteritems, iterkeys, viewkeys, viewitems, itervalues, viewvalues
from builtins import input as input3
import sys,os
from starfish.SolverLib.class1DflowSolver import FlowSolver
from starfish.UtilityLib import moduleXML
import progressbarsimple as cPB
import gc,time
import multiprocessing

## module to run simulations as a batch jobs

def runBatchAsSingleProcess(batchDataList, quiet = False):
    '''
    Run a set of simulations on one core without multiprocessing
    
    Args:
        batchDataList (list) :
            with data for each batch job [batchData1, batchData .. ]
        batchData (dict) :
             dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
        
    '''
    timeStartBatch = time.time()
    print('=====================================')
    print('------Single Process Batch Job-------')
    print('numberOfEval.:   {}'.format(len(batchDataList)))
    #progressBar = cPB.ProgressBar(35, len(batchDataList))
    for completed,batchData in enumerate(batchDataList):
        minutesSolve,secsSolve = runSingleBatchSimulation(batchData)
        if quiet == False:
            print('____________Batch   {:5} ___________'.format(batchDataList.index(batchData))) 
            print('Runtime:        {} min {} sec'.format(minutesSolve,secsSolve))
        #progressBar.progress()
    timeBatchJob= time.time()-timeStartBatch
    minutesBatch = int(timeBatchJob/60.)
    secsBatch = timeBatchJob-minutesBatch*60.
    print('=====================================')
    print('total runtime:  {} min {} sec'.format(minutesBatch,secsBatch))
    print('=====================================')
  
def runSingleBatchSimulation(batchData):
    '''
    Function which runs a single simulation the defined with batchData
    
    Args:
        batchData (dict):
             dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
    '''
    
    networkName              = batchData['networkName']
    dataNumber               = batchData['dataNumber']
    networkXmlFileLoad       = batchData['networkXmlFileLoad']
    networkXmlFileSave       = batchData['networkXmlFileSave']
    pathSolutionDataFilename = batchData['pathSolutionDataFilename']
    try:
        timeStart = time.clock()
        #simulationIndex          = batchData['simulationIndex']
        
        vascularNetworkTemp = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = networkXmlFileLoad, pathSolutionDataFilename = pathSolutionDataFilename)
        vascularNetworkTemp.quiet = True
        flowSolver = FlowSolver(vascularNetworkTemp, quiet=True)
        flowSolver.solve()
        vascularNetworkTemp.saveSolutionData()
        moduleXML.writeNetworkToXML(vascularNetworkTemp, dataNumber, networkXmlFileSave)
        del flowSolver
        gc.collect()
        timeSolverSolve = time.clock()-timeStart
        minutesSolve = int(timeSolverSolve/60.)
        secsSolve = timeSolverSolve-minutesSolve*60.
    except Exception as e: #TODO handle exceptions with more information
        minutesSolve= 0
        secsSolve = 0
        print("Error in running {}".format(networkXmlFileLoad))
        print("The exception was: ")
        print(e) #TODO add a traceback
    
    return minutesSolve,secsSolve

def runBatchAsMultiprocessing(batchDataList, numberWorkers = None, quiet = False):
    '''
    Run a set of simulations on one core without multiprocessing
    
    Args:
        batchDataList (list) :
            with data for each batch job [batchData1, batchData .. ]
        batchData (dict) :
             dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
    '''
    if numberWorkers == None: numberWorkers = multiprocessing.cpu_count()
    
    timeStartBatch = time.time()
    print('=====================================')
    print('------Multiprocessing Batch Job------')
    print('numberWorkers:   {}'.format(numberWorkers))
    print('numberOfEval.:   {}'.format(len(batchDataList)))
    #progressBar = cPB.ProgressBar(35, len(batchDataList))
    pool = multiprocessing.Pool(numberWorkers, maxtasksperchild = None)
    results = pool.imap(runSingleBatchSimulation,batchDataList)
    pool.close() 
    last_update  = 0
    while (True):
        completed = results._index
        while last_update < completed:
            #progressBar.progress()
            last_update +=1
        if (completed == len(batchDataList)): break
    pool.join()
    if quiet == False:
        print('=====================================')
        for batchJobIndex,[minutesSolve,secsSolve] in enumerate(results):
            print('____________Batch   {:5} ___________'.format(batchJobIndex+1)) 
            print('Runtime:        {} min {} sec'.format(minutesSolve,secsSolve))
        print('=====================================')
    timeBatchJob= time.time()-timeStartBatch
    minutesBatch = int(timeBatchJob/60.)
    secsBatch = timeBatchJob-minutesBatch*60.
    print("total runtime {:d} min {:.2f} sec".format(minutesBatch, secsBatch)) 
    print('=====================================')
