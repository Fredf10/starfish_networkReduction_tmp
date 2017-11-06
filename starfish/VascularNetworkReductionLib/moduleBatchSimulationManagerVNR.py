'''
Created on Jul 5, 2017

@author: fredrik
'''


import os

from starfish.SolverLib.class1DflowSolver import FlowSolver

from starfish.UtilityLib import moduleXML
import starfish.UtilityLib.progressBar as cPB

import starfish.VascularNetworkReductionLib.classPostprocessReduction as cVNRpostProcess

import gc,time

import multiprocessing

## module to run simulations as a batch jobs

def runBatchAsSingleProcess(batchDataList, quiet = False):
    '''
    Run a set of simulations on one core without multiprocessing
    
    Input:
        batchDataList <list> := with data for each batch job [batchData1, batchData .. ]
            batchData <dict> := dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
        
    '''
    timeStartBatch = time.time()
    print '====================================='
    print '------Single Process Batch Job-------'
    print 'numberOfEval.:   {}'.format(len(batchDataList))
    progressBar = cPB.ProgressBar(35, len(batchDataList))
    for completed,batchData in enumerate(batchDataList):
        progressBar.progress(completed)
        minutesSolve,secsSolve = runSingleBatchSimulation(batchData)
        if quiet == False:
            print '____________Batch   {:5} ___________'.format(batchDataList.index(batchData)) 
            print 'Runtime:        {} min {} sec'.format(minutesSolve,secsSolve)
    timeBatchJob= time.time()-timeStartBatch
    minutesBatch = int(timeBatchJob/60.)
    secsBatch = timeBatchJob-minutesBatch*60.
    print '====================================='
    print 'total runtime:  {} min {} sec'.format(minutesBatch,secsBatch)
    print '====================================='
  
def runSingleBatchSimulation(batchData):
    '''
    Function which runs a single simulation the defined with batchData
    
    Input:
        batchData <dict> := dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
    '''
    timeStart = time.clock()
    #simulationIndex          = batchData['simulationIndex']
    networkName              = batchData['networkName']
    dataNumber               = batchData['dataNumber']
    networkXmlFileLoad       = batchData['networkXmlFileLoad']
    networkXmlFileSave       = batchData['networkXmlFileSave']
    pathSolutionDataFilename = batchData['pathSolutionDataFilename']
    vascularNetworkTemp = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = networkXmlFileLoad, pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkTemp.quiet = True
    flowSolver = FlowSolver(vascularNetworkTemp, quiet=True)

    flowSolver.solve()

    vascularNetworkTemp.saveSolutionData()
    moduleXML.writeNetworkToXML(vascularNetworkTemp, dataNumber, networkXmlFileSave)
    del flowSolver
    del vascularNetworkTemp
    t_tmp = time.time()
    print "###### starting to postprocess"
    cVnRpost = cVNRpostProcess.PostprocessReduction()
    cVnRpost.postprocessSingle(batchData)
    print "postprocessing time: ", time.time() - t_tmp
    gc.collect()
    timeSolverSolve = time.clock()-timeStart
    minutesSolve = int(timeSolverSolve/60.)
    secsSolve = timeSolverSolve-minutesSolve*60.
    return minutesSolve,secsSolve

def runBatchAsMultiprocessing(batchDataList, numberWorkers=None, quiet = False, postProcess=True, CPUTimeFile=None):
    '''
    Run a set of simulations on one core without multiprocessing
    
    Input:
        batchDataList <list> := with data for each batch job [batchData1, batchData .. ]
            batchData <dict> := dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
        
    '''
    if numberWorkers == None: numberWorkers = multiprocessing.cpu_count()
    print "\n"
    print numberWorkers
    timeStartBatch = time.time()
    print '====================================='
    print '------Multiprocessing Batch Job------'
    print 'numberWorkers:   {}'.format(numberWorkers)
    print 'numberOfEval.:   {}'.format(len(batchDataList))
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


    
    if CPUTimeFile is not None:
        
        writeCPUTimeFile(CPUTimeFile, len(batchDataList), minutesBatch, secsBatch)

def writeCPUTimeFile(CPUTimeFile, Nsolve, minutesBatch, secsBatch):
    
    f = open(CPUTimeFile, 'w')
    f.write("Solved {0} cases in {1} minutes and {2} sec".format(Nsolve, minutesBatch, secsBatch))
    f.close()