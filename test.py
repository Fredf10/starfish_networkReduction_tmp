import os
import sys
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import starfish
import starfish.NetworkLib as NetworkLib
import starfish.SolverLib as SolverLib
import starfish.UtilityLib as UtilityLib
import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath('__file__') )
import logging
logger = logging.getLogger('starfish')
logger.setLevel(logging.DEBUG)

import starfish.SolverLib.class1DflowSolver as c1DFlowSolv 
import starfish.UtilityLib.moduleXML  as mXML
import starfish.UtilityLib.moduleStartUp as mStartUp #import parseOptions
import starfish.UtilityLib.moduleFilePathHandler as mFPH 

import matplotlib.pyplot as plt

import gc

import subprocess

def main():
    optionsDict = mStartUp.parseOptions(['f','n','d','s','v','r','w','p'])
    
    networkName           = optionsDict['networkName']
    save                  = optionsDict['save']
    dataNumber            = optionsDict['dataNumber']
    simulationDescription = optionsDict['simulationDescription']
    vizOutput             = optionsDict['vizOutput']
    resimulate            = optionsDict['resimulate']
    
    filename = str(networkName+'.xml')
        
    logger.info('____________Simulation_______________')
    logger.info('%-20s %s' % ('Network name',networkName))
    logger.info('%-20s %s' % ('Data number', dataNumber))
    logger.info('%-20s %s' % ('Save simulation', save))
    logger.info('%-20s %s' % ('Case description', simulationDescription))
    logger.info('%-20s %s' % ('Resimulate', resimulate))
    logger.info('%-20s %s' % ('Visualisationmode', vizOutput))
    
    ## check if template
    if '_template' in networkName:
        networkName = mFPH.createWorkingCopyOfTemplateNetwork(networkName)
    
    # load network from the path!
    if resimulate == False:
        vascularNetwork = mXML.loadNetworkFromXML(networkName) # moved to vascularNetowrk constror
    else:
        # resimulate network
        vascularNetwork = mXML.loadNetworkFromXML(networkName, dataNumber = dataNumber)        
        if simulationDescription == '':
            simulationDescription = vascularNetwork.description
    
    if vascularNetwork == None: exit()
    
    
    vascularNetwork.update({'description':simulationDescription,
                            'dataNumber' :dataNumber})
    
    #initialize Solver
    flowSolver = c1DFlowSolv.FlowSolver(vascularNetwork)
    #solve the system
    flowSolver.solve()
    
    
    vascularNetwork.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetwork, dataNumber = dataNumber) # needs to be moved to vascularNetwork
    
    
    del flowSolver
    gc.collect()
    
    mFPH.updateSimulationDescriptions(networkName, dataNumber, simulationDescription)
    
    gc.collect()
            
if __name__ == '__main__':
    main()

