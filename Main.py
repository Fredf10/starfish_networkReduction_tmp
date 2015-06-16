########################################################################################
#                            Vascular1dFlow v0.2
########################################################################################
## 
# created by Vinzenz Eck vinzenz.eck@mytum.de
# uses polynomial Chaos toolbox from Jonathan Feinberg, Simula Center Oslo
##

#---------------------------------------------------------------------------------------#
import time 
import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath('__file__') )

#import NetworkLib.classVascularNetwork as cVascNw
### Unnecessary import, never used later.

import SolverLib.class1DflowSolver as c1DFlowSolv

import UtilityLib.moduleXML as mXML
import UtilityLib.moduleStartUp as mStartUp #import parseOptions
import UtilityLib.moduleFilePathHandler as mFPH

import pprint
import matplotlib.pyplot as plt

import gc

import subprocess


def main():
    print ""
    print '====================================='
    print '#     STARFiSh_v0.3_development     #'
    print '====================================='
    
    optionsDict = mStartUp.parseOptions(['f','n','d','s','v','r','w'])
    
    networkName           = optionsDict['networkName']
    save                  = optionsDict['save']
    dataNumber            = optionsDict['dataNumber']
    simulationDescription = optionsDict['simulationDescription']
    vizOutput             = optionsDict['vizOutput']
    resimulate            = optionsDict['resimulate']
    
    filename = str(networkName+'.xml')
        
    print '____________Simulation_______________'
    print '%-20s %s' % ('Network name',networkName)
    print '%-20s %s' % ('Data number', dataNumber)
    print '%-20s %s' % ('Save simulation', save)
    print '%-20s %s' % ('Case description', simulationDescription)
    print '%-20s %s' % ('Resimulate', resimulate)
    print '%-20s %s' % ('Visualisationmode', vizOutput)
    
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
    
    timeSolverInitStart = time.clock()
    #initialize Solver
    flowSolver = c1DFlowSolv.FlowSolver(vascularNetwork)
    timeSolverInit = time.clock()-timeSolverInitStart
    timeSolverSolveStart = time.clock()
    #solve the system
    flowSolver.solve()
    timeSolverSolve = time.clock()-timeSolverSolveStart
    
    minutesInit = int(timeSolverInit/60.)
    secsInit = timeSolverInit-minutesInit*60.
    minutesSolve = int(timeSolverSolve/60.)
    secsSolve = timeSolverSolve-minutesSolve*60.
    
    
    #print '====================================='
    print '____________ Solver time _____________'
    print 'Initialisation: {} min {} sec'.format(minutesInit,secsInit)
    print 'Solving:        {} min {} sec'.format(minutesSolve,secsSolve)
    print '====================================='
    
    vascularNetwork.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetwork, dataNumber = dataNumber) # needs to be moved to vascularNetwork
    
    
    del flowSolver
    gc.collect()
    
    mFPH.updateSimulationDescriptions(networkName, dataNumber, simulationDescription)
    
    gc.collect()
    
    
    if vizOutput == "2D":
        string = ' '.join(['python',cur+'/Visualisation/class2dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber)])                
        subprocess.Popen(string, shell=True)
        
    
    if vizOutput == "3D":
        string = ' '.join(['python',cur+'/Visualisation/class3dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        subprocess.Popen(string, shell=True)
        
        
    if vizOutput == "2D+3D":
           
        string1 = ' '.join(['python',cur+'/Visualisation/class2dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        string2 = ' '.join(['python',cur+'/Visualisation/class3dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        
        viz2d = subprocess.Popen(string1, shell = True )
        viz3d = subprocess.Popen(string2, shell = True )
        
        while True:
            
            if viz2d.poll() != None:
                viz3d.terminate()
                exit()
                
            if viz3d.poll() != None:
                viz2d.terminate()
                exit()
        
if __name__ == '__main__':
    
    profiling = False
    callGraph = False
    
    if profiling == True:
        import cProfile
        cProfile.run('main()')
        
    elif callGraph == True:
        import pycallgraph
        pycallgraph.start_trace()
        main()
        pycallgraph.make_dot_graph('STARFiSh-CallGraph.png')
    else:
        main()
