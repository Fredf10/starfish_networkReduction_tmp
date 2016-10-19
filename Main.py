########################################################################################
#                           STARFiSh v0.4 
########################################################################################
## 
# http://www.ntnu.no/starfish
#
# created by Vinzenz Eck TODO:vinzenz.eck@mytum.de
# Developed by:
# TODO: Hallvard M Nydal, Yappi, ... Fredrik Fossan, Yvan Gugler, Jacob Sturdy, Einar ..., 
# TODO: ADD LICENSE (MIT) and COPYRIGHT
#
#Copyright (c) <year> <copyright holders>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this 
#software and associated documentation files (the "Software"), to deal in the Software 
#without restriction, including without limitation the rights to use, copy, modify, 
#merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
#permit persons to whom the Software is furnished to do so, subject to the following 
# conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or 
#substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
#CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##

#---------------------------------------------------------------------------------------#
from __future__ import division
import time 
import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath('__file__') )
import logging
logger = logging.getLogger('starfish')
logger.setLevel(logging.DEBUG)


import SolverLib.class1DflowSolver as c1DFlowSolv

import UtilityLib.moduleXML as mXML
import UtilityLib.moduleStartUp as mStartUp #import parseOptions
import UtilityLib.moduleFilePathHandler as mFPH

import matplotlib.pyplot as plt

import gc

import subprocess

def main():
    file_h = logging.FileHandler('starfish.log', mode='w')
    file_h.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_h.setFormatter(file_formatter)
    logger.addHandler(file_h)

    console_h = logging.StreamHandler(sys.stdout)
    console_h.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(message)s')
    console_h.setFormatter(console_formatter)
    logger.addHandler(console_h)

    logger.info("")
    logger.info('=====================================')
    logger.info('#     STARFiSh_v0.4.2016.10.19      #')
    logger.info('=====================================')
    
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
    
    
    logger.info('____________ Solver time _____________')
    logger.info('Initialisation: {} min {} sec'.format(minutesInit,secsInit))
    logger.info('Solving:        {} min {} sec'.format(minutesSolve,secsSolve))
    logger.info('=====================================')
    
    vascularNetwork.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetwork, dataNumber = dataNumber) # needs to be moved to vascularNetwork
    
    
    del flowSolver
    gc.collect()
    
    mFPH.updateSimulationDescriptions(networkName, dataNumber, simulationDescription)
    
    gc.collect()
    
    
    if vizOutput == "2D":
        string = ' '.join(['python',cur+'/VisualisationLib/class2dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber)])                
        subprocess.Popen(string, shell=True)
        
    
    if vizOutput == "3D":
        string = ' '.join(['python',cur+'/VisualisationLib/class3dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        subprocess.Popen(string, shell=True)
        
        
    if vizOutput == "2D+3D":
           
        string1 = ' '.join(['python',cur+'/VisualisationLib/class2dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        string2 = ' '.join(['python',cur+'/VisualisationLib/class3dVisualisation.py','-f',vascularNetwork.name, '-n',str(dataNumber), '-c True']) 
        
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
    main()
