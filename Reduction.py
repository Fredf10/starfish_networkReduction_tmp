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
import shutil
import subprocess
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath('__file__') )

#import NetworkLib.classVascularNetwork as cVascNw
### Unnecessary import, never used later.


import NetworkLib.classNetworkReduction as cNred

import UtilityLib.moduleXML as mXML
import UtilityLib.moduleLogFile as mLOG
import UtilityLib.moduleStartUp as mStartUp #import parseOptions
import UtilityLib.moduleFilePathHandler as mFPH



def main():
    print ""
    print '====================================='
    print '#     STARFiSh_v0.3_development     #'
    print '====================================='
    
    optionsDict = mStartUp.parseOptions(['f', 'e', 'n','d','s','v','r','w','p'])
    
    networkName           = optionsDict['networkName']
    newNetworkName        = optionsDict['NewNetworkName']
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
    
    oldNetworkDirectory = mFPH.getDirectory('networkXmlFileXXXDirectory', networkName, "xxx", 'read')
    
    vascularNetwork.update({'description':simulationDescription,
                            'dataNumber' :dataNumber})
    
    New_network = cNred.NetworkReduction(vascularNetwork)
    truncateFile =  ''.join([oldNetworkDirectory,'/','truncate.txt'])
    New_network.reduceNetwork(truncateFile)

        
    New_network.name = newNetworkName
    newNetworkXmlFile =  mFPH.getFilePath('networkXmlFile', newNetworkName, "xxx", 'write')
    
    
    mXML.writeNetworkToXML(New_network, dataNumber = dataNumber, networkXmlFile=newNetworkXmlFile)

    if New_network.initialsationMethod == 'FromSolution':
        oldInitialValuePath = mFPH.getDirectory('initialValueFileDirectory', networkName, dataNumber, 'write')
        newInitialValuePath = mFPH.getDirectory('initialValueFileDirectory', newNetworkName, dataNumber, 'write')
        if os.path.isdir(newInitialValuePath):
            shutil.rmtree(newInitialValuePath)

        shutil.copytree(oldInitialValuePath, newInitialValuePath)
    
    copyFlowFile = True
    if copyFlowFile:
            oldInflowFile = ''.join([oldNetworkDirectory,'/','inflow.csv'])
            newNetworkDirectory = mFPH.getDirectory('networkXmlFileXXXDirectory', newNetworkName, "xxx", 'read')
            newInflowFile = ''.join([newNetworkDirectory,'/','inflow.csv'])
            
            shutil.copyfile(oldInflowFile, newInflowFile)

    copyTruncateFile = True
    if copyTruncateFile:
            
            newNetworkDirectory = mFPH.getDirectory('networkXmlFileXXXDirectory', newNetworkName, "xxx", 'read')
            newTruncateFile = ''.join([newNetworkDirectory,'/','truncate.txt'])
            
            shutil.copyfile(truncateFile, newTruncateFile)

    
    
#     if save:
#         string1 = ' '.join(['/usr/bin/python',cur+'/Main.py','-f',newNetworkName, '-n',dataNumber, '-d', simulationDescription])
#         print "gonna run new network"
#         print string1
#         subprocess.Popen(string1, shell = True )
        
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
