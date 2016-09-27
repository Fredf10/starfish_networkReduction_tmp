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
    
    optionsDict = mStartUp.parseOptions(['f','n','d','s','v','r','w','p'])
    
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
    
    New_network = cNred.NetworkReduction(vascularNetwork)
    
    mXML.writeNetworkToXML(New_network, dataNumber = dataNumber)
    #mLog2 = mLOG.NetworkLogFile(vascularNetwork, dataNumber = dataNumber)
    #mLog2.writeNetworkLogfile()

        
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
