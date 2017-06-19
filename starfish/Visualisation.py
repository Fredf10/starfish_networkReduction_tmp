#!/usr/bin/env python
from __future__ import print_function, absolute_import
from builtins import input
import subprocess

from optparse import OptionParser
import os,sys
cur = os.path.dirname( os.path.realpath( __file__ ) )

#sys.path.append(cur+'/UtilityLib')
import starfish
import starfish.UtilityLib.moduleStartUp as mStartUp

def main():    
    print("")
    print('=====================================')
    print('#    STARFiSh_v0.4 Visualisation    #')
    print('=====================================')
    print("")
    
    optionsDict = mStartUp.parseOptions(['f','n','v'], visualisationOnly = True)
    
    networkName    = optionsDict['networkName']
    dataNumber     = optionsDict['dataNumber']
    vizOutput      = optionsDict['vizOutput']
    if vizOutput == 'non':
        vizOutput = "2D+3D"
        
    string2d = ' '.join([sys.executable, '-c' 
        '"import starfish.VisualisationLib.class2dVisualisation as viz;', 
        "viz.main()", '"', 
        '-f', networkName, 
        '-n',str(dataNumber)])

    string3d = ' '.join([sys.executable, '-c' 
        '"import starfish.VisualisationLib.class3dVisualisation as viz;', 
        "viz.main()", '"', 
        '-f', networkName, 
        '-n',str(dataNumber)])
    
    if vizOutput == "2D":
        subprocess.Popen(string2d, shell=True)

        
    if vizOutput == "3D":
        viz3d = subprocess.Popen(string3d, shell = True )
        
    if vizOutput == "2D+3D":
        viz2d = subprocess.Popen(string2d, shell = True )
        viz3d = subprocess.Popen(string3d, shell = True )
        while True:
            if viz2d.poll() is not None:
                viz2d.terminate()
                exit()
            if viz3d.poll() is not None:
                viz3d.terminate()
                exit()
                
if __name__ == '__main__':
    main()
