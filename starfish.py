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

import subprocess,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath('__file__') )

import Main as simulator
import Visualisation as visualisationToolBox
import vnc 
import VascularPolynomialChaos as uqsaToolBox

def main():
        
    print("")
    print('=====================================')
    print('#     STARFiSh_v0.4.19.10.2016      #')
    print('=====================================')
    
    mainMenuInput = ""
    while mainMenuInput not in ['q']:
        print('\n Main menu:\n')
        print " [1] - run simulation | Main.py"
        print " [2] - run vnc (vascular network creator) | vnc.py"
        print " [3] - run visualisation | Visualisation.py"
        print " [4] - run uncertainty quantification tool box | VascularPolynomialChaos.py"
        print " [q] - quit \n"
        while  mainMenuInput not in ('1','2','3','4','q'):
            mainMenuInput = raw_input("what to do? ")
        
        if mainMenuInput == '1':
            print "\n .. running simulation \n"
            simulator.main()
            mainMenuInput = ""
            
        if mainMenuInput == '2':
            print "\n .. running vnc \n"
            vnc.main()
            mainMenuInput = ""
            
        if mainMenuInput == '3':
            print "\n .. running visualisation \n"
            visualisationToolBox.main()
            mainMenuInput = ""
            
        if mainMenuInput == '4':
            print "\n .. running uncertainty quantification tool box \n"
            uqsaToolBox.uncertaintyPropagation()
            mainMenuInput = ""
            
    print "bye bye .." 
    
if __name__ == '__main__':
    main()