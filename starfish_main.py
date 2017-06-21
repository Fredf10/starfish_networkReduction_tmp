#!/usr/bin/env python
# -*- coding: utf-8 -*- 
########################################################################################
#                           STARFiSh v0.4 
########################################################################################
## 
# http://www.ntnu.no/starfish
#
# Contributors:
# Leif Rune Hellevik, Vinzenz Gregor Eck, Jacob Sturdy, Fredrik Eikeland Fossan, 
# Einar Nyberg Karlsen, Yvan Gugler, Yapi Donatien Achou, Hallvard Moian Nydal, 
# Knut Petter Maråk #TODO: THIS MAY CAUSE UNICODE ISSUES, Paul Roger Leinan 
#
# TODO: ADD LICENSE (MIT) and COPYRIGHT
#
#Copyright (c) <2012-> <NTNU>
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
from __future__ import print_function, absolute_import
from future.utils import iteritems, iterkeys, viewkeys, viewitems, itervalues, viewvalues
from builtins import input as input3
import starfish.Simulator as simulator
import starfish.Visualisation as visualisationToolBox
import starfish.vnc  as vnc
import starfish.VascularPolynomialChaos as uqsaToolBox

def main():
        
    print("")
    print('=====================================')
    print('#     STARFiSh_v0.4.19.10.2016      #')
    print('=====================================')
    
    mainMenuInput = ""
    while mainMenuInput not in ['q']:
        print('\n Main menu:\n')
        print(" [1] - run simulation | Main.py")
        print(" [2] - run vnc (vascular network creator) | vnc.py")
        print(" [3] - run visualisation | Visualisation.py")
        print(" [4] - run uncertainty quantification tool box | VascularPolynomialChaos.py")
        print(" [q] - quit \n")
        while  mainMenuInput not in ('1','2','3','4','q'):
            mainMenuInput =input3("what to do? ")
        
        if mainMenuInput == '1':
            print("\n .. running simulation \n")
            simulator.main()
            mainMenuInput = ""
            
        if mainMenuInput == '2':
            print("\n .. running vnc \n")
            vnc.main()
            mainMenuInput = ""
            
        if mainMenuInput == '3':
            print("\n .. running visualisation \n")
            visualisationToolBox.main()
            mainMenuInput = ""
            
        if mainMenuInput == '4':
            print("\n .. running uncertainty quantification tool box \n")
            uqsaToolBox.uncertaintyPropagation()
            mainMenuInput = ""
            
    print("bye bye ..") 
    
if __name__ == '__main__':
    main()
