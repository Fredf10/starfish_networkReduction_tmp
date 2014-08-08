
# import dependencies

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )

# syspaths and functions for vascular1DFlow_v0.2
sys.path.append(cur+'/UtilityLib')
sys.path.append(cur+'/NetworkLib')
sys.path.append(cur+'/UtilityLib')
sys.path.append(cur+'/VncLib')
from classVascularNetwork import VascularNetwork 
from classBoundaryConditions import *
from moduleXML import writeNetworkToXML 
from moduleXML import loadNetworkFromXML 

from moduleCSV import readVesselDataFromCSV 
from moduleCSV import writeVesselDataToCSV 
from moduleCSV import writeBCToCSV 
from moduleCSV import readBCFromCSV

import networkXml041 as nxml

### import units of all variables in the SI system
from constants import variableUnitsSI as variableUnits
from constants import unitsDictSI as unitsDict

from modulePickle import loadSolutionDataFile
from moduleStartUp import chooseSolutionDataCase

### import units of all variales in the Medical System
#from constants import variableUnitsMed as variableUnits
#from constants import unitsDictMed as unitsDict

import pydot
import xdot
import gtk

from vnc_classes import *

import cPickle
import pprint as pprint
import numpy as np
import thread
from copy import deepcopy

def main():
    # set graphs directory
    graphPath = str(cur+'/NetworkFiles/')
    
    # create a new window
    window = MyDotWindow()    
    window.connect('destroy',gtk.main_quit)
    window.show()
    
    #create the main graph instance of a graph class
    mainGraph = Graph()
    
    #START THE MAIN LOOP
    menuInput = ""
    subMenuInput = ''
    
    # create vascularNetwork instance
    vascularNetwork = VascularNetwork()
    filename = None
    k = None
    while menuInput != "q":
        menuInput = ""
        print ""
        print '====================================='
        print '#    VascularNetworkCreator_v2.1    #'
        print '====================================='
        print " [a] - add vessel to network"
        print " [d] - delete vessel in network"
        print " [n] - new network"
        print " [b] - set boundary conditions"
        print " [f] - set global fluid properties"
        print " [l] - load network"
        print " [s] - save network"
        print " [u] - update XML form CSV file(s)"
        print " [g] - print network graph"
        print " [p] - print network informations"
        print " [q] - quit"
        print ""
        print '  current network: ', filename
        while  menuInput not in ("l","b","q","a","s","g","f","d","u",'n','p'):
            menuInput = raw_input("what to do? ")
        
        if menuInput == "a": 
            print "Add new vessel"
            
            existing = False
            vesselId = raw_input(" enter the vessel id:  ")
            while True:
                try:
                    vesselId = int(vesselId)
                    if vesselId not in vascularNetwork.vessels:
                        break
                    else:
                        existing = True
                except ValueError:
                    print "TYPE-ERROR: vessel id must be type(int) not type(string)"
                    vesselId = raw_input(" enter non existing id: ")
                if existing == True:
                    print " the vessel id exists already enter a new one"
                    vesselId = raw_input(" enter non existing id: ")
                    existing = False
            
            
            if vascularNetwork.vessels != {}:
                                    
                existing = False
                mother = raw_input(" enter existing mother id:  ")
                while True:
                    try:
                        mother = int(mother)
                        if mother in vascularNetwork.vessels.keys() and vascularNetwork.vessels[mother].rightDaughter == None:
                            break
                        else:
                            existing = True
                    except ValueError:
                        print "TYPE-ERROR: mother id must be type(int) not type(string)"
                        mother = raw_input(" enter existing mother id:  ")
                    if existing == True:
                        if mother not in vascularNetwork.vessels: print " there exists no vessel with this id"
                        else: print "   only bifurcations possible!"
                        mother = raw_input(" enter existing mother id:  ")
                        existing = False
                    
                if vascularNetwork.vessels[mother].leftDaughter == None: vascularNetwork.vessels[mother].leftDaughter = vesselId
                else: vascularNetwork.vessels[mother].rightDaughter = vesselId
                                
            vascularNetwork.addVessel(vesselId)
            print " define vessel compliance!"
            
            
            inputType = '0'
            print "     available compliance types:"
            print ""
            # get all defined boundaryConditions from constants-dict save as bcTypes
            complianceTypes = nxml.vesselComplianceElements.keys()
            compTypes = ['default (Hayashi)']
            for compType in complianceTypes: 
                compTypes.append(compType)
            # show all compliance types in the compTypes
            index = 0
            for key in compTypes:
                print "       [",str(index).rjust(2),"]    ",key
                index = index+1
            # get user input and check if it was correct to define the bcType 
            existing = False
            inputType = raw_input ("      choose type ")
            while True:
                # check if right input
                try:
                    inputType = int(inputType)
                    if inputType in np.linspace(0,len(compTypes)-1,len(compTypes)):
                        break
                    else:
                        existing = True
                # if not int but string
                except ValueError:
                    print "      TYPE-ERROR: vessel id must be type(int) not type(string)"
                    inputType = (raw_input ("      choose type "))
                # if int but to low or high
                if existing == True:
                    print "       the type does not exist"
                    inputType = (raw_input ("      choose type "))
                    existing = False
            
            compType = compTypes[int(inputType)] 
            
            if compType != 'default (Hayashi)':
                vesselData = {'complianceType':compType}
                nxml.vesselComplianceElements[compType]
                                    
                print ""
                print "      set values for the Compliance: ", compType
                question = True
                for arg in nxml.vesselComplianceElements[compType]:
                    if arg != 'complianceType':
                        currValue = raw_input (str("            set value for "+str(arg)+' '))
                        test = True
                        try: float(currValue)
                        except:
                            print '            VALUE or TYPE ERROR, set to None'
                            test = False
                        if test == True: vesselData[arg] = (float(currValue))
                        else: vesselData[arg] = None
                vascularNetwork.updateNetwork({'vesselData':  { vesselId : vesselData}})
            
            mainGraph.update_graph(vascularNetwork, window)
                    
        if menuInput == "d":
            print "Delete a vessel and all its daugthers"
            if vascularNetwork.vessels.keys() != []:
                
                existing = False
                vesselId = raw_input(" enter existing vessel id: ")
                while True:
                    try:
                        vesselId = int(vesselId)
                        if vesselId in vascularNetwork.vessels:
                            break
                        else: 
                            existing = True
                    except ValueError:
                        print "TYPE-ERROR: vessel id must be type(int) not type(string)"
                        vesselId = raw_input(" enter existing vessel id: ")
                    if existing == True:
                        print " the vessel does not exist"
                        vesselId = raw_input(" enter existing vessel id: ")
                        existing = False
                
                #travers the tree starting with the vessel and collect all ids
                toDelete = findAllDaughters(vascularNetwork,vesselId)
                
                toDelete.append(vesselId)
                
                for vesselToDelete in toDelete:
                    vascularNetwork.deleteVessel(vesselToDelete)
                
                # empty the graph to redraw it
                if vascularNetwork.vessels.keys() == []:
                    mainGraph.update_graph(None, window)
                    vascularNetwork = VascularNetwork()
                else:
                    mainGraph.update_graph(vascularNetwork, window)
                
            else:
                print " there are no vessels to delete"
                
        elif menuInput == "n":
            print "new network"
            question = raw_input(" are u sure to delete all current data? [y] - yes: ")
            if question == 'y':
                # delete vascularNetwork
                del vascularNetwork
                # create vascularNetwork instance
                mainGraph.update_graph(None, window)
                vascularNetwork = VascularNetwork()
                    
        elif menuInput == "p":
            vascularNetwork.showVessels()
            vascularNetwork.showNetwork()
                        
        elif menuInput == "g":
            print mainGraph.getGraph()
            
        elif menuInput == "b":
            subMenuInput = ''
            
            while  subMenuInput not in ['1','2','3','4','5','b']:
                if vascularNetwork.getVariableValue('vessels') == {}:
                    print " there are no vessels defined and thus no boundarys available";break
                else:
                    # evaluate boundarys in Network
                    boundarys = []
                    notDefinedBoundarys = []
                    
                    boundarys.extend(vascularNetwork.boundaryVessels)
                    if vascularNetwork.root != None and vascularNetwork.root not in boundarys: 
                        boundarys.append(vascularNetwork.root)
                        
                    boundarysSaved = vascularNetwork.boundaryConditions.keys()
                    
                    # update saved boundary conditions
                    for boundarysCurrent in boundarys:
                        if boundarysCurrent not in boundarysSaved:
                            print " boundary added to vascularNetwork"
                            vascularNetwork.boundaryConditions[boundarysCurrent] = []
                            boundarysSaved.append(boundarysCurrent)
                            
                        if vascularNetwork.boundaryConditions[boundarysCurrent] == []:
                            notDefinedBoundarys.append(boundarysCurrent)
                        
                    nonBoundarys = list(set(boundarys).symmetric_difference(set(boundarysSaved)))
                    for nonBoundary in nonBoundarys:
                        print " boundary removed from vacularNetwork"
                        del(vascularNetwork.boundaryConditions[nonBoundary])
                        if nonBoundary in notDefinedBoundarys: notDefinedBoundarys.remove(nonBoundary)
                        
                    vascularNetwork.evaluateConnections()
                    print ""
                    print "    sub menu: boundary conditions"
                    print ""
                    print "     [1] - show  boundary conditions"
                    print "     [2] - add   boundary condition "
                    print "     [3] - del   boundary condition "
                    print "     [4] - load  boundary conditions from CSV"
                    print "     [5] - write boundary conditions to CSV"
                    print "     [b] - back to the main menu"
                    print ""     
                    subMenuInput = raw_input("     what to do? ") 
                    
                    if subMenuInput == '1':
                        print "     boundary conditions"
                        pprint.pprint(vascularNetwork.boundaryConditions)
                        subMenuInput = ''
                        
                    elif subMenuInput == '2' and vascularNetwork.root != []:
                        
                        print "     add   boundary condition"
                        print ""
                        
                        definedBoundarys = list(set(notDefinedBoundarys).symmetric_difference(set(vascularNetwork.boundaryConditions.keys())))
                        print "     vessels with defined boundary condition:"
                        print "       ",'  '.join(str(i) for i in definedBoundarys)
                        
                        print "     vessels with undefined boundary condition:"
                        print "       ",'  '.join(str(i) for i in notDefinedBoundarys)
                        print ""
                                            
                        existing = False
                        vesselId = raw_input(" enter existing vessel id: ")
                        while True:
                            try:
                                vesselId = int(vesselId)
                                if vesselId in vascularNetwork.vessels:
                                    break
                                else:
                                    existing = True
                            except ValueError:
                                print " TYPE-ERROR: vessel id must be type(int) not type(string)"
                                vesselId = raw_input(" enter existing vessel id: ")
                            if existing == True:
                                print " the vessel does not exist"
                                vesselId = raw_input(" enter existing vessel id: ")
                                existing = False
                                                
                        inputType = '0'
                        print "     add boundary condition type:"
                        print ""
                        # get all defined boundaryConditions from constants-dict save as bcTypes
                        bcTypesAll = nxml.bcTagsClassReferences.keys()
                        bcTypes = []
                        for bcType in bcTypesAll: 
                            if "_" is not bcType[0]:
                                bcTypes.append(bcType)
                        bcTypes.sort()
                        # show all boundaryConditions in the bcTypes
                        index = 0
                        for key in bcTypes:
                            print "       [",str(index).rjust(2),"]    ",key
                            index = index+1
                        # get user input and check if it was correct to define the bcType 
                        existing = False
                        inputType = raw_input ("      choose type ")
                        while True:
                            # check if right input
                            try:
                                inputType = int(inputType)
                                if inputType in np.linspace(0,len(bcTypes)-1,len(bcTypes)):
                                    break
                                else:
                                    existing = True
                            # if not int but string
                            except ValueError:
                                print "      TYPE-ERROR: vessel id must be type(int) not type(string)"
                                inputType = (raw_input ("      choose type "))
                            # if int but to low or high
                            if existing == True:
                                print "       the type does not exist"
                                inputType = (raw_input ("      choose type "))
                                existing = False
                        
                        bcType = bcTypes[int(inputType)] 
                        
                        boundaryInstance = eval(nxml.bcTagsClassReferences[bcType])()
                        boundaryDataDict = {}
                        boundaryDataDict['name']= bcType
                        
                        print ""
                        print "      set values for the BC condition: ", bcType
                        print "          enter 'b' for the first value to skip this procedure"
                        question = True
                        for arg in nxml.bcTags[bcType]:
                            if question == True: 
                                currValue = raw_input (str("            set value for "+str(arg)+' '))
                                if currValue == 'b': question=False
                                test = True
                                try: float(currValue)
                                except:
                                    print '            VALUE or TYPE ERROR, set to None'
                                    test = False
                                if test == True: boundaryDataDict[arg] = (float(currValue))
                                else: boundaryDataDict[arg] = None
                            else: boundaryDataDict[arg] = None
                        if len(vascularNetwork.boundaryConditions.keys()) == 1:
                            print "      set position of the BC condition"
                            position = '2'
                            while position not in ['1','0']:
                                position = raw_input ("          enter '0' for the start or '1' for the end of the vessel ")
                            if position == '1':
                                #bcType = ''.join(['_',bcType])
                                boundaryDataDict['name']= ''.join(['_',bcType])
                        
                        print bcType,boundaryDataDict
                        
                        boundaryInstances = []
                        boundaryInstance.update(boundaryDataDict)
                        boundaryInstances.append(boundaryInstance)
                        
                        if vesselId not in vascularNetwork.boundaryConditions.keys():
                            vascularNetwork.boundaryConditions[vesselId] = boundaryInstances
                        else:
                            vascularNetwork.boundaryConditions[vesselId].extend(boundaryInstances)
                        
                        #if vascularNetwork.getVariableValue('root') != None:
                        mainGraph.update_graph(vascularNetwork, window)
                        subMenuInput = ''
                            
                    
                    elif subMenuInput == '3' and vascularNetwork.root != []:
                        print "     delete boundary condition"
                        print ""
                        pprint.pprint(vascularNetwork.boundaryConditions)
                        
                        vesselId = -1
                        while vesselId not in vascularNetwork.boundaryConditions.keys():
                            vesselId = int(raw_input ("      choose vessel id "))
                        
                        bcs = vascularNetwork.boundaryConditions[vesselId]
                        if bcs != []:                 
                            print ""
                            index = 0
                            for bc in bcs:
                                print "       [",str(index).rjust(2),"]    ",bc.name
                                index = index+1
                            print ""
                            inType = '0'
                            while inType not in np.linspace(0,len(bcs)-1,len(bcs)):
                                inType = int(raw_input ("      choose condition to delete "))
                             
                            print ""
                            print "     boundary condition ",bcs[inType]," removed!"
                            print ""
                            vascularNetwork.boundaryConditions[vesselId].remove(bcs[inType])
                            
                            mainGraph.update_graph(vascularNetwork, window)
                        else:
                            print "     nothing to delete!"
                        subMenuInput = ''
                    
                    elif subMenuInput == '4' and vascularNetwork.root != []:
                        print "     load  boundary conditions from CSV"
                        print ""
                        filename = enterFilename(filename,'')
                        boundaryConditions,boundaryConditionPolyChaos = readBCFromCSV(filename = ''.join([filename,'BC','.csv']))
                        vascularNetwork.update({'boundaryConditions':boundaryConditions,
                                                'boundaryConditionPolyChaos':boundaryConditionPolyChaos})
                        
                        mainGraph.update_graph(vascularNetwork, window)
                                                                
                    elif subMenuInput == '5' and vascularNetwork.root != []:
                        print "     write boundary conditions to CSV"
                        print ""
                        filename = enterFilename(filename,'')
                        boundaryConditions = vascularNetwork.getVariableValue('boundaryConditions')
                        boundaryConditionPolyChaos = deepcopy(vascularNetwork.getVariableValue('boundaryConditionPolyChaos'))
                        writeBCToCSV(boundaryConditions,boundaryConditionPolyChaos,filename = ''.join([filename,'BC','.csv']))         
                        
                    elif subMenuInput == 'b':
                        break
        
        elif menuInput == "f":
            
            subMenuInput = ''
            while  subMenuInput not in ["1","2","b"]:
                print ""
                print "    sub menu: set global fluid properties"
                print ""
                print "     [1] - set all"
                print "     [2] - set individual"
                print "     [b] - back to the main menu"
                print ""
                print "    current fluid properties:"
                for key,value in vascularNetwork.globalFluid.iteritems():
                    print "     {0:20}     {1:10}  {2:10}".format(key,value,variableUnits[key])
                print ""
                subMenuInput = raw_input("    what to do? ")
                
                if subMenuInput == '1':
                    print "     set all fluid properties"
                    
                    for key in vascularNetwork.globalFluid.keys():
                        inputType = "1"
                        typeFalse = False
                        while typeFalse == False:
                            try: 
                                inputType = raw_input("      type value for property: {} with unit {} : ".format(key,variableUnits[key]))
                                inputType = float(inputType)
                                typeFalse = True
                                vascularNetwork.globalFluid[key] = inputType
                            except: pass
                    subMenuInput = ''
                    
                elif subMenuInput == '2':
                    print "     set individual fluid property:"
                    i = 0
                    properties = vascularNetwork.globalFluid.keys()
                    for property_i in properties:
                        print "          [",i,'] - ',property_i
                        i = 1+i
                    inputType = 0
                    while inputType not in [str(i) for i in range(0,len(properties))]:
                        inputType = raw_input("      choose property to set ")
                    
                    inputProperty = ""   
                    inputFalse = False
                    while inputFalse == False:
                        try: 
                            inputProperty = raw_input("      type value for property: {} with unit {} : ".format(properties[int(inputType)],variableUnits[properties[int(inputType)]]))
                            inputProperty = float(inputProperty)
                            vascularNetwork.globalFluid[properties[int(inputType)]] = inputProperty
                            inputFalse = True
                        except: pass
                    subMenuInput = ''
                     
                elif subMenuInput == 'b':
                    break
                
            
        elif menuInput == "l":
            try:
                FILE = open('.recentFilenames.pickle',"rb")
                # store pickle
                recentFilenames = cPickle.load(FILE)
                FILE.close()
            except:
                recentFilenames = []
                
            subMenuInput = ''
            print ""
            print "    sub menu: load data"
            print ""
            print "     [1] - load network from XML"
            print "     [2] - load vessel data from CSV"
            print "     [3] - load vessel data and boundary conditions from CSV"
            print "     [4] - load network from SolutionData"
            print "     [b] - back to the main menu"
            print ""
            while  subMenuInput not in ["1","2","3","b"]:
                subMenuInput = raw_input("what to do? ")
            
                print ""
                print "         resent used networks"
                i = 1
                for name in recentFilenames:
                    print "          [",i,'] - ',name
                    i = 1+i
                print ""
                if subMenuInput == '1':
                    print "     load from XML"
                    
                    filename = enterFilename(filename,'.xml',recentFilenames = recentFilenames)
                    if filename == None:break
                    # delete the old network
                    del vascularNetwork
                    #load the new network
                    vascularNetwork = loadNetworkFromXML(filename= filename)
                    if vascularNetwork == None:
                        mainGraph.update_graph(None, window)
                        vascularNetwork = VascularNetwork()
                        filename = None
                        break
                    mainGraph.update_graph(vascularNetwork, window)
                    filename, value = filename.split(".",1)
                    break
                
                elif subMenuInput == '2':
                    print "     load vessel data from CSV - non existing vessels are added automatically"
                    print ""
                    filename = enterFilename(filename,'.csv',recentFilenames = recentFilenames)
                    if filename == None:break
                    vascularNetwork.updateNetwork(readVesselDataFromCSV(filename=filename))
                    
                    mainGraph.update_graph(vascularNetwork, window)
                    filename, value = filename.split(".",1)
                    break
                
                elif subMenuInput == '3':
                    print "     load vessel data and boundary conditions from CSV"
                    filename = enterFilename(filename,'.csv',recentFilenames = recentFilenames)
                    if filename == None:break
                    vascularNetwork.updateNetwork(readVesselDataFromCSV(filename=filename))
                    mainGraph.update_graph(vascularNetwork, window)
                    
                    
                    filename, value = filename.split(".",1)
                    boundaryConditions,boundaryConditionPolyChaos = readBCFromCSV(filename = ''.join([filename,'BC','.csv']))
                    vascularNetwork.update({'boundaryConditions':boundaryConditions,
                                            'boundaryConditionPolyChaos':boundaryConditionPolyChaos})
                    
                    mainGraph.update_graph(vascularNetwork, window)
                    break
                
                elif subMenuInput == '4':
                    print "     load network from SolutionData"
                    try:
                        networkName,dataSetNumber = chooseSolutionDataCase()
                        vascularNetwork = loadSolutionDataFile(networkName, dataSetNumber)
                        mainGraph.update_graph(vascularNetwork, window)
                        filename, value = filename.split(".",1)
                    except: print "\n ERROR occured could not open requested network, file does not exist or is out dated"
                    break
                
                elif subMenuInput == 'b':
                    break
            
            if filename != None:    
                if filename not in recentFilenames: recentFilenames.insert(0,filename)
                else: 
                    recentFilenames.remove(filename)
                    recentFilenames.insert(0,filename)
                if len(recentFilenames) > 9: recentFilenames.pop(-1)
                            
                FILE = open('.recentFilenames.pickle',"w")
                # store pickle
                cPickle.dump(recentFilenames, FILE, protocol=2)
                FILE.close()
            
        elif menuInput == "s":
            subMenuInput = ''
            print ""
            print "    sub menu: save data"
            print ""
            print "     [1] - write to XML"
            print "     [2] - write vessel data to CSV"
            print "     [3] - write vessel data and boundary conditions to CSV"
            print "     [4] - write graph to .png and .dot files"
            print "     [b] - back to the main menu"
            print ""
            while subMenuInput not in ["1","2","3","b"]:
                subMenuInput = raw_input("what to do? ")
                     
                if subMenuInput == '1':
                    print "     write to XML"
                    filename = enterFilename(filename,'.xml')
                    if filename == None:break
                    writeNetworkToXML(vascularNetwork,filename = filename)
                    filename, value = filename.split(".",1)
                    break
                    
                elif subMenuInput == '2':
                    print "     write vessel data to CSV"
                    filename = enterFilename(filename,'.csv') 
                    if filename == None:break
                    writeVesselDataToCSV(vessels = vascularNetwork.vessels, filename = filename)
                    filename, value = filename.split(".",1)
                    break
                
                elif subMenuInput == '3':
                    print "     write vessel data and boundary conditions to CSV"
                    filename = enterFilename(filename,'.csv') 
                    if filename == None:break
                    writeVesselDataToCSV(vessels = vascularNetwork.vessels, filename = filename)
                    filename, value = filename.split(".",1)
                    boundaryConditions = vascularNetwork.getVariableValue('boundaryConditions')
                    boundaryConditionPolyChaos = deepcopy(vascularNetwork.getVariableValue('boundaryConditionPolyChaos'))
                    writeBCToCSV(boundaryConditions,boundaryConditionPolyChaos,filename = ''.join([filename,'BC','.csv']))  
                    break
                
                elif subMenuInput == '4':
                    print "     write graph to .png and .dot files"
                    pictureName = str(raw_input("     enter pictureName (only!):"))
                    if pictureName == "":
                        pictureName = 'pydotTest'
                    mainGraph.graph.write(graphPath+filename+'/'+pictureName+'.dot')
                    mainGraph.graph.write_png(graphPath+filename+'/'+pictureName+'.png')
                    break
                    
                if subMenuInput == 'b':
                    break
        
        elif menuInput == "u":
            subMenuInput = ''
            print ""
            print "    sub menu: update XML from CSV"
            print ""
            print "     load from XML"
            
            filename = enterFilename(filename,'.xml')
            if filename == None:break
            # delete the old network
            del vascularNetwork
            #load the new network
            vascularNetwork = loadNetworkFromXML(filename= filename)
            if vascularNetwork == None:
                vascularNetwork = VascularNetwork()
                break
            
            mainGraph.update_graph(vascularNetwork, window)
            
            filename, value = filename.split(".",1)
            
            print "     load vessel data from CSV - non existing vessels are added automatically"
            filenameCSV = ''.join([filename,'.csv'])
            vascularNetwork.updateNetwork(readVesselDataFromCSV(filename=filenameCSV))
            print "     load boundaryData from csv as well? press [u]" 
            subMenuInput = raw_input("yes [u]? ")
            if subMenuInput == 'u':
                boundaryConditions,boundaryConditionPolyChaos = readBCFromCSV(filename = ''.join([filename,'BC','.csv']))
                vascularNetwork.update({'boundaryConditions':boundaryConditions,
                                        'boundaryConditionPolyChaos':boundaryConditionPolyChaos})          
            
            try: mainGraph.update_graph(vascularNetwork, window)
            except:mainGraph.update_graph(None, window)
            
            print "     write to XML"
            filenameXML = ''.join([filename,'.xml'])
            writeNetworkToXML(vascularNetwork,filename = filenameXML)
            
            
    print "bye bye .."

if __name__ == '__main__':
    main()
