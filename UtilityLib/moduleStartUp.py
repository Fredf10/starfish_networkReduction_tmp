#############################################################################
#
# moduleStartUp
#
# includes option parser function and 
# provides functions to ask user for need data such as NetworkName, dataNumber Visualisation cases
# if not given with options while staring Main or Visualistion or VascularPolynomialChaos
#
# created by Vinzenz Eck // vinzenz.g.eck@ntnu.no
##

import os
cur = os.path.dirname( os.path.realpath( __file__ ) )
import inspect
from pprint import pprint as pp

from moduleXML import savePolyChaosXML



from optparse import OptionParser

def parseOptions(activeOptions, visualisationOnly = False):
    '''
    parse options for the code
    
    enable options to parse with str in activeOptions
    input:     activeOptions := list with str 
    
            'f' : parse filename
            'n' : dataNumber
            'd' : simulation description
            'r' : resimulated case
            's' : save
            'v' : visualisation type
            'c' : connect visualisations       
    
            visualisationOnly := bool if True proposal of visualisation cases are made if non is given
            
    return: dict with options and arguments
    
    Usage e.g., :
    
        optionsDict = parseOptions(['f','n','d','s','v','r'])
        
        networkName           = optionsDict['networkName']
        save                  = optionsDict['save']
        dataNumber            = optionsDict['dataNumber']
        simulationDescription = optionsDict['simulationDescription']
        vizOutput             = optionsDict['vizOutput']
        resimulate            = optionsDict['resimulate']
    '''
    
    parser = OptionParser()
    
    for activeOpotion in activeOptions:
        if activeOpotion == 'f':
            parser.add_option("-f", "--file", dest='networkName',
                              help = "open file with given network name")
        elif activeOpotion == 'n':
                parser.add_option("-n", "--dataNumber", dest='dataNumber',
                                  help = "number of the solution data (last number in filename), default = 999, max 999;")
        elif activeOpotion == 's':
            parser.add_option("-s", "--save", action="store_true", dest="save", 
                              help = "if set solution data is saved")
        elif activeOpotion == 'd':
                parser.add_option('-d', '--description', dest='description',
                                  help = "simulation case description; NB: no space subported")
        elif activeOpotion == 'v':
                parser.add_option('-v', '--vizBool', dest='vizBool',
                                  help = "choose visualisation mode, 0: no visualisation, 1: 2d and 3d, 2: 2d plots, 3: 3d visualisation default = 1")
        elif activeOpotion == 'c':       
                parser.add_option("-c", "--connect", action="store_true",  dest='connect',
                                  help="connect to 3dViz (True) or not (False); currently not working")
        elif activeOpotion == 'r':
            parser.add_option("-r", "--resimulate", action="store_true", dest="resimulate", 
                              help = "resimulate case with same network saved in datanumber file, 0 = False, 1 = True")

    (options, args) = parser.parse_args()
    optionsDict = options.__dict__
        
    networkName = None
    save        = False
    dataNumber  = False
    dataSetNumber = False
    vizOutput   = 'non'
    simulationDescription = None
    connect     = False
    resimulate  = False
    
    
    for option,optionArgument in optionsDict.iteritems():
        # -f network name
        if option == 'networkName':
            if optionArgument != None:
                networkName = optionArgument
        # -n dataNumber
        elif option == 'dataNumber':
            dataNumber,dataSetNumber = evaluateDataNumber(optionArgument)  
            if dataNumber == None:
                dataNumber = '999'
            else: save = True     
        # -s save
        elif option == 'save':
            if optionArgument != None:
                save = optionArgument
            #save solution data and the vascularNetwork in c pickle, if no save take temporary slot 999
            if save == False:
                dataNumber = '999'
        # -d simulation Description
        elif option == 'description':
            if optionArgument != None:
                simulationDescription = optionArgument
        # -v visialisation type        
        elif option == 'vizBool':
            if optionArgument != None:
                if optionArgument == '1':
                    vizOutput = "2D+3D"
                elif optionArgument == '2':
                    vizOutput = "2D"
                elif optionArgument == '3':
                    vizOutput = "3D"
        # -c connect
        elif option == 'connect':
            if optionArgument != None:
               connect = optionArgument
        # -r resimulate
        elif option == 'resimulate':
            if optionArgument != None:
                    resimulate = optionArgument
                            
    # catch up non given but necessary options
    if networkName == None and visualisationOnly == False:
        networkName = chooseNetworkName()
    if visualisationOnly == True and (dataSetNumber == None or networkName == None): 
            print "\n  No networkName passed, choose between all available networks:"
            print "  (NB: use -f networkName to define a specific file you want to open) \n"
            networkName,dataNumber = chooseSolutionDataCase()
    if simulationDescription == None and visualisationOnly == False:
        simulationDescription = defineSimulationDescription()
        
    del parser
    
    return {'networkName'           : networkName,
            'save'                  : save,
            'dataNumber'            : dataNumber,
            'dataSetNumber'         : dataSetNumber,
            'vizOutput'             : vizOutput,
            'simulationDescription' : simulationDescription,
            'connect'               : connect,
            'resimulate'            : resimulate}


def chooseNetworkName():
    '''
    console Interface to choose a VascularNetwork for simulation / vascularPolynomial Chaos
    '''
    path = ''.join([cur,'/../','NetworkFiles/'])
    dirNames = os.listdir(path)
    i=0
    for dirName in dirNames:
        print "   [",dirNames.index(dirName),'] - ',dirName
        i=i+1
    filenameT = str(raw_input("\n  Choose simulation case you want to open according to its number: "))
    try:
        filenameT = int(filenameT)
    except:
        print "no integer given, system exit"
        exit()
    try:
        networkName = dirNames[filenameT]
    except:
        print "given number not in list, system exit"
        exit()
    print ""
    print '====================================='
    return networkName

def evaluateDataNumber(dataNumberString):
    '''
    Function to evaluate DataNumbers given as a string
    Max lenght of dataNumber = 3
    Input:
        dataNumber String sparate datanumbers with ,
        (e.g. '12,23' or '4')
    output:
        dataNumber = String (of the first one)
        dataSetNumber = [ String, ... ,String] 
        (e.g. dataNumber = '012' dataSetNumber = ['012','023']
              dataNumber = '004' dataSetNumber = ['004'] 
    '''
    dataSetNumber = None
    dataNumber = None
    if dataNumberString != None:
        dataNumber = dataNumberString
        if ',' in dataNumber:
            dataNumbers = dataNumber.split(',')
            dataSetNumber = []
            for dataNum in dataNumbers:
                if len(dataNum) < 4:
                    dataSetNumber.append(dataNum.zfill(3))
                else:
                    print 'ERROR: Datanumer to high! system exit'
                    exit()
            dataNumber = dataSetNumber[0]
        else:
            dataSetNumber = [dataNumber.zfill(3)]
            dataNumber = dataNumber.zfill(3)
            if len(dataNumber) > 3:
                print "ERROR: dataNumber is to high! system exit"
                exit()      
    return dataNumber,dataSetNumber
  
  
def defineSimulationDescription():
    simulationDescription = str(raw_input("\n  Type in description of the simulation case: "))
    try:
        simulationDescription = str(simulationDescription)
    except:
        print 'ERROR: no String Convertable input given, system exit'
        exit()
    if simulationDescription in ['',' ']: simulationDescription = '-'
    return simulationDescription
    
    
def chooseVPCconfigFile(networkName):
    '''
    console Interface to choose a vascularPolynomialChaos Config File
     including the possility to create a template Config File
    Input:
        networkName of VascularNetwork
    Output:
        networkName,dataNumber of the Config File (networkName should be the same)
    '''
    path = ''.join([cur,'/../','NetworkFiles/',networkName,'/'])
    allFilenames = os.listdir(path)
    filenames = []
    for filename in allFilenames:
        if ".xml" in filename and "vpcConfig" in filename:
            filenames.append(filename)
    if filenames != []:
        print "\n  No dataNumber for Config-File passed, choose between all available Config Files:"
        print "  (NB: use -n dataNumber to define a specific Config-file you want to open)\n"
        print '   [ 0 ] - Create template Config-File and exit'
        print '   [ 1 ] - Create template Config-File and run vascularPolynomialChaos'
        index = 2
        for filename in filenames:
                print "   [",str(index),'] -',filename
                index = index+1
        filenameT = str(raw_input("\n  Choose Config-File you want to open according to its number: "))
        try:
            filenameT = int(filenameT)
        except:
            print "no integer given, system exit"
        if filenameT == 0 or filenameT == 1:
            dataNumber = '999'
            vpcConfigFilenameTemplate = '_'.join([networkName,'vpcConfig',dataNumber,'template'])
            savePolyChaosXML({},vpcConfigFilenameTemplate)
            if filenameT == 0: exit()
        else:
            networkName = filenames[filenameT-2]
            dataNumber = networkName.split('.')[0].split('_')[2]
            networkName = filenames[filenameT-2].split('_')[0]
    else:
        print "\n  For this network no vascularPolynomialChaos Config-Files is available!!"
        createTemplate = str(raw_input("  Do you want to create template file? (n=exit) y/n: "))
        if createTemplate == "y":
            dataNumber = '999'
            vpcConfigFilenameTemplate = '_'.join([networkName,'vpcConfig',dataNumber,'template'])
            savePolyChaosXML({},vpcConfigFilenameTemplate)
            print ""
            runPolyChaos = str(raw_input("  Do you want to run vascularPolynomialChaos now? (n=exit) y/n: "))
            if runPolyChaos != 'y':
                print "  .. exiting!"
                exit()
        else: 
            print "  .. exiting!"
            exit()
        
    return networkName,dataNumber   

def chooseSolutionDataCase():
    '''
    console Interface to choose a vascular1DFlow simualtion case for e.g. Visualisation
    Output:
        networkName <string>
        dataNumber  <string>
        
    '''
    path = ''.join([cur,'/../','NetworkFiles/'])
    index = 0
    fileNameDataNumber = []
    for dirName, dirNames, fileNames in os.walk(path):
        # print path to all filenames.
        for fileName in fileNames:
            if ".pickle" in fileName and "polyChaos" not in fileName and "pure" not in fileName:
                fileName = fileName.split('.')[0]
                splitName = fileName.split('_SolutionData_')
                fileName = splitName[0]
                dataNumber = splitName[-1] 
                if len(dataNumber) == 3:
                    print "   [",str(index).rjust(2),'] - ', fileName.ljust(25) ,' : ','DataSet ',dataNumber
                    fileNameDataNumber.append([fileName,dataNumber])
                    index = index+1
    print ""
    inputInt = str(raw_input("  Choose simulation case you want to open according to its number: "))
    try:
        inputInt = int(inputInt)
        networkName = fileNameDataNumber[inputInt][0]
        dataNumbers = fileNameDataNumber[inputInt][1]
    except:
        if inputInt > len(fileNameDataNumber):
            print " given number to high, system exit"
            
        else:    
            print "no integer given, system exit"
        exit()
    
    return networkName,dataNumbers



    
    