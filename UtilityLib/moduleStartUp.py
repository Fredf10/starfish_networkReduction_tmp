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

import os,shutil,sys
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))

import inspect
from pprint import pprint as pp

import moduleXML as mXML
import moduleFilePathHandler as mFPH

# sys.path.append([cur,'/../VascularPolynomialChaosLib/'])
import VascularPolynomialChaosLib.moduleFilePathHandlerVPC as mFPH_VPC

from optparse import OptionParser

def parseOptions(activeOptions, visualisationOnly = False, vascularPolynomialChaos = False):
    """
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
            'w' : set working directory

            visualisationOnly := bool if True proposal of visualisation cases are made if non is given

    return: dict with options and arguments

    Usage e.g., :

        optionsDict = parseOptions(['f','n','d','s','v','r','w'])

        networkName           = optionsDict['networkName']
        save                  = optionsDict['save']
        dataNumber            = optionsDict['dataNumber']
        simulationDescription = optionsDict['simulationDescription']
        vizOutput             = optionsDict['vizOutput']
        resimulate            = optionsDict['resimulate']
    """
    
    parser = OptionParser()
    
    for activeOption in activeOptions:
        if activeOption == 'f':
            parser.add_option("-f", "--file", dest='networkName',
                              help = "open file with given network name")
        elif activeOption == 'n':
                parser.add_option("-n", "--dataNumber", dest='dataNumber',
                                  help = "number of the solution data (last number in filename), default = 999, max 999;")
        elif activeOption == 's':
            parser.add_option("-s", "--save", action="store_true", dest="save", 
                              help = "if set solution data is saved")
        elif activeOption == 'd':
                parser.add_option('-d', '--description', dest='description',
                                  help = "simulation case description; NB: no space subported")
        elif activeOption == 'v':
                parser.add_option('-v', '--vizBool', dest='vizBool',
                                  help = "choose visualisation mode, 0: no visualisation, 1: 2d and 3d, 2: 2d plots, 3: 3d visualisation default = 1")
        elif activeOption == 'c':       
                parser.add_option("-c", "--connect", action="store_true",  dest='connect',
                                  help="connect to 3dViz (True) or not (False); currently not working")
        elif activeOption == 'r':
            parser.add_option("-r", "--resimulate", action="store_true", dest="resimulate", 
                              help = "resimulate case with same network saved in datanumber file, 0 = False, 1 = True")
            
        elif activeOption == 'w':
            parser.add_option("-w", "--workingDirectory", dest="workingDirectory", 
                              help = "set the absolute path of your working Directory where you the networkfiles are stored")

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
        elif option == 'workingDirectory':
            if optionArgument != None:
                print "Setting new working directory"
                if os.path.isdir(optionArgument):
                    mFPH.saveConfigFile({'WorkingDirectory':optionArgument})
                    print "   working directory set!"
                else:
                    print "  working directory does not exist! try to create folder"
                    try:
                        os.mkdir(optionArgument)
                        mFPH.saveConfigFile({'WorkingDirectory':optionArgument})
                        print "   created working directory folder sucessfully"
                        print "   working directory set!"
                    except:
                        print "  WARNING: moduleStartUp.parseOptions() could not set WorkingDirectory {} directory does not exists!".format(optionArgument)
                exit()
                    
    # catch up non given but necessary options
    ## simulation and visualisation
    if visualisationOnly == False and vascularPolynomialChaos == False:
        if networkName == None:
            networkName = chooseNetwork()
        if simulationDescription == None:
            simulationDescription = defineSimulationDescription()
    ## visualisation only
    if visualisationOnly == True and (dataSetNumber == None or networkName == None) and vascularPolynomialChaos == False: 
            print "\n  No networkName passed, choose between all available networks:"
            print "  (NB: use -f networkName to define a specific file you want to open) \n"
            networkName,dataNumber = chooseSolutionDataCase()
    ## polynomial chaos
    if vascularPolynomialChaos == True:
        if networkName == None:
            networkName = chooseNetwork(showTemplates = False)
        if dataSetNumber == None:
            dataNumber = chooseVPCconfigFile(networkName)
        
    del parser
    
    return {'networkName'           : networkName,
            'save'                  : save,
            'dataNumber'            : dataNumber,
            'dataSetNumber'         : dataSetNumber,
            'vizOutput'             : vizOutput,
            'simulationDescription' : simulationDescription,
            'connect'               : connect,
            'resimulate'            : resimulate}


def chooseNetwork(showTemplates = True):
    """
    console Interface to choose a VascularNetwork for simulation / vascularPolynomial Chaos
    """
    dirNamesTemplate = []
    if showTemplates:
        # network templates
        templatePath = mFPH.getDirectory('networkXmlFileTemplateDirectory','','','read')
        dirNamesTemplate = [d for d in os.listdir(templatePath) if '.' not in  d]
        
        nTemplates = len(dirNamesTemplate)
    
        print "\n Template Networks: \n"
        for dirName in dirNamesTemplate:
            print "   [ {:3} ] - {}".format(dirNamesTemplate.index(dirName),dirName)
    
    # working directory
    workingDirectoryPath = mFPH.getDirectory('workingDirectory','','','read')
    dirWorkingDirectory = [d for d in os.listdir(workingDirectoryPath) if '.' not in  d]
    
    dirNames = dirNamesTemplate+dirWorkingDirectory
    
    
    print "\n WorkingDirectory Networks: \n"
    for dirName in dirWorkingDirectory:
        print "   [ {:3} ] - {}".format(dirNames.index(dirName),dirName)
        
    filenameT = str(raw_input("\n  Choose network you want to open according to its number: "))
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

# TODO: (einar) fix exception variable
def evaluateDataNumber(dataNumberString, exception = "Error"):
    """
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
    """
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
                    if exception == "Error":
                        raise ValueError('moduleStartUp.evaluateDataNumber. Datanumer {} to high! system exit'.format(dataSetNumber))
                        exit()
                    elif exception == 'Warning':
                        print 'moduleStartUp.evaluateDataNumber. Datanumer {} to high'.format(dataSetNumber)
                        return False,False
                    else:
                        raise Exception
                    
            dataNumber = dataSetNumber[0]
        else:
            dataSetNumber = [dataNumber.zfill(3)]
            dataNumber = dataNumber.zfill(3)
            if len(dataNumber) > 3:
                if exception == "Error":
                    raise ValueError('moduleStartUp.evaluateDataNumber. Datanumer {} to high! system exit'.format(dataSetNumber))
                    exit()
                elif exception == 'Warning':
                    print 'moduleStartUp.evaluateDataNumber. Datanumer {} to high'.format(dataSetNumber)
                    return False,False
                else:
                    raise Exception   
                
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
    
    
def chooseSolutionDataCase():
    """
    console Interface to choose a vascular1DFlow simualtion case for e.g. Visualisation
    Output:
        networkName <string>
        dataNumber  <string>
    """
    workingDirectory = mFPH.getDirectory('workingDirectory','','','read')
    networkCases = [d for d in os.listdir(workingDirectory) if '.' not in  d]
        
    enumerationIndex = 0
    fileNameDataNumber = []
    
    for networkName in networkCases:
        simulationCaseDict = mFPH.getSimulationCaseDescriptions(networkName )#, exception = 'No')
        networkDirectory = mFPH.getDirectory('networkXmlFileXXXDirectory',networkName,'xxx','read')
        
        if simulationCaseDict != None:
            first = True
            for root, dirs, files in os.walk(networkDirectory):
                for file in files:
                    if ".hdf5" in file and "polyChaos" not in file:
                        if first: 
                            print "\n        {}".format(networkName)
                            first = False
                        file = file.split('.')[0]
                        dataNumber = file.split('_SolutionData_')[-1]
                        
                        if len(dataNumber) == 3:
                            if dataNumber not in simulationCaseDict:
                                description =  "'{}' not listed in simulation descriptions of network '{}'.".format(dataNumber,networkName)
                                # update descriptions
                            else:
                                description = simulationCaseDict[dataNumber]
                            print "[ {:3} ]     {} : {}".format(enumerationIndex, dataNumber, description)
                            fileNameDataNumber.append([networkName,dataNumber])
                            enumerationIndex = enumerationIndex+1
                        
                       
    inputInt = str(raw_input("  Choose simulation case you want to open according to its number: "))
    try:
        inputInt = int(inputInt)
        networkName = fileNameDataNumber[inputInt][0]
        dataNumber = fileNameDataNumber[inputInt][1]
    except:
        if inputInt > len(fileNameDataNumber):
            print " given number to high, system exit"
            
        else:    
            print "no integer given, system exit"
        exit()
    
    return networkName,dataNumber


def chooseVPCconfigFile(networkName):
    """
    console Interface to choose a vascularPolynomialChaos Config File
     including the possility to create a template Config File
    Input:
        networkName of VascularNetwork
    Output:
        networkName,dataNumber of the Config File (networkName should be the same)
    """
    workingDirectory = mFPH.getDirectory('workingDirectory','','','read')
    networkDirectory = '/'.join([workingDirectory,networkName])
    
    filesNetworkDir = os.listdir(networkDirectory)
    filenames = []
    
    for fileNetworkDir in filesNetworkDir:
        # check if polychaos directory:
        if 'vascularPolynomialChaos' in fileNetworkDir:
            allFilenames = os.listdir('/'.join([workingDirectory,networkName,fileNetworkDir]))
            for filename in allFilenames:
                if ".xml" in filename and "vpcConfig" in filename:
                    filenames.append(filename)
            
    print "\n  No dataNumber for Config-File passed, choose between all available Config Files:"
    print "  (NB: use -n dataNumber to define a specific Config-file you want to open)\n"
    print "   [   0 ] - Create new template Config-File"
    #print "   [ 1 ] - Create template Config-File and run vascularPolynomialChaos."
    index = 1
    if filenames != []:
        for i,filename in enumerate(filenames):
                print "   [ {:3} ]     {}".format(i+1, filename)
        index = 1+len(filenames)
    print index, range(index)
    
    userInput = 'a'
    while userInput not in [str(i) for i in xrange(index)]:
        userInput = str(raw_input("\n  Choose Option or Config-File you want to open according to its number: "))
    try:
        userInput = int(userInput)
    except:
        print "no integer given, system exit"
    if userInput == 0 :#or userInput == 1:
        
        userInputDataNumber = 'xxxx'
        dataNumber = False
        while dataNumber == False:
            userInputDataNumber = str(raw_input("\n  Insert dataNumber for polynomial Chaos case (overwrites if existing): "))
            dataNumber = evaluateDataNumber(userInputDataNumber, exception = "Warning")[0]
        # save polychaos config file
        vpcConfigXmlFile =  mFPH_VPC.getFilePath('vpcConfigXmlFile', networkName, dataNumber, 'write')
        mXML.savePolyChaosXML(vpcConfigXmlFile,networkName,dataNumber)
        # copy network file
        toCopyFile = mFPH.getFilePath('networkXmlFile', networkName, 'xxx','write')
        destinationFile = mFPH_VPC.getFilePath('vpcNetworkXmlFile', networkName, dataNumber,'write')
        shutil.copy(toCopyFile, destinationFile)
    else:
        networkName = filenames[userInput-2]
        dataNumber = networkName.split('.')[0].split('_')[2]
        
    return dataNumber   
    
    
