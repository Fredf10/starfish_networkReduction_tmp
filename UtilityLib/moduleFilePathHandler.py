#####
# all functions to save, load, maipulate, parse etc. pickle files of simulation cases
# .v1d
####

#############################################################################
#
# moduleFilePathHandler
#
# provide functions to save, load, maipulate, parse etc. pickle files of simulation cases
#
# loadSolutionDataFile(networkName,dataNumbers)
# parseDirectoryForSimulationCases(networkName)
# updateSimulationDescriptions(networkName)
#
#
# created by Vinzenz Eck // vinzenz.g.eck@ntnu.no
##

import cPickle
import os,sys,shutil
cur = os.path.dirname( os.path.realpath( __file__ ) )

from copy import copy as copy 

from pprint import pprint as pp


def getFilePath(fileType, networkName, dataNumber, mode, exception = 'Error'):
    '''
    Function return a requested file path, if this file exists
    
    fileType:
        'vesselCSVFile',
        'boundaryCSVFile',
        'networkXmlFileTemplate',
        'networkXmlFile',
        'solutionFile',
        'configFile',
        'simulationDescriptionFile',
        'vncRescentNetworksFile',
        'vncNetworkGraphFile,
        'vpcConfigXmlFile'
        'vpcNetworkXmlFile'
    
    networkName:
    
        name of the network file
    
    dataNumber:
    
        data number of solution file or xml file
    
    mode:
    
        read or write
        
    exception: (for read mode)
    
        Error (default): raise error and exit if file is not exiting
        Warning: just raise Warning and return with error string      
    '''
    existingFileTypes = ['vesselCSVFile',
                         'boundaryCSVFile',
                         'networkXmlFileTemplate',
                         'networkXmlFile',
                         'solutionFile',
                         'configFile',
                         'simulationDescriptionFile',
                         'vncRescentNetworksFile',
                         'vncNetworkGraphFile',
                         'vpcConfigXmlFile',
                         'vpcNetworkXmlFile']
    
    if fileType not in existingFileTypes:
        raise ValueError("ERROR: getFilePath, requested file type {}\
                          is not in existingFileTypess {}".format(fileType, existingFileTypes))
    
    if '_template' in networkName:
        if fileType == 'networkXmlFile':
            fileType = 'networkXmlFileTemplate'
    if '_vpc' in networkName:
        if fileType == 'networkXmlFile':
            fileType = 'vpcNetworkXmlFile'
            networkName = networkName.split('_vpc')[0]
    
    # file names
    filenames = {
                 'configFile'                : 'STARFiSh.config',
                 'vesselCSVFile'             : ''.join([networkName,'.csv']),
                 'boundaryCSVFile'           : ''.join([networkName,'BC.csv']),
                 'networkXmlFileTemplate'    : ''.join([networkName,'.xml']),
                 'networkXmlFileXXX'         : ''.join([networkName,'.xml']),
                 'networkXmlFileSim'         : ''.join([networkName,'_SolutionData_',dataNumber,'.xml']),
                 'solutionFile'              : ''.join([networkName,'_SolutionData_',dataNumber,'.hdf5']),
                 'simulationDescriptionFile' : ''.join(['simulationCaseDescriptions.txt']),
                 'vncRescentNetworksFile'    : '.recentNetworkNames.pickle',
                 'vncNetworkGraphFile'       : ''.join([networkName,'.png']),
                 ## vascularPolynomialChaos
                 'vpcConfigXmlFile'          : ''.join([networkName,'_vpcConfig_',dataNumber,'.xml']),
                 'vpcNetworkXmlFile'         : ''.join([networkName,'_vpc_',dataNumber,'.xml'])
                 }    
                
    ## if fileType=networkXmlFile check if master (dn = XXX) or simulated is meant
    if fileType == 'networkXmlFile':
        if dataNumber == 'xxx':
            fileType = ''.join([fileType,'XXX'])
        else:
            fileType = ''.join([fileType,'Sim'])
    
    ## find requested file name
    requestedFilename  = filenames[''.join([fileType])]
    ## find directory    
    requestedDirectory = getDirectory(''.join([fileType,'Directory']), networkName, dataNumber, mode, exception = exception)
    if requestedDirectory == None:
        if exception == "Warning":
            print "WARNING: moduleFilePathHandler.getFileAndPaths() directory of file '{}' does not exits. Exit()".format(requestedFilename)
            return None
        elif exception == "No":
            pass
        else:
            raise ValueError("ERROR: moduleFilePathHandler.getFileAndPaths() directory of file '{}' does not exits. Exit()".format(requestedFilename))
              
    ## requested file path 
    requestedFilePath = ''.join([requestedDirectory,'/',requestedFilename])
                
    if mode == 'read':    
        # if mode read
        ## ensure that the file exists
        if not os.path.isfile(requestedFilePath):
            if exception == "Warning":
                print "WARNING: moduleFilePathHandler.getFileAndPaths() file '{}' does not exits. Exit()".format(requestedFilePath)
                return None
            elif exception == "No":
                print "raise no exception"
                return None
            else:
                raise ValueError("ERROR: moduleFilePathHandler.getFileAndPaths() file '{}' does not exits. Exit()".format(requestedFilePath))
              
    return requestedFilePath
    
def getDirectory(directoryType, networkName, dataNumber, mode, exception = 'Error'):
    '''
    Function returns a requested directory path, if this directory does not exists
    it is created.
    
    directoryType:
        'workingDirectory',
        'configFileDirectory',
        'vesselCSVFileDirectory',
        'boundaryCSVFileDirectory',
        'networkXmlFileTemplateDirectory',
        'networkXmlFileXXXDirectory',
        'networkXmlFileSimDirectory',
        'solutionFileDirectory',
        'screenshotDirectory',
        'movieDirectory',
        'simulationDescriptionFileDirectory',
        'vncRescentNetworksFileDirectory',
        'vncNetworkGraphFileDirectory',
        'vpcConfigXmlFileDirectory',
        'vpcNetworkXmlFileDirectory',
    
    networkName:
    
        name of the network file
    
    dataNumber:
    
        data number of solution file or xml file
    
    mode:
    
        read or write
        
    exception: (for read mode)
    
        Error (default): raise error and exit if file is not exiting
        Warning: just raise Warning and return with error string      
    '''
    
    existingDirectoryTypes = {'workingDirectory',
                              'configFileDirectory',
                              'vesselCSVFileDirectory',
                              'boundaryCSVFileDirectory',
                              'networkXmlFileTemplateDirectory',
                              'networkXmlFileXXXDirectory',
                              'networkXmlFileSimDirectory',
                              'solutionFileDirectory',
                              'screenshotDirectory',
                              'movieDirectory',
                              'simulationDescriptionFileDirectory',
                              'vncRescentNetworksFileDirectory',
                              'vncNetworkGraphFileDirectory',
                              'vpcConfigXmlFileDirectory',
                              'vpcNetworkXmlFileDirectory'} 
    
    if directoryType not in existingDirectoryTypes:
        raise ValueError("ERROR: getDirectory, requested directoryType {}\
                          is not in existingDirectoryTypes{}".format(directoryType, existingDirectoryTypes))
    ##definitions
    starfishHomeDirectory = ''.join([cur,'/..'])
    
    if directoryType != 'configFileDirectory': 
        # load working directory from config file
        workingDirectory = readConfigFile(['WorkingDirectory'])['WorkingDirectory']
    else:
        workingDirectory = starfishHomeDirectory
           
    networkXmlFileTemplateDirectory = ''.join([starfishHomeDirectory,'/TemplateNetworks/',networkName])
    networkXmlFileDirectory         = ''.join([workingDirectory,'/',networkName])
    solutionFileDirectory           = ''.join([networkXmlFileDirectory,'/SolutionData_',str(dataNumber)])
    movieDirectory              = ''.join([solutionFileDirectory,'/Movies'])
    screenshotDirectory         = ''.join([solutionFileDirectory,'/Screenshots'])
    polynomialChaosCaseDirectory = ''.join([networkXmlFileDirectory,'/vascularPolynomialChaos_',str(dataNumber)]) 
    ## look up tables
    # directories
    directories = {
                   'workingDirectory'                   : workingDirectory,
                   'configFileDirectory'                : starfishHomeDirectory,
                   'vesselCSVFileDirectory'             : networkXmlFileDirectory,
                   'boundaryCSVFileDirectory'           : networkXmlFileDirectory,
                   'networkXmlFileTemplateDirectory'    : networkXmlFileTemplateDirectory,
                   'networkXmlFileXXXDirectory'         : networkXmlFileDirectory,
                   'networkXmlFileSimDirectory'         : solutionFileDirectory,
                   'solutionFileDirectory'              : solutionFileDirectory,
                   'simulationDescriptionFileDirectory' : networkXmlFileDirectory,
                   # 3d viz
                   'screenshotDirectory'                : screenshotDirectory,
                   'movieDirectory'                     : movieDirectory,
                   # vnc
                   'vncRescentNetworksFileDirectory'    : workingDirectory,
                   'vncNetworkGraphFileDirectory'       : networkXmlFileDirectory,
                   # vascular polynomial chaos
                   'vpcConfigXmlFileDirectory'          : polynomialChaosCaseDirectory,
                   'vpcNetworkXmlFileDirectory'         : polynomialChaosCaseDirectory
                   }
    
    requestedDirectory = os.path.normpath(directories[directoryType])
    
    # if mode write
    if mode == 'write':
        ## ensure that the directory exists
        if not os.path.exists(requestedDirectory):
            os.makedirs(requestedDirectory)  
    if mode == 'read':
        if not os.path.exists(requestedDirectory):
            requestedDirectory = None
    
    return requestedDirectory

def createWorkingCopyOfTemplateNetwork(templateNetworkName, destinationNetworkName = None):
    '''
    Function which copys all data from a template network wtih template network name into
    the working directory.
    It uses the same name of the template network if destinationNetworkName is not defined
    '''
    if destinationNetworkName == None:
        destinationNetworkName = templateNetworkName.split('_template')[0]
    
    pathTemplateNetwork     = getDirectory('networkXmlFileTemplateDirectory', templateNetworkName, 'xxx', 'read')
    pathDestinationNetwork  = getDirectory('networkXmlFileXXXDirectory', destinationNetworkName, 'xxx', 'write')
    
    #loop through files
    for file in os.listdir(pathTemplateNetwork):
        # remove _template from name
        renamedFile = ''.join(file.split('_template'))
        # check if new name needs to be applied
        oldName = templateNetworkName.split('_template')[0]
        if oldName in renamedFile:
            renamedFile = ''.join([destinationNetworkName,renamedFile.split(oldName)[-1]])
            
        shutil.copy('/'.join([pathTemplateNetwork,file]), '/'.join([pathDestinationNetwork,renamedFile]))
        
    return destinationNetworkName


def readConfigFile(options):
    '''
    Function to read from options from STARFiSh.ini config file
    input:
        options = list with options
                existing options: 'WorkingDirectory'
    output:
        configurations = dict with {option: configuration from file}
    ''' 
    import ConfigParser
    config = ConfigParser.ConfigParser()
    config.read(getFilePath('configFile', "", '', 'read'))
        
    
    configurations = {}
    
    for option in options:    
        if option == 'WorkingDirectory':
            try:
                workingDirectory = config.get('Directory Paths', option)
            except:
                workingDirectory = None
                raise ValueError("ERROR pathAndFilenameHandler.readConfigFile reading WorkingDirectory failed ini file corrupted, exit()")
            if workingDirectory == '':
                raise ValueError("ERROR pathAndFilenameHandler.readConfigFile reading WorkingDirectory failed: no path defined, exit()")
            configurations['WorkingDirectory'] = workingDirectory
            
    return configurations    

def saveConfigFile(configurations):
    '''
    Function to save configurations to options in the STARFiSh.ini config file
    The file will update the existing config file, it is not neccessary to pass all 
    configurations which exist.
    input:
        configurations = dict with {option: configuration from file}
        
    ''' 
    # open config to get current states
    existingOptions = ['WorkingDirectory']
    
    import ConfigParser
    Config = ConfigParser.ConfigParser()
    
    configFilePath = getFilePath('configFile', '','',  'read',exception = 'No')
    
    
    
    if configFilePath is not None:  #  file exists
        Config.read(configFilePath)
    else: #  file does not exist
        Config.add_section('Directory Paths')
        
    for option,config in configurations.iteritems(): 
            if option in existingOptions:
                Config.set('Directory Paths', option, config)   
            
    with open(getFilePath('configFile', '','',  'write'), 'wb') as configfile:
        Config.write(configfile)
    
def updateSimulationDescriptions(networkName, currentDataNumber, currentDescription):
    '''
    Function to update the text-file with the simulation description for the given network:
    Input:
        networkName <String> (= name of the network to process)
    
    workflow:
        1. open all pickle files and write out the simulation descriptions and datanumbers
        2. write information into file
    '''
    # open File
    #simCaseDescFilePath = getFilePath('simulationDescriptionFile', networkName, currentDataNumber, 'read')#, exception = 'No')
    try:
        simCaseDescFilePath = getFilePath('simulationDescriptionFile', networkName, currentDataNumber, 'read', exception = 'No')
        simCaseDescFile = open(simCaseDescFilePath, 'r')
    except:
        simCaseDescFilePath = getFilePath('simulationDescriptionFile', networkName, currentDataNumber, 'write')#, exception = 'No')
    
        simCaseDescFile = file(simCaseDescFilePath, 'w+')
        simCaseDescFile.write("DataNumber   Description \n")
        simCaseDescFile.close()
        simCaseDescFile = open(simCaseDescFilePath, 'r')
            
    alreadyWritten = False
    # check if datanumber in file:
    simCaseDescriptionFileLines = simCaseDescFile.readlines()
    simCaseDescFile.close() 
    
    for Line in simCaseDescriptionFileLines:
        if currentDataNumber in Line: 
            if currentDescription not in Line or currentDescription == '':
                index = simCaseDescriptionFileLines.index(Line)
                simCaseDescriptionFileLines.remove(Line)
                simCaseDescriptionFileLines.insert(index,"  {} {} \n".format(currentDataNumber.ljust(10),currentDescription))
            alreadyWritten = True
            
    if alreadyWritten == False:
        simCaseDescriptionFileLines.append("  {} {} \n".format(currentDataNumber.ljust(10),currentDescription))
    
    #print simCaseDescriptionFileLines
    
    simCaseDescFile = open(simCaseDescFilePath, 'w')
    simCaseDescFile.writelines(simCaseDescriptionFileLines)
    simCaseDescFile.close()
    
def getSimulationCaseDescriptions(networkName, exception = 'Warning'):
    '''
    
    '''
    simCaseDescFilePath = getFilePath('simulationDescriptionFile', networkName, 'xxx', 'read', exception = 'No')
    try:
        simCaseDescFile = open(simCaseDescFilePath, 'r')
    except:
        if exception == 'Warning':
            print "WARNING getSimulationCaseDescriptions() simulation description file of network {} does not exist!".format(networkName)
        elif exception == 'No':
            pass
        else: raise ValueError(exception)
        return None
    simCaseDescriptionFileLines = simCaseDescFile.readlines()
    simCaseDescFile.close() 
    
    dataNumberDescriptionDict = {}
    for Line in simCaseDescriptionFileLines:
        sol = Line.split('{:8}'.format(' '))
        if len(sol)>1:
            #dataNumber
            dataNumberDescriptionDict[sol[0].split(' ')[-1]] = sol[1].split(' \n')[0] 
    
    return dataNumberDescriptionDict
    
   
def regenerateSimulationCaseDescription(networkName):
    '''
    Function to regenerate the simulation case descriptions of a network
    '''
    
        
def loadExternalDataSet(fileName):
    '''
    Function to open external (preprocessed) DataSets with the file ending *.v1dfExD
    Input:
        fileName <string> (with the total path)
    Output:
        extData <dict> 
    '''
    
    externalData = {'PressureTime':None,'Pressure':None,'PressureUnit':'',
                    'FlowTime':None,'Flow':None,'FlowUnit':'',
                    'AreaTime':None,'Area':None,
                    'Description': ''}
    try:
        externalDataFile = open(fileName,'rb')
        externalData = cPickle.load(externalDataFile)
    except:
        print "Error: no or corrupted external data-file found"
    
    return externalData


if __name__ == '__main__':

    print getFileAndPaths('configFile', "networkName", '485', 'read')
    print getFileAndPaths('networkXmlFile', "haiojaHsd", 'xxx', 'read', exception = 'Warning')
    print getFileAndPaths('networkXmlFile', "haiojaHsd", '124','write')
    print getFileAndPaths('solutionFile', "haiojaHsd", '124', 'write')
    print getFileAndPaths('networkXmlFileTemplate', "haiojaHsd", '124', 'write')
    

## Old unneeded functions replaced by functions in classVascularNetwork due to new data handling

# def parseDirectoryForSimulationCases(networkName):
#     '''
#     Function to search for all simulationCases of a given network
#     Input:
#         networkName <String> (= name of the network to process) 
#     Output:
#         simulationCases <Dict> = { dataNumber1<String:SimCase1<String>, ...
#                                    dataNumberN<String:SimCaseN<String>}
#     '''
#     path = ''.join([cur,'/../','NetworkFiles/',networkName,'/','SolutionData','/'])
#     simulationCases = {}
#     
#     for dirName, dirNames, fileNames in os.walk(path):
#         for fileName in fileNames:
#             if ".hdf5" in fileName and "polyChaos" not in fileName and "pure" not in fileName:
#                 dataNumber = fileName.split('.')[0].split('_SolutionData_')[-1]
#                 simulationCases[dataNumber] = ''.join([path,fileName])
#                 
#     return simulationCases