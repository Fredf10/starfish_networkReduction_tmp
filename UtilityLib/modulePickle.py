#####
# all functions to save, load, maipulate, parse etc. pickle files of simulation cases
# .v1d
####

#############################################################################
#
# modulePickle
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
import os,sys
cur = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append(cur+'/../NetworkLib')
from classVascularNetwork import VascularNetwork 

from moduleXML import writeNetworkToXML
from moduleXML import loadNetworkFromXML

from processing import memoryUsagePsutil

from copy import copy as copy 

from pprint import pprint as pp


def getDirectoryPath(networkName, type):
    '''
    Function which returns the path for "type" relative to the Main directory!
        'networkDirectory'
        'solutionDirectory'
    '''
    
    networkDirectory  = ''.join(['/NetworkFiles/',networkName])
    solutionDirectory = ''.join(['/NetworkFiles/',networkName,'/SolutionData'])
    
    pathDicts = {'networkDirectory' : networkDirectory,
                 'solutionDirectory': solutionDirectory}
    
    requestedPath = pathDicts[type]
    
    if not os.path.exists(''.join([cur,requestedPath])):
        os.makedirs(''.join([cur,requestedPath]))  
    
    return requestedPath    
    

def saveSolutionDataFile(vascularNetwork, dataNumber):
    '''
    Function to save solutionData after a simulation is run or the solutionData was manipulated
    
    Input: 
        - vascularNetwork <VascularNetwork> (= VascularNetwork instance to save)
        - dataNumbers <String> (= dataNumber with three digits as string)
        
    Output: 
        None
    '''
    print "modulePickle: saveSolutionDataFile "
    networkName = vascularNetwork.name
    
    simulationCaseDescription = vascularNetwork.description
    solutionData = vascularNetwork.prepareToSave()
            
    #standart solution path
    solutionDirectory = ''.join([cur,'/../','NetworkFiles/',networkName,'/SolutionData'])
    if not os.path.exists(solutionDirectory):
        os.makedirs(solutionDirectory)  
    pathSolutionData = ''.join([solutionDirectory,'/',networkName,'_SolutionData_',dataNumber,'.pickle'])
                
    # save solution data
    FILE = open(pathSolutionData,"w")
    cPickle.dump(solutionData, FILE, protocol=2)
    
    FILE.close()
        
    # save to network xml 
    writeNetworkToXML(vascularNetwork, dataNumber = dataNumber, filename = networkName, networkPath = solutionDirectory)


def loadSolutionDataFile(networkName, dataNumber):
    '''
    Function to open SolutionDataFiles for the given Network with the given dataNumbers
    
    Input:
        - networkName <String> (= name of the network to process)
        - dataNumber  <String> (= data number with 3 digits)
    
    Output:
        - vascularNetwork <vascularNetwork> # inclusive solution data stored in the vessel class
    '''
    
    print "load Solution File"
    
    vascularNetwork = None
    
    # 
    simulationCasesAvailable = parseDirectoryForSimulationCases(networkName)    
    try: simCasePathPickle = simulationCasesAvailable[dataNumber]    
    except: 
        print """ERROR: modulePickle.loadSolutionDataFile:
        cannot load solution case of network "{}"
        with data number "{}", it does not exist!,
        system exit()""".format(networkName, dataNumber)
        exit()
    simCasePathXml = '/'.join(simCasePathPickle.split('/')[0:-1])
        
    # 1. load network xml and create vascularNetwork instance
    vascularNetwork = loadNetworkFromXML(filename = networkName, dataNumber = dataNumber, networkPath = simCasePathXml)
    vascularNetwork.initialize()
    # 2. try load solution data dictionary
    solutioDataFile = open(simCasePathPickle,"rb")
    solutionData = cPickle.load(solutioDataFile)
    solutioDataFile.close()    
    # 3. update vascular network with solution data
    vascularNetwork.prepareSolutionDataAfterLoad(solutionData)
    
    # 4. return vascular network  
    return vascularNetwork


def parseDirectoryForSimulationCases(networkName):
    '''
    Function to search for all simulationCases of a given network
    Input:
        networkName <String> (= name of the network to process) 
    Output:
        simulationCases <Dict> = { dataNumber1<String:SimCase1<String>, ...
                                   dataNumberN<String:SimCaseN<String>}
    '''
    path = ''.join([cur,'/../','NetworkFiles/',networkName,'/','SolutionData','/'])
    simulationCases = {}
    
    for dirName, dirNames, fileNames in os.walk(path):
        for fileName in fileNames:
            if ".pickle" in fileName and "polyChaos" not in fileName and "pure" not in fileName:
                dataNumber = fileName.split('.')[0].split('_SolutionData_')[-1]
                simulationCases[dataNumber] = ''.join([path,fileName])
                
    return simulationCases

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
    simCaseDescFilePath = ''.join([cur,'/../','NetworkFiles/',networkName,'/simulationCaseDescriptions.txt'])
    try:
        simCaseDescFile = open(simCaseDescFilePath, 'r')
    except:
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
