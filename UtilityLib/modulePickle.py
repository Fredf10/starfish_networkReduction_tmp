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

from copy import copy as copy 

from pprint import pprint as pp

def saveSolutionDataFile(networkName,dataNumber,vascularNetwork,solutionData,simulationCaseDescription):
    '''
    Function to save solutionData after a simulation is run or the solutionData was manipulated
    
    Input: 
        - networkName <String> (= name of the network to process)
        - dataNumbers <String> (= dataNumber with three digits as string)
        - solutionData <List> with SolutionData <Dict>
            [ {'Pressure': P, 'Flow': Q, 'Area': A, 'Name': name}, ... ]
        - simulationCaseDescription <String>
        
    Output: 
        None
    '''
    
    vascularNetwork.prepareToSave()
        
    #standart solution path
    solutionDirectory = ''.join([cur,'/../','NetworkFiles/',networkName,'/SolutionData/'])
    if not os.path.exists(solutionDirectory):
        os.makedirs(solutionDirectory)  
    pathSolutionData = ''.join([solutionDirectory,networkName,'_SolutionData_',dataNumber,'.pickle'])
    
    FILE = open(pathSolutionData,"w")
    dataToSave = [vascularNetwork,solutionData,simulationCaseDescription]
    cPickle.dump(dataToSave, FILE, protocol=2)
    FILE.close()
    
    pathSolutionData = ''.join([solutionDirectory,networkName,'_SolutionData_',dataNumber,'_pure.pickle'])
    FILE = open(pathSolutionData,"w")
    dataToSave = solutionData
    cPickle.dump(dataToSave, FILE, protocol=2)
    FILE.close()

def loadSolutionDataFile(networkName, dataNumbersInput):
    '''
    Function to open SolutionDataFiles for the given Network with the given dataNumbers
    
    Input:
        - networkName <String> (= name of the network to process)
        - dataNumbers <List> ['dataNumber', ... , 'dataNumber']
    
    if dataNumbers == ['all'] all solutionData is loaded and all dataSets are opened        
    
    Output:
        - vascularNetwork <vascularNetwork>
        - solutionDataSets <List> with SolutionData <Dict>
            [ {'Pressure': P, 'Flow': Q, 'Area': A, 'Name': name}, ... ]
        - simulationCaseDescriptions <nested list> with [[Description1<String>,dataNumber1<String], ...
                                                       [DescriptionN<String>,dataNumberN<String]]
    '''
    
    simulationCasesAvailable = parseDirectoryForSimulationCases(networkName)
    
    
    vascularNetwork = None
    solutionDataSets = []
    simulationCaseDescriptions =  []
    
    vascularNetworks = []
    
    dataNumbers = copy(dataNumbersInput)
    
    if dataNumbers == ['all']: dataNumbers = simulationCasesAvailable.keys()
    
    for dataNumber in dataNumbers:
        simCasePathAv = simulationCasesAvailable[dataNumber]
        #try:
        FILE = open(simCasePathAv,"rb")
        loadedData = cPickle.load(FILE)
        if vascularNetwork == None: 
            
            vascularNetwork = loadedData[0]
            vascularNetwork.quiet = True
            try:
                vascularNetwork.initialize()
            except:
                try:
                    ##
                    print "WARNING: some <instances> of the vascularNetwork are out of date.. try to update"
                    ## update vessel variables
                    vesselDataDict = {} 
                    
                    # get vessel data
                    for Id_i,vessel_i in vascularNetwork.vessels.iteritems():
                        vesselDataDict[Id_i] = copy(vessel_i.__dict__)
                                            
                    # create new data
                    for Id_i,data_i in vesselDataDict.iteritems():
                        vascularNetwork.deleteVessel(Id_i)
                        #print type(data_i['As'])
                        if type(data_i['As']) is not float and data_i['As'] is not None:
                            data_i['As'] = data_i['As'][0]
                        if type(data_i['beta']) is not float and data_i['beta'] is not None:
                            data_i['beta'] = data_i['beta'][0]
                        
                        vascularNetwork.addVessel(vesselId = Id_i, dataDict = data_i)
                    vascularNetwork.quiet = False
                except:
                    print "ERROR: could not update <instances> -> modulePickle could initialize network, exit()"
                    exit()
            
            vascularNetworks.append(vascularNetwork)
        else: loadedData[0].quiet = True
        solutionDataLoad = loadedData[1]
        for sol in solutionDataLoad:
            solutionDataSets.append(sol)
        try:
            simulationCaseDescriptions.append([loadedData[2],dataNumber])
        except:
            print "ERROR modulePickle: no Simulation-Case description for solution file {} - {} description to <<not available>>".format(dataNumber,networkName)
            simulationCaseDescriptions.append(['not available',dataNumber])
        FILE.close() 
#         except:
#             print "ERROR modulePickle: no solution file with given data number {} and networkName {} or it is corrupted! system exit".format(dataNumber,networkName)
#             exit()
    
    return vascularNetwork, solutionDataSets, simulationCaseDescriptions


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

def updateSimulationDescriptions(networkName, dataNumberInput='all'):
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
        
    vascularNetwork, solutionDataSets, simulationCaseDescriptions = loadSolutionDataFile(networkName,[dataNumberInput])
    
    alreadyWritten = False
    # check if datanumber in file:
    simCaseDescriptionFileLines = simCaseDescFile.readlines()
    simCaseDescFile.close() 
    
    currentDescription = simulationCaseDescriptions[0][0]
    currentDataNumber  = simulationCaseDescriptions[0][1]
    for Line in simCaseDescriptionFileLines:
        if currentDataNumber in Line: 
            if currentDescription not in Line or currentDescription == '':
                index = simCaseDescriptionFileLines.index(Line)
                simCaseDescriptionFileLines.remove(Line)
                simCaseDescriptionFileLines.insert(index,"  {} {} \n".format(dataNumberInput.ljust(10),simulationCaseDescriptions[0][0]))
            alreadyWritten = True
            
    if alreadyWritten == False:
        simCaseDescriptionFileLines.append("  {} {} \n".format(dataNumberInput.ljust(10),simulationCaseDescriptions[0][0]))
    
    #print simCaseDescriptionFileLines
    
    simCaseDescFile = open(simCaseDescFilePath, 'w')
    simCaseDescFile.writelines(simCaseDescriptionFileLines)
    simCaseDescFile.close()
    
    if vascularNetwork: vascularNetwork.quiet = True
        
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
