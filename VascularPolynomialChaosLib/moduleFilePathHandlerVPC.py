#!/usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
#
# moduleFilePathHandler VPC
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
sys.path.append(cur+'/../')

from copy import copy as copy 
from pprint import pprint as pp
import UtilityLib.moduleFilePathHandler as mFPH

from classConfigurationUQSA import ConfigurationUQSA

def getFilePath(fileType, networkName, dataNumber, mode, gPCEmethod = "None", gPCEorder = "None", evaluationNumber = "None", exception = 'Error'):
    '''
    Function return a requested file path, if this file exists
    
    Args:
        fileType (str): 'vpcConfigXmlFile',
                        'vpcConfigTemplateFile',
                        'vpcNetworkXmlFile',
                        'vpcSampleFile',
                        'vpcEvaluationNetworkXmlFile',
                        'vpcEvaluationSolutionDataFile',
                        'vpcProcessedSolutionDataFile',
                        'evaluationLogFile',
                        'vpcSolutionDataFile'
        
        networkName (str):name of the network file
        dataNumber (str):data number of solution file or xml file
        mode (str): read or write
        exception (str): (for read mode),
                    Error (default): raise error and exit if file is not exiting,
                    Warning: just raise Warning and return with error string,  
    '''
    existingFileTypes = ['vpcConfigXmlFile',
                         'vpcConfigTemplateFile',
                         'vpcNetworkXmlFile',
                         'vpcSampleFile',
                         'vpcEvaluationNetworkXmlFile',
                         'vpcEvaluationSolutionDataFile',
                         'vpcProcessedSolutionDataFile',
                         'evaluationLogFile',
                         'vpcSolutionDataFile']
    
    if fileType not in existingFileTypes:
        raise ValueError("ERROR: getFilePath, requested file type {}\
                          is not in existingFileTypess {}".format(fileType, existingFileTypes))
        
    # file names
    filenames = {
                 'vpcConfigTemplateFile'        : 'template_vpcConfig_xxx.xml',
                 'vpcConfigXmlFile'             : ''.join([networkName,'_vpcConfig_',dataNumber,'.xml']),
                 'vpcNetworkXmlFile'            : ''.join([networkName,'_vpc_',dataNumber,'.xml']),
                 'vpcSampleFile'                : ''.join(['samples.hdf5']),
                 'vpcEvaluationNetworkXmlFile'  : ''.join([networkName,'_evaluation_',str(evaluationNumber).zfill(7),'.xml']),
                 'vpcEvaluationSolutionDataFile': ''.join([networkName,'_evaluation_',str(evaluationNumber).zfill(7),'.hdf5']),
                 'evaluationLogFile'            : ''.join(['evaluationLogFile.txt']),
                 'vpcSolutionDataFile'          : ''.join([networkName,'_vpc-SolutionData_',dataNumber,'.hdf5'])
                 }    
        
    ## find requested file name
    requestedFilename  = filenames[''.join([fileType])]
    ## find directory    
    requestedDirectory = getDirectory(''.join([fileType,'Directory']), 
                                      networkName, 
                                      dataNumber,
                                      mode,
                                      gPCEmethod = gPCEmethod, 
                                      gPCEorder= gPCEorder,
                                      exception = exception)
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
    
def getDirectory(directoryType, networkName, dataNumber, mode, exception = 'Error', gPCEmethod = 'None', gPCEorder = 'None'):
    '''
    Function returns a requested directory path, if this directory does not exists
    it is created.
    
    Args:
        directoryType (str):  'workingDirectory',
                              'vpcConfigXmlFileDirectory',
                              'vpcNetworkXmlFileDirectory',
                              'vpcEvaluationNetworkXmlFileDirectory',
                              'vpcEvaluationSolutionDataFileDirectory',
                              'vpcSampleFileDirectory',
                              'evaluationLogFileDirectory',
                              'vpcSolutionDataFileDirectory',
                              'vpcConfigTemplateFileDirectory'
        networkName (str): name of the network file
        dataNumber (str): data number of solution file or xml file
        mode (str): read or write
        exception (str): (for read mode) 
                    Error (default): raise error and exit if file is not exiting-
                    Warning: just raise Warning and return with error string      
    '''
    
    existingDirectoryTypes = {'workingDirectory',
                              'vpcConfigXmlFileDirectory',
                              'vpcNetworkXmlFileDirectory',
                              'vpcEvaluationNetworkXmlFileDirectory',
                              'vpcEvaluationSolutionDataFileDirectory',
                              'vpcSampleFileDirectory',
                              'evaluationLogFileDirectory',
                              'vpcSolutionDataFileDirectory',
                              'vpcConfigTemplateFileDirectory'} 
    
    if directoryType not in existingDirectoryTypes:
        raise ValueError("ERROR: getDirectory, requested directoryType {}\
                          is not in existingDirectoryTypes{}".format(directoryType, existingDirectoryTypes))
    ##definitions
    starfishHomeDirectory       = ''.join([cur,'/..'])
    workingDirectory            = mFPH.readConfigFile(['WorkingDirectory'])['WorkingDirectory']
    networkXmlFileDirectory     = ''.join([workingDirectory,'/',networkName])
    vpcCaseDirectory            = ''.join([networkXmlFileDirectory,'/vascularPolynomialChaos_',str(dataNumber)]) 
    vpcOrderMethodDirectory     = ''.join([vpcCaseDirectory,'/','method_',gPCEmethod,'_order_',str(gPCEorder).zfill(2)])
    vpcEvaluationNetDirectory   = ''.join([vpcOrderMethodDirectory,'/evaluationNetworkFiles'])
    vpcEvaluationSolDirectory   = ''.join([vpcOrderMethodDirectory,'/evaluationSolutionData'])
    vpcConficTemplateDicrectory = ''.join([starfishHomeDirectory,'/TemplateNetworks','/vascularPolynomialChaos'])
    
    ## look up tables
    # directories
    directories = {
                   'workingDirectory'                      : workingDirectory,
                   # vascular polynomial chaos
                   'vpcConfigXmlFileDirectory'             : vpcCaseDirectory,
                   'vpcNetworkXmlFileDirectory'            : vpcCaseDirectory,
                   'vpcSampleFileDirectory'                : vpcOrderMethodDirectory,
                   'vpcEvaluationNetworkXmlFileDirectory'  : vpcEvaluationNetDirectory,
                   'vpcEvaluationSolutionDataFileDirectory': vpcEvaluationSolDirectory,
                   'evaluationLogFileDirectory'            : vpcOrderMethodDirectory,
                   'vpcSolutionDataFileDirectory'          : vpcOrderMethodDirectory,
                   'vpcConfigTemplateFileDirectory'        : vpcConficTemplateDicrectory,
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



def createConfigurationUQSAXMLfromTemplate(networkName, dataNumber):
    '''
    Function creates configuration file for UQSA for the case defined by given networkName and dataNumber from a template file
    
    Args:
        networkName (str): name of the network
        dataNumber (str): data number of case
    '''
    configurationFilePathTemplate = getFilePath('vpcConfigTemplateFile', networkName, dataNumber, 'read')
    configurationUQSA = ConfigurationUQSA()
    configurationUQSA.loadXMLFile(configurationFilePathTemplate)
    configurationFilePath = getFilePath('vpcConfigXmlFile', networkName, dataNumber, 'write')
    configurationUQSA.writeXMLFile(configurationFilePath)
        
def loadConfigurationUQSAFromXMLFile(networkName, dataNumber):
    '''
    Function loads configuration file for UQSA for the case defined by given networkName and dataNumber from a template file and returns a ConfigurationUQSA class instance
    
    Args:
        networkName (str): name of the network
        dataNumber (str): data number of case
    
    Returns: 
        configurationUQSA (object): ConfigurationUQSA class instance with loaded data
    '''
    configurationFilePath = getFilePath('vpcConfigXmlFile', networkName, dataNumber, 'read')
    configurationUQSA = ConfigurationUQSA()
    configurationUQSA.loadXMLFile(configurationFilePath)
    return configurationUQSA
    
def writeConfigutationUQSAToXMLFile(networkName, dataNumber, configurationUQSA):
    '''
    Function writes a configuration file for a given configurationUQSA class for the case defined by given networkName and dataNumber from a template file
    
    Args:
        networkName (str): name of the network
        dataNumber (str): data number of case
        configurationUQSA (object): ConfigurationUQSA class instance with data to write
    '''
    configurationFilePath = getFilePath('vpcConfigXmlFile', networkName, dataNumber, 'write')
    configurationUQSA.writeXMLFile(configurationFilePath)
