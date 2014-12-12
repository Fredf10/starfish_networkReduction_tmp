
try:
    from lxml import etree
except:
    from xml.etree import ElementTree as etree

import os,sys

#from pprint import pprint as pp


# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../'+'/NetworkLib')

from classVascularNetwork import VascularNetwork
from classBoundaryConditions import *

from constants import variablesDict 
from constants import unitsDictSI as unitsDict
from constants import newestNetworkXmlVersion
### import units of all variales in the Medical System
#from constants import variableUnitsMed as variableUnits
#from constants import unitsDictMed as unitsDict


def writeXMLsetUnit(xmlElement, variable, unit = 'unitSI'):
    '''
    Checks if element has a unit and adds it in the XML file
    '''
    try: 
        if variablesDict[variable][unit]: xmlElement.set('unit', variablesDict[variable][unit])
    except: print """ERROR: moduleXML.writeXML():
            variable {} of element {} is not proper defined
            in variablesDict, system exit()""".format(variable,xmlElement),exit()

def writeXMLsaveValues(xmlElement, variable, variableValues, polychaos = False):
    '''
    Writes the variable values to the xml file with xmlElement.text() function
    if variable can have multiple values they are saved with space as delimiter
    '''
    if variablesDict[variable]['multiVar'] or polychaos:
        xmlElement.text = ' '.join(str(i) for i in variableValues)
    else:
        xmlElement.text = str(variableValues)


def writeNetworkToXML(vascularNetwork, dataNumber = "xxx", filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    '''
    This function creates an XML file and writes all variable data of a vascularNetwork into it (except solution)
    The forma of the XML and all variable data are defined in constants.py
    '''
        
    #print filename
    if filename == None:
        filename = "newNetwork.xml"
        
    if dataNumber is not 'xxx':
        filename = filename.split('.')[0]
        networkDirectory = filename
        filename = ''.join([filename,'_SolutionData_',dataNumber,'.xml'])
    else: 
        networkDirectory = filename.split('.')[0]
    
    if networkDirectory not in networkPath.split('/'):
        folderPath = ''.join([networkPath,networkDirectory])
        if not os.path.exists(folderPath):
            os.makedirs(folderPath)  
    else:
        folderPath = networkPath
    
    networkXmlFile = ''.join([folderPath,'/',filename])
        
    try:
        root = etree.Element(filename, id = dataNumber, version = newestNetworkXmlVersion)
    except:
        print " Error: path / file does not exist"
        return
        
    ## import current network xml description as nxmlW(rite) to avoid version clash
    from constants import newestNetworkXml as nxmlW
        
    xmlFile = etree.ElementTree(root)
    
    for xmlElementName in nxmlW.xmlElements:
        xmlFileElement = etree.SubElement(root, xmlElementName)
        xmlElement = nxmlW.xmlElementsReference[xmlElementName]
        
        if xmlElementName == 'boundaryConditions':
            for vesselId,boundaryConditions in vascularNetwork.boundaryConditions.iteritems():
                subElement = etree.SubElement(xmlFileElement, 'boundaryCondition', vesselId = str(vesselId))
                try: polyChaosList = vascularNetwork.boundaryConditionPolyChaos[vesselId]
                except: pass
                # loop saved condition instances
                for boundaryCondition in boundaryConditions:
                    boundaryType = boundaryCondition.getVariableValue('name')
                    typeElement = etree.SubElement(subElement, boundaryType)
                    try: polyChaosDict = [polyDict for polyDict in polyChaosList if polyDict['name'] == boundaryType][0]
                    except: polyChaosDict = {}
                    # loop variables of the instance to be saved in xml
                    for variable in nxmlW.boundaryConditionElements[boundaryType]:
                        variableElement = etree.SubElement(typeElement, variable)
                        writeXMLsetUnit(variableElement,variable)
                        writeXMLsaveValues(variableElement,variable,boundaryCondition.getVariableValue(variable))  
                        if variable in polyChaosDict.keys():
                            subElementIntervalValues = polyChaosDict[variable]
                            intervalElement = etree.SubElement(typeElement, "-".join([variable,'polyChaos']))
                            writeXMLsetUnit(intervalElement,variable)
                            writeXMLsaveValues(intervalElement,variable,subElementIntervalValues,polychaos = True)  
                            
        elif xmlElementName == 'vessels':
            for vessel in vascularNetwork.vessels.itervalues():
                attributes = {}
                for attribute in nxmlW.vesselAttributes:
                    attributes[attribute] = str(vessel.getVariableValue(attribute))
                vesselSubElement = etree.SubElement(xmlFileElement,'vessel',attributes)
                for vesselElement in xmlElement:
                    # add vesselElement
                    subElement = etree.SubElement(vesselSubElement, vesselElement)
                    # check if compliance and adjust variables
                    if vesselElement == 'compliance':
                        complianceType = vessel.getVariableValue('complianceType')
                        variables = nxmlW.vesselElementReference[vesselElement][complianceType]
                    else:
                        variables = nxmlW.vesselElementReference[vesselElement]
                    # save variables
                    for variable in variables:
                        subsubElement = etree.SubElement(subElement, variable)
                        variableValues =  vessel.getVariableValue(variable)
                        writeXMLsetUnit(subsubElement,variable)
                        writeXMLsaveValues(subsubElement,variable,variableValues)  
                        # check for polychaos data of vessel
                        try: 
                            subElementIntervalValues = vessel.polyChaos[variable]
                            intervalElement = etree.SubElement(subElement, "-".join([variable,'polyChaos']))
                            writeXMLsetUnit(intervalElement,variable)
                            writeXMLsaveValues(intervalElement,variable,subElementIntervalValues,polychaos = True) 
                        except : pass
        
        elif xmlElementName == 'communicators':
            for comId,comData in vascularNetwork.communicators.iteritems():
                comType = comData['comType']
                subElement = etree.SubElement(xmlFileElement, comType)
                for variable in nxmlW.communicatorReference[comData['comType']]:
                    subsubElement = etree.SubElement(subElement, variable)
                    writeXMLsetUnit(subsubElement,variable)
                    writeXMLsaveValues(subsubElement,variable,comData[variable])  
                    
        elif xmlElementName == 'baroreceptors':
            for baroId,baroData in vascularNetwork.baroreceptors.iteritems():
                baroType = baroData['baroType']
                subElement = etree.SubElement(xmlFileElement, baroType)
                for variable in nxmlW.baroreceptorReference[baroData['baroType']]:
                    subsubElement = etree.SubElement(subElement, variable)
                    writeXMLsetUnit(subsubElement,variable)
                    writeXMLsaveValues(subsubElement,variable,baroData[variable])  
                
        else: # vascularNetwork
            for variable in xmlElement:
                subElement = etree.SubElement(xmlFileElement, variable) # add subElement
                writeXMLsetUnit(subElement,variable)  # check unit
                if xmlElementName == 'globalFluid': # get variable values from varscularNetwork     
                    variableValues = vascularNetwork.globalFluid[variable]
                else: 
                    variableValues = vascularNetwork.getVariableValue(variable)
                writeXMLsaveValues(subElement,variable,variableValues)  
                # check for interaval data of global fluid
                try: 
                    subElementIntervalValues = vascularNetwork.globalFluidPolyChaos[variable]
                    intervalElement = etree.SubElement(xmlFileElement, "-".join([variable,'polyChaos']))
                    writeXMLsetUnit(intervalElement,variable)
                    intervalElement.text = ' '.join(str(i) for i in subElementIntervalValues)
                except : pass   
                                
    
    xmlFile.write(networkXmlFile,encoding='iso-8859-1',pretty_print = True)
    

def loadVariablesConversion(variable, variableValueStr, variableUnit, unit = 'unitSI', polychaos = False):
    '''
    checks the element.text string and
    evaluates it corresponding to the definition
    of the variable in the constants variablesDict 
    
    return converted evaluated value of variable
    '''
    
    polynomialChaosDistributions = ['Uniform','Normal'] # need to be moved and changed to networkXML files
    
    multiVariable = False
    variableValue = 'notConvertable'
    convertError = []
    variableTypes = variablesDict[variable]['type']
    # check if variable is a multiple variable (means list of variables)
    if variablesDict[variable]['multiVar'] or polychaos == True:
        variableValueStrings = variableValueStr.split(' ')
        multiVariable = True
        variableValues = []
        if polychaos == True: variableTypes = ' '.join(['str',variablesDict[variable]['type']])    
    else:
        variableValueStrings = [variableValueStr]
    
    # check if variable can have multiple types
    if ' ' in variableTypes:
        variableTypes = variableTypes.split(' ')
    else: variableTypes = [variableTypes]
    variableTypes.sort()
    
    # start conversion loop over variable and types
    for variableValueString in variableValueStrings:
        for variableType in variableTypes:
                        
            if variableType in ['float','int']:
                try: variableValue = float(eval(variableValueString))
                except: convertError.append('float') 
                
                if variablesDict[variable][unit]:
                    try:
                        if ' ' in variableUnit:
                            variableUnits = variableUnit.split(' ')
                            for variableUnit in variableUnits: variableValue = variableValue*unitsDict[variableUnit]
                        else: variableValue = variableValue*unitsDict[variableUnit]
                    except: pass 
                    
                if variableType == 'int': 
                    try: variableValue = int(variableValueString)
                    except: convertError.append('int') 
                    
            elif variableType == 'bool': 
                try: 
                    if variableValueString == 'False': variableValue = False
                    elif variableValueString == 'True': variableValue = True
                    else: variableValue = eval(variableValueString)
                except: 
                    convertError.append('bool') 
                
            elif variableType in ['str']:
                if polychaos: 
                    if variableValueString in polynomialChaosDistributions: variableValue = variableValueString
                elif variableValueString in variablesDict[variable]['strCases']: variableValue = variableValueString
                elif variablesDict[variable]['strCases'][0] == 'anything': variableValue = variableValueString
                else: convertError.append(''.join(['str == ',str(variablesDict[variable]['strCases'])]))
            
            elif variableType in ['None']:
                if variableValueString == 'None' or variableValueString == '' or variableValueString == None:
                    variableValue = None
                else: convertError.append('None')
                    
        if variableValue == 'notConvertable':
            print """ERROR: moduleXML.loadVariablesConversion():
                  Cannot convert given value "{}" of variable "{}"
                  to {}!
                  Check if it is of type {}, system exit!""".format(variableValueString,variable,convertError,variableTypes); exit()
        if multiVariable == False: return variableValue
        else: variableValues.append(variableValue)
    
    return variableValues
   
def loadingErrorMessageValueError(variableName, element, elementName):
    variableDict =   variablesDict[variableName]
    print """ERROR loadNetworkFromXML():
          value for variable <<{}>> of {} {} is not defined.
          (Hint:{}) , system exit!""".format(variableName, element, elementName, variableDict)
    exit()
   
def loadingErrorMessageVariableError(variableName, element, elementName):
    variableDict =   variablesDict[variableName]
    print """ERROR loadNetworkFromXML():
          variable "{}" of {} {} is not defined.
          (Hint:{}) , system exit!""".format(variableName, element, elementName, variableDict)
    exit()
      
   
def loadNetworkFromXML(filename = None, dataNumber = "xxx", networkPath = str(cur+"/../NetworkFiles/")):
    '''
    Function laods network from XML-file
    
    version of XML files supported: 4.0, 4.1
    '''
    currentVersions = ['4.0','4.1']
    
    # read from file
    if filename == None:
        print 'ERROR: moduleXML.loadNetworkFromXML() : load XML - no filename passed'
        return None
    
    networkDirectory = filename.split('.')[0]
    
    if dataNumber is not 'xxx':        
        filename = ''.join([filename.split('.')[0],'_SolutionData_',dataNumber,'.xml'])
    
    if '.xml' not in filename:
        filename = ''.join([filename,'.xml'])  
        
    if networkDirectory not in networkPath.split('/'):
        folderPath = ''.join([networkPath,networkDirectory])
    else:
        folderPath = networkPath
    
    
      
    if not os.path.exists(folderPath):
        print 'ERROR: moduleXML.loadNetworkFromXML(): directory and file does not exists'
        return None
    
    networkXmlFile = ''.join([folderPath,'/',filename])
    
    # create vascularNetwork instance
    vascularNetwork = VascularNetwork()
    # set name
    vascularNetwork.name = networkDirectory
    
    try:
        parser = etree.XMLParser(encoding='iso-8859-1')
        tree = etree.parse(''.join([networkXmlFile]), parser)
    except (etree.ParseError, ImportError) as e:
        if isinstance(e, etree.ParseError):
            print " ERROR moduleXML.loadNetworkFromXML() on line {} {}: ".format(e.position[0], e)
            exit()
      
    # create root
    root = tree.getroot()
    xmlFileVersion = root.attrib['version']
    if xmlFileVersion not in currentVersions:
        print "ERROR moduleXML.loadNetworkFromXML(): XML file is outdated file-version {} " \
        "current supported version {}, could not parse file! system exit".format(root.attrib['version'],currentVersions); exit()
    
    
    if xmlFileVersion == newestNetworkXmlVersion:
        from constants import newestNetworkXml as nxml
    
    elif xmlFileVersion == '4.0':
        import networkXml040 as nxml  
    
    if xmlFileVersion != newestNetworkXmlVersion:
        print " WARNING the version of the network xml file you try to load is outdated it may occure some problems!"
    
    for xmlElementName in nxml.xmlElements:
        for xmlElement in root.findall(''.join([".//",xmlElementName])):
                      
            if xmlElementName == 'boundaryConditions':
                # loop through all boundaryCondition
                for boundaryConditionElement in xmlElement.findall(''.join(['.//','boundaryCondition'])):
                    try: vesselId = int(boundaryConditionElement.attrib['vesselId'])
                    except: loadingErrorMessageVariableError('vesselId', 'one boundaryCondition', '')
                       
                    boundaryInstances = []
                    boundaryIntervals = []
                    # loop through possible communicator class types
                    for boundaryType in nxml.xmlElementsReference[xmlElementName]:
                        # find all bcs of this type
                        for bcElements in boundaryConditionElement.findall(''.join(['.//',boundaryType])):
                            boundaryInstance = eval(nxml.bcTagsClassReferences[boundaryType])()
                            boundaryDataDict = {}
                            boundaryDataDict['name'] = boundaryType
                            boundaryIntervalDataDict = {}
                            # loop through all variables of this type, convert and save values of these
                            for variable in nxml.xmlElementsReference[xmlElementName][boundaryType]: 
                                # find normal variables
                                try: element = bcElements.findall(''.join(['.//',variable]))[0]
                                except: loadingErrorMessageVariableError(variable, 'boundaryCondition', boundaryType)
                                # get variable value                        
                                try: variableValueStr = element.text
                                except: loadingErrorMessageValueError(variable, 'boundaryCondition', boundaryType)
                                # get unit
                                try: variableUnit = element.attrib['unit']
                                except: variableUnit = None 
                                # save converted XML-value
                                boundaryDataDict[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                                # find polynomial chaos variable
                                try:                                                       
                                    # find polychaos variables
                                    element = bcElements.findall(''.join(['.//',variable,'-polyChaos']))[0]
                                    # get variable value                        
                                    try: variableValueStr = element.text
                                    except: loadingErrorMessageValueError(variable, 'boundaryCondition', boundaryType)
                                    # get unit
                                    try: variableUnit = element.attrib['unit']
                                    except: variableUnit = None 
                                    # save converted XML-value                      
                                    boundaryIntervalDataDict[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit, polychaos = True)
                                    boundaryIntervalDataDict['name'] = boundaryType  
                                except: pass
                                # adjust path to boundary condition file
                                if variable == 'filePathName':
                                    path = ''.join([networkDirectory,'/'])
                                    if path not in variableValueStr:  variableValueStr = variableValueStr.join([path,''])
                                    boundaryDataDict['filePathName'] = variableValueStr
                                    
                            boundaryInstance.update(boundaryDataDict)
                            boundaryInstances.append(boundaryInstance)             
                            if boundaryIntervalDataDict != {}:
                                boundaryIntervals.append(boundaryIntervalDataDict)
                                
                    # apply read data to vascularNetwork
                    if vesselId not in vascularNetwork.boundaryConditions.keys():
                        vascularNetwork.boundaryConditions[vesselId] = boundaryInstances
                    else:
                        vascularNetwork.boundaryConditions[vesselId].extend(boundaryInstances)
                    if boundaryIntervals != []: 
                        vascularNetwork.boundaryConditionPolyChaos[vesselId] = boundaryIntervals
                        
            elif xmlElementName == 'vessels':
                for vesselXMLnode in xmlElement.findall(''.join(['.//','vessel'])):
                    vesselData = {}
                    vesselPolychaosData = {}
                    # load vessel attributes
                    for attribute in nxml.vesselAttributes:
                        try: vesselData[attribute] = loadVariablesConversion(attribute, vesselXMLnode.attrib[attribute], '')
                        except: 
                            try:    loadingErrorMessageVariableError(attribute, 'vessel', vesselData['Id'])
                            except: loadingErrorMessageVariableError(attribute, 'one vessel', '')
                        
                        
                    for vesselElement in nxml.xmlElementsReference[xmlElementName]: 
                        # check if compliance and adjust variables
                        if vesselElement == 'compliance':
                            try: complianceTypeElement = vesselXMLnode.findall(''.join(['.//','complianceType']))[0]
                            except: loadingErrorMessageVariableError('complianceType', 'vessel', vesselData['Id'])
                            complianceType = loadVariablesConversion('complianceType', complianceTypeElement.text, '')
                            variables = nxml.vesselElementReference[vesselElement][complianceType]
                        else:
                            variables = nxml.vesselElementReference[vesselElement]
                        # load variables
                        for variable in variables:
                            try: element = vesselXMLnode.findall(''.join(['.//',variable]))[0]
                            except: loadingErrorMessageVariableError(variable, 'vessel', vesselData['Id'])
                            # get variable value                        
                            try: variableValueStr = element.text
                            except: loadingErrorMessageValueError(variable, 'vessel', vesselData['Id'])
                            # get unit 
                            try: variableUnit = element.attrib['unit']
                            except: variableUnit = None 
                            # save converted XML-value
                            vesselData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                            # find polynomial chaos variable
                            try:                        
                                # find polychaos variables
                                element = vesselXMLnode.findall(''.join(['.//',variable,'-polyChaos']))[0]
                                # get variable value                        
                                try: variableValueStr = element.text
                                except: loadingErrorMessageValueError(variable, 'vessel', vesselData['Id'])
                                # get unit
                                try: variableUnit = element.attrib['unit']
                                except: variableUnit = None 
                                # save converted XML-value                      
                                vesselPolychaosData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit, polychaos = True)
                            except: pass
                        vesselData['polyChaos'] = vesselPolychaosData
                    vascularNetwork.updateNetwork({'vesselData':{vesselData['Id']:vesselData}})
            
            elif xmlElementName == 'communicators':
                # loop through possible communicator class types
                for comunicatorType in nxml.xmlElementsReference[xmlElementName]:
                    # find all communicator of this type
                    for comElements in xmlElement.findall(''.join(['.//',comunicatorType])):
                        # loop through all variables of this type, convert and save values of these
                        communicatorData = {}
                        for variable in nxml.xmlElementsReference[xmlElementName][comunicatorType]: 
                            try: element = comElements.findall(''.join(['.//',variable]))[0]
                            except: loadingErrorMessageVariableError(variable, 'communicator', comunicatorType)
                            # get variable value                        
                            try: variableValueStr = element.text
                            except: loadingErrorMessageValueError(variable, 'communicator', comunicatorType)
                            # get unit
                            try: variableUnit = element.attrib['unit']
                            except: variableUnit = None 
                            communicatorData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                            if variable == 'comId': comId = communicatorData[variable]
                        vascularNetwork.updateNetwork({'communicators':{comId:communicatorData}})
                        
            elif xmlElementName == 'baroreceptors':
                # loop through possible communicator class types
                for baroreceptorType in nxml.xmlElementsReference[xmlElementName]:
                    # find all communicator of this type
                    for baroElements in xmlElement.findall(''.join(['.//',baroreceptorType])):
                        # loop through all variables of this type, convert and save values of these
                        baroreceptorData = {}
                        for variable in nxml.xmlElementsReference[xmlElementName][baroreceptorType]: 
                            try: element = baroElements.findall(''.join(['.//',variable]))[0]
                            except: loadingErrorMessageVariableError(variable, 'baroreceptors', baroreceptorType)
                            # get variable value                        
                            try: variableValueStr = element.text
                            except: loadingErrorMessageValueError(variable, 'baroreceptors', baroreceptorType)
                            # get unit
                            try: variableUnit = element.attrib['unit']
                            except: variableUnit = None 
                            baroreceptorData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                            if variable == 'baroId': baroId = baroreceptorData[variable]
                        baroreceptorData['baroType'] = baroreceptorType
                        vascularNetwork.updateNetwork({'baroreceptors':{baroId:baroreceptorData}})
                        
            elif xmlElementName == 'globalFluid':
                    
                globalFluidData = {}
                globalFluidPolychaosData = {}
                
                for variable in nxml.xmlElementsReference[xmlElementName]: 
                    # find normal variables
                    try: element = xmlElement.findall(''.join(['.//',variable]))[0]
                    except: loadingErrorMessageVariableError(variable, 'global fluid', '')
                    # get variable value                        
                    try: variableValueStr = element.text
                    except: loadingErrorMessageValueError(variable, 'global fluid', '')
                    # get unit
                    try: variableUnit = element.attrib['unit']
                    except: variableUnit = None 
                    # save converted XML-value
                    globalFluidData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit) 
                    
                    # find polynomial chaos variable
                    try:                        
                        # find polychaos variables
                        element = xmlElement.findall(''.join(['.//',variable,'-polyChaos']))[0]
                        # get variable value                        
                        try: variableValueStr = element.text
                        except: loadingErrorMessageValueError(variable, 'global fluid', '')
                        # get unit
                        try: variableUnit = element.attrib['unit']
                        except: variableUnit = None 
                        # save converted XML-value                      
                        globalFluidPolychaosData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit, polychaos = True)
                    except: pass
                vascularNetwork.updateNetwork({'globalFluid':globalFluidData,'globalFluidPolyChaos':globalFluidPolychaosData})
                    
            elif xmlElementName in nxml.vascularNetworkElements: # vascularNetwork
                vascularNetworkData = {}
                for variable in nxml.xmlElementsReference[xmlElementName]: 
                    try: element = xmlElement.findall(''.join(['.//',variable]))[0]
                    except: loadingErrorMessageVariableError(variable, xmlElementName, '')
                    # get variable value                        
                    try: variableValueStr = element.text
                    except: loadingErrorMessageValueError(variable, xmlElementName, '')
                    # get 
                    try: variableUnit = element.attrib['unit']
                    except: variableUnit = None 
                    # save converted XML-value
                    vascularNetworkData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                vascularNetwork.update(vascularNetworkData)
        
    return vascularNetwork

def savePolyChaosXML(vPCconfiguration, filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    '''
    Function to write a xml file with vascularPolynomialChaos Configurations
    '''
    
    if '_template' in filename:
        vPCconfiguration = vPCconfigurationTemplate
        filename = filename.rsplit('_template')[0]
    
    if '.xml' not in filename:
        filename = ''.join([filename,'.xml']) 
    
    networkDirectory = filename.split('_vpcConfig_')[0]
    if not os.path.exists(str(networkPath+networkDirectory)):
        os.makedirs(str(networkPath+networkDirectory))  
    
    pathAndFileName = ''.join([networkPath,networkDirectory,'/',filename])

    try:
        root = etree.Element(filename, id= "1.0", version="1.0")
    except:
        print " Error: path / file does not exist"
        return
    
    
    xmlFile = etree.ElementTree(root)
    #----------------------------------------------------------------------------------------# 
    ### Write Control Variables
    controlVariables = etree.SubElement(root, "controlVariables")
    
    preProcessing = etree.SubElement(controlVariables, "preProcessing")
    
    
    createOrthoPoly = etree.SubElement(preProcessing, "createDistributions")
    createOrthoPoly.text = str(vPCconfiguration['createDistributions'])
    createOrthoPoly = etree.SubElement(preProcessing, "createOrthoPoly")
    createOrthoPoly.text = str(vPCconfiguration['createOrthoPoly'])
    createSample = etree.SubElement(preProcessing, "createSample")
    createSample.text = str(vPCconfiguration['createSample'])
       
    runSimulations = etree.SubElement(controlVariables, "runSimulations")
    runSimulations.text = str(vPCconfiguration['runSimulations'])
    
    calculations = etree.SubElement(controlVariables, "calculations")
    
    calculateGPCE = etree.SubElement(calculations, "calculateGPCE")
    calculateGPCE.text = str(vPCconfiguration['calculateGPCE'])
    preProcessData = etree.SubElement(calculations, "preProcessData")
    preProcessData.text = str(vPCconfiguration['preProcessData'])
    plotMinMaxPoints = etree.SubElement(calculations, "plotMinMaxPoints")
    plotMinMaxPoints.text = str(vPCconfiguration['plotMinMaxPoints'])
    
    postProcessing = etree.SubElement(controlVariables, "postProcessing")
    #postProcessing.text = str(vPCconfiguration['postProcessing'])    
    plotMeanSTD = etree.SubElement(postProcessing, "plotMeanSTD")
    plotMeanSTD.text = str(vPCconfiguration['plotMeanSTD'])
    plotPeaks = etree.SubElement(postProcessing, "plotPeaks")
    plotPeaks.text = str(vPCconfiguration['plotPeaks'])

    #----------------------------------------------------------------------------------------#    
    #### POLYNOMIAL CHAOS DEFINITIONS    
    polyChaosDefinition = etree.SubElement(root, "polyChaosConfig")
        
    polynomialOrders = etree.SubElement(polyChaosDefinition, "polynomialOrders")
    polynomialOrders.text = ' '.join(str(i) for i in vPCconfiguration['polynomialOrders'])
    sampleMethod = etree.SubElement(polyChaosDefinition, "sampleMethod")
    sampleMethod.text = vPCconfiguration['sampleMethod'] #.replace('><',' ').strip("><")
    
    #----------------------------------------------------------------------------------------#
    #### WAVE SPLITTING
    waveSplitting = etree.SubElement(root, "waveSplitting")
    
    linearWaveSplit = etree.SubElement(waveSplitting, "linearWaveSplit")
    linearWaveSplit.text = str(vPCconfiguration['linearWaveSplit'])
    velocityProfileCoefficient = etree.SubElement(waveSplitting, "velocityProfileCoefficient")
    velocityProfileCoefficient.text = str(vPCconfiguration['velocityProfileCoefficient'])
    
    #---------------------------------------------------------------------------------------#
    ### POSTPROCESSING: PLOTTING
    
    postprocessing = etree.SubElement(root, "postProcessingConfig")
    
    generalPlotting = etree.SubElement(postprocessing, "generalPlotting")
    ##########################################################################################
    plotDirectory = etree.SubElement(generalPlotting, "plotDirectory")
    plotDirectory.text = str(vPCconfiguration['plotDirectory'])
    deterministicDataSetNumbers = etree.SubElement(generalPlotting, "deterministicDataSetNumbers")
    deterministicDataSetNumbers.text = ' '.join(str(i) for i in vPCconfiguration['deterministicDataSetNumbers'])
    polynomsToPlotOrder = etree.SubElement(generalPlotting, "polynomsToPlotOrder")
    polynomsToPlotOrder.text = ' '.join(str(i) for i in vPCconfiguration['polynomsToPlotOrder'])
    
    meanSTDplots = etree.SubElement(postprocessing, "plotsMeanSTD")
    ##########################################################################################
    plotMeanConfidenceInterval = etree.SubElement(meanSTDplots, "plotConfidenceInterval")
    plotMeanConfidenceInterval.text = str(vPCconfiguration['plotMeanConfidenceInterval'])
    plotMeanConfidenceAlpha = etree.SubElement(meanSTDplots, "confidenceAlpha")
    plotMeanConfidenceAlpha.text = str(vPCconfiguration['plotMeanConfidenceAlpha'])
    plotMeanSigmaInterval = etree.SubElement(meanSTDplots, "sigmaInterval")
    plotMeanSigmaInterval.text = str(vPCconfiguration['plotMeanSigmaInterval'])

    plotsPeaks = etree.SubElement(postprocessing, "plotsPeaks")
    ##########################################################################################
    peakAnalysis = etree.SubElement(plotsPeaks, "peakAnalysis")
    peakAnalysis.text = str(vPCconfiguration['peakAnalysis'])
    plotPeaksConfidenceAlpha = etree.SubElement(plotsPeaks, "plotPeaksConfidenceAlpha")
    plotPeaksConfidenceAlpha.text = str(vPCconfiguration['plotPeaksConfidenceAlpha'])
    plotPeaksAnalyticSensitivity = etree.SubElement(plotsPeaks, "plotPeaksAnalyticSensitivity")
    plotPeaksAnalyticSensitivity.text = str(vPCconfiguration['plotPeaksAnalyticSensitivity'])
    plotPeaksMeanSTDBoxPlotsSingle = etree.SubElement(plotsPeaks, "plotPeaksMeanSTDBoxPlotsSingle")
    plotPeaksMeanSTDBoxPlotsSingle.text = str(vPCconfiguration['plotPeaksMeanSTDBoxPlotsSingle'])

    
    #----------------------------------------------------------------------------------------#
    # INVESTIGATION POINTS
    evaluationPoints = etree.SubElement(root, "evaluationPoints")
    
    evalPointCount = 0
    for idNodeTuple in vPCconfiguration['polynomsToCalculate']:
        try:
            name = str(vPCconfiguration['names'][evalPointCount])
        except: "Error: no name for id,node {} tuple defined".format(idNodeTuple)
        evaluationPoint = etree.SubElement(evaluationPoints, "evaluationPoint",  id = str(idNodeTuple[0]), name = str(name), gridNode = str(idNodeTuple[1]))        
        try:
            deltas = etree.SubElement(evaluationPoint, "deltasMinMaxFunction")
            deltaPressure = etree.SubElement(deltas, "Pressure", unit = variableUnits['Pressure'])
            deltaPressure.text = str(vPCconfiguration['delta'][name]['Pressure'])
            deltaPressureForward = etree.SubElement(deltas, "Pressure_f", unit = variableUnits['Pressure'])
            deltaPressureForward.text = str(vPCconfiguration['delta'][name]['Pressure_f'])
            deltaPressureBackward = etree.SubElement(deltas, "Pressure_b", unit = variableUnits['Pressure'])
            deltaPressureBackward.text = str(vPCconfiguration['delta'][name]['Pressure_b'])
            deltaFlow = etree.SubElement(deltas, "Flow", unit = variableUnits['Flow'])
            deltaFlow.text = str(vPCconfiguration['delta'][name]['Flow'])
            deltaFlowForward = etree.SubElement(deltas, "Flow_f", unit = variableUnits['Flow'])
            deltaFlowForward.text = str(vPCconfiguration['delta'][name]['Flow_f'])
            deltaFlowBackward = etree.SubElement(deltas, "Flow_b", unit = variableUnits['Flow'])
            deltaFlowBackward.text = str(vPCconfiguration['delta'][name]['Flow_b'])
        except: "Error: in deltas for id,node {} tuple defined".format(idNodeTuple)
        try:
            peaksToEvaluate = etree.SubElement(evaluationPoint, "peaksToEvaluate")
            
            extremaPressure = etree.SubElement(peaksToEvaluate, "Pressure", unit = variableUnits['Pressure'])
            extremaPressure.text = ' '.join(str(i) for i in vPCconfiguration['peaksToEvaluate'][name]['extremaPressure'])
            extremaPressureBackward = etree.SubElement(peaksToEvaluate, "Pressure_f", unit = variableUnits['Pressure'])
            extremaPressureBackward.text = ' '.join(str(i) for i in vPCconfiguration['peaksToEvaluate'][name]['extremaPressure_f'])
            deltaPressureBackward = etree.SubElement(peaksToEvaluate, "Pressure_b", unit = variableUnits['Pressure'])
            deltaPressureBackward.text = ' '.join(str(i) for i in vPCconfiguration['peaksToEvaluate'][name]['extremaPressure_b'])
            deltaFlow = etree.SubElement(peaksToEvaluate, "Flow", unit = variableUnits['Flow'])
            deltaFlow.text = ' '.join(str(i) for i in vPCconfiguration['peaksToEvaluate'][name]['extremaFlow'])
            deltaFlowForward = etree.SubElement(peaksToEvaluate, "Flow_f", unit = variableUnits['Flow'])
            deltaFlowForward.text = ' '.join(str(i) for i in vPCconfiguration['peaksToEvaluate'][name]['extremaFlow_f'])
            deltaFlowBackward = etree.SubElement(peaksToEvaluate, "Flow_b", unit = variableUnits['Flow'])
            deltaFlowBackward.text = ' '.join(str(i) for i in vPCconfiguration['peaksToEvaluate'][name]['extremaFlow_b']) 
            
        except: "Error: in deltas for id,node {} tuple defined".format(idNodeTuple)
        evalPointCount = evalPointCount+1

    xmlFile.write(pathAndFileName,encoding='iso-8859-1',pretty_print = True)
    print " ... file saved"
  
  
###----------------------------------------------------------------------------------------
### Polynomial chaos
# Conversion
def unitConversion(unitsDict,value,unit):
    try: 
        value = float(value)
    except:
        if value is None: return value
        else: print 'ERROR: unitConversion not feasible: value {} not convertable to float'.format(value)
    if ' ' in unit:
        unit = unit.split(' ')
        for un in range (0,len(unit),1): value = value*unitsDict[unit[un]]
    else: 
        value = value*unitsDict[unit]
    return value
  
    
def loadPolyChaosXML(filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    
    
    ### import units of all variables in the SI system outdated used of only few functions -> to be changed
    ## just used for polynomial chaos config.
    from constants import variableUnitsSI as variableUnits
    from constants import vPCconfigurationTemplate
    
    if '.xml' not in filename:
        filename = ''.join([filename,'.xml'])
    
    vPCconfiguration =  {}
    
    networkDirectory = filename.split('_vpcConfig_')[0]
    pathAndFileName = ''.join([networkPath,networkDirectory,'/',filename])
    
    # load the data!
    try:
        parser = etree.XMLParser(encoding='iso-8859-1')
        tree = etree.parse(pathAndFileName, parser)
    except (etree.ParseError, ImportError) as e:
        if isinstance(e, etree.ParseError):
            print " Error in XML file on line: "
            print e.position[0]
            exit()     
            
    # create root
    root = tree.getroot()
    
    #----------------------------------------------------------------------------------------# 
    ### Read Control Variables
    controlVariables = ['createDistributions','createOrthoPoly','createSample','runSimulations','calculateGPCE','preProcessData','plotMinMaxPoints','plotMeanSTD','plotPeaks']
    for control in root.findall(".//controlVariables"):
        for data in control:
            if data.tag in controlVariables:
                vPCconfiguration[data.tag] = eval(data.text)
    for preProcessing in root.findall(".//preProcessing"):
        for data in preProcessing:
            if data.tag in controlVariables:
                vPCconfiguration[data.tag] = eval(data.text)
    for calculations in root.findall(".//calculations"):
        for data in calculations:
            if data.tag in controlVariables:
                vPCconfiguration[data.tag] = eval(data.text)
    for preProcessing in root.findall(".//postProcessing"):
        for data in preProcessing:
            if data.tag in controlVariables:
                vPCconfiguration[data.tag] = eval(data.text)
    
    if vPCconfiguration['plotPeaks'] or vPCconfiguration['plotMeanSTD'] is True: vPCconfiguration['postProcessing'] = True
    else: vPCconfiguration['postProcessing'] = False
    
    #----------------------------------------------------------------------------------------#    
    #### POLYNOMIAL CHAOS DEFINITIONS    
    for preProcessing in root.findall(".//polyChaosConfig"):
        for data in preProcessing:
            if data.tag == "polynomialOrders":
                vPCconfiguration['polynomialOrders'] = map(int, data.text.split(' '))
            if data.tag == 'sampleMethod':
                vPCconfiguration['sampleMethod'] = data.text#''.join(['<',data.text.replace(' ','><'),'>'])
                        
    
    #----------------------------------------------------------------------------------------#
    #### WAVE SPLITTING
    for waveSplitting in root.findall(".//waveSplitting"):
        for data in waveSplitting:
            if data.tag == 'linearWaveSplit':
                vPCconfiguration[data.tag] = eval(data.text)
            if data.tag == 'velocityProfileCoefficient':
                if data.text != 'None' and data.text !=  None: 
                    vPCconfiguration['velocityProfileCoefficient'] = float(data.text)
                    print "asd"

    #---------------------------------------------------------------------------------------#
    ### POSTPROCESSING: PLOTTING
    for preProcessing in root.findall(".//generalPlotting"):
        for data in preProcessing:
            if data.tag == 'plotDirectory':
                vPCconfiguration[data.tag] = data.text
            if data.tag == "deterministicDataSetNumbers":
                if data.text: vPCconfiguration[data.tag] = map(int, data.text.split(' '))
                else: vPCconfiguration[data.tag] = []
            if data.tag == "polynomsToPlotOrder":
                vPCconfiguration[data.tag] = map(int, data.text.split(' '))
    
    for plotsMeanSTD in root.findall(".//plotsMeanSTD"):
        for data in plotsMeanSTD:
            if data.tag == 'plotConfidenceInterval':
                vPCconfiguration['plotMeanConfidenceInterval'] = eval(data.text)
            if data.tag == 'confidenceAlpha':
                if data.text != 'None' and data.text !=  None: 
                    vPCconfiguration['plotMeanConfidenceAlpha'] = float(data.text)
            if data.tag == 'sigmaInterval':
                vPCconfiguration['plotMeanSigmaInterval'] = eval(data.text)
    
    for plotsPeaks in root.findall(".//plotsPeaks"):
        for data in plotsPeaks:
            if data.tag in ['peakAnalysis','plotPeaksMeanSTDBoxPlotsSingle','plotPeaksAnalyticSensitivity']:
                vPCconfiguration[data.tag] = eval(data.text)
                
            if data.tag == 'plotPeaksConfidenceAlpha':
                if data.text != 'None' and data.text !=  None: 
                    vPCconfiguration['plotPeaksConfidenceAlpha'] = float(data.text)
    
    #----------------------------------------------------------------------------------------#
    # INVESTIGATION POINTS
    polynomsToCalculate = []
    names = []
    delta = {}
    peaksToEvaluate = {}
    
    for evaluationPoint in root.findall(".//evaluationPoint"):
        evaluationPointsDict = evaluationPoint.attrib
        names.append(evaluationPointsDict['name'])
        polynomsToCalculate.append([int(evaluationPointsDict['id']),int(evaluationPointsDict['gridNode'])])
        
        deltaDict = {}
        peaksToEvaluateDict = {}
        
        for data in evaluationPoint:
            if data.tag == "deltasMinMaxFunction":
                for value in data.findall(".//"):
                    deltaDict[value.tag] = unitConversion(unitsDict, float(value.text), value.attrib['unit'])
            if data.tag == "peaksToEvaluate":
                for value in data.findall(".//"):
                    if value.text: peaksToEvaluateDict[''.join(['extrema',value.tag])] = map(int, value.text.split(' '))
                    else: peaksToEvaluateDict[''.join(['extrema',value.tag])] = []
        
        delta[evaluationPointsDict['name']] = deltaDict
        peaksToEvaluate[evaluationPointsDict['name']] = peaksToEvaluateDict
        
    vPCconfiguration['polynomsToCalculate'] = polynomsToCalculate
    vPCconfiguration['names'] = names
    vPCconfiguration['delta'] = delta
    vPCconfiguration['peaksToEvaluate'] = peaksToEvaluate
                
    return vPCconfiguration
