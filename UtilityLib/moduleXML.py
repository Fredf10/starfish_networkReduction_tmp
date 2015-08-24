
try:
    from lxml import etree
except:
    from xml.etree import ElementTree as etree

import os,sys

#from pprint import pprint as pp


# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))
#sys.path.append(cur+'/../'+'/NetworkLib')

from NetworkLib.classVascularNetwork import VascularNetwork
from NetworkLib.classBoundaryConditions import *

from constants import variablesDict 
from constants import unitsDictSI as unitsDict
from constants import newestNetworkXmlVersion
import moduleFilePathHandler as mFPH

#sys.path.append(cur + '/../VascularPolynomialChaosLib')
from VascularPolynomialChaosLib.classRandomInputManager import RandomInputManager

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

def writeRandomInputElement(subElement, variable, randomInputManager,randomInputLocation):
    '''
    writes a random variable in xml file
    
    input: 
        subElement := mother xml element where the randomInputElement should be added
        variable   := name of the deterministic random variable
        randomInputManager := reference to the randomInputManager of the case
        randomInputNameLocation := specific location of the randomInput in the vascular network    
    '''
    ## import current network xml description as nxmlW(rite) to avoid version clash
    from constants import newestNetworkXml as nxmlW
    
    randomInput = randomInputManager(randomInputManager.map[randomInputLocation])   
    
    if 'generalRandomInput' in randomInputLocation:
        attributes = {}
        for attribute in nxmlW.generalRandomInputsAttributes:
            attributes[attribute] = str(randomInput.getVariableValue(attribute))
        subsubElement = etree.SubElement(subElement,'generalRandomInput',attributes)
    else:
        subsubElement = etree.SubElement(subElement, "-".join([variable,'randomInput']))
    
    for randomInputElement in nxmlW.randomInputDistributionElements:
        subsubsubElement = etree.SubElement(subsubElement, randomInputElement) 
        writeXMLsetUnit(subsubsubElement,randomInputElement)
        value = randomInput.getVariableValue(randomInputElement)
        writeXMLsaveValues(subsubsubElement,randomInputElement,value)  


def writeNetworkToXML(vascularNetwork, dataNumber = "xxx", networkXmlFile = None):
    '''
    This function creates an XML file and writes all variable data of a vascularNetwork into it (except solution)
    The forma of the XML and all variable data are defined in constants.py
    '''
    networkName = vascularNetwork.getVariableValue('name')
        
    if networkXmlFile == None:
        networkXmlFile =  mFPH.getFilePath('networkXmlFile', networkName, dataNumber, 'write')
        
    try:
        root = etree.Element(networkName, id = dataNumber, version = newestNetworkXmlVersion)
    except:
        print " Error: path / file does not exist"
        return
        
    ## import current network xml description as nxmlW(rite) to avoid version clash
    from constants import newestNetworkXml as nxmlW
        
    xmlFile = etree.ElementTree(root)
    
    ### polynomial chaos distributionstuff
    randomInputManager = vascularNetwork.randomInputManager
    
    for xmlElementName in nxmlW.xmlElements:
        xmlFileElement = etree.SubElement(root, xmlElementName)
        xmlElement = nxmlW.xmlElementsReference[xmlElementName]
        
        if xmlElementName == 'boundaryConditions':
            for vesselId,boundaryConditions in vascularNetwork.boundaryConditions.iteritems():
                subElement = etree.SubElement(xmlFileElement, 'boundaryCondition', vesselId = str(vesselId))
                # loop saved condition instances
                for boundaryCondition in boundaryConditions:
                    boundaryType = boundaryCondition.getVariableValue('name')
                    subsubElement = etree.SubElement(subElement, boundaryType)
                    # loop variables of the instance to be saved in xml
                    for variable in nxmlW.boundaryConditionElements[boundaryType]:
                        variableElement = etree.SubElement(subsubElement, variable)
                        writeXMLsetUnit(variableElement,variable)
                        writeXMLsaveValues(variableElement,variable,boundaryCondition.getVariableValue(variable))  
                        # polynomial chaos
                        randomInputLocation = '_'.join(['boundaryCondition',boundaryType,str(vesselId),variable])
                        if randomInputLocation in randomInputManager.map:
                            writeRandomInputElement(subsubElement, variable, randomInputManager, randomInputLocation)
                            
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
                        # polynomial chaos
                        randomInputLocation = '_'.join(['vessel',str(vessel.Id),variable])
                        if randomInputLocation in randomInputManager.map:
                            writeRandomInputElement(subElement, variable, randomInputManager, randomInputLocation)
                            
        
        elif xmlElementName == 'communicators':
            for comId,comData in vascularNetwork.communicators.iteritems():
                comType = comData['comType']
                subElement = etree.SubElement(xmlFileElement, comType)
                for variable in nxmlW.communicatorReference[comData['comType']]:
                    subsubElement = etree.SubElement(subElement, variable)
                    writeXMLsetUnit(subsubElement,variable)
                    writeXMLsaveValues(subsubElement,variable,comData[variable])  
                    
        elif xmlElementName == 'baroreceptors':
            for baroId, baro in vascularNetwork.baroreceptors.iteritems():
                # add baroElement
                subElement = etree.SubElement(xmlFileElement, 'baroreceptor', Id = str(baroId))
                for subsubElementTag in nxmlW.xmlElementsReference[xmlElementName]:
                    # check for sub element model and add attribute
                    attributes = {}
                    if subsubElementTag == "model":
                        classType = baro.getVariableValue('modelName')
                        attributes['type'] = classType
                        variables = nxmlW.baroreceptorElementReference[subsubElementTag][classType]
                    else: variables = nxmlW.baroreceptorElementReference[subsubElementTag]
                    subsubElement = etree.SubElement(subElement, subsubElementTag,attributes)
                    # save variables
                    for variable in variables:
                        subsubsubElement = etree.SubElement(subsubElement, variable)
                        writeXMLsetUnit(subsubsubElement,variable)
                        writeXMLsaveValues(subsubsubElement,variable,baro.getVariableValue(variable))  
                        # polynomial chaos
                        randomInputLocation = '_'.join(['baroreceptor',str(baroId),variable])
                        if randomInputLocation in randomInputManager.map:
                            writeRandomInputElement(subsubElement, variable, randomInputManager, randomInputLocation)
                    
        elif xmlElementName == 'generalRandomInputs':
            for randomInput in randomInputManager():
                if randomInput.randomInputType == 'generalRandomInput':
                    writeRandomInputElement(xmlFileElement, None, randomInputManager, randomInput.location)
                            
                                
        else: # vascularNetwork
            for variable in xmlElement:
                subElement = etree.SubElement(xmlFileElement, variable) # add subElement
                writeXMLsetUnit(subElement,variable)  # check unit
                if xmlElementName == 'globalFluid': # get variable values from varscularNetwork     
                    variableValues = vascularNetwork.globalFluid[variable]
                else: 
                    variableValues = vascularNetwork.getVariableValue(variable)
                writeXMLsaveValues(subElement,variable,variableValues)  
                ## check for interaval data of global fluid TODO: reimplement polyomial chaos for global fluid
#                 try: 
#                     subElementIntervalValues = vascularNetwork.globalFluidPolyChaos[variable]
#                     intervalElement = etree.SubElement(xmlFileElement, "-".join([variable,'polyChaos']))
#                     writeXMLsetUnit(intervalElement,variable)
#                     intervalElement.text = ' '.join(str(i) for i in subElementIntervalValues)
#                 except : pass   
                                
    
    xmlFile.write(networkXmlFile,encoding='iso-8859-1',pretty_print = True)
    

def loadVariablesConversion(variable, variableValueStr, variableUnit, unit = 'unitSI'):
    '''
    checks the element.text string and
    evaluates it corresponding to the definition
    of the variable in the constants variablesDict 
    
    return converted evaluated value of variable
    '''
    
    multiVariable = False
    variableValue = 'notConvertable'
    convertError = []
    variableTypes = variablesDict[variable]['type']
    # check if variable is a multiple variable (means list of variables)
    if variablesDict[variable]['multiVar']:
        # TODO: Check with Vinz if there's a reason to use ' ', as sep=None seems better
        variableValueStrings = variableValueStr.split() #variableValueStr.split(' ')
        multiVariable = True
        variableValues = []
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
                    try: variableValue = int(variableValue)
                    except: convertError.append('int') 
                    
            elif variableType == 'bool': 
                try: 
                    if variableValueString == 'False': variableValue = False
                    elif variableValueString == 'True': variableValue = True
                    else: variableValue = eval(variableValueString)
                except: 
                    convertError.append('bool') 
                
            elif variableType in ['str']:
                if variableValueString in variablesDict[variable]['strCases']: variableValue = variableValueString
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
    try: variableDict =   variablesDict[variableName]
    except : variableDict = "No entry defined for This element"
    print """ERROR loadNetworkFromXML():
          variable "{}" of {} {} is not defined.
          (Hint:{}) , system exit!""".format(variableName, element, elementName, variableDict)
    exit()


def loadRandomInputElement(xmlElement,nxml,variableName,randomInputManager, randomInputLocation):
    '''
    Function to load random variable xml element
    
    input:
        xmlElement
        xmlElementReferences
        variableName := name of the random variable
        randomInputManager
        randomInputLocation    
    '''    
    
    dataDict = {'location'        : randomInputLocation,
                'variableName'    : [variableName],
                'randomInputType' : 'parametricRandomInput'}
    
    if variableName == None: dataDict['randomInputType'] = 'generalRandomInput'
    
    if variableName == None:
        for attribute in nxml.generalRandomInputsAttributes:
            try: dataDict[attribute]  = loadVariablesConversion(attribute, xmlElement.attrib[attribute], '')
            except: loadingErrorMessageVariableError(attribute, 'generalRandomInput', '')
            dataDict['location'] = '_'.join(['generalRandomInput',dataDict['name']])
                                                                 
    for variable in nxml.randomInputDistributionElements:
        try: element = xmlElement.findall(''.join(['.//',variable]))[0]
        except: loadingErrorMessageVariableError(variable, 'randomInput', variableName)
        # get variable value                        
        try: variableValueStr = element.text
        except: loadingErrorMessageValueError(variable, 'randomInput', variableName)
        # get unit 
        try: variableUnit = element.attrib['unit']
        except: variableUnit = None 
        # save converted XML-value
        dataDict[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
    randomInputManager.addRandomInput(dataDict)


def loadNetworkFromXML(networkName , dataNumber = "xxx", exception = 'Error', networkXmlFile = None, pathSolutionDataFilename = None):
    '''
    Function laods network from XML-file
    
    version of XML files supported: 4.0, 4.1, 4.2
    '''
    currentVersions = ['4.0','4.1','4.2']
    
    # read from file
    if networkName == None:
        print 'ERROR: moduleXML.loadNetworkFromXML() : load XML - no networkName passed'
        return None
    
    if networkXmlFile == None:
        networkXmlFile = mFPH.getFilePath('networkXmlFile', networkName, dataNumber, 'read', exception = exception)
    
    # create vascularNetwork instance
    vascularNetwork = VascularNetwork()
    # set name
    vascularNetwork.update({'name': networkName,
                            'dataNumber':dataNumber,
                            'pathSolutionDataFilename': pathSolutionDataFilename})
    ## create random vector
    randomInputManager = RandomInputManager()
    vascularNetwork.randomInputManager = randomInputManager
        
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
    elif xmlFileVersion == '4.1':
        import networkXml041 as nxml  
    
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
                    # loop through possible communicator class types
                    for boundaryType in nxml.xmlElementsReference[xmlElementName]:
                        # find all bcs of this type
                        for bcElements in boundaryConditionElement.findall(''.join(['.//',boundaryType])):
                            boundaryInstance = eval(nxml.bcTagsClassReferences[boundaryType])()
                            boundaryDataDict = {}
                            boundaryDataDict['name'] = boundaryType
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
                                randomInputName = ''.join([variable,'-randomInput'])
                                elementRandomInput = bcElements.findall(''.join(['.//',randomInputName]))
                                if len(elementRandomInput) == 1:  
                                    randomInputLocation = '_'.join(['boundaryCondition',boundaryType,str(vesselId),variable])
                                    loadRandomInputElement(elementRandomInput[0],
                                                              nxml,
                                                              variable,
                                                              randomInputManager,
                                                              randomInputLocation
                                                              )
                                
                                # adjust path to boundary condition file
                                if variable == 'filePathName':
                                    ## TODO: fix problem when loading an absolute path
                                    # networkDirectory = '/'.join(networkXmlFile.split('/')[0:-1])
                                    # variableValueStr = '/'.join([networkDirectory,variableValueStr])
                                    
                                    boundaryDataDict['filePathName'] = variableValueStr
                                    
                            boundaryInstance.update(boundaryDataDict)
                            boundaryInstances.append(boundaryInstance)             
                                
                    # apply read data to vascularNetwork
                    if vesselId not in vascularNetwork.boundaryConditions.keys():
                        vascularNetwork.boundaryConditions[vesselId] = boundaryInstances
                    else:
                        vascularNetwork.boundaryConditions[vesselId].extend(boundaryInstances)
                        
            elif xmlElementName == 'vessels':
                for vesselXMLnode in xmlElement.findall(''.join(['.//','vessel'])):
                    vesselData = {}
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
                            randomInputName = ''.join([variable,'-randomInput'])
                            elementRandomInput = vesselXMLnode.findall(''.join(['.//',randomInputName]))
                            if len(elementRandomInput) == 1:  
                                randomInputLocation = '_'.join(['vessel',str(vesselData['Id']),variable])
                                loadRandomInputElement(elementRandomInput[0],
                                                          nxml,
                                                          variable,
                                                          randomInputManager,
                                                          randomInputLocation)
                                
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
                for baroreceptorElement in xmlElement.findall(''.join(['.//','baroreceptor'])):
                    try: baroId = int(baroreceptorElement.attrib['Id'])
                    except: loadingErrorMessageVariableError('baroId', 'one baroreceptor', '')
                    
                    baroreceptorData = {'baroId':baroId}
                    # loop through top level elements:
                    for baroreceptorElementTag in nxml.xmlElementsReference[xmlElementName]:
                        baroreceptorTopLevelElement = baroreceptorElement.findall(''.join(['.//',baroreceptorElementTag]))[0]
                        
                        if baroreceptorElementTag == "model":
                            try: baroType = baroreceptorTopLevelElement.attrib['type']
                            except: loadingErrorMessageVariableError('type', 'baroreceptor', str(baroId))
                            variables = nxml.baroreceptorElementReference[baroreceptorElementTag][baroType]
                            baroreceptorData['modelName'] = baroType
                        else: variables = nxml.baroreceptorElementReference[baroreceptorElementTag]
                        # read all variables
                        for variable in variables: 
                            try: 
                                element = baroreceptorTopLevelElement.findall(''.join(['.//',variable]))[0]
                            except: loadingErrorMessageVariableError(variable, 'baroreceptors', baroId)
                            # get variable value                        
                            try: variableValueStr = element.text
                            except: loadingErrorMessageValueError(variable, 'baroreceptors', baroId)
                            # get unit
                            try: variableUnit = element.attrib['unit']
                            except: variableUnit = None 
                            baroreceptorData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                            if variable == 'baroId': baroId = baroreceptorData[variable]
                                                        
                            # find polynomial chaos variable
                            randomInputName = ''.join([variable,'-randomInput'])
                            elementRandomInput = baroreceptorTopLevelElement.findall(''.join(['.//',randomInputName]))
                            if len(elementRandomInput) == 1:  
                                randomInputLocation = '_'.join(['baroreceptor',str(baroId),variable])
                                loadRandomInputElement(elementRandomInput[0],
                                                       nxml,
                                                       variable,
                                                       randomInputManager,
                                                       randomInputLocation)
                            
                    vascularNetwork.updateNetwork({'baroreceptors': {baroId:baroreceptorData}})
                
            elif xmlElementName == 'globalFluid':
                    
                globalFluidData = {}
                ## TODO: reimplement global fluid polynomial chaos
                ##globalFluidPolychaosData = {}
                
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
                    
                    ## TODO: reimplement global fluid polynomial chaos
#                     # find polynomial chaos variable
#                     try:                        
#                         # find polychaos variables
#                         element = xmlElement.findall(''.join(['.//',variable,'-polyChaos']))[0]
#                         # get variable value                        
#                         try: variableValueStr = element.text
#                         except: loadingErrorMessageValueError(variable, 'global fluid', '')
#                         # get unit
#                         try: variableUnit = element.attrib['unit']
#                         except: variableUnit = None 
#                         # save converted XML-value                      
#                         globalFluidPolychaosData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
#                     except: pass
                vascularNetwork.updateNetwork({'globalFluid':globalFluidData}) ##,'globalFluidPolyChaos':globalFluidPolychaosData})
                    
            elif xmlElementName == 'generalRandomInputs':
                # loop through possible communicator class types
                for generalRandomInputXMLnode in xmlElement.findall(''.join(['.//','generalRandomInput'])):
                    loadRandomInputElement(generalRandomInputXMLnode,
                                           nxml,
                                           None,
                                           randomInputManager,
                                           None)
                    
                    
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
                
    # link random variables
    randomInputManager.linkRandomInputUpdateFunctions(vascularNetwork)
                           
    return vascularNetwork


###----------------------------------------------------------------------------------------
### Polynomial chaos

def savePolyChaosXML(vpcConfigXmlFile,networkName,dataNumber, vPCconfiguration = None):
    '''
    Function to write a xml file with vascularPolynomialChaos Configurations
    '''
    from constants import variableUnitsSI as variableUnits
    from constants import vPCconfigurationTemplate
        
    root = etree.Element(networkName, id= dataNumber, version="1.0")
        
    xmlFile = etree.ElementTree(root)
    
    if vPCconfiguration == None:
        vPCconfiguration = vPCconfigurationTemplate
    
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
    for idNodeTuple in vPCconfiguration['locationsToEvaluate']:
        try:
            name = str(vPCconfiguration['locationNames'][evalPointCount])
        except: "Error: no name for id,node {} tuple defined".format(idNodeTuple)
        evaluationPoint = etree.SubElement(evaluationPoints, "evaluationPoint",  vesselId = str(idNodeTuple[0]), name = str(name), gridNode = str(idNodeTuple[1]))        
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

    xmlFile.write(vpcConfigXmlFile,encoding='iso-8859-1',pretty_print = True)
    print " ... vpcConfig template saved"
  
  

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
  
    
def loadPolyChaosXML(vpcConfigXmlFile):
    ### import units of all variables in the SI system outdated used of only few functions -> to be changed
    ## just used for polynomial chaos config.
    from constants import variableUnitsSI as variableUnits
    from constants import vPCconfigurationTemplate
    
    vPCconfiguration =  {}
    
    # load the data!
    try:
        parser = etree.XMLParser(encoding='iso-8859-1')
        tree = etree.parse(vpcConfigXmlFile, parser)
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
        polynomsToCalculate.append([int(evaluationPointsDict['vesselId']),int(evaluationPointsDict['gridNode'])])
        
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
        
    vPCconfiguration['locationsToEvaluate'] = polynomsToCalculate
    vPCconfiguration['locationNames'] = names
    vPCconfiguration['delta'] = delta
    vPCconfiguration['peaksToEvaluate'] = peaksToEvaluate
                
    return vPCconfiguration
