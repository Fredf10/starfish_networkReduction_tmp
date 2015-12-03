from NetworkLib import classVenousPool
try:
    from lxml import etree
except:
    # TODO: This produces an error!! with the pretty_print argument
    from xml.etree import ElementTree as etree

import os,sys
import numbers

# from pprint import pprint as pp


# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))
# sys.path.append(cur+'/../'+'/NetworkLib')

import NetworkLib.classVascularNetwork as cVascNw
from NetworkLib.classBoundaryConditions import *
# XML data with names referencing ccBC's functions
# won't work if this import changes

from constants import variablesDict
from constants import unitsDictSI as unitsDict
from constants import newestNetworkXmlVersion
import moduleFilePathHandler as mFPH

#sys.path.append(cur + '/../VascularPolynomialChaosLib')
from VascularPolynomialChaosLib.classRandomInputManager import RandomInputManager

### import units of all variales in the Medical System
#from constants import variableUnitsMed as variableUnits
#from constants import unitsDictMed as unitsDict


def writeXMLsetUnit(xmlElement, variable, value, unit = 'unitSI'):
    """
    Checks if element has a unit and adds it in the XML file
    """
    try:
        # TODO:Should we add a field to write out non SI to xml?
        unitName = variablesDict[variable][unit]
        if unitName:
            xmlElement.set('unit', unitName)
            if ' ' in unitName:
                unitNames = unitName.split(' ')
                conversionFactor = 1.
                for unitName in unitNames:
                    conversionFactor = conversionFactor*unitsDict[unitName]
            else:
                conversionFactor = unitsDict[unitName]
            if conversionFactor and value and isinstance(value,numbers.Number):
                value = value/conversionFactor
    except (KeyError,TypeError) as e:
        print """ERROR: moduleXML.writeXML():
            variable {} of element {} is not properly defined
            in variablesDict""".format(variable,xmlElement)
        raise e
    return value

def writeXMLsaveValues(xmlElement, variable, variableValues, polychaos = False):
    """
    Writes the variable values to the xml file with xmlElement.text() function
    if variable can have multiple values they are saved with space as delimiter
    """
    if variablesDict[variable]['multiVar'] or polychaos:
        xmlElement.text = ' '.join(str(i) for i in variableValues)
    else:
        xmlElement.text = str(variableValues)

# def writeRandomInputElement(subElement, variable, randomInputManager,randomInputLocation):
#     """
#     writes a random variable in xml file
#
#     input:
#         subElement := mother xml element where the randomInputElement should be added
#         variable   := name of the deterministic random variable
#         randomInputManager := reference to the randomInputManager of the case
#         randomInputNameLocation := specific location of the randomInput in the vascular network
#     """
#     ## import current network xml description as nxmlW(rite) to avoid version clash
#     from constants import newestNetworkXml as nxmlW
#
#     randomInput = randomInputManager(randomInputManager.map[randomInputLocation])
#
#     if 'generalRandomInput' in randomInputLocation:
#         attributes = {}
#         for attribute in nxmlW.generalRandomInputsAttributes:
#             attributes[attribute] = str(randomInput.getVariableValue(attribute))
#         subsubElement = etree.SubElement(subElement,'generalRandomInput',attributes)
#     else:
#         subsubElement = etree.SubElement(subElement, "-".join([variable,'randomInput']))
#
#     for randomInputElement in nxmlW.randomInputDistributionElements:
#         subsubsubElement = etree.SubElement(subsubElement, randomInputElement)
#         value = randomInput.getVariableValue(randomInputElement)
#         writeXMLsetUnit(subsubsubElement,randomInputElement,value)
#         writeXMLsaveValues(subsubsubElement,randomInputElement,value)


def writeNetworkToXML(vascularNetwork, dataNumber = "xxx", networkXmlFile = None):
    """
    This function creates an XML file and writes all variable data of a vascularNetwork into it (except solution)
    The forma of the XML and all variable data are defined in constants.py
    """
    networkName = vascularNetwork.getVariableValue('name')

    if networkXmlFile == None:
        networkXmlFile =  mFPH.getFilePath('networkXmlFile', networkName, dataNumber, 'write')

    try:
        root = etree.Element(networkName, id = dataNumber, version = newestNetworkXmlVersion)
    except: #TODO: find out what errors etree.Element raises and fix this except statement.
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
                        value = boundaryCondition.getVariableValue(variable)
                        value = writeXMLsetUnit(variableElement,variable,value)
                        writeXMLsaveValues(variableElement,variable,value)
#                         # polynomial chaos
#                         randomInputLocation = '_'.join(['boundaryCondition',boundaryType,str(vesselId),variable])
#                         if randomInputLocation in randomInputManager.map:
#                             writeRandomInputElement(subsubElement, variable, randomInputManager, randomInputLocation)

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
                        variableValues = vessel.getVariableValue(variable)
                        variableValues = writeXMLsetUnit(subsubElement,variable,variableValues)
                        writeXMLsaveValues(subsubElement,variable,variableValues)
#                         # check for polychaos data of vessel
#                         # polynomial chaos
#                         randomInputLocation = '_'.join(['vessel',str(vessel.Id),variable])
#                         if randomInputLocation in randomInputManager.map:
#                             writeRandomInputElement(subElement, variable, randomInputManager, randomInputLocation)


        elif xmlElementName == 'communicators':
            for comId,comData in vascularNetwork.communicators.iteritems():
                comType = comData['comType']
                subElement = etree.SubElement(xmlFileElement, comType)
                for variable in nxmlW.communicatorReference[comData['comType']]:
                    subsubElement = etree.SubElement(subElement, variable)
                    value = comData[variable]
                    value = writeXMLsetUnit(subsubElement,variable, value)
                    writeXMLsaveValues(subsubElement,variable,value)

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
                        value = baro.getVariableValue(variable)
                        value = writeXMLsetUnit(subsubsubElement,variable, value)
                        writeXMLsaveValues(subsubsubElement,variable,value)
##                        polynomial chaos
#                         randomInputLocation = '_'.join(['baroreceptor',str(baroId),variable])
#                         if randomInputLocation in randomInputManager.map:
#                             writeRandomInputElement(subsubElement, variable, randomInputManager, randomInputLocation)

#         elif xmlElementName == 'generalRandomInputs':
#             for randomInput in randomInputManager():
#                 if randomInput.randomInputType == 'generalRandomInput':
#                     writeRandomInputElement(xmlFileElement, None, randomInputManager, randomInput.location)

        elif xmlElementName == 'venousPool':
            if vascularNetwork.venousPool is not None:
                venousPoolWrapper = classVenousPool.VenousPoolXMLWrapper()
                venousPoolWrapper.venousPoolContent = vascularNetwork.venousPool
                venousPoolWrapper.writeDataToXmlNode(xmlFileElement)

        elif xmlElementName == 'randomInputManager':
            if vascularNetwork.randomInputManager is not None:
                vascularNetwork.randomInputManager.writeDataToXmlNode(xmlFileElement)

        elif xmlElementName == "externalStimuli":
            for stimulusId, stimulus in vascularNetwork.externalStimuli.iteritems():
                stimulusType = stimulus['type']
                # add stimuliElement
                subElement = etree.SubElement(xmlFileElement, 'externalStimulus', Id=str(stimulusId), type=stimulusType)
                variables = nxmlW.xmlElementsReference[xmlElementName][stimulusType]

                for variable in variables:
                    subsubElement = etree.SubElement(subElement, variable)
                    value = stimulus[variable]
                    value = writeXMLsetUnit(subsubElement, variable, value)
                    writeXMLsaveValues(subsubElement, variable, value)
#                     # polynomial chaos
#                     randomInputLocation = '_'.join(['externalStimulus', str(stimulusId), variable])
#                     if randomInputLocation in randomInputManager.map:
#                         writeRandomInputElement(subsubElement, variable, randomInputManager, randomInputLocation)


        else: # vascularNetwork
            for variable in xmlElement:
                subElement = etree.SubElement(xmlFileElement, variable) # add subElement
                if xmlElementName == 'globalFluid': # get variable values from varscularNetwork
                    variableValues = vascularNetwork.globalFluid[variable]
                else:
                    variableValues = vascularNetwork.getVariableValue(variable)
                variableValues = writeXMLsetUnit(subElement,variable,variableValues)  # check unit

                writeXMLsaveValues(subElement,variable,variableValues)
                ## check for interval data of global fluid
                ## TODO: reimplement polyomial chaos for global fluid
#                 try:
#                     subElementIntervalValues = vascularNetwork.globalFluidPolyChaos[variable]
#                     intervalElement = etree.SubElement(xmlFileElement, "-".join([variable,'polyChaos']))
#                     writeXMLsetUnit(intervalElement,variable)
#                     intervalElement.text = ' '.join(str(i) for i in subElementIntervalValues)
#                 except : pass


    xmlFile.write(networkXmlFile,encoding='iso-8859-1',pretty_print = True)


def loadVariablesConversion(variable, variableValueStr, variableUnit, unit = 'unitSI'):
    """
    checks the element.text string and
    evaluates it corresponding to the definition
    of the variable in the constants variablesDict

    return converted evaluated value of variable
    """

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
                continueConversion = False
                try:
                    variableValue = float(eval(variableValueString))
                    continueConversion = True
                except (ValueError, TypeError, NameError): convertError.append('float')

                if continueConversion == True:
                    if variablesDict[variable][unit] != None and variableUnit != None:
                        if ' ' in variableUnit:
                            variableUnits = variableUnit.split(' ')
                            for variableUnit in variableUnits: variableValue = variableValue*unitsDict[variableUnit]
                        else: variableValue = variableValue*unitsDict[variableUnit]

                    if variableType == 'int':
                        try: variableValue = int(variableValue)
                        except ValueError: convertError.append('int')

            elif variableType == 'bool':
                try:
                    variableValue = eval(variableValueString)
                except TypeError:
                    print "Warning: loadVariablesConversion() bool: variableValueString is not a valid expression for eval()"
                    convertError.append('bool')
                except NameError:
                    print "Warning: loadVariablesConversion() bool: variableValueString is not defined"
                    convertError.append('bool')

            elif variableType in ['str']:
                if variableValueString in variablesDict[variable]['strCases']: variableValue = variableValueString
                elif variablesDict[variable]['strCases'][0] == 'anything': variableValue = variableValueString
                else: convertError.append(''.join(['str == ',str(variablesDict[variable]['strCases'])]))
                #TODO: fix exception handling here

            elif variableType in ['None']:
                if variableValueString == 'None' or variableValueString == '' or variableValueString == None:
                    variableValue = None
                else: convertError.append('None')
                #TODO: fix exception handling here

        if variableValue == 'notConvertable':
            raise TypeError("""ERROR: moduleXML.loadVariablesConversion():
                  Cannot convert given value "{}" of variable "{}"
                  to {}!
                  Check if it is of type {}, system exit!
                  """.format(variableValueString,variable,convertError,variableTypes))

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
    except KeyError: variableDict = "No entry defined for This element"
    print """ERROR loadNetworkFromXML():
          variable "{}" of {} {} is not defined.
          (Hint:{}) , system exit!""".format(variableName, element, elementName, variableDict)
    exit()


# def loadRandomInputElement(xmlElement,nxml,variableName,randomInputManager, randomInputLocation):
#     """
#     Function to load random variable xml element
#
#     Args:
#         xmlElement
#         xmlElementReferences
#         variableName : name of the random variable
#         randomInputManager
#         randomInputLocation
#     """
#
#     dataDict = {'location'        : randomInputLocation,
#                 'variableName'    : [variableName],
#                 'randomInputType' : 'parametricRandomInput'}
#
#     if variableName == None: dataDict['randomInputType'] = 'generalRandomInput'
#
#     if variableName == None:
#         for attribute in nxml.generalRandomInputsAttributes:
#             try: dataDict[attribute]  = loadVariablesConversion(attribute, xmlElement.attrib[attribute], '')
#             except: loadingErrorMessageVariableError(attribute, 'generalRandomInput', '')
#             dataDict['location'] = '_'.join(['generalRandomInput',dataDict['name']])
#
#     for variable in nxml.randomInputDistributionElements:
#         try: element = xmlElement.findall(''.join(['.//',variable]))[0]
#         except: loadingErrorMessageVariableError(variable, 'randomInput', variableName)
#         # get variable value
#         try: variableValueStr = element.text
#         except: loadingErrorMessageValueError(variable, 'randomInput', variableName)
#         # get unit
#         try: variableUnit = element.attrib['unit']
#         except: variableUnit = None
#         # save converted XML-value
#         dataDict[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
#     randomInputManager.addRandomInput(dataDict)

# TODO: (einar) fix exception variable, also exceptionhandler should be at end to make it consistent with other functions.
# TODO: (einar) function should pass an exception out of scope if variables passed lack necessary info.
# TODO: (einar) variable = "default" style of defining should only be used for optional arguments
def loadNetworkFromXML(networkName ,
                       dataNumber = "xxx",
                       exception = 'Error',
                       networkXmlFile = None,
                       pathSolutionDataFilename = None):
    """
    Function loads network from XML-file

    version of XML files supported: 4.0, 4.1, 4.2
    """
    currentVersions = ['4.0','4.1','4.2']

    # read from file
    if networkName == None:
        print 'ERROR: moduleXML.loadNetworkFromXML() : load XML - no networkName passed'
        return None

    if networkXmlFile == None:
        networkXmlFile = mFPH.getFilePath('networkXmlFile', networkName, dataNumber, 'read', exception = exception)

    # create vascularNetwork instance
    vascularNetwork = cVascNw.VascularNetwork()
    # set name
    vascularNetwork.update({'name': networkName,
                            'dataNumber':dataNumber,
                            'pathSolutionDataFilename': pathSolutionDataFilename})

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
        print " WARNING the version of the network xml file you try to load is outdated it may cause some problems!"

    # TODO: (einar) This for loop is too complicated, should be simplified or split up.
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
                            boundaryDataDict = {"vesselId":vesselId}
                            boundaryDataDict['name'] = boundaryType
                            # loop through all variables of this type, convert and save values of these
                            for variable in nxml.xmlElementsReference[xmlElementName][boundaryType]:
                                # find normal variables
                                try:
                                    element = bcElements.findall(''.join(['.//',variable]))[0]
                                    # get variable value
                                    try: variableValueStr = element.text
                                    except: loadingErrorMessageValueError(variable, 'boundaryCondition', boundaryType)
                                    # get unit
                                    try: variableUnit = element.attrib['unit']
                                    except: variableUnit = None
                                except: loadingErrorMessageVariableError(variable, 'boundaryCondition', boundaryType)
                                # save converted XML-value
                                boundaryDataDict[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)

#                                 # find polynomial chaos variable
#                                 randomInputName = ''.join([variable,'-randomInput'])
#                                 elementRandomInput = bcElements.findall(''.join(['.//',randomInputName]))
#                                 if len(elementRandomInput) == 1:
#                                     randomInputLocation = '_'.join(['boundaryCondition',boundaryType,str(vesselId),variable])
#                                     loadRandomInputElement(elementRandomInput[0],
#                                                               nxml,
#                                                               variable,
#                                                               randomInputManager,
#                                                               randomInputLocation
#                                                               )

                                # adjust path to boundary condition file
                                if variable == 'filePathName':
                                    ## TODO: fix problem when loading an absolute path
                                    networkDirectory = '/'.join(networkXmlFile.split('/')[0:-1])
                                    # variableValueStr = '/'.join([networkDirectory,variableValueStr])
                                    boundaryDataDict['networkDirectory'] = networkDirectory
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
                            try:
                                element = vesselXMLnode.findall(''.join(['.//',variable]))[0]
                                # get variable value
                                try: variableValueStr = element.text
                                except: loadingErrorMessageValueError(variable, 'vessel', vesselData['Id'])
                                # get unit
                                try: variableUnit = element.attrib['unit']
                                except KeyError: variableUnit = None
                            except: loadingErrorMessageVariableError(variable, 'vessel', vesselData['Id'])

                            # save converted XML-value
                            vesselData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)

#                             # find polynomial chaos variable
#                             randomInputName = ''.join([variable,'-randomInput'])
#                             elementRandomInput = vesselXMLnode.findall(''.join(['.//',randomInputName]))
#                             if len(elementRandomInput) == 1:
#                                 randomInputLocation = '_'.join(['vessel',str(vesselData['Id']),variable])
#                                 loadRandomInputElement(elementRandomInput[0],
#                                                           nxml,
#                                                           variable,
#                                                           randomInputManager,
#                                                           randomInputLocation)

                    vascularNetwork.updateNetwork({'vesselData':{vesselData['Id']:vesselData}})


            elif xmlElementName == 'communicators':
                # loop through possible communicator class types
                for comunicatorType in nxml.xmlElementsReference[xmlElementName]:
                    # find all communicator of this type
                    for comElements in xmlElement.findall(''.join(['.//',comunicatorType])):
                        # loop through all variables of this type, convert and save values of these
                        communicatorData = {}
                        for variable in nxml.xmlElementsReference[xmlElementName][comunicatorType]:
                            try:
                                element = comElements.findall(''.join(['.//',variable]))[0]
                                try: variableValueStr = element.text
                                except: loadingErrorMessageValueError(variable, 'communicator', comunicatorType)
                                # get unit
                                try: variableUnit = element.attrib['unit']
                                except: variableUnit = None
                            except: loadingErrorMessageVariableError(variable, 'communicator', comunicatorType)
                            # get variable value

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
                                try: variableValueStr = element.text
                                except: loadingErrorMessageValueError(variable, 'baroreceptors', baroId)
                                # get unit
                                try: variableUnit = element.attrib['unit']
                                except: variableUnit = None
                            except: loadingErrorMessageVariableError(variable, 'baroreceptors', baroId)
                            # get variable value

                            baroreceptorData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                            if variable == 'baroId': baroId = baroreceptorData[variable]

#                             # find polynomial chaos variable
#                             randomInputName = ''.join([variable,'-randomInput'])
#                             elementRandomInput = baroreceptorTopLevelElement.findall(''.join(['.//',randomInputName]))
#                             if len(elementRandomInput) == 1:
#                                 randomInputLocation = '_'.join(['baroreceptor',str(baroId),variable])
#                                 loadRandomInputElement(elementRandomInput[0],
#                                                        nxml,
#                                                        variable,
#                                                        randomInputManager,
#                                                        randomInputLocation)

                    vascularNetwork.updateNetwork({'baroreceptors': {baroId:baroreceptorData}})


            elif xmlElementName == 'externalStimuli':
                for externalStimulusElement in xmlElement.findall(''.join(['.//','externalStimulus'])):
                    try: stimulusId = int(externalStimulusElement.attrib['Id'])
                    except: loadingErrorMessageVariableError('stimulusId', 'one stimulus', '')
                    try: stimulusType = externalStimulusElement.attrib['type']
                    except: loadingErrorMessageVariableError('type', 'one stimulus', '')

                    stimulusData = {'stimulusId':stimulusId}
                    stimulusData['type'] = stimulusType
                    externalStimulusTopLevelElement = externalStimulusElement
                    variables = nxml.externalStimulusElements[stimulusType]
                    for variable in variables:
                        try:
                            element = externalStimulusTopLevelElement.findall(''.join(['.//',variable]))[0]
                            try: variableValueStr = element.text
                            except: loadingErrorMessageValueError(variable, 'externalStimuli', stimulusId)
                            # get unit
                            try: variableUnit = element.attrib['unit']
                            except: variableUnit = None
                        except: loadingErrorMessageVariableError(variable, 'externalStimuli', stimulusId)
                        # get variable value

                        stimulusData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                        if variable == 'stimulusId': stimulusId = stimulusData[variable]

#                         # find polynomial chaos variable
#                         randomInputName = ''.join([variable,'-randomInput'])
#                         elementRandomInput = externalStimulusTopLevelElement.findall(''.join(['.//',randomInputName]))
#                         if len(elementRandomInput) == 1:
#                             randomInputLocation = '_'.join(['externalStimulus',str(stimulusId),variable])
#                             loadRandomInputElement(elementRandomInput[0],
#                                                    nxml,
#                                                    variable,
#                                                    randomInputManager,
#                                                    randomInputLocation)

                    vascularNetwork.updateNetwork({'externalStimuli': {stimulusId:stimulusData}})

            elif xmlElementName == 'globalFluid':

                globalFluidData = {}
                ## TODO: reimplement global fluid polynomial chaos
                ##globalFluidPolychaosData = {}

                for variable in nxml.xmlElementsReference[xmlElementName]:
                    # find normal variables
                    try:
                        element = xmlElement.findall(''.join(['.//',variable]))[0]
                        # get variable value
                        try: variableValueStr = element.text
                        except: loadingErrorMessageValueError(variable, 'global fluid', '')
                        # get unit
                        try: variableUnit = element.attrib['unit']
                        except: variableUnit = None
                    # save converted XML-value
                    except: loadingErrorMessageVariableError(variable, 'global fluid', '')

                    globalFluidData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)

                    ## TODO: reimplement global fluid polynomial chaos
                vascularNetwork.updateNetwork({'globalFluid':globalFluidData}) ##,'globalFluidPolyChaos':globalFluidPolychaosData})

#             elif xmlElementName == 'generalRandomInputs':
#                 # loop through possible communicator class types
#                 for generalRandomInputXMLnode in xmlElement.findall(''.join(['.//','generalRandomInput'])):
#                     loadRandomInputElement(generalRandomInputXMLnode,
#                                            nxml,
#                                            None,
#                                            randomInputManager,
#                                            None)

            elif xmlElementName == "venousPool":
                if xmlElement:
                    venousPoolWrapper = classVenousPool.VenousPoolXMLWrapper()
                    venousPoolWrapper.loadDataFromXmlNode(xmlElement)
                    vascularNetwork.venousPool = venousPoolWrapper.venousPoolContent

            elif xmlElementName == 'randomInputManager':
                ## create random vector
                if xmlElement:
                    randomInputManager = RandomInputManager()
                    vascularNetwork.randomInputManager = randomInputManager
                    vascularNetwork.randomInputManager.loadDataFromXmlNode(xmlElement)

            elif xmlElementName in nxml.vascularNetworkElements: # vascularNetwork
                vascularNetworkData = {}
                for variable in nxml.xmlElementsReference[xmlElementName]:
                    try:
                        element = xmlElement.findall(''.join(['.//',variable]))[0]
                        # get variable value
                        try: variableValueStr = element.text
                        except: loadingErrorMessageValueError(variable, xmlElementName, '')
                        # get
                        try: variableUnit = element.attrib['unit']
                        except: variableUnit = None
                    except: loadingErrorMessageVariableError(variable, xmlElementName, '')

                    # save converted XML-value
                    vascularNetworkData[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                vascularNetwork.update(vascularNetworkData)

    return vascularNetwork







