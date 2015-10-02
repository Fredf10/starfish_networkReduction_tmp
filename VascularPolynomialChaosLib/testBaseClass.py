

# TODO: Fix imports
import os,sys
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append('/'.join([cur,'..','UtilityLib']))
from constants import unitsDictSI as unitsDict
        
try:
    from lxml import etree
except:
    from xml.etree import ElementTree as etree
        
class TestBaseClass(object):
     
    externVariables      = {}
    externXmlAttributes  = []
    externXmlElements    = []
    
    class ExtValue(object):
        
        def __init__(self, variableType, unit = None, strCases = None, multiVar = False):
            '''
            
            Note:
                        
            Args:
            
            Attributes:
            
            Examples:
            
            Raises:
            
            Returns: 
            '''
            
            if type(variableType) is not list:
                variableType = [variableType]
                
                #raise ValueError("ERROR: extValue in {}, variableType is not a list:  <<{}>>.".format(self.__class__.__name__,variableType))
                  
            for varType in variableType:
                if varType not in [bool,int,float,str,None]:
                    raise ValueError("ERROR: extValue in {} has non supported type <<{}>>.".format(self.__class__.__name__,varType))
            
            self.variableType = variableType
            self.unit     = unit
            self.multiVar = multiVar
            
            if str in self.variableType:
                if strCases is not None:
                    self.strCases = strCases
                else: 
                    raise ValueError("ERROR: extValue of type <str> in {} has no <strCases> -list defined.".format(self.__class__.__name__,variableType))
                
    class ExtDict(object):
        
        def __init__(self, dictObjName, dictObjType):
            '''
            
            Note:
                        
            Args:
            
            Attributes:
            
            Examples:
            
            Raises:
            
            Returns: 
            '''
            
            self.variableType = ['dict']
            self.dictObjName = dictObjName
            
            self.dictObjType = dictObjType
            
            #TODO: approptiate testing for dictObjType  
#             print type(dictObjType)
#             if type(dictObjType) is type(object):
#             else: raise ValueError("ERROR: ExtDict in {}, dictObjType is not a class instance: <<{}>>".format(self.__class__.__name__,dictObjType))
#             
    class ExtObject(object):
        
        def __init__(self, classCases):
            '''
            
            Note:
                        
            Args:
            
            Attributes:
            
            Examples:
            
            Raises:
            
            Returns: 
            '''
            
            self.variableType = ['object']
            
            if type(classCases) is dict:            
                self.classCases = classCases
            else: raise ValueError("ERROR: ExtObject in {}, classCases is not a dictionary: <<{}>>".format(self.__class__.__name__, classCases))
        
        
    def writeDataToXmlNode(self,xmlNode):
        '''
        function writes data to xml-node
        
        
        Note:
                    
        Args:
        
        Attributes:
        
        Examples:
        
        Raises:
        
        Returns: 
        '''        
        # 1. write attributes
        for attribute in self.externXmlAttributes:
            xmlNode.set(attribute, str(self.getVariable(attribute)))
            
        # 2. loop through data to save
        for externXmlElement in self.externXmlElements: 
            # check if lists are proper defined
            if externXmlElement in self.externVariables:
                # get the extern variable definition class
                externVariable  = self.externVariables[externXmlElement]
                # create xml node
                externXmlNode  = etree.SubElement(xmlNode, externXmlElement)  
                ## find out what type the externXmlElement variable is:
                # 2.1 if dict variable -> writeExtValueXml
                if isinstance(externVariable,self.ExtDict): 
                    self.writeExtDictXml(externXmlNode, self.getVariable(externXmlElement), externVariable)
                # 2.2 if object variable -> writeExtDictXml       
                elif isinstance(externVariable,self.ExtObject):
                    self.writeExtObjectXml(externXmlNode, self.getVariable(externXmlElement), externVariable)
                # 2.3 if value variable -> writeExtObjectXml
                elif isinstance(externVariable,self.ExtValue):
                    self.writeExtValueXml(externXmlNode, self.getVariable(externXmlElement), externVariable)
                    
            else: raise KeyError("""ERROR: try to write <<{}>> to xml-node {},
               however <<{}>> is not defined in self.externVariables""".format(externXmlElement, xmlNode, externXmlElement))
            
    def writeExtDictXml(self,externXmlNode, dictToWrite, externVariable):
        '''
        Writes dictionary xml note and envokes the dictionary data to be written
        
        Note:
                    
        Args:
            externXmlNode (etree element)  : node in the xml file
            externXmlElement (dict)        : dict to be written
            externVariable (ExtObject)     : instance of ExtObject defining the variable properties
        '''
        # iterate through dictionary
        for key,value in dictToWrite.iteritems():
            # find out the name of the variable and create a xml node
            externXmlDictNode  = etree.SubElement(externXmlNode, externVariable.dictObjName) 
            # write id
            externXmlDictNode.set('Id', str(key))
            # if object variable -> writeExtDictXml       
            if isinstance(externVariable.dictObjType,self.ExtObject):
                self.writeExtObjectXml(externXmlDictNode, value, externVariable.dictObjType)
            # if value variable -> writeExtObjectXml
            elif isinstance(externVariable.dictObjType,self.ExtValue):
                self.writeExtValueXml(externXmlDictNode, value, externVariable.dictObjType)
        
    def writeExtObjectXml(self,externXmlNode, classToWrite, externVariable):
        '''
        externXmlNode := node in the xml file
        classToWrite (object) : class-instance to be written
        externVariable := instance of ExtObject defining the variable properties
        '''
        classNameDefined = False
        for definedClassName, definedClasses in externVariable.classCases.iteritems():
            if isinstance(classToWrite,definedClasses):
                externXmlNode.set('class',definedClassName)
                classNameDefined = True
        if classNameDefined == False:
            className = classToWrite.__class__.__name__
            raise KeyError("""ERROR: try to write class <<{}>>,
       however <<{}>> is not defined in externVariable.classCases: {}""".format(classToWrite, className,externVariable.classCases ))
        classToWrite.writeDataToXmlNode(externXmlNode)
        
    def writeExtValueXml(self,externXmlNode, variableValues, externVariable):
        '''
        Args:
            externXmlNode : node in the xml file
            variableValues : str of the variable values
            externVariable : instance of ExtObject defining the variable properties
        '''
        # set unit if existing
        if externVariable.unit != None: externXmlNode.set('unit', externVariable.unit)
        # write variable value
        if externVariable.multiVar == True:
            externXmlNode.text = ' '.join(str(i) for i in variableValues)
        else:
            externXmlNode.text = str(variableValues)
        
    def readDataFromXmlNode(self, xmlNode):
        '''
        parser which parses the xml-node given and search for data defined in
        
        externAttributes
        externElements
        externVariables
        
        '''
        newData = {}
        # load object attributes
        for attribute in self.externXmlAttributes:            
            # check if lists are proper defined
            if attribute in self.externVariables:
                if isinstance(self.externVariables[attribute], self.ExtValue):
                    newData[attribute] = self.loadVariableConversion(xmlNode.attrib[attribute], '', self.externVariables[attribute])
                else:raise KeyError("""ERROR: try to read attribute <<{}>> of xml-node {},
                 however <<{}>> is not defined as class-instance <self.ExtValue>  in self.externVariables""".format(attribute, xmlNode, attribute))
                    
            else: raise KeyError("""ERROR: try to read attribute <<{}>> of xml-node {},
               however <<{}>> is not defined in self.externVariables""".format(attribute, xmlNode, attribute))
            
        ## TODO: check if there is not-needed information in the xml file
        for externXmlElement in self.externXmlElements: 
            # check if lists are proper defined
            if externXmlElement in self.externVariables:
                # get the extern variable definition class
                externVariable  = self.externVariables[externXmlElement]
                                
                # try to get the corresponding xml element:
                try: externXmlNode = xmlNode.findall(''.join(['.//',externXmlElement]))[0]
                except IndexError: self.loadingErrorMessageVariableError(externXmlElement, xmlNode)
                
                ## find out what type the externXmlElement variable is:
                if isinstance(externVariable,self.ExtDict): 
                    newData[externXmlElement] = self.loadExtDict(externXmlNode, externXmlElement, externVariable)
                        
                elif isinstance(externVariable,self.ExtObject):
                    newData[externXmlElement] = self.loadExtObject(externXmlNode, externXmlElement, externVariable)
                    
                elif isinstance(externVariable,self.ExtValue):
                    newData[externXmlElement] = self.loadExtValue(externXmlNode, externXmlElement, externVariable)
                    
            else: raise KeyError("""ERROR: try to read <<{}>> of xml-node {},
               however <<{}>> is not defined in self.externVariables""".format(externXmlElement, xmlNode, externXmlElement))
            
        self.setVariablesDict(newData)
     
    
    def loadExtDict(self,externXmlNode, externXmlElement, externVariable):
        # TODO: rename loadExtDictXml
        '''
        
        externXmlNode := node in the xml file
        externXmlElement := str of the variable name
        externVariable := instance of ExtDict defining the variable properties
        '''
        elementDictData= {}
                
        for dictXmlNode in externXmlNode.getchildren():
        
            dictXmlElement = dictXmlNode.tag
            
            if dictXmlElement == externVariable.dictObjName:
                dictVariable = externVariable.dictObjType
                # check if Id == dict.key is defined and unique
                if 'Id' in dictXmlNode.attrib:
                    dictXmlNodeId = dictXmlNode.attrib['Id']
                    if dictXmlNodeId  not in elementDictData:
                        
                        # check if ExtObject is expected or if ExtValue is expected
                        if isinstance(dictVariable,self.ExtObject):
                            elementDictData[dictXmlNodeId] = self.loadExtObject(dictXmlNode,dictXmlElement,dictVariable)
                            
                        elif isinstance(dictVariable,self.ExtValue):
                            elementDictData[dictXmlNodeId] = self.loadExtValue(dictXmlNode,dictXmlElement,dictVariable)
                            
                        else: raise ValueError("""ERROR try to read <<{}>> of xml-node {},
                   however <{}> is not defined as class-instance <self.ExtValue> or <self.ExtObject> in self.externVariables""".format(dictXmlElement, externXmlElement, dictXmlElement))
                        
                    else: raise KeyError("""ERROR try to read <<{}>> of xml-node {}, however dict-key Id=<<{}>> is defined multiple times""".format(dictXmlElement, externXmlElement,dictXmlNodeId))
                    
                else: raise ValueError("""ERROR try to read <<{}>> of xml-node {},
                   however attribute <<Id>> is not defined in the XML-tag""".format(dictXmlElement, externXmlElement))
                
            else: print """WARNING: try to read xml-node <<{}>> as dict element for <<{}>>,
        however this sub-type is not defined as variable-type of <<{}>>. Skipping xml-node""".format(dictXmlElement,externXmlElement,externXmlElement)
                
                
        return elementDictData
            
            
    def loadExtObject(self,externXmlNode, externXmlElement, externVariable):
        #TODO: rename to loadExtObjectXml
        '''
        externXmlNode := node in the xml file
        externXmlElement := str of the variable name
        externVariable := instance of ExtObject defining the variable properties
        '''
        # check if class in attribute
        if 'class' in externXmlNode.attrib:
            classXmlNode = externXmlNode.attrib['class']
            if classXmlNode in externVariable.classCases:
                return externVariable.classCases[classXmlNode](externXmlNode)
            else:
                raise ValueError("""ERROR try to read <<{}>> of xml-node {},
               however class-tpye <<{}>> is not defined ExtObject.classCases""".format(externXmlElement, externXmlNode, classXmlNode))
        
        else: raise ValueError("""ERROR try to read <<{}>> of xml-node {},
               however attribute <<class>> is not defined in the XML-tag""".format(externXmlElement, externXmlNode))
            
        
    def loadExtValue(self,externXmlNode, externXmlElement, externVariable):
        #TODO: rename to loadExtValueXML
        '''
        Load the value of externXmlElement:
        
        Checks the externXmlNode.text string and
        evaluates it corresponding to the properties in externVariable
        
        externXmlNode := node in the xml file
        externXmlElement := str of the variable name
        externVariable := instance of ExtValue defining the variable properties
        '''
        ## retrieve data from xml load
        # get variable value str from xml node
        variableValueStr = externXmlNode.text
        # get unit
        if 'unit' in externXmlNode.attrib: 
            variableUnit = externXmlNode.attrib['unit']
        else: variableUnit = None 
        
        return self.loadVariableConversion(variableValueStr, variableUnit, externVariable)
        
    def loadVariableConversion(self, variableValueStr, variableUnit, externVariable):
        
        # start conversion process
        multiVariable = False
        variableValue = 'notConvertable'
        convertError = []
        variableTypes = externVariable.variableType
        
        # check if variable is a multiple variable (means list of variables)
        if externVariable.multiVar:
            variableValueStrings = (variableValueStr.replace(',',' ')).split() 
            multiVariable = True
            variableValues = []
        else:
            variableValueStrings = [variableValueStr]
                
        # start conversion loop over variable and types
        for variableValueString in variableValueStrings:
            for variableType in variableTypes:
                
                if variableType in [float,int]:
                    try: variableValue = float(eval(variableValueString))
                    except ValueError: convertError.append('float') 
                    
                    if externVariable.unit != None and variableUnit != None:
                        if ' ' in variableUnit:
                            variableUnits = variableUnit.split(' ')
                            for variableUnit in variableUnits: variableValue = variableValue*unitsDict[variableUnit]
                        else: variableValue = variableValue*unitsDict[variableUnit]
                        
                    if variableType == int: 
                        try: variableValue = int(variableValue)
                        except ValueError: convertError.append('int') 
                        
                elif variableType == bool: 
                    if variableValueString == 'False': variableValue = False
                    elif variableValueString == 'True': variableValue = True
                    else: convertError.append('bool') 
                    
                elif variableType == str:
                    if variableValueString in externVariable.strCases: variableValue = variableValueString
                    elif externVariable.strCases[0] == 'anything': variableValue = variableValueString
                    else: convertError.append(''.join(['str == ',str(externVariable.strCases)]))
                
                elif variableType == None:
                    if variableValueString == 'None' or variableValueString == '' or variableValueString == None:
                        variableValue = None
                    else: convertError.append('None')
                        
            if variableValue == 'notConvertable':
                raise ValueError("""ERROR: {}.loadExtValue():
                      Cannot convert given value "{}" of variable "{}"
                      to {}!
                      Check if it is of type {}, system exit!""".format(self.__class__.__name__,
                                                                        variableValueString,
                                                                        variableValueStr,
                                                                        convertError,
                                                                        variableTypes))
            if multiVariable == False: return variableValue
            else: variableValues.append(variableValue)
        
        return variableValues   
    
    def loadingErrorMessageVariableError(self, variableName, element):
        # TODO: is this really needed?
        if variableName in self.externVariables:
            variableType = self.externVariables[variableName].variableType
        else: variableType = "No entry defined for This element"
        raise ValueError("""ERROR loadNetworkFromXML():
              variable "{}" of {} is not defined.
              (Hint:{}) , system exit!""".format(variableName, element, variableType))
        
        
    def setVariablesDict(self, Dict):
        '''
        updates the class data using a dictionary in from of 
        dataDict = {'variableName': value}
        '''
        for key,value in Dict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except KeyError:
                print "WARNING {}.updateData - Wrong key: {}, could not update varibale".format(self.__class__.__name__, key)
                
    def getVariable(self,variableName):
        '''
        Returns value of variable with name : variableName
        States Error if not such variable
        '''
        try:
            return self.__getattribute__(variableName)
        except: 
            print "ERROR Vessel.getVariable() : vessel has no variable {}".format(variableName)