#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys,os
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur+'/../')


import moduleFilePathHandlerVPC as mFPH_VPC

class VpcConfiguration(object):
    '''
    Configuration class of a vascular polynomial chaos class
    '''
    def __init__(self, xmlNode):
        '''
        define all variables with comments here
        
        try to:
        update variables from vpcConfig-file (networkName,dataNumber)
        '''
        ### network name and datanumber
        #self.networkName = networkName
        #self.dataNumber  = dataNumber
        
        #control variables
        ##  0.2 collocation method ( TRUE == create and save, FALSE == load existing)
        self.createSample     = True
        
        ### 1.step genrealized polynomial chaos evaluations + data storing
        self.createEvaluationXmlFiles = True
        
        self.simulateEvaluations    = True
        self.local                  = True #TODO: add functions for server
        self.multiprocessing        = True
        self.numberOfProcessors     = 12
        self.evaluationNumbers      = []
        
        # 2.1 pre process data for all locations of interest
        self.preProcessData   = True
        
        ### 3.step post processing - orthogonal polynomials, gpce expansion, uncertainty quantification and sensitivity analysis
        self.postProcessing  = True
        
        #----POLYNOMIAL CHAOS DEFINITIONS -------------------------------------------------------------#
        self.runPolynomialChaos = True
        #polynomialOrders of the polynomial chaos expansion || if more then one given they are processed consecutevely
        self.polynomialOrders = [2,3]
        # method of the spares grid collocation 
        self.sampleMethod = 'M'
        #Parameters
        #----------
        #sample : str
        #         Normal sampling schemes
        #         Key     Name                Nested
        #         ----    ----------------    ------
        #         "K"     Korobov             no
        #         "R"     (Pseudo-)Random     no
        #         "L"     Latin hypercube     no
        #         "S"     Sobol               yes
        #         "H"     Halton              yes
        #         "M"     Hammersley          yes
        #     
        #         Grided sampling schemes
        #         Key     Name                Nested
        #         ----    ----------------    ------
        #         "C"     Chebyshev nodes     no
        #         "NC"    Nested Chebyshev    yes
        #         "G"     Gaussian quadrTrueature no
        #         "RG"    Regular grid        no
        # 
        #  
        
        #----MONTE CARLO DEFINITIONS -------------------------------------------------------------#
        self.runMonteCarlo = False
        
        #----External Variables -------------------------------------------------------------#
        
        
        
        
        self.externVariables = { 'createSample'             : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},
                                 'createEvaluationXmlFiles' : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},
                                 'simulateEvaluations'      : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},
                                 'preProcessData'           : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},
                                 'postProcessing'           : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},
                                 'localEvaluation'          : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},               
                                 'multiprocessing'          : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},  
                                 'numberOfProcessors'       : {'type':'int',    'unitSI': None,   'strCases': None, 'multiVar': False},
                                 'evaluationNumbers'        : {'type':'int',    'unitSI': None,   'strCases': None, 'multiVar': True},
                                 'runPolynomialChaos'       : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},
                                 'polynomialOrders'         : {'type':'int',    'unitSI': None,   'strCases': None, 'multiVar': True}, 
                                 'sampleMethod'             : {'type':'str',    'unitSI': None,   'strCases': ['K','R','L','S','H','M','C','NC','G','RG'], 'multiVar': False},
                                 'runMonteCarlo'            : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False}
                                 }
        
        self.externXmlAttributes  = []
        self.externXmlElements    = self.externVariables.keys()
        # optional set order to elements for germans only :)
        self.externXmlElements    = ['createSample',
                                    'createEvaluationXmlFiles',
                                    'simulateEvaluations',
                                    'preProcessData' ,
                                    'postProcessing',
                                    'localEvaluation',                 
                                    'multiprocessing',     
                                    'numberOfProcessors',   
                                    'evaluationNumbers',
                                    'runPolynomialChaos',
                                    'polynomialOrders', 
                                    'sampleMethod',
                                    'runMonteCarlo']
        
        
        self.readDataFromXmlNode(xmlNode)
        
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
            newData[attribute] = self.loadVariablesConversion(attribute, xmlNode.attrib[attribute], '')
                                            
        for externXmlElement in self.externXmlElements: 
            print externXmlElement
            # check if externXmlElement exists in the current node
            try: element = xmlNode.findall(''.join(['.//',externXmlElement]))[0]
            except: self.loadingErrorMessageVariableError(externXmlElement, 'vessel', 'Id')
            # get variable value     
            try: variableValueStr = element.text
            except: self.loadingErrorMessageValueError(externXmlElement, 'vessel', 'Id')
            # get unit
            if 'unit' in element.attrib: variableUnit = element.attrib['unit']
            else: variableUnit = None 
            # save converted XML-value
            newData[externXmlElement] = self.loadVariablesConversion(externXmlElement, variableValueStr, variableUnit)
     
        print newData
     
    def loadVariablesConversion(self, variable, variableValueStr, variableUnit, unit = 'unitSI'):
        '''
        checks the element.text string and
        evaluates it corresponding to the definition
        of the variable in the constants variablesDict 
        
        return converted evaluated value of variable
        '''
        
        sys.path.append('/'.join([cur,'..','UtilityLib']))
        from constants import unitsDictSI as unitsDict
        
        multiVariable = False
        variableValue = 'notConvertable'
        convertError = []
        variableTypes = self.externVariables[variable]['type']
        # check if variable is a multiple variable (means list of variables)
        if self.externVariables[variable]['multiVar']:
            variableValueStrings = variableValueStr.replace(',',' ').split() 
            multiVariable = True
            variableValues = []
        else:
            variableValueStrings = [variableValueStr]
        
        # check if variable can have multiple types
        if ' ' in variableTypes:
            variableTypes = variableTypes.replace(',',' ').split()
        else: variableTypes = [variableTypes]
        variableTypes.sort()
        
        # start conversion loop over variable and types
        for variableValueString in variableValueStrings:
            for variableType in variableTypes:
                            
                if variableType in ['float','int']:
                    try: variableValue = float(eval(variableValueString))
                    except: convertError.append('float') 
                    
                    if self.externVariables[variable][unit]:
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
                    if variableValueString in self.externVariables[variable]['strCases']: variableValue = variableValueString
                    elif self.externVariables[variable]['strCases'][0] == 'anything': variableValue = variableValueString
                    else: convertError.append(''.join(['str == ',str(self.externVariables[variable]['strCases'])]))
                
                elif variableType in ['None']:
                    if variableValueString == 'None' or variableValueString == '' or variableValueString == None:
                        variableValue = None
                    else: convertError.append('None')
                        
            if variableValue == 'notConvertable':
                raise ValueError("""ERROR: moduleXML.loadVariablesConversion():
                      Cannot convert given value "{}" of variable "{}"
                      to {}!
                      Check if it is of type {}, system exit!""".format(variableValueString,variable,convertError,variableTypes))
            if multiVariable == False: return variableValue
            else: variableValues.append(variableValue)
        
        return variableValues   
        
    def loadingErrorMessageValueError(self, variableName, element, elementName):
        variableDict =  self.externVariables[variableName]
        raise ValueError("""ERROR loadNetworkFromXML():
              value for variable <<{}>> of {} {} is not defined.
              (Hint:{}) , system exit!""".format(variableName, element, elementName, variableDict))

    
    def loadingErrorMessageVariableError(self, variableName, element, elementName):
        if variableName in self.externVariables:
            variableDict = self.externVariables[variableName]
        else: variableDict = "No entry defined for This element"
        raise ValueError("""ERROR loadNetworkFromXML():
              variable "{}" of {} {} is not defined.
              (Hint:{}) , system exit!""".format(variableName, element, elementName, variableDict))
        
        
    def update(self, Dict):
        '''
        updates the class data using a dictionary in from of 
        dataDict = {'variableName': value}
        '''
        for key,value in Dict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except:
                print "ERROR VpcConfiguration.updateData - Wrong key: {}, could not update varibale".format(key)
    
