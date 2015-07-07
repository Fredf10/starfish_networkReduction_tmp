import csv
from numpy import sqrt,pi

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(''.join([cur,'/../']))

#sys.path.append(''.join([cur,'/../NetworkLib/']))
from NetworkLib.classBoundaryConditions import *

from constants import newestNetworkXml as nxml
from constants import variablesDict

import moduleXML as mXML
import moduleFilePathHandler as mFPH

"""
This is moduleCSV, should be imported as mCSV

Some other information about the module 
"""

def writeVesselDataToCSV(networkName, vessels, delimiter=';'):
    '''
    Functions writes vessel data to *.csv file inclusive polynomial chaos definitions
    
    input:
        networkName <string>
        vessels     <dict>    := vessels dict of class vascularNetwork {vesselId : vesselInstance} 
        delimiter   <string>  (default = ';')
    
    '''
        
    # find all tags needed # TODO: read write polynomial chaos variables
#     polyChaosTags = {}
#     for vessel in vessels.itervalues():
#         pcTags = vessel.getVariableValue('polyChaos').keys()
#         for pcTag in pcTags: 
#             if pcTag not in polyChaosTags.keys(): polyChaosTags[pcTag] = len(vessel.getVariableValue('polyChaos')[pcTag])     
    tags = []
    for tag in nxml.vesselAttributes: tags.append(tag)
    for vesselElement in nxml.vesselElements:
        if vesselElement == 'compliance':
            for specificCompElements in nxml.vesselElementReference[vesselElement].values():
                for tag in specificCompElements:
                    if tag not in tags:
                        tags.append(tag)
#                         if tag in polyChaosTags.keys():
#                             for count in range(polyChaosTags[tag]): tags.append(''.join([tag,'-pC',str(int(count)+1)]))
        else: 
            for tag in nxml.vesselElementReference[vesselElement]:
                tags.append(tag)
#                 if tag in polyChaosTags.keys():
#                     for count in range(polyChaosTags[tag]): tags.append(''.join([tag,'-pC',str(int(count)+1)]))
    
    ## openFile and create writer
    vesselCSVFile = mFPH.getFilePath('vesselCSVFile', networkName, 'xxx', 'write')
    writer = csv.DictWriter(open(vesselCSVFile,'wb'),tags,delimiter=delimiter)
    
    # write first row == tags
    firstRow = {}
    for item in tags: firstRow[item] = item
    writer.writerow(firstRow)
    
    # write unit row
    unitRow = {}
    for tag in tags:
        try:
#             if '-pC' in tag: 
#                 tagUnit = tag.split('-pC')[0]
#                 unitRow[tag] = ''.join(['#',variablesDict[tagUnit]['unitSI']])
#             else:
                unitRow[tag] = ''.join(['#',variablesDict[tag]['unitSI']])
        except: unitRow[tag] = ''
    unitRow['Id'] = 'unit'
    
    # write all data
    writer.writerow(unitRow)
    data = [] 
    for vessel in vessels.itervalues():
        vesselDict = {}
        for tag in tags:
#             if '-pC' in tag:
#                 try:
#                     variable,number = tag.split('-pC')
#                     vesselDict[tag] = vessel.getVariableValue('polyChaos')[variable][int(number)-1]
#                 except: pass            
#             else: 
                vesselDict[tag] = vessel.getVariableValue(tag)
        data.append(vesselDict)
    writer.writerows(data)

def readVesselDataFromCSV(networkName, delimiter=';'):
    '''
    Functions loads vessel data from *.csv file inclusive polynomial chaos definitions
    
    input:
        networkName 
        delimiter   (default = ';')
    
    return:
            dict := {'vesselData': vesselData} which is used by vascular network to
                    update its vessel data with the function vascularNetwork.updateNetwork(dataDict)
    '''
        
    vesselCSVFile = mFPH.getFilePath('vesselCSVFile', networkName, 'xxx', 'read', exception = 'Warning')
    
    # load data    
    reader = csv.DictReader(open(vesselCSVFile,'rb'),delimiter=delimiter)
    # hash data with in dictionary and separate units
    columUnits = {}
    vesselData = {}
    for row in reader:
        Id = row.pop('Id')
        if Id == 'unit': columUnits = row
        else:
            Id = int(Id)
            vesselData[Id] = row
        
    # TODO: read write polynomial chaos variables
    
    variablesToDiscard = []
    for Id,data in vesselData.iteritems():
#         polyChaos = {}
        for variable,variableValueStr in data.iteritems():
            # check if value is defined
            if variableValueStr not in ['', None]:
                #find units 
                if '#' in columUnits[variable]: #  '#' in variable or 
                    nothing,variableUnit = columUnits[variable].split('#',1)
#                 # check for polyChaos variables
#                 if '-pC' not in variable:
                # convert variables to corret unit and type
                data[variable] = mXML.loadVariablesConversion(variable, variableValueStr, variableUnit)
#                 else:
#                     variable,number = variable.split('-pC')
#                     if variable not in polyChaos.keys():
#                         polyChaos[variable] = variableValueStr
#                     else:
#                         polyChaos[variable] = ' '.join([polyChaos[variable],variableValueStr])   
                                             
            else: variablesToDiscard.append([Id,variable]) # find out variables which have no values
        # convert polynomial chaos variables to corret unit and type
#         for variable,variableValueStr in polyChaos.iteritems():
#             variableUnit = columUnits[variable].split('#',1)
#             polyChaos[variable] = mXML.loadVariablesConversion(variable, variableValueStr, variableUnit, polychaos = True)
#         data['polyChaos'] = polyChaos
#         for variable in data.iterkeys():
#             if '-pC' in variable: variablesToDiscard.append([Id,variable])
            
    # remove variables which have no values 
    for Id,variableToDiscard in variablesToDiscard:
        del vesselData[Id][variableToDiscard]
    
    return {'vesselData':vesselData}    



def writeBCToCSV(networkName, boundaryConditionDict, boundaryConditionPolyChaos, delimiter=';'):
    '''
    Functions writes boundaryCondition data to *.csv file inclusive polynomial chaos definitions
    
    input:
        networkName <string>
        boundaryConditionDict      <dict>  := boundaryConditionDict dict of class VascularNetwork {vesselId : boundaryConditionInstance} 
        boundaryConditionPolyChaos <dict>  := boundaryConditionPolyChaos dict of class VascularNetwork 
        delimiter   <string>  (default = ';')
    
    '''
    # find all polychaos tags
    # TODO: read write polynomial chaos variables
    polyChaosTags = {}
    for id,bcPolyChaosList in boundaryConditionPolyChaos.iteritems():
        for bcPolyChaosDict in bcPolyChaosList:
            for variable,interval in bcPolyChaosDict.iteritems():
                if variable != 'name': polyChaosTags[variable] = len(interval)     
                
    # find all tags which are known for all boundary conditions
    tagsBCType1 = ['Id','boundaryType']
    tagsBCType2 = [] 
    for boundaryCondition,elementTags in nxml.boundaryConditionElements.iteritems():
        # None - bc for tag evaluation
        if 'None' not in boundaryCondition and '_' not in boundaryCondition:
            #check out type of BoundaryCondition:
            bcType = eval(nxml.bcTagsClassReferences[boundaryCondition])().getVariableValue('type')
            for elementTag in elementTags:
                if bcType == 1:
                    if elementTag not in tagsBCType1:
                        tagsBCType1.append(elementTag)
                        if elementTag in polyChaosTags.keys():
                            for count in range(polyChaosTags[elementTag]):
                                tagsBCType1.append(''.join([elementTag,'-pC',str(int(count)+1)]))
                elif bcType == 2:
                    if elementTag not in tagsBCType2:
                        tagsBCType2.append(elementTag)
                        if elementTag in polyChaosTags.keys():
                            for count in range(polyChaosTags[elementTag]):
                                tagsBCType2.append(''.join([elementTag,'-pC',str(int(count)+1)]))
        
    tagsBCType1.extend(tagsBCType2)
    tags = tagsBCType1
    
    boundaryCSVFile = mFPH.getFilePath('boundaryCSVFile', networkName, 'xxx', 'write')
    writer = csv.DictWriter(open(boundaryCSVFile,'wb'),tags,delimiter=delimiter)
    
    # write Tag row
    firstRow = {}
    for item in tags:
        firstRow[item] = item
    writer.writerow(firstRow)
    
    # write unit row
    unitRow = {}
    for tag in tags:
        try:
            if '-pC' in tag: 
                tagUnit = tag.split('-pC')[0]
                unitRow[tag] = ''.join(['#',variablesDict[tagUnit]['unitSI']])
            else:
                unitRow[tag] = ''.join(['#',variablesDict[tag]['unitSI']])
        except: unitRow[tag] = ''
    unitRow['Id'] = 'unit'
    writer.writerow(unitRow)
    
    ## fill out data of defined boundary conditions
    for Id,boundaryConditions in boundaryConditionDict.iteritems():
        for boundaryCondition in boundaryConditions:
            boundaryType  = boundaryCondition.getVariableValue('name')
            dataRow       = {}
            dataRow['Id'] = Id
            dataRow['boundaryType'] = boundaryType
            for variable in nxml.boundaryConditionElements[boundaryType]:
                dataRow[variable] = boundaryCondition.getVariableValue(variable)
            try:
                for bcPolyChaosDict in boundaryConditionPolyChaos[Id]:
                    bcPolyChaosDict.pop('name')
                    for variable,interval in bcPolyChaosDict.iteritems():
                        for count,value in enumerate(interval): 
                            dataRow[''.join([variable,'-pC',str(int(count)+1)])] = value
            except: pass               
                    
            writer.writerow(dataRow) 
    
    
def readBCFromCSV(networkName, delimiter=';'):
    '''
    Functions loads boundaryCondition data from \*.csv file inclusive polynomial chaos definitions

    Args:
        networkName (str): The name of the network
        delimiter (str): Delimiter (default = ';')

    Returns:
        VascularNetwork.boundaryConditionDict 
            A description of the VascNw.BCD instance returned
        
        VascularNetwork.boundaryConditionPolyChaos
            A description of the VascNw.BCPC instance returned
       
    '''
    
    boundaryCSVFile = mFPH.getFilePath('boundaryCSVFile', networkName, 'xxx', 'read', exception = 'Warning')
    reader = csv.DictReader(open(boundaryCSVFile,'rb'),delimiter=delimiter)
    
    # hash data with in dictionary and separate units
    columUnits = {}
    boundaryData = []
    for row in reader:
        if row['Id'] == 'unit': columUnits = row
        else:  boundaryData.append(row)
    
    boundaryConditionPolyChaos = {}       
    BCconditionData = {}
    for bcData in boundaryData:
        Id = int(bcData['Id'])
        # TODO: read write polynomial chaos variables
        # create class instance
        boundaryType = bcData['boundaryType']
        try: boundaryInstance = eval(nxml.bcTagsClassReferences[boundaryType])()
        except: 'ERROR moduleCSV.readBCFromCSV: boundaryType <<{}>> does not exist'.format(boundaryType)
        boundaryDataDict = {'name':boundaryType}
        polyChaos = {}
        for variable,variableValueStr in bcData.iteritems():
            if variable in nxml.boundaryConditionElements[boundaryType]: 
                try: variableUnit = columUnits[variable]
                except: variableUnit = None 
                # save converted XML-value
                if variable == 'filePathName':
                    path = ''.join([boundaryCSVFile])
                    if path not in variableValueStr:  variableValueStr = variableValueStr.join([path,''])
                print variableValueStr
                if variableValueStr != '':
                    try:
                        boundaryDataDict[variable] = mXML.loadVariablesConversion(variable, variableValueStr, variableUnit)
                    except:
                        pass
            if '-pC' in variable and variableValueStr != '':
                polyChaos['name'] = boundaryType
                variable,number = variable.split('-pC')
                if variable in polyChaos.keys():
                    polyChaos[variable] = ' '.join([polyChaos[variable],variableValueStr])
                else: polyChaos[variable] = variableValueStr
        # convert polynomial chaos variables to corret unit and type
        if polyChaos != {}:
            for variable,variableValueStr in polyChaos.iteritems():
                try:
                    variableUnit = columUnits[variable].split('#',1)
                    polyChaos[variable] = mXML.loadVariablesConversion(variable, variableValueStr, variableUnit, polychaos = True)
                except: pass          
            if Id not in boundaryConditionPolyChaos.keys(): 
                boundaryConditionPolyChaos[Id] =[polyChaos]
            else: boundaryConditionPolyChaos[Id].append(polyChaos)           
            
        boundaryInstance.update(boundaryDataDict)
        
        if Id not in BCconditionData.keys(): BCconditionData[Id] = [boundaryInstance]
        else: BCconditionData[Id].append(boundaryInstance)
    
    return BCconditionData, boundaryConditionPolyChaos
  
