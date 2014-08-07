import csv
from numpy import sqrt,pi

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append(''.join([cur,'/../NetworkLib/']))
from classBoundaryConditions import *

from constants import newestNetworkXml as nxml
from constants import variablesDict

from moduleXML import loadVariablesConversion


def writeVesselDataToCSV(vessels,delimiter=';',filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    '''
    
    '''
    if filename == None:
        print "ERROR: moduleCSV.writeVesselDataToCSV(): no filename passed to csv-writer!"
        return
        
    networkDirectory = filename.split('.')[0]
    if not os.path.exists(str(networkPath+networkDirectory)):
        os.makedirs(str(networkPath+networkDirectory))
    
    # find all tags needed
    polyChaosTags = {}
    for vessel in vessels.itervalues():
        pcTags = vessel.getVariableValue('polyChaos').keys()
        for pcTag in pcTags: 
            if pcTag not in polyChaosTags.keys(): polyChaosTags[pcTag] = len(vessel.getVariableValue('polyChaos')[pcTag])     
    tags = []
    for tag in nxml.vesselAttributes: tags.append(tag)
    for vesselElement in nxml.vesselElements:
        if vesselElement == 'compliance':
            for specificCompElements in nxml.vesselElementReference[vesselElement].values():
                for tag in specificCompElements:
                    if tag not in tags:
                        tags.append(tag)
                        if tag in polyChaosTags.keys():
                            for count in range(polyChaosTags[tag]): tags.append(''.join([tag,'-pC',str(int(count)+1)]))
        else: 
            for tag in nxml.vesselElementReference[vesselElement]:
                tags.append(tag)
                if tag in polyChaosTags.keys():
                    for count in range(polyChaosTags[tag]): tags.append(''.join([tag,'-pC',str(int(count)+1)]))
    
    writer = csv.DictWriter(open(networkPath+networkDirectory+'/'+filename,'wb'),tags,delimiter=delimiter)
    
    # write first row == tags
    firstRow = {}
    for item in tags: firstRow[item] = item
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
    
    # write all data
    writer.writerow(unitRow)
    data = [] 
    for vessel in vessels.itervalues():
        vesselDict = {}
        for tag in tags:
            if '-pC' in tag:
                try:
                    variable,number = tag.split('-pC')
                    vesselDict[tag] = vessel.getVariableValue('polyChaos')[variable][int(number)-1]
                except: pass            
            else: 
                vesselDict[tag] = vessel.getVariableValue(tag)
        data.append(vesselDict)
    writer.writerows(data)

def readVesselDataFromCSV(filename= None,delimiter=';',networkPath = str(cur+"/../NetworkFiles/")):
    '''
    
    '''
    if filename == None:
        print "ERROR: moduleCSV.readVesselDataFromCSV(): no filename passed to csv-writer!"
        return {}
               
    networkDirectory = filename.split('.')[0]
    if not os.path.exists(str(networkPath+networkDirectory)):
        print ""
        print 'ERROR: moduleCSV.readVesselDataFromCSV(): file not in the correct directory / directory does not exists'
        return {}
    
    # load data    
    reader = csv.DictReader(open(networkPath+networkDirectory+'/'+filename,'rb'),delimiter=delimiter)
    # hash data with in dictionary and separate units
    columUnits = {}
    vesselData = {}
    for row in reader:
        Id = row.pop('Id')
        if Id == 'unit': columUnits = row
        else:
            Id = int(Id)
            vesselData[Id] = row
        
    variablesToDiscard = []
    for Id,data in vesselData.iteritems():
        polyChaos = {}
        for variable,variableValueStr in data.iteritems():
            # check if value is defined
            if variableValueStr not in ['', None]:
                #find units 
                if '#' in columUnits[variable]: #  '#' in variable or 
                    nothing,variableUnit = columUnits[variable].split('#',1)
                # check for polyChaos variables
                if '-pC' not in variable:
                    # convert variables to corret unit and type
                    data[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
                else:
                    variable,number = variable.split('-pC')
                    if variable not in polyChaos.keys():
                        polyChaos[variable] = variableValueStr
                    else:
                        polyChaos[variable] = ' '.join([polyChaos[variable],variableValueStr])   
                                             
            else: variablesToDiscard.append([Id,variable]) # find out variables which have no values
        # convert polynomial chaos variables to corret unit and type
        for variable,variableValueStr in polyChaos.iteritems():
            variableUnit = columUnits[variable].split('#',1)
            polyChaos[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit, polychaos = True)
        data['polyChaos'] = polyChaos
        for variable in data.iterkeys():
            if '-pC' in variable: variablesToDiscard.append([Id,variable])
            
    # remove variables which have no values 
    for Id,variableToDiscard in variablesToDiscard:
        del vesselData[Id][variableToDiscard]
    
    return {'vesselData':vesselData}    



def writeBCToCSV(boundaryConditionDict,boundaryConditionPolyChaos,delimiter=';',filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    '''
    
    '''
    if filename == None:
        print "ERROR: moduleCSV.writeBCToCSV(): no filename passed to csv-writer!"
        return {}
    
    networkDirectory = filename.split('BC')[0]
    if not os.path.exists(str(networkPath+networkDirectory)):
        os.makedirs(str(networkPath+networkDirectory))
        
    # find all polychaos tags
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
    
    writer = csv.DictWriter(open(networkPath+networkDirectory+'/'+filename,'wb'),tags,delimiter=delimiter)
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
    
    
def readBCFromCSV(delimiter=';',filename = None, networkPath = str(cur+"/../NetworkFiles/")):
    '''
    
    
    '''
    if filename == None:
        print "ERROR: moduleCSV.readBCFromCSV() no filename passed to csv-writer!"
        return {}
        
    networkDirectory = filename.split('BC')[0]
    
    if os.path.isfile(networkPath+networkDirectory+'/'+filename) == False:
        print "ERROR: moduleCSV.readBCFromCSV() file does not exists!"
        return {}
    
    reader = csv.DictReader(open(networkPath+networkDirectory+'/'+filename,'rb'),delimiter=delimiter)
    
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
                    path = ''.join([networkDirectory,'/'])
                    if path not in variableValueStr:  variableValueStr = variableValueStr.join([path,''])
                boundaryDataDict[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit)
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
                    polyChaos[variable] = loadVariablesConversion(variable, variableValueStr, variableUnit, polychaos = True)
                except: pass          
            if Id not in boundaryConditionPolyChaos.keys(): 
                boundaryConditionPolyChaos[Id] =[polyChaos]
            else: boundaryConditionPolyChaos[Id].append(polyChaos)           
            
        boundaryInstance.update(boundaryDataDict)
        
        if Id not in BCconditionData.keys(): BCconditionData[Id] = [boundaryInstance]
        else: BCconditionData[Id].append(boundaryInstance)
    
    return BCconditionData,boundaryConditionPolyChaos
  