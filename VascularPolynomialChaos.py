########################################################################################
#                            Vascular Polynomial Chaos 0.3
########################################################################################
## 
# created by Vinzenz Eck vinzenz.eck@mytum.de
# uses polynomial Chaos toolbox from Jonathan Feinberg, Simula Center Oslo
##

import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))

sys.path.append(cur+'/NetworkLib')
from classVascularNetwork import VascularNetwork

sys.path.append('/'.join([cur,'Solver']))
from class1DflowSolver import FlowSolver

sys.path.append('/'.join([cur,'VascularPolynomialChaosLib']))
from classVpcConfiguration import VpcConfiguration
from classDistributionManager import DistributionManagerChaospy
import moduleFilePathHandlerVPC as mFPH_VPC

sys.path.append('/'.join([cur,'UtilityLib']))
import moduleStartUp as mSU
import moduleXML

import chaospy as cp
import pprint

import numpy as np

import cPickle

def vascularPolyChaos():
    '''
    Perform vascular polynomial chaos 
    or MonteCarlo analysis for STARFiSh simulation case
    # steps
    # 1. load vpc case and configuration 
    
    # 2. create distributions
    
    # 3. add correlation if existent
    
    # 4. create samples
        
    # 5. evaluate model / on local machine or on server
    
    # 6. postprocess evaluated data, peak finding etc
    
    # 7. create Orthogonal polynomials
    
    # 8. calculate polynomial chaos expansion
    
    # 9. uncertainty quantfication, sensitivity analysis
    '''
    print ""
    print '=============================================='
    print '#        VascularPolynomialChaos_v0.3        #'
    print '=============================================='
    # steps
    # 1. load vpc case and configuration 
    optionsDict = mSU.parseOptions(['f','n'],vascularPolynomialChaos=True)
    networkName           = optionsDict['networkName']
    dataNumber            = optionsDict['dataNumber']
    # 1.1 load configuration  
    vpcConfiguration = VpcConfiguration(networkName,dataNumber)
    # 1.2 load vascular network file polynomial chaos
    vpcNetworkXmlFile = mFPH_VPC.getFilePath('vpcNetworkXmlFile', networkName, dataNumber, 'write')
    vascularNetwork = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = vpcNetworkXmlFile)
    # 1.3 print distributions
    vascularNetwork.randomInputManager.printOutInfo()
    if len(vascularNetwork.randomInputManager.randomInputs) == 0:
        print "VascularPolynomialChaos_v0.3: no random inputs defined!"
        exit()
        
    # 2. create distributions    
    distributionManager = DistributionManagerChaospy(vascularNetwork.randomInputManager.randomInputVector)
    distributionManager.createRandomVariables()
    
    # 3. add correlation if existent
        
    ## do the analysis for all defined polynomial orders:
    for polynomialOrder in vpcConfiguration.polynomialOrders:
        # 4. create samples
        if vpcConfiguration.createSample == True:      
            distributionManager.createSamples(vpcConfiguration.sampleMethod, expansionOrder = polynomialOrder)
        else:
            distributionManager.loadSamples()
            
        # 5. evaluate model / on local machine or on server
        # 5.1 save/create xml files
        simulationCaseFiles = {}
        # loop through number of samples
        for sampleIndex in xrange(distributionManager.samplesSize):
            distributionManager.passRealisation(sampleIndex)
            vpcNetworkXmlEvaluationFile = mFPH_VPC.getFilePath('vpcEvaluationNetworkXmlFile', networkName, dataNumber, 'write',
                                                               gPCEmethod=vpcConfiguration.sampleMethod, gPCEorder= polynomialOrder, evaluationNumber=sampleIndex)
            moduleXML.writeNetworkToXML(vascularNetwork,  dataNumber = dataNumber, networkXmlFile= vpcNetworkXmlEvaluationFile)
        
        # 5.2 run evaluations
        startIndex = 0
        endIndex = int(distributionManager.samplesSize)
        for simulationIndex in np.arange(startIndex,endIndex):
            vpcNetworkXmlEvaluationFile = mFPH_VPC.getFilePath('vpcEvaluationNetworkXmlFile', networkName, dataNumber, 'read',
                                                               gPCEmethod=vpcConfiguration.sampleMethod, gPCEorder= polynomialOrder, evaluationNumber=simulationIndex)
            vpcEvaluationSolutionDataFile = mFPH_VPC.getFilePath('vpcEvaluationSolutionDataFile', networkName, dataNumber, 'write',
                                                               gPCEmethod=vpcConfiguration.sampleMethod, gPCEorder= polynomialOrder, evaluationNumber=simulationIndex)
            moduleXML.writeNetworkToXML(vascularNetwork,  dataNumber = dataNumber, networkXmlFile= vpcNetworkXmlEvaluationFile)
            vascularNetworkTemp = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = vpcNetworkXmlEvaluationFile, pathSolutionDataFilename = vpcEvaluationSolutionDataFile)
            vascularNetworkTemp.quiet = True
            flowSolver = FlowSolver(vascularNetworkTemp)
            flowSolver.solve()
            vascularNetwork.saveSolutionData()
            del flowSolver
            gc.collect()
            
        # 6. postprocess evaluated data, peak finding etc
        
        # 7. create Orthogonal polynomials
        # ( TRUE == create and save, FALSE == load existing)#
#         if createOrthoPoly == True:
#             if not os.path.exists(polynomPathOrtho):
#                 os.makedirs(polynomPathOrtho)
#             
#             print " create and save orthogonal polynoms "
#             #create orthogonal polynom  
#             #orthoPoly = pc.orth_gs(order,distributions)
#             orthoPoly = pc.orth_ttr(order,distributions)
#             #orthoPoly = pc.orth_ttr(order,distributions)
#             #orthoPoly = pc.orth_chol(order,distributions)
#             #orthoPoly = pc.orth_svd(order,distributions)
#             
#             #save file
#             saveFile = open(orthoFile,"wb")       
#             cPickle.dump(orthoPoly,saveFile,protocol=2)
#             saveFile.close()
#             print ".. done"
#         else:
#             try:
#                 print " load orthogonal polynoms "
#                 loadFile = open(orthoFile,"rb")
#                 # load pickle
#                 orthoPoly = cPickle.load(loadFile)
#                 loadFile.close()
#                 print ".. done"
#             except:
#                 print 'File does not exits:'
#                 print orthoFile
#                 exit()
        # 8. calculate polynomial chaos expansion    
        
        # 9. uncertainty quantfication, sensitivity analysis


if __name__ == '__main__':
    vascularPolyChaos()