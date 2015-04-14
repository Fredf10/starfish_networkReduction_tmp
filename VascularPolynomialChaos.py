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
import moduleBatchSimulationManager as mBSM
from classLocationOfInterestManager import LocationOfInterestManager

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
    # 3. add correlation if existent TODO:
    
    ## do the analysis for all defined polynomial orders:
    for polynomialOrder in vpcConfiguration.polynomialOrders:
        # 4. create samples
        if vpcConfiguration.createSample == True:      
            distributionManager.createSamples(networkName, dataNumber, vpcConfiguration.sampleMethod, expansionOrder = polynomialOrder)
        else:
            distributionManager.loadSamples(networkName, dataNumber, vpcConfiguration.sampleMethod, polynomialOrder)
        # 5. evaluate model / on local machine or on server
        # 5.1 create evaluation case file list
        evaluationCaseFiles = [] # list of [ [networkName,dataNumber,xml-filePath, hdf-filePath] for each evaluation
        for simulationIndex in xrange(distributionManager.samplesSize):
            vpcNetworkXmlEvaluationFile = mFPH_VPC.getFilePath('vpcEvaluationNetworkXmlFile', networkName, dataNumber, 'write',
                                                               gPCEmethod=vpcConfiguration.sampleMethod, gPCEorder= polynomialOrder, evaluationNumber=simulationIndex)
            vpcEvaluationSolutionDataFile = mFPH_VPC.getFilePath('vpcEvaluationSolutionDataFile', networkName, dataNumber, 'write',
                                                               gPCEmethod=vpcConfiguration.sampleMethod, gPCEorder= polynomialOrder, evaluationNumber=simulationIndex)
            evaluationCaseFiles.append([networkName,dataNumber,vpcNetworkXmlEvaluationFile,vpcEvaluationSolutionDataFile])        
        # 5.2 save/create xml files
        if vpcConfiguration.createEvaluationXmlFiles == True:
            for sampleIndex in xrange(distributionManager.samplesSize):
                distributionManager.passRealisation(sampleIndex)
                vpcNetworkXmlEvaluationFile = evaluationCaseFiles[sampleIndex][2]
                moduleXML.writeNetworkToXML(vascularNetwork,  dataNumber = dataNumber, networkXmlFile= vpcNetworkXmlEvaluationFile)
            vascularNetwork.randomInputManager.saveRealisationLog(networkName, dataNumber, vpcConfiguration.sampleMethod, polynomialOrder)
        # 5.3 run evaluation simulations
        if vpcConfiguration.simulateEvaluations == True:
            if vpcConfiguration.local == True:
                startIndex = 0
                endIndex   = int(distributionManager.samplesSize)
                newRange = vpcConfiguration.evaluationNumbers
                if len(newRange) == 2:
                    # check if the indices are avaliable
                    if all([i in xrange(distributionManager.samplesSize) for i in newRange]):
                        if newRange[0] < newRange[1]:
                            startIndex = newRange[0]
                            endIndex   = newRange[1]
                batchFileList = evaluationCaseFiles[startIndex:endIndex+1]
                if vpcConfiguration.multiprocessing == False:
                    mBSM.runBatchAsSingleProcess(batchFileList, quiet = True)
                else:
                    mBSM.runBatchAsMultiprocessing(batchFileList, vpcConfiguration.numberOfProcessors , quiet = True)
            else: print "server simulations not implemented yet";exit() # TODO: server simulations not implemented yet
        # 6. process quantity of interest
        
        ## defined query location and quantities to process
        quantitiesOfInterestToProcess = ['ForwardPressure', 'Pressure']
        queryLocation = 'vessel_1'
        xVals = 0.25
        # create locationOfInterestManager
        locationOfInterestManager = LocationOfInterestManager(distributionManager.samplesSize)
        # add location of interest
        locationOfInterestManager.addLocationOfInterest(queryLocation, quantitiesOfInterestToProcess, xVals)
        
        if vpcConfiguration.preProcessData == True:
            ## process the data of interest
            locationOfInterestManager.preprocessSolutionData(evaluationCaseFiles)
        
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
#def calculatePolynomialChaosExpansion(self):
        
                # 8. calculate polynomial chaos expansion    
#         print "    starting the polychaos polynomial calculation from polychaos simulation result!!"
#                 startTime = time.clock()
#                 
#                 polynomsT = []
#                 count = 0
#                 for sol in solutionInterpolated:
#                     print "      Polynomial calculation for point ", names[count]     
#                     count = count+1
#                     polyDict = {}
#                     
#                     for tag,data in sol.iteritems():
#                         # polynoms for the total pressure signal
#                         print "        polynoms for ",str(tag)
#                         if 'extrema' not in tag:
#                             
#                             polynomial = pc.fitter_lr(orthoPoly, sample.T, data)     
#                             polyDict[tag]= polynomial
#                             
#                         else:
#                             print sample.shape
#                             
#                             #print orthoPoly.dim
#                             #print orthoPoly.shape
#                             
#                             #print sample.ravel().shape
#                             #print len(data[0])
#                             polynomialTime = pc.fitter_lr(orthoPoly, sample.T, data[0])    
#                             polynomialAmp  = pc.fitter_lr(orthoPoly, sample.T, data[1])    
#                             
#                             extremaDict = {'Time':polynomialTime,'Amp':polynomialAmp}
#                             polyDict[tag]= extremaDict
#                             
#                     endTime = time.clock()
#                     polynomsT.append(polyDict)
#                     
#                 print "    It took %1.2f seconds to create the polynoms with order %d\n" % (endTime-startTime,order)
#                 
#                 if not os.path.exists(polynomPath):
#                     os.makedirs(polynomPath)    
#                 
#                 #save the polynom
#                 #create file with the network name in solution path directory
#                 FILE = open(polynomFile,"w")
#                 # store pickle and close file
#                 cPickle.dump(polynomsT,FILE,protocol=2)
#                 FILE.close()
#                 polynoms = [polynomsT]
        # 9. uncertainty quantfication, sensitivity analysis


if __name__ == '__main__':
    vascularPolyChaos()