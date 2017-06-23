########################################################################################
#                            Vascular Polynomial Chaos 0.3
########################################################################################
##
# created by Vinzenz Eck vinzenz.eck@mytum.de
# uses polynomial Chaos toolbox from Jonathan Feinberg, Simula Center Oslo
##

from __future__ import print_function, absolute_import
from future.utils import iteritems, iterkeys, viewkeys, viewitems, itervalues, viewvalues
from builtins import range
from builtins import input as input3
import os

import starfish.VascularPolynomialChaosLib.classDistributionManager as cDistMng
import starfish.VascularPolynomialChaosLib.moduleFilePathHandlerVPC as mFPH_VPC
import starfish.VascularPolynomialChaosLib.moduleBatchSimulationManager as mBSM
import starfish.VascularPolynomialChaosLib.classUqsaCase as cUqsaCase

import starfish.UtilityLib.moduleStartUp as mStartUp
import starfish.UtilityLib.moduleXML as mXML

def uncertaintyPropagation():
    '''
    Perform vascular polynomial chaos
    or MonteCarlo analysis for STARFiSh simulation case
    # steps
    # 1. load vpc case and configuration

    # 2. create distributions

    # 3. add dependentCase if existent

    # 4. create samples

    # 5. evaluate model / on local machine or on server

    # 6. postprocess evaluated data, peak finding etc

    # 7. create Orthogonal polynomials

    # 8. calculate polynomial chaos expansion

    # 9. uncertainty quantfication, sensitivity analysis
    '''
    print("")
    print('==============================================')
    print('#        VascularPolynomialChaos_v0.3        #')
    print('==============================================')
    # steps
    # 1. load vpc case and configuration
    optionsDict = mStartUp.parseOptions(['f','n'],vascularPolynomialChaos=True)
    networkName = optionsDict['networkName']
    dataNumber  = optionsDict['dataNumber']
    
    # 1.1 load configuration and locations of interest      
    uqsaCase = cUqsaCase.UqsaCase() 
    uqsaCaseFile = mFPH_VPC.getFilePath('uqsaCaseXmlFile', networkName, dataNumber, 'read')
    uqsaCase.loadXMLFile(uqsaCaseFile)
    uqsaCase.initialize(networkName,dataNumber)
    
    # 1.2 load vascular network file polynomial chaos
    vpcNetworkXmlFile = mFPH_VPC.getFilePath('vpcNetworkXmlFile', networkName, dataNumber, 'read')
    vascularNetwork = mXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = vpcNetworkXmlFile)
    run_uqsa_case(uqsaCase, vascularNetwork)
    
def run_uqsa_case(uqsaCase, vascularNetwork):

    # 1.3 initialized defined random inputs
    vascularNetwork.randomInputManager.initialize(vascularNetwork)
    assert len(vascularNetwork.randomInputManager.randomInputs.keys()) != 0, "VascularPolynomialChaos_v0.3: no random inputs defined!"
    vascularNetwork.randomInputManager.printOutInfo()
 
    # 2. create distributions
    distributionManager = cDistMng.DistributionManagerChaospy(vascularNetwork.randomInputManager.randomInputsExtDist)
    distributionManager.createRandomVariables()

    # 3. add dependentCase if existent        
    if uqsaCase.sampleManager.dependentCase == True:
        # this enables dependentCase in Distribution Manager
        distributionManager.createDependentDistribution(vascularNetwork.randomInputManager.correlationMatrix)
        
    # 4. create or load samples
    uqsaCase.aquireSamples(distributionManager, vascularNetwork.randomInputManager.randomInputsExtDist)
    # 5. evaluate model / on local machine or on server
    # 5.1 create evaluation case file list
    uqsaCase.createEvaluationCaseFiles()
    
    # 5.2 save/create simulation xml files
    if uqsaCase.createEvaluationXmlFiles == True:
        
        for sampleIndex in range(uqsaCase.sampleManager.currentSampleSize):
            # update network with current evaluation number
            # get sample
            sample = uqsaCase.sampleManager.getSample(sampleIndex)            
            # pass realisation of Z to network
            distributionManager.passRealisation(sample,sampleIndex)
            # save evaluation file
            networkXmlFileSave = uqsaCase.evaluationCaseFiles[sampleIndex]['networkXmlFileSave']
            mXML.writeNetworkToXML(vascularNetwork,  dataNumber = vascularNetwork.dataNumber, networkXmlFile= networkXmlFileSave)
        
        # save evaluation to log file
        evaluationLogFile = mFPH_VPC.getFilePath('evaluationLogFile', vascularNetwork.name, vascularNetwork.dataNumber, mode = "write", caseName = uqsaCase.sampleManager.samplingMethod)
        vascularNetwork.randomInputManager.saveRealisationLog(evaluationLogFile, vascularNetwork.name, vascularNetwork.dataNumber, caseName = uqsaCase.sampleManager.samplingMethod)
        
    # 5.3 run evaluation simulations
    if uqsaCase.simulateEvaluations == True:
        batchFileList = uqsaCase.getSimulationBatchFileList()
        if uqsaCase.localEvaluation == True:
            if uqsaCase.multiprocessing == False:
                mBSM.runBatchAsSingleProcess(batchFileList, quiet = True)
            else:
                mBSM.runBatchAsMultiprocessing(batchFileList, uqsaCase.numberOfProcessors , quiet = True)
        else:
            # TODO: server simulations not implemented yet
            raise NotImplementedError("server simulations not implemented yet")
    
    # 6. process quantity of interest
    uqsaCase.preprocessSolutionData()
    
    # 7. uncertainty quantification and sensitivty analysis
    uqsaCase.quantifyUncertaintyAndAnalyseSensitivtiy(distributionManager)
    # 8. plotting of variables

if __name__ == '__main__':
    uncertaintyPropagation()
