# In the real world
# 1. run singleVesselUQSA.xml with standard simulator?
# 2. run VascularPolynomialChaos.py  on the singleVesselUQSA.xml to produce vpcCase
# 3. run VascularPolynomialChaos.py on the resulting vpc.xml to evaluate all

from __future__ import print_function
import sys, os
import shutil
import numpy as np
import math
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import starfish.VascularPolynomialChaosLib.classDistributionManager as cDistMng
import starfish.VascularPolynomialChaosLib.moduleFilePathHandlerVPC as mFPH_VPC
import starfish.VascularPolynomialChaosLib.moduleBatchSimulationManager as mBSM
import starfish.VascularPolynomialChaosLib.classUqsaCase as cUqsaCase
import starfish.UtilityLib.moduleXML as mXML


def test_singleVessel():

    # Set which values should be tested, and the threshold for mean square error
    testDict = {"Pressure": 10.0,
                "Area": 10.0,
                "Flow": 10.0}

    # load and run a new simulation of reference
    networkName = "singleVesselUQSA"
    dataNumber = "tst"
    dataNumberSol = "tst"


   # 1.1 load configuration and locations of interest      
    uqsaCase = cUqsaCase.UqsaCase() #cConfigUQSA.ConfigurationUQSA()
    uqsaCaseFileSol = os.path.join(*[cur, "singleVesselUQSA", "vascularPolynomialChaos_tst",  "singleVesselUQSA_uqsaCase_tst.xml"])
    uqsaCaseFile = mFPH_VPC.getFilePath('uqsaCaseXmlFile', networkName, dataNumber, 'write')
    path_to_delete = os.path.split(uqsaCaseFile)[0]
    if not os.path.exists(os.path.split(uqsaCaseFile)[0]):
        os.makedirs(os.path.split(uqsaCaseFile)[0])
    shutil.copy(uqsaCaseFileSol, uqsaCaseFile)
    uqsaCase.loadXMLFile(uqsaCaseFile)
    uqsaCase.initialize(networkName,dataNumber)

    
    # 1.2 load vascular network file polynomial chaos
    vpcNetworkXmlFileSol = os.path.join(*[cur, "singleVesselUQSA", "vascularPolynomialChaos_tst",  "singleVesselUQSA_vpc_tst.xml"])
    vpcNetworkXmlFile = mFPH_VPC.getFilePath('vpcNetworkXmlFile', networkName, dataNumber, 'write')
    if not os.path.exists(os.path.split(vpcNetworkXmlFileSol)[0]):
        os.makedirs(os.path.split(vpcNetworkXmlFileSol)[0])
    shutil.copy(vpcNetworkXmlFileSol, vpcNetworkXmlFile)

    vascularNetwork = mXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = vpcNetworkXmlFile)
    # TODO: steps 1.1 and 1.2 should be entirely external to the UQSA core

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
        
        for sampleIndex in xrange(uqsaCase.sampleManager.currentSampleSize):
            # update network with current evaluation number
            # get sample
            sample = uqsaCase.sampleManager.getSample(sampleIndex)            
            # pass realisation of Z to network
            distributionManager.passRealisation(sample,sampleIndex)
            # save evaluation file
            networkXmlFileSave = uqsaCase.evaluationCaseFiles[sampleIndex]['networkXmlFileSave']
            mXML.writeNetworkToXML(vascularNetwork,  dataNumber = dataNumber, networkXmlFile= networkXmlFileSave)
        
        # save evaluation to log file
        evaluationLogFile = mFPH_VPC.getFilePath('evaluationLogFile', networkName, dataNumber, mode = "write", caseName = uqsaCase.sampleManager.samplingMethod)

        vascularNetwork.randomInputManager.saveRealisationLog(evaluationLogFile, networkName, dataNumber, caseName = uqsaCase.sampleManager.samplingMethod)
        
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

    #TODO: add tests to compare to saved solution dataNumber
    shutil.rmtree(path_to_delete)

  

if __name__ == "__main__":
    test_singleVessel()

