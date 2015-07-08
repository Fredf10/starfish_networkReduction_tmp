import unittest
import sys, os, gc
import matplotlib.pyplot as plt
import numpy as np

cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur + '/../')

import UtilityLib.moduleXML as mXML
import SolverLib.class1DflowSolver as c1dFS


def test_singleVessel():

    # Set which values should be tested, and the threshold for mean square error
    testDict = {"pressure": 10.0,
                "area": 10.0,
                "flow": 10.0}

    ## load and run a new simulation of reference
    networkName = "singleVessel"
    dataNumber = "999"
    dataNumberSol = "012"
    networkXmlFileLoad = cur + "/singleVessel/singleVessel.xml"
    networkXmlSolutionLoad = cur + "/singleVessel/singleVessel_SolutionData_012.xml"
    networkXmlFileSave = cur + "/tmp.xml"
    pathSolutionDataFilename = cur + "/tmpSol.hdf5"


    vascularNetworkNew = mXML.loadNetworkFromXML(networkName,
                                                 dataNumber,
                                                 networkXmlFile = networkXmlFileLoad,
                                                 pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkNew.quiet = True
    flowSolver = c1dFS.FlowSolver(vascularNetworkNew, quiet=True)
    flowSolver.solve()
    vascularNetworkNew.saveSolutionData()
    mXML.writeNetworkToXML(vascularNetworkNew, dataNumber, networkXmlFileSave)
    del flowSolver
    gc.collect()

    # link simulation data
    vascularNetworkNew.linkSolutionData()

    # load reference data and link it
    vascularNetworkRef = mXML.loadNetworkFromXML(networkName,
                                                      dataNumber = dataNumberSol,
                                                      networkXmlFile = networkXmlSolutionLoad,
                                                      pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkRef.linkSolutionData()

    for vesselId in vascularNetworkNew.vessels:
        print "testing vessel nr. {} \n".format(vesselId)

        dataDictNew = vascularNetworkNew.getSolutionData( vesselId, ['Pressure'], vascularNetworkNew.simulationTime, [0.0, 0.2])
        dataDictRef = vascularNetworkRef.getSolutionData( vesselId, ['Pressure'], vascularNetworkRef.simulationTime, [0.0, 0.2])

#   lag dict = {"pressure": tolerance}
#   "pressure", "area", "flow"
    # root mean square error

    print dataDictNew
    print dataDictRef


if __name__ == "__main__":
    test_singleVessel()