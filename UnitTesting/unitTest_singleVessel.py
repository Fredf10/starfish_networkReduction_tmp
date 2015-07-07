


#TODO: 1. fix imports

def test_singleVessel():
    
    
    ## load and run a new simulation of reference
    
    #TODO: 2. define paths and names
    networkName,dataNumber,networkXmlFileLoad,networkXmlFileSave,pathSolutionDataFilename 
    
    vascularNetworkNew = moduleXML.loadNetworkFromXML(networkName, dataNumber, networkXmlFile = networkXmlFileLoad, pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkNew.quiet = True
    flowSolver = FlowSolver(vascularNetworkNew, quiet=True)
    flowSolver.solve()
    #vascularNetworkRef.saveSolutionData()
    #moduleXML.writeNetworkToXML(vascularNetworkNew, dataNumber, networkXmlFileSave)
    del flowSolver
    gc.collect()
        
    # load reference solution
    #TODO: define paths and names
    vascularNetworkRef = moduleXML.loadNetworkFromXML(networkName, dataNumber = dataNumber, networkXmlFile = networkXmlFile, pathSolutionDataFilename = pathSolutionDataFilename)
    vascularNetworkRef.linkSolutionData()
    
    
    # compare 
    #TODO: 3. write comparison as in Masterthesis of Fredirk for pressure and flow in all vessels
    
    # uses assert to check solution within tolerances