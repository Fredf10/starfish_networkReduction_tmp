########################################################################################
#                            Vascular Polynomial Chaos 0.2
########################################################################################
## 
# created by Vinzenz Eck vinzenz.eck@mytum.de
# uses polynomial Chaos toolbox from Jonathan Feinberg, Simula Center Oslo
##

#---------------------------------------------------------------------------------------#
#imports:

import pprint
import time 

import cPickle

import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append(cur+'/NetworkLib')
from classVascularNetwork import VascularNetwork

sys.path.append(cur+'/UtilityLib')
from moduleXML import loadNetworkFromXML
from moduleXML import writeNetworkToXML
from moduleXML import savePolyChaosXML
from moduleXML import loadPolyChaosXML
from processing import linearWaveSplitting
from processing import nonLinearWaveSplitting
from processing import minMaxFunction
from processing import calculateWaveShoulderPoint
from processing import calculateWaveShoulderPointTangent
from moduleStartUp import parseOptions
from moduleStartUp import chooseVPCconfigFile

sys.path.append(cur+'/Solver')
from class1DflowSolver import FlowSolver

import matplotlib.pyplot as plt   
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d 
from pprint import pprint as pp

from optparse import OptionParser

import gc
gc.enable()

# path to polynomialChaos toolbox
polynomialChaosPath = "/../../PolynomialChaos"
sys.path.append(cur+polynomialChaosPath)
try :import polychaos as pc
except: pass
polynomialChaosPath2 = "/../../"
sys.path.append(cur+polynomialChaosPath2)
try :import polychaos as pc
except: pass

import polychaos as pc

def main():
    print ""
    print '=============================================='
    print '#        VascularPolynomialChaos_v0.2        #'
    print '=============================================='
       
    optionsDict = parseOptions(['f','n'])
    
    networkName           = optionsDict['networkName']
    dataNumber            = optionsDict['dataNumber']
        
    if networkName == None or dataNumber == '999': 
        networkName,dataNumber = chooseVPCconfigFile(networkName)
     
    vpcConfigFilename = '_'.join([networkName,'vpcConfig',dataNumber])

    ##########################################################################################  
    ########################################################################################## 
    #----LOADING CONFIGURATION OF POLYNOMIAL CHAOS VARIABLES-----------------------------#    
    # network Name
          
    vPCconfiguration = loadPolyChaosXML(vpcConfigFilename)

    createDistributions = vPCconfiguration['createDistributions']
    createOrthoPoly  = vPCconfiguration['createOrthoPoly']
    createSample     = vPCconfiguration['createSample']
    runSimulations   = vPCconfiguration['runSimulations']
    calculateGPCE    = vPCconfiguration['plotMinMaxPoints']
    preProcessData   = vPCconfiguration['preProcessData']
    plotMinMaxPoints = vPCconfiguration['plotMinMaxPoints']
    postProcessing   = vPCconfiguration['postProcessing']
    plotMeanSTD      = vPCconfiguration['plotMeanSTD']
    plotPeaks        = vPCconfiguration['plotPeaks']
    
    polynomialOrders = vPCconfiguration['polynomialOrders']
    sampleMethod     = vPCconfiguration['sampleMethod']
    
    linearWaveSplit            = vPCconfiguration['linearWaveSplit']
    velocityProfileCoefficient = vPCconfiguration['velocityProfileCoefficient']
    
    polynomsToCalculate = vPCconfiguration['polynomsToCalculate'] 
    names               = vPCconfiguration['names']
    delta               = vPCconfiguration['delta']
    peaksToEvaluate     = vPCconfiguration['peaksToEvaluate']
        
    plotDirectory               = vPCconfiguration['plotDirectory']
    deterministicDataSetNumbers = vPCconfiguration['deterministicDataSetNumbers']
    polynomsToPlotOrder         = vPCconfiguration['polynomsToPlotOrder']
    
    plotMeanConfidenceInterval = vPCconfiguration['plotMeanConfidenceInterval'] 
    plotMeanConfidenceAlpha    = vPCconfiguration['plotMeanConfidenceAlpha']
    plotMeanSigmaInterval      = vPCconfiguration['plotMeanSigmaInterval']
    
    peakAnalysis                   = vPCconfiguration['peakAnalysis']
    plotPeaksConfidenceAlpha       = vPCconfiguration['plotPeaksConfidenceAlpha']
    plotPeaksAnalyticSensitivity   = vPCconfiguration['plotPeaksAnalyticSensitivity']
    plotPeaksMeanSTDBoxPlotsSingle = vPCconfiguration['plotPeaksMeanSTDBoxPlotsSingle']
        
        
    manualControl = False
    if manualControl:
        #----CONTROL VARIABLES: WHAT SHOULD BE DONE------------------------------------------------------#    
        ### 0. set up polychaos variables
        createDistributions = False  #(TRUE == create and save, FALSE == load existing)
        ##  0.1 orthogonal polynoms ( TRUE == create and save, FALSE == load existing)
        createOrthoPoly  = False
        ##  0.2 collocation method ( TRUE == create and save, FALSE == load existing)
        createSample     = False
        
        ### 1.step genrealized polynomial chaos evaluations + data storing
        runSimulations     = False
        
        ### 2.step Construct generalized polynomial chaos expansion and  pre process data
        calculateGPCE    = False
        #control variables
        # 2.1 is the data already preprocessed?
        preProcessData   = True
        # 2.2 create plots for min max points
        plotMinMaxPoints = True
            
        ### 3.step post processing - sensitivity analysis
        postProcessing  = True
        plotMeanSTD     = True
        plotPeaks       = True
        #plotSensitiviy  = False
        
        #----POLYNOMIAL CHAOS DEFINITIONS -------------------------------------------------------------#
        ####    
        
        #polynomialOrders of the polynomial chaos expansion || if more then one given they are processed consecutevely
        polynomialOrders = [3]
        # method of the spares grid collocation 
        sampleMethod = 'M'
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
        
        
        #----CHOICE OF INVESTIGATION POINTS FOR WHICH THE POLYNOMIAL CHAOS IS CALCULATED-----------#

        
        # [[vesselId, node],[vesselId, node],[vesselId, node]]
        polynomsToCalculate = [[0,0] ]#   , [0,99]  , [0,199],[1,99],[1,199]]
        names = ['A' ]#,'B','C','D','E']
        
        #### WAVE SPLITTING
        # linear or non linear wave splitting ( non linear =  default)
        linearWaveSplit = True
        # velocity profile coefficient (1 = linear, 2= peuiseulle ..)
        velocityProfileCoefficient = 2.0
        
        #### PEACK EVALUATION
        # fine tune min max function delta to find all the peaks
        delta ={'A': 
                {'Pressure' :1.0,'Flow':1.0e-8,
                'Pressure_f':1.0,'Flow_f':1.0e-8,
                'Pressure_b':1.0,'Flow_b':1.0e-8} ,
                'B':
                {'Pressure' :2.7,'Flow':0.5e-7,
                'Pressure_f':2.7,'Flow_f':0.5e-7,
                'Pressure_b':2.7,'Flow_b':0.5e-7} , 
                'C':
                {'Pressure' :2.7,'Flow':0.5e-7,
                'Pressure_f':2.7,'Flow_f':0.5e-7,
                'Pressure_b':2.7,'Flow_b':0.5e-7} ,
                'D':
                {'Pressure' :2.7,'Flow':0.5e-7,
                'Pressure_f':2.7,'Flow_f':0.5e-7,
                'Pressure_b':2.7,'Flow_b':0.5e-7} ,
                'E':
                {'Pressure' :2.7,'Flow':0.5e-7,
                'Pressure_f':2.7,'Flow_f':0.5e-7,
                'Pressure_b':2.7,'Flow_b':0.5e-7} }
        
        #check extrema in minMaxPointsand write number of peaks you want compare starting with 0 in the brackets
        peaksToEvaluate=    {'A': 
                            {'extremaPressure':[0],'extremaFlow':[0],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[0],'extremaFlow_b':[1]} ,
                            'B':
                            {'extremaPressure':[0],'extremaFlow':[0,1],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[0],'extremaFlow_b':[0]},  
                            'C':
                            {'extremaPressure':[0],'extremaFlow':[0],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[0],'extremaFlow_b':[0]} , 
                            'D':
                            {'extremaPressure':[0],'extremaFlow':[0],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[],'extremaFlow_b':[]} , 
                            'E':
                            {'extremaPressure':[0],'extremaFlow':[0],
                            'extremaPressure_f':[0],'extremaFlow_f':[0],
                            'extremaPressure_b':[],'extremaFlow_b':[]}  }
        ########################################################################################## 
        #--POSTPROCESSING: PLOTTING-------------------------------------------------------------------------------------#
        ### 
        
        plotDirectory = 'PolynomialChaos'
        ## load additional data for postprocessing
        # load additional simulation results from deterministic simulations
        # (stored in the same directory as the xml file
        deterministicDataSetNumbers = [999] # [200,199,201] # not in use
        # load additional polynomials with a different order then defined above 
        # if 0 then the current order is respected only
        polynomsToPlotOrder = [0]
        
        ##########################################################################################
        ## MEAN and STD plots
        # bool for confidence interval area for 100-plotMeanConfidenceAlpha % confidence
        plotMeanConfidenceInterval = False 
        plotMeanConfidenceAlpha    = 1.0 
        # bool for sigma interval area
        plotMeanSigmaInterval = False
            
        ##########################################################################################
        ## Plot of Peak Analysis
        
        # do peak analysis and save it if True // else: load peak file 
        peakAnalysis = False 
        # Confidence interval with 100-plotPeaksConfidenceAlpha % confidence
        plotPeaksConfidenceAlpha = 1.0
        
        ### analytic bar plots
        # compare occ.time and amplitude to analytic values (e.g. Sherwin Bifurcation) in bar plots
        plotPeaksAnalyticSensitivity = False
        
        ### Mean and STD Box plots 
        # save plots for each parameters in a seperate file if True // else: save one countious
        plotPeaksMeanSTDBoxPlotsSingle = False
    
    
    #---Additional hidden fine tuning parameters---
    # set start and end value for the time-line, if plotxTimeEnd = 'totalTime' simulation-endtime is taken
    totalTime = 'vascularNetwork.totalTime'
    
    plotxTimeStart = 0.0 
    plotxTimeEnd   = 'totalTime'
    plotXLablesTime = ['0','0.2','0.4','0.6','0.8']
    plotYLabelNumber = 5
        
    startEvaluation = 0 ## start the polychaos evaluations from number xxx
    
    plotFileType = '.pdf'#'.png'
    ## MEAN and STD plots
    latexUnits = False
    if latexUnits: plotMeanStdQuantities = {'Pressure':r'$[mmHg]$', 'Pressure_f':r'$[mmHg]$', 'Pressure_b':r'$[mmHg]$',
                                            'Flow':r'$[\frac{ml}{s}]$','Flow_f':r'$[\frac{ml}{s}]$','Flow_b':r'$[\frac{ml}{s}]$'}
    else: plotMeanStdQuantities = {'Pressure':'[mmHg]'} #, 'Pressure_f':'[mmHg]', 'Pressure_b':r'[mmHg]',
                                   #'Flow':'[ml/s]','Flow_f':r'[ml/s]','Flow_b':r'[ml/s]'}
    
    #fine tune y-axis for Mean  and STD plots if [0,0] use automatic-limit-caluclation
    limits= {'Pressure':  [ 0,0],'Flow':  [ 0,0],
             'Pressure_b':[ 0,0],'Flow_b':[ 0,0],
             'Pressure_f':[ 0,0],'Flow_f':[ 0,0]}
    
    ##########################################################################################  
    ########################################################################################## 
    #### INITIALIZE OF THE UNCERTAIN VARIABLES
    
    method = sampleMethod
    filename = str(networkName+'.xml')
    
    # load network from the path!
    vascularNetwork = loadNetworkFromXML(filename=filename)
    # update the connections of the network, as they prob weren't saved correct
    vascularNetwork.evaluateConnections()
    
    
    ########################################################################################## 
    ########################################################################################## 
    #---MAIN CONTROL LOOP FOR ALL ORDERS #
    
    ageTest = False
       
    for order in polynomialOrders:
        ### define the paths for saving data and plots
        print 'current order is ', order, ' choosen method is', method
        
        ## Polynomial Chaos
        saveDirectory = ''.join([cur,"/NetworkFiles/",vascularNetwork.name,'/',"VascularPolynomialChaos",'/','polyChaos_solution_',dataNumber,'_order_',str(order).zfill(2),'_method_',method,'/'])
        
        saveDirectory2 = ''.join([cur,'/../../../Daten/NetworkFiles/',vascularNetwork.name,'/',"VascularPolynomialChaos",'/','polyChaos_solution_',dataNumber,'_order_',str(order).zfill(2),'_method_',method,'/'])
        saveDirectoryPC = ''.join([saveDirectory2,'rawPolySimResults'])
        SolutionFileNamePath1 = ''.join([saveDirectoryPC,'/polyChaos_solution_'])
        numberOfEvaluationsPath = ''.join([saveDirectoryPC,'/polyChaos_solution_NumberOfEvaluations.pickel'])
        solutionInterpolatedPath = ''.join([saveDirectory,'preProcessedPolySimResults/'])
        solutionInterpolatedFile = ''.join([solutionInterpolatedPath,'polyChaos_solution_',dataNumber,'_order_',str(order).zfill(2),'_method_',method,'_preProcessed','.pickle'])
        
        distributionsDictPath = ''.join([saveDirectory,'distributions/']) 
        distributionsDictFile = ''.join([distributionsDictPath,'polyChaos_solution_',dataNumber,'_distributions','.pickle'])
        
        polynomPath = ''.join([saveDirectory,'polynoms/'])
        polynomFile = ''.join([polynomPath,'polyChaos_solution_',dataNumber,'_order_',str(order).zfill(2),'_method_',method,'_polynom','.pickle'])
        
        polynomPathOrtho = ''.join([saveDirectory,'orthoPolyAndSample/'])
        orthoFile = ''.join([polynomPathOrtho,'polyChaos_solution_',dataNumber,'_order_',str(order).zfill(2),'_method_',method,'_orthogonalPolynoms','.pickle'])
        sampleFile = ''.join([polynomPathOrtho,'polyChaos_solution_',dataNumber,'_order_',str(order).zfill(2),'_method_',method,'_sampling','.pickle'])
        
        ## PLOTS
        plotDirectoryPath         = ''.join([cur,"/NetworkFiles/",vascularNetwork.name,'/2dTimePlots/',plotDirectory])
        saveDirectoryPolySolPlots = ''.join([plotDirectoryPath,'/PQextrema/'])
        plotDirectoryMeanVar      = ''.join([plotDirectoryPath,'/Mean-STD/'])
        plotDirectorySensitivity  = ''.join([plotDirectoryPath,'/Sensitivity/'])
        plotDirectoryPeakP_bar    = ''.join([plotDirectoryPath,'/Peaks/AnalyticBarPlot/'])
        plotDirectoryPeakP_box    = ''.join([plotDirectoryPath,'/Peaks/MeanStdBox/'])
        plotDirectoryPeakP_Conf   = ''.join([plotDirectoryPath,'/Peaks/ConfidenceInterval/'])
        plotDirectoryPeakPLogs    = ''.join([plotDirectoryPath,'/Peaks/Logs/'])
        peakAnalysisFile          = ''.join([plotDirectoryPeakPLogs,'polyChaos_peakAnlaysis_',dataNumber,'_order_',str(order).zfill(2),'_method_',method,'.pickle'])
        
        
           
        ### create distribution for the random variables given with intervals
        if createDistributions == True:
            uncertainVarNames = []
            #---VESSEL DATA------------------------------------------------------------------------------------#
            
            ## global Variables 
            number = []
            numberVessel = []
            vesselIds = []
            varNames = []
            intervals = []
            distributions = pc.J()
            
            count = 0 
            ## random variables for vessel data
            for Id,vessel in vascularNetwork.vessels.iteritems():
                for varName,interval in vessel.polyChaos.iteritems():
                    # save data of the intervals, to know the distributions
                    if interval[0] != interval[1]:
                        vesselIds.append(Id)
                        varNames.append(varName)
                        uncertainVarNames.append(' '.join([varName,str(Id)]))
                        intervals.append(interval)
                        numberVessel.append(count)
                        number.append(count)
                        count += 1
                        # create uniform distribution with the interval 
                        distributions = pc.J(distributions,pc.Uniform(interval[0],interval[1]))  
            #----BOUNDARY CONDITIONS-----------------------------------------------------------------------------------#            
            ## global Variables 
            BCtype = []
            BCvesselID = []
            numberBC = []
            BCCount = 0
            ## random variables for boundaryConditions
            
            for vesselId,boundaryConditionList in vascularNetwork.boundaryConditionIntervals.iteritems():
                for boundaryCondition in boundaryConditionList:
                    for variableName,interval in boundaryCondition.iteritems():
                        if variableName != "name":
                            BCvesselID.append(vesselId)
                            BCtype.append(variableName)
                            uncertainVarNames.append(' '.join([boundaryCondition["name"],variableName,str(vesselId)]))
                            intervals.append(interval)         
                            distributions = pc.J(distributions,pc.Uniform(interval[0],interval[1]))         
                            numberBC.append(BCCount)
                            number.append(count)
                            count += 1
                            BCCount += 1           
                           
            #----FLUID GLOBAL PROPERTIES-----------------------------------------------------------------------------------#
            #### FLUID GLOBAL PROPERTIES
            numberFluid = []
            FluidType = []
            FluidCount = 0
            for fluidProperty,interval in vascularNetwork.globalFluidPolyChaos.iteritems():
                print fluidProperty,interval
                if interval[0] != interval[1]:
                    FluidType.append(fluidProperty)
                    uncertainVarNames.append(str(fluidProperty))
                    intervals.append(interval)         
                    distributions = pc.J(distributions,pc.Uniform(interval[0],interval[1]))         
                    numberFluid.append(FluidCount)
                    number.append(count)
                    count += 1
                    FluidCount += 1
                    
            #----NETWORK SOLVER PROPERTIES-----------------------------------------------------------------------------------#
            #### NETWORK SOLVER PROPERTIES
            numberNetworkSolver = []
            networkSolverCount = 0
            
            intervalsNetworkSolver = [[20.,60.]] #[[1.15,2.0]] #
            uncertainParameterNetworkSolver = ['Age']# ['totalComplianceFactor'] #['totalComplianceFactor']
                        
            ageTest = False
            
        #     for interval in intervalsNetworkSolver:
        #         intervals.append(interval)
        #         distributions = pc.J(distributions,pc.Uniform(interval[0],interval[1])) 
        #         uncertainVarNames.append(uncertainParameterNetworkSolver[networkSolverCount])
        #         numberNetworkSolver.append(networkSolverCount)
        #         number.append(count)
        #         count += 1
        #         networkSolverCount += 1
            
            #----MULTI VARIABLE-----------------------------------------------------------------------------------#
            #### MULTI VARIABLE (one distibution for several variables changed at the same time)
            multiTrue = False
            ## global Variables 
            intervalsMulti = []
            multiVesselIDS = []
            if multiTrue == True:
                intervalsMulti = [[1.0, 0.3]] #[0.707383288834, 1.34743195311] #[[0.75,2.0],[0.75,2.0],[0.75,2.0],[0.75,2.0],[0.75,2.0],[0.75,2.0],[0.75,2.0],[0.75,2.0]]#[] #[[0.75,2.0]] # #[[0.9,1.5],[0.9,1.5] ]#,[1,1.75],[1,1.75],[1,1.75],[1,1.75],[1,1.75],[1,1.75]]#[[0.75,1.75] ,[0.75,1.75],[0.75,1.75],[0.75,1.75]]
                multiVesselIDS = [[1,2,14,18,27,28,35,37,39,41,3,4,5,6,15,19,20,26,29,30,31,32,33,34,36,38,40,7,8,9,10,11,21,22,23,24,25,42,50,51,52,53,54,55,43,44,45,46,47,48,49]]
        
                
                #[[1,2,14,18,27,28,35,37,39,41]] # polychaos aeging just aorta
                
        #                         [[1],
        #                         [2,14,18,27,28,35,37,39,41],
        #                         [3,4,5,6,15,19,20],
        #                         [26,29,30,31,32,33,34,36,38,40],
        #                         [7,8,9,10,11],
        #                         [21,22,23,24,25],
        #                         [42,50,51,52,53,54,55],
        #                         [43,44,45,46,47,48,49] ] # polychaos aeging compartments
                                #[]
                                #[[1,2,3,4,5,6,14,15,18,19,20,27,28,35,37,39,41]]
                                #[[1],
                                #  [2,14,18,27,28,35,37,39,41],
                                #  [3,4,5,6,15,19,20],
                                #  [26,29,30,31,32,33,34,36,38,40],
                                #  [7,8,9,10,11],
                                #  [21,22,23,24,25],
                                #  [42,50,51,52,53,54,55],
                                #  [43,44,45,46,47,48,49] ]
                                ###
                                #[2,14,18,27,28,35,37,39,41],
                                #[26,29,30,31,32,33,34,36,38,40],
                                #[3,4,5,6,15,19,20],
                                #[21,22,23,24,25],
                                #[7,8,9,10,11],
                                #[43,44,45,46,47,48,49],
                                #[42,50,51,52,53,54,55]] #
                                ##
                                #[] #[[0,1,2,3,13,17,18,25,26,27,28,29,30,31,32,33,34]
                                    #,[4,11,12,14,15,16]
                                    #,[5,6,7,8,9,10,19,20,21,22,23,24]
                                    #,[35,36,37,38,39,40,41,42,43,44,45,46,47,48]]
                      
            uncertainParameter = ['beta']#['beta','beta','beta','beta','beta','beta','beta','beta'] #['beta']# #,'beta','beta','beta','beta','beta','beta'] #['beta','beta','beta','beta','beta','beta','beta','beta']
            numberMulti = []
            countMulti = 0
            groupNames = ['aortic arteries'] #['ascending aorta','aortic arteries',r'cerebral \& shoulder','organ arteries','right arm','left arm','right leg', 'left leg'] # ['Age']# #'Aortic Arch','Aortic Ateries','Head+shoulder','Left Arm','Right Arm', 'Left Leg','Right Leg'] #'Aortic Arch','Head Arteries','Arm Arteries', 'Leg Arteries']
            ## random variables for multiple parametric uncertainity in vessel Data
            for interval in intervalsMulti:
                intervals.append(interval)
                distributions = pc.J(distributions,pc.Normal(interval[0],interval[1])) 
                uncertainVarNames.append(groupNames[countMulti])
                numberMulti.append(countMulti)
                number.append(count)
                count += 1
                countMulti += 1
                
            
            distributionsDict = {   'uncertainVarNames':uncertainVarNames,
                                    'number':number,
                                    'numberVessel':numberVessel,
                                    'vesselIds':vesselIds,
                                    'intervals':intervals,
                                    'varNames':varNames,
                                    'BCtype':BCtype,
                                    'BCvesselID':BCvesselID,
                                    'numberBC':numberBC,
                                    'numberFluid':numberFluid,
                                    'FluidType':FluidType,
                                    'FluidCount':FluidCount,
                                    'multiVesselIDS':multiVesselIDS,
                                    'numberMulti':numberMulti,
                                    'intervalsMulti':intervalsMulti,
                                    'uncertainParameter':uncertainParameter,
                                    'uncertainParameterNetworkSolver':uncertainParameterNetworkSolver,
                                    'numberNetworkSolver':numberNetworkSolver}
            
            if not os.path.exists(distributionsDictPath):
                os.makedirs(distributionsDictPath)     
        
            saveFile = open(distributionsDictFile,"wb")       
            cPickle.dump(distributionsDict,saveFile,protocol=2)
            saveFile.close()
            print ".. done"
            
        else:
            try:
                print " load and create distributions"
                loadFile = open(distributionsDictFile,"rb")
                # load pickle
                distributionsDict = cPickle.load(loadFile)
                loadFile.close()
                print ".. done"
            except:
                print 'File does not exits:'
                print distributionsDictFile
                exit()
                
            uncertainVarNames   = distributionsDict['uncertainVarNames']
            number              = distributionsDict['number']
            numberVessel        = distributionsDict['numberVessel']
            varNames            = distributionsDict['varNames']
            vesselIds           = distributionsDict['vesselIds']
            intervals           = distributionsDict['intervals']
            BCtype              = distributionsDict['BCtype']
            BCvesselID          = distributionsDict['BCvesselID']
            numberBC            = distributionsDict['numberBC']
            numberFluid         = distributionsDict['numberFluid']
            FluidType           = distributionsDict['FluidType']
            FluidCount          = distributionsDict['FluidCount']
            multiVesselIDS      = distributionsDict['multiVesselIDS']
            numberMulti         = distributionsDict['numberMulti']
            intervalsMulti      = distributionsDict['intervalsMulti']
            uncertainParameter  = distributionsDict['uncertainParameter']
            uncertainParameterNetworkSolver       = distributionsDict['uncertainParameterNetworkSolver']
            numberNetworkSolver = distributionsDict['numberNetworkSolver']
            
            distributions = pc.J()     
            for interval in intervals:
                distributions = pc.J(distributions,pc.Uniform(interval[0],interval[1])) 
                        
        if varNames != []: nameLength = int(len(max(varNames))+2)
        elif uncertainParameter != [] : nameLength = int(len(max(uncertainParameter))+2)
        else: nameLength = 15
        
        #---CREATE THE UPDATE DICTIONARY FOR VESSEL DATA-------------------------------------------#
        updateDict = {'networkData':{'globalFluid': {} },'vesselData': {}}
        # single vessel uncertanity to updateDict
        for num in numberVessel:
            if vesselIds[num] not in (updateDict['vesselData']): (updateDict['vesselData'])[vesselIds[num]] = {}
            ((updateDict['vesselData'])[vesselIds[num]])[varNames[num]] = intervals[num][0]
        # fluid properties
        for num in numberFluid:
            updateDict['networkData']['globalFluid'][FluidType[num]] = intervals[num+len(numberBC)+len(numberVessel)][0]
        # multi uncertainity to updateDict
        for num in numberMulti:
            for vesselId in multiVesselIDS[num]:
                if vesselId not in (updateDict['vesselData']): (updateDict['vesselData'])[vesselId] = {}
                ((updateDict['vesselData'])[vesselId])[uncertainParameter[num]] = intervalsMulti[num][0]
        
        ### age correlation stuff
        #for vesselId in vascularNetwork.vessels.keys():
        #    if vesselId not in (updateDict['vesselData']): (updateDict['vesselData'])[vesselId] = {}
        #    ((updateDict['vesselData'])[vesselId])['beta'] = 20.0
        
        initialGridUpdateDict = {'networkData':{'globalFluid': {} },'vesselData': {}}
        for vesselId,vessel_i in vascularNetwork.vessels.iteritems():
            if vesselId not in (initialGridUpdateDict['vesselData']): (initialGridUpdateDict['vesselData'])[vesselId] = {}
            ((initialGridUpdateDict['vesselData'])[vesselId])['N'] = vessel_i.N
        
        
        ### print the intervals and distributions
        printArguments = [number,nameLength,numberVessel,vesselIds,varNames,numberBC,BCvesselID,BCtype,numberFluid,FluidType,numberNetworkSolver,uncertainParameterNetworkSolver,numberMulti,multiVesselIDS,uncertainParameter,uncertainVarNames]
        printIntervalData(intervals, 'Intervals', printArguments)
        printIntervalData(distributions, 'Distribution', printArguments)
        
        if totalTime == 'vascularNetwork.totalTime': totalTime = vascularNetwork.totalTime
        if plotxTimeEnd == 'totalTime': plotxTimeEnd = totalTime
    
        ##########################################################################################    
        #---0.1 orthogonal polynoms ( TRUE == create and save, FALSE == load existing)--#
        if createOrthoPoly == True:
            if not os.path.exists(polynomPathOrtho):
                os.makedirs(polynomPathOrtho)
            
            print " create and save orthogonal polynoms "
            #create orthogonal polynom  
            #orthoPoly = pc.orth_gs(order,distributions)
            orthoPoly = pc.orth_ttr(order,distributions)
            orthoPoly = pc.orth_ttr(order,distributions)
            #orthoPoly = pc.orth_chol(order,distributions)
            #orthoPoly = pc.orth_svd(order,distributions)
            
            #save file
            saveFile = open(orthoFile,"wb")       
            cPickle.dump(orthoPoly,saveFile,protocol=2)
            saveFile.close()
            print ".. done"
        else:
            try:
                print " load orthogonal polynoms "
                loadFile = open(orthoFile,"rb")
                # load pickle
                orthoPoly = cPickle.load(loadFile)
                loadFile.close()
                print ".. done"
            except:
                print 'File does not exits:'
                print orthoFile
                exit()
                
        #--- 0.2 collocation method ( TRUE == create and save, FALSE == load existing)
        if createSample == True:
            if not os.path.exists(polynomPathOrtho):
                os.makedirs(polynomPathOrtho)
            print " create and save collocation methods sample space "
    
            
            sampleOrder = 2*pc.terms(order,len(distributions))
            sample = distributions.sample(sampleOrder,sampleMethod)
            
            if len(distributions) == 1: sample = np.array([sample])
            sample = sample.transpose()
            
            print "number of simulations",len(sample)
            #save file
            saveFile = open(sampleFile,"wb")       
            cPickle.dump(sample,saveFile,protocol=2)
            saveFile.close()
            print ".. done"
        else:
            try:
                print " load collocation methods sample space "
                loadFile = open(sampleFile,"rb")
                # load pickle
                sample = cPickle.load(loadFile)
                loadFile.close()
                print ".. done"
            except:
                print 'File does not exits:'
                print sampleFile
                exit()
        
        ################################################################################
        
        #---1.step genrealized polynomial chaos evaluations + data storing ------------------------
        
        if runSimulations == True:
        #standart solution path
            
            if not os.path.exists(saveDirectoryPC):
                os.makedirs(saveDirectoryPC)  
            #--set up log file------------------------------------------------------------------  
            logfile = open(''.join([SolutionFileNamePath1, 'log.txt']), "wb")
            logfile.write(''.join(['Vascular polyChaos Simulation ',str(networkName),'\n','\n']))
            logfile.write(''.join(['Order :',str(order),'  grid method: ',method,'\n','\n']))
            logfile.write(printIntervalData(intervals, 'Intervals', printArguments ,printTrue=False))
            logfile.write('  '.join(['evalNr', '        '.join( ' '.join([varNames[i],str(vesselIds[i])]) for i in numberVessel),'        '.join( ' '.join([BCtype[i],str(BCvesselID[i])]) for i in numberBC),'        '.join( ' '.join(['Multiple ',str(numberMulti[i])]) for i in numberMulti),'\n']))
            for polyCount in np.linspace(0, len(sample)-1, len(sample)):
                logfile.write('  '.join([str(int(polyCount)).center(5), '   '.join( str('%.4f') %i for i in sample[polyCount,:]   ),'\n']))
            logfile.close()
            
            numberOfSmulations = len(sample)
            
            print " starting the calculation of the samples! %.0f Simulations will be run" %(numberOfSmulations)
            startTime = time.clock()
            #---start the simulations
            for polyCount in np.linspace(0, len(sample)-1, len(sample)):
                
                if polyCount >= startEvaluation:
                
                    print '\n polyChaos run', polyCount+1  
                    #--- GET POlYVARABLES OF CURRENT RUN --------------------------------------------------
                    ####               
                    polyVariable = sample[polyCount,:]
                    if  type(polyVariable) == np.float64:
                        polyVariable=np.array([polyVariable])
                    #print polyVariable
                    printIntervalData(polyVariable, 'polyValue', printArguments, number)
                    #--- Set grid Nodes to the initial values to enhance better sim results with cfl meshing etc.
                    vascularNetwork.updateNetwork(initialGridUpdateDict)
                    
                    #### UPDATE VASCULAR NETWORK WITH POLYVARIABLES                 
                    #  update Vessel Data
                    for num in numberVessel:
                        ((updateDict['vesselData'])[vesselIds[num]])[varNames[num]] = polyVariable[num]
                    
                    # update Boundary Conditions
                    splitNumber = len(numberVessel)
                    #for num in numberBC:
                    for num in numberBC:
                        boundaryConditionName, boundaryParameter, vesselId = uncertainVarNames[num+splitNumber].split(' ')
                        for boundaryCondition in vascularNetwork.boundaryConditions[int(vesselId)]:
                            if boundaryCondition.name == boundaryConditionName:
                                boundaryCondition.update({boundaryParameter: polyVariable[num+splitNumber]})
                    splitNumber = splitNumber+ len(numberBC)
                    
                    # update global Fluid Properties
                    for num in numberFluid:
                        updateDict['networkData']['globalFluid'][FluidType[num]] = polyVariable[num+splitNumber]                  
                    splitNumber = splitNumber+ len(numberFluid)            
                    
                    if ageTest == False: 
                        # update network Solver properties
                        for num in numberNetworkSolver:
                            vascularNetwork.networkSolver[uncertainParameterNetworkSolver[num]] = polyVariable[num+splitNumber]
                        splitNumber = splitNumber+ len(numberNetworkSolver)
                    
                    # update Dictionary for multi uncertainty in vessel data
                    for num in numberMulti:
                        for vesselId in multiVesselIDS[num]:
                            ((updateDict['vesselData'])[vesselId])[uncertainParameter[num]] = vascularNetwork.vessels[vesselId].getIntervalDictionary()[uncertainParameter[num]][1] *polyVariable[num+splitNumber]
                            
                            #print vesselId, polyVariable[num+splitNumber]
                            #print vascularNetwork.vessels[vesselId].getIntervalDictionary()[uncertainParameter[num]][1]
                            #print ((updateDict['vesselData'])[vesselId])[uncertainParameter[num]]
                            #print '\n'
                            
                    if ageTest == True:       
                        ## update age test
                        for num in numberNetworkSolver:
                            age = polyVariable[num+splitNumber]
                        
                        ## totalComplianceFactor transver function
                        tCF = -0.00875*age+1.675
                        vascularNetwork.networkSolver['totalComplianceFactor'] = tCF
                        splitNumber = splitNumber+ len(numberNetworkSolver)
                        
                        
                        ## beta transver function
                        beta = 0.025*age+0.5
                        for vesselId in vascularNetwork.vessels.keys():
                            ((updateDict['vesselData'])[vesselId])['beta'] = vascularNetwork.vessels[vesselId].getIntervalDictionary()['beta'][0]*beta  
    
                    # update Network with dictionary
                    vascularNetwork.updateNetwork(updateDict)
                    # ---RUN THE DETERMINISTIC SIMULATIONS-------------------------------------------------------------------
                    timeSolverSolveStart = time.clock()
                    flowSolver = FlowSolver(vascularNetwork,quiet = True)
                    P,Q,A = flowSolver.solve()             
                    timeSolverSolve = time.clock()-timeSolverSolveStart
                    minutesSolve = int(timeSolverSolve/60.)
                    secsSolve = timeSolverSolve-minutesSolve*60.
                    print '      Solver time {} min {} sec'.format(minutesSolve,secsSolve)
                    
                    # ---POSTPROCESS AND SAVE THE DATA---------------------------------------------------------------
                    #calculate wave speed array
                    c = {}
                    for Id,vessel in vascularNetwork.vessels.iteritems():
                        csol = vessel.waveSpeed(A[Id],vessel.C(P[Id]))
                        c[Id]=csol
                    # save the result
                    currentData= {'Pressure': P, 'Flow': Q, 'Area': A,'WaveSpeed':c,'Name': str('polychaos'+str(int(polyCount)).zfill(2)) }
                    #currentData= {'Pressure': P, 'Flow': Q, 'Name': str('polychaos'+str(polyCount[0]).zfill(2)) }
                    saveFile = ''.join([SolutionFileNamePath1,(str(int(polyCount)).zfill(4)),'.pickel'])
                    outFile = open(saveFile,"wb")       
                    cPickle.dump(currentData,outFile,protocol=2)
                    outFile.close()
                  
            
            endTime = time.clock()
            print "It took %1.2f seconds to run polychaos simulation evaluations\n" % (endTime-startTime)
            
            # --saving numberOfEvaluations-File-----------------------------------------------------------
            log = [int(polyCount+1),order,numberOfSmulations]
            logFILE = open(numberOfEvaluationsPath,"w")
            cPickle.dump(log,logFILE)
            logFILE.close()
        
        ################################################################################   
        #---2.step pre process data (if preProcessed == true) and construct generalized polynomial chaos expansion
        if calculateGPCE == True:
            polynoms = []
            solutionInterpolated = None
            if preProcessData == True:
                
                
                print '   loading polychaos solution data'
                startTime = time.clock()
                #load log file
                logFILE = open(numberOfEvaluationsPath,"rb")
                log = cPickle.load(logFILE) 
                logFILE.close()
                numberOfPoly = log[0]
                
                # initialize solution Data
                solutionDataAdapted = []
                for pTC in polynomsToCalculate:
                    solutionDataAdapted.append({'Pressure':[],'Flow':[],'WaveSpeed':[],'Area':[],'Tsteps':999999999})
                
                
                
                # load soutionData from file
                for polyCount in range(0,numberOfPoly,1):
                    loadSolutionFile = ''.join([SolutionFileNamePath1,str(polyCount).zfill(4),'.pickel'])
                    #Load solution data 
                    loadFile = open(loadSolutionFile,"rb")
                    # load pickle
                    solDataCurrent = cPickle.load(loadFile)
                    loadFile.close()
                    # hash the data:
                    for pTC in polynomsToCalculate: 
                        index = polynomsToCalculate.index(pTC, )
                        p = np.ravel(solDataCurrent['Pressure'][pTC[0]][:,[pTC[1]]])
                        solutionDataAdapted[index]['Pressure'].append(p)
                        solutionDataAdapted[index]['Flow'].append(np.ravel(solDataCurrent['Flow'][pTC[0]][:,[pTC[1]]]))        
                        solutionDataAdapted[index]['WaveSpeed'].append(np.ravel(solDataCurrent['WaveSpeed'][pTC[0]][:,[pTC[1]]]))        
                        solutionDataAdapted[index]['Area'].append(np.ravel(solDataCurrent['Area'][pTC[0]][:,[pTC[1]]]))        
                        if len(p) < solutionDataAdapted[index]['Tsteps']:
                            solutionDataAdapted[index]['Tsteps'] = len(p)
                                             
                     
                endTime = time.clock()
                print "    It took %1.2f seconds to load the solution \n" % (endTime-startTime)
                # solution pre processing
                
                # interpolate solution
                solutionInterpolated = []
                rho = vascularNetwork.globalFluid['rho']
                
                #--- 2.1 pre process the data
                count = 0
                for dictionary in solutionDataAdapted:
                    if dictionary['Tsteps'] == 999999999:
                        print 'ERROR in data hashing no/bad data found polyCount is probably zero!'
                        break
                    name = names[count]
                    count = count+1
                    timeValues = np.linspace(0.0,totalTime,dictionary['Tsteps'])
                    solutionInterpolatedCurrentPressure = np.ones((dictionary['Tsteps']))
                    solutionInterpolatedCurrentP_f = np.ones((dictionary['Tsteps']-1))
                    solutionInterpolatedCurrentP_b = np.ones((dictionary['Tsteps']-1))
                    solutionInterpolatedCurrentFlow = np.ones((dictionary['Tsteps']))
                    solutionInterpolatedCurrentQ_f = np.ones((dictionary['Tsteps']-1))
                    solutionInterpolatedCurrentQ_b = np.ones((dictionary['Tsteps']-1))
                    
                    extremaP   = [[],[]] # [[points in time],[amplitudes]]
                    extremaP_f = [[],[]] # [[points in time],[amplitudes]]
                    extremaP_b = [[],[]] # [[points in time],[amplitudes]]
                    extremaQ   = [[],[]] # [[points in time],[amplitudes]]
                    extremaQ_f = [[],[]] # [[points in time],[amplitudes]]
                    extremaQ_b = [[],[]] # [[points in time],[amplitudes]]
    
                    for arrayNumber in range(0,numberOfPoly,1):        
                        # if array has the right amount of timesteps, append it directly
                        pArray = dictionary['Pressure'][arrayNumber]
                        qArray = dictionary['Flow'][arrayNumber]
                        aArray = dictionary['Area'][arrayNumber]
                        cArray = dictionary['WaveSpeed'][arrayNumber]
                        
                        if len(pArray) != (dictionary['Tsteps']):   
                        # if the array has more time steps as required, interpolate to gain an array with right length (/timesteps)
                            timeInt = np.linspace(0.0,totalTime,len(pArray))
                            # correct pressure
                            intFunctionP = interp1d(timeInt,pArray)
                            pArray = intFunctionP(timeValues)
                            # correct flow
                            intFunctionQ = interp1d(timeInt,qArray)
                            qArray = intFunctionQ(timeValues)
                            # correct area
                            intFunctionA = interp1d(timeInt,aArray)
                            aArray = intFunctionA(timeValues)
                            # correct wavespeed
                            intFunctionC = interp1d(timeInt,cArray)
                            cArray = intFunctionC(timeValues)
                        # set sol data    
                        solutionInterpolatedCurrentPressure = np.vstack((solutionInterpolatedCurrentPressure,pArray))
                        solutionInterpolatedCurrentFlow = np.vstack((solutionInterpolatedCurrentFlow,qArray))
                        
                        # linear wave splitting of flow and pressure
                        if linearWaveSplit == True:
                            pArray_f,pArray_b,qArray_f,qArray_b = linearWaveSplitting(pArray,qArray,aArray,cArray,rho)
                        # non linear wave splitting
                        else:
                            pArray_f,pArray_b,qArray_f,qArray_b = nonLinearWaveSplitting(pArray,qArray,aArray,cArray,rho,velocityProfileCoefficient)
                        
                        solutionInterpolatedCurrentP_f = np.vstack((solutionInterpolatedCurrentP_f,pArray_f))
                        solutionInterpolatedCurrentP_b = np.vstack((solutionInterpolatedCurrentP_b,pArray_b))
                        solutionInterpolatedCurrentQ_f = np.vstack((solutionInterpolatedCurrentQ_f,qArray_f))
                        solutionInterpolatedCurrentQ_b = np.vstack((solutionInterpolatedCurrentQ_b,qArray_b))
                        
                        #calculate pressure and flow extrem points
                        minMaxAmpP,minMaxTimeP = minMaxFunction(pArray,timeValues=timeValues,delta=delta[name]['Pressure'], seperateMinMax = False )
                        minMaxAmpQ,minMaxTimeQ = minMaxFunction(qArray,timeValues=timeValues,delta=delta[name]['Flow'], seperateMinMax = False)
                        extremaP[0].append(minMaxTimeP)
                        extremaP[1].append(minMaxAmpP) 
                        extremaQ[0].append(minMaxTimeQ)
                        extremaQ[1].append(minMaxAmpQ) 
                                                
                        #calculate pressure forward/backward and flowforward/backward extrem points
                        minMaxAmpP_f,minMaxTimeP_f = minMaxFunction(pArray_f,timeValues=timeValues,delta=delta[name]['Pressure_f'], seperateMinMax = False )                    
                        minMaxAmpP_b,minMaxTimeP_b = minMaxFunction(pArray_b,timeValues=timeValues,delta=delta[name]['Pressure_b'], seperateMinMax = False )
                        minMaxAmpQ_f,minMaxTimeQ_f = minMaxFunction(qArray_f,timeValues=timeValues,delta=delta[name]['Flow_f'], seperateMinMax = False)
                        minMaxAmpQ_b,minMaxTimeQ_b = minMaxFunction(qArray_b,timeValues=timeValues,delta=delta[name]['Flow_b'], seperateMinMax = False)
                        
                                               
                        # calculate tanget and shoulder points using diastolic point
                        # point [waveData, timeData]
                        shoulderPoint = calculateWaveShoulderPoint(pArray_b,timeValues,minMaxTimeP_b[0],minMaxTimeP_b[1])
                        tangentPoint  = calculateWaveShoulderPointTangent(pArray_b,timeValues,minMaxTimeP_b[0],minMaxTimeP_b[1])
                        
                        # extremaP_b[0] == time, extremaP_b[1] == waveData
                        minMaxTimeP_b.append(shoulderPoint[1])
                        minMaxAmpP_b.append(shoulderPoint[0]) 
                        
                        minMaxTimeP_b.append(tangentPoint[1])
                        minMaxAmpP_b.append(tangentPoint[0]) 
                        
                        extremaP_f[0].append(minMaxTimeP_f)
                        extremaP_f[1].append(minMaxAmpP_f)
                        extremaP_b[0].append(minMaxTimeP_b)
                        extremaP_b[1].append(minMaxAmpP_b) 
                        
                        extremaQ_f[0].append(minMaxTimeQ_f)
                        extremaQ_f[1].append(minMaxAmpQ_f)
                        extremaQ_b[0].append(minMaxTimeQ_b)
                        extremaQ_b[1].append(minMaxAmpQ_b) 
                        
                        
                    #delete the ones in the solution
                    solutionInterpolatedCurrentPressure = np.delete(solutionInterpolatedCurrentPressure, 0,0) 
                    solutionInterpolatedCurrentFlow = np.delete(solutionInterpolatedCurrentFlow, 0,0)
                    solutionInterpolatedCurrentP_f = np.delete(solutionInterpolatedCurrentP_f, 0,0) 
                    solutionInterpolatedCurrentP_b = np.delete(solutionInterpolatedCurrentP_b, 0,0)
                    solutionInterpolatedCurrentQ_f = np.delete(solutionInterpolatedCurrentQ_f, 0,0) 
                    solutionInterpolatedCurrentQ_b = np.delete(solutionInterpolatedCurrentQ_b, 0,0)
                    
                    
                    ## cheat an delete if more then 2 pressure solutions
                    for solNumber in range(len(solutionInterpolatedCurrentPressure)):
                        if len(extremaP[0][solNumber]) == 4:
                            extremaP[0][solNumber].pop()
                            extremaP[0][solNumber].pop()
                            extremaP[1][solNumber].pop()
                            extremaP[1][solNumber].pop()
                    
                    ## Check if extremas are found in each simulation
                    countExtremaP = []
                    countExtremaP_f = []
                    countExtremaP_b = []
                    countExtremaQ = []
                    countExtremaQ_f = []
                    countExtremaQ_b = []
                    notSameExtremaInAllSolData = False
                    for solNumber in range(len(solutionInterpolatedCurrentPressure)):
                        if len(extremaP[0][solNumber]) == len(extremaP[1][solNumber]):
                            countExtremaP.append(len(extremaP[0][solNumber]))
                        else: print "ERROR of extremaP"
                        if len(extremaP_f[0][solNumber]) == len(extremaP_f[1][solNumber]):
                            countExtremaP_f.append(len(extremaP_f[0][solNumber]))
                        else: print "ERROR of extremaP_f"
                        if len(extremaP_b[0][solNumber]) == len(extremaP_b[1][solNumber]):
                            countExtremaP_b.append(len(extremaP_b[0][solNumber]))
                        else: print "ERROR of extremaP_b"
                        
                        if len(extremaQ[0][solNumber]) == len(extremaQ[1][solNumber]):
                            countExtremaQ.append(len(extremaQ[0][solNumber]))
                        else: print "ERROR of extremaQ"
                        if len(extremaQ_f[0][solNumber]) == len(extremaQ_f[1][solNumber]):
                            countExtremaQ_f.append(len(extremaQ_f[0][solNumber]))
                        else: print "ERROR of extremaQ_f"
                        if len(extremaQ_b[0][solNumber]) == len(extremaQ_b[1][solNumber]):
                            countExtremaQ_b.append(len(extremaQ_b[0][solNumber]))
                        else: print "ERROR of extremaQ_b"
                    
                            
                    if sum(np.asarray(countExtremaP) == countExtremaP[0]) != len(countExtremaP):
                        print "ERROR: vPC: not same minMaxPoints extremaP are found for all solutions"
                        print countExtremaP
                        notSameExtremaInAllSolData = True
                    if sum(np.asarray(countExtremaP_f) == countExtremaP_f[0]) != len(countExtremaP_f):
                        print "ERROR: vPC: not same minMaxPoints countExtremaP_f are found for all solutions"
                        print countExtremaP_f
                        notSameExtremaInAllSolData = True
                    if sum(np.asarray(countExtremaP_b) == countExtremaP_b[0]) != len(countExtremaP_b):
                        print "ERROR: vPC: not same minMaxPoints countExtremaP_b are found for all solutions"
                        print countExtremaP_b
                        notSameExtremaInAllSolData = True
                    if sum(np.asarray(countExtremaQ) == countExtremaQ[0]) != len(countExtremaQ):
                        print "ERROR: vPC: not same minMaxPoints countExtremaQ are found for all solutions"
                        print countExtremaQ
                        notSameExtremaInAllSolData = True
                    if sum(np.asarray(countExtremaQ_f) == countExtremaQ_f[0]) != len(countExtremaQ_f):
                        print "ERROR: vPC: not same minMaxPoints countExtremaQ_f are found for all solutions"
                        print countExtremaQ_f
                        notSameExtremaInAllSolData = True
                    
                    if sum(np.asarray(countExtremaQ_b) == countExtremaQ_b[0]) != len(countExtremaQ_b):
                        print "ERROR: vPC: not same minMaxPoints countExtremaQ_b are found for all solutions"
                        print countExtremaQ_b
                        notSameExtremaInAllSolData = True
                    
                    solutionInterpolated.append({'Pressure':solutionInterpolatedCurrentPressure,'Pressure_f':solutionInterpolatedCurrentP_f,'Pressure_b':solutionInterpolatedCurrentP_b,
                                                                      'Flow':solutionInterpolatedCurrentFlow,'Flow_f':solutionInterpolatedCurrentQ_f,'Flow_b':solutionInterpolatedCurrentQ_b,
                                                                      'extremaPressure': extremaP,'extremaFlow': extremaQ,
                                                                      'extremaPressure_f': extremaP_f,'extremaFlow_f': extremaQ_f ,
                                                                      'extremaPressure_b': extremaP_b,'extremaFlow_b': extremaQ_b  })
                
                if not os.path.exists(solutionInterpolatedPath):
                            os.makedirs(solutionInterpolatedPath)              
                
                FILE = open(solutionInterpolatedFile,"w")
                # store pickle and close file
                cPickle.dump(solutionInterpolated,FILE,protocol=2)
                FILE.close()
            else:
                notSameExtremaInAllSolData = False
                print '   loading interpolated and adapted polychaos solution data'
                startTime = time.clock()
                #Load solution data 
                loadFile = open(solutionInterpolatedFile,"rb")
                # load pickle
                solutionInterpolated = cPickle.load(loadFile)
                loadFile.close()
                endTime = time.clock()
                print "    It took %1.2f seconds to load the solution \n" % (endTime-startTime)
    
            ################################################################################
            #--- 2.2 plot min max points
            if plotMinMaxPoints == True:
                print '   plotting: P,Q extrema'
                if not os.path.exists(saveDirectoryPolySolPlots):
                    os.makedirs(saveDirectoryPolySolPlots)  
                
                fig = plt.figure(1, edgecolor='k',dpi=100)
                fig.subplots_adjust(hspace=0.5)   
                fig.set_figwidth(8.27*3)
                fig.set_figheight((11.69/3)*2)
                nIndex = 0 
                extremPairs = {'extremaPressure':'extremaFlow','extremaPressure_f':'extremaFlow_f','extremaPressure_b':'extremaFlow_b'}
                extremNormalP = {'extremaPressure':'Pressure','extremaPressure_f':'Pressure_f','extremaPressure_b':'Pressure_b'}
                extremNormalQ = {'extremaFlow':'Flow','extremaFlow_f':'Flow_f','extremaFlow_b':'Flow_b'}
                
                for sol in solutionInterpolated:
                    count = 0
                    plt.rc('text', usetex=1)
                    for pressurePeaks,flowPeaks in extremPairs.iteritems():    
                        count = count +1   
                        for timeP in  sol[pressurePeaks][0]:
                            plt.subplot(2,3,count)
                            index = sol[pressurePeaks][0].index(timeP, )
                            if index !=0:
                                timeInt = np.linspace(0.0,totalTime,len(sol[extremNormalP[pressurePeaks]][index]))
                                plt.plot(timeP, sol[pressurePeaks][1][index],color=(0.6,0.6,1) ,marker = 'x',linestyle = '', linewidth = 1)
                                plt.plot(timeInt, sol[extremNormalP[pressurePeaks]][index],color=(0.6,1 , 0.5),linestyle = ':', linewidth = 1)
                            plt.xlabel('Time [s]')
                            plt.ylabel('Pressure [Pa]')
                            pressurePeaksTitle = pressurePeaks
                            if '_' in pressurePeaksTitle: pressurePeaksTitle = pressurePeaksTitle.replace('_','-')
                            plt.title (pressurePeaksTitle)
                            #plt.grid(True,linestyle='-',color='0.1',linewidth = 0.4)    
                            ax = fig.gca() 
                            ax.yaxis.major.formatter.set_powerlimits((0,0))
                            
                            ax.spines['top'].set_visible(False)
                            ax.tick_params(axis='x',top='off')
                            ax.spines['right'].set_visible(False)
                            ax.tick_params(axis='y',right='off')
                            
                            
                            #x1,x2,y1,y2 = plt.axis()
                            #plt.axis((0,totalTime,y1+0.1*y1,y2+0.1*y2))
                        timeInt = np.linspace(0.0,totalTime,len(sol[extremNormalP[pressurePeaks]][0]))
                        plt.plot(sol[pressurePeaks][0][0], sol[pressurePeaks][1][0],color=(0.0,0.0,1) ,marker = 'o',linestyle= '' ,linewidth = 1)
                        plt.plot(timeInt, sol[extremNormalP[pressurePeaks]][0],color='g' ,linestyle = '-', linewidth = 1)
                         
                        for timeQ in  sol[flowPeaks][0]:
                            plt.subplot(2,3,count+3)
                            index = sol[flowPeaks][0].index(timeQ, )
                            if index != 0: 
                                timeInt = np.linspace(0.0,totalTime,len(sol[extremNormalQ[flowPeaks]][index]))
                                plt.plot(timeQ, sol[flowPeaks][1][index],color=(1.0 ,0.6 , 0.6) ,marker = 'x',linestyle = '', linewidth = 1)
                                plt.plot(timeInt, sol[extremNormalQ[flowPeaks]][index],color=(0.6,1 , 0.6) ,linestyle = ':', linewidth = 1)
                            plt.xlabel('Time [s]')
                            #plt.ylabel('Flow [$m^3/s$]')
                            plt.ylabel('Flow')
                            flowPeaksTitle = flowPeaks
                            if '_' in flowPeaksTitle: flowPeaksTitle = flowPeaksTitle.replace('_','-')
                            plt.title (flowPeaksTitle)
                            #plt.grid(True,linestyle='-',color='0.1',linewidth = 0.4)    
                            ax = fig.gca() 
                            ax.yaxis.major.formatter.set_powerlimits((0,0))                       
                            #x1,x2,y1,y2 = plt.axis()
                            #plt.axis((0,totalTime,y1+0.1*y1,y2+0.1*y2))
                            ax.spines['top'].set_visible(False)
                            ax.tick_params(axis='x',top='off')
                            ax.spines['right'].set_visible(False)
                            ax.tick_params(axis='y',right='off')
                           
                        timeInt = np.linspace(0.0,totalTime,len(sol[extremNormalQ[flowPeaks]][0]))
                        plt.plot(sol[flowPeaks][0][0], sol[flowPeaks][1][0],color= (1 ,0 , 0) ,marker = 'o',linestyle = '', linewidth = 1)
                        plt.plot(timeInt, sol[extremNormalQ[flowPeaks]][0],color='g' ,linestyle = '-', linewidth = 1)
                    savePath = ''.join([saveDirectoryPolySolPlots,'Polychaos_point_',names[nIndex],'_N.',str(order),'_method_',method,'_Extrema',plotFileType])
                    plt.savefig(savePath)
                    plt.clf()
                    nIndex = nIndex+1
            
            if notSameExtremaInAllSolData == True: exit()
            
            ################################################################################
              
            if solutionInterpolated != None:
            #--- 2.3 claculate gPCE
                print "    starting the polychaos polynomial calculation from polychaos simulation result!!"
                startTime = time.clock()
                
                polynomsT = []
                count = 0
                for sol in solutionInterpolated:
                    print "      Polynomial calculation for point ", names[count]     
                    count = count+1
                    polyDict = {}
                    
                    for tag,data in sol.iteritems():
                        # polynoms for the total pressure signal
                        print "        polynoms for ",str(tag)
                        if 'extrema' not in tag:
                            
                            polynomial = pc.fitter_lr(orthoPoly, sample.T, data)     
                            polyDict[tag]= polynomial
                            
                        else:
                            print sample.shape
                            
                            #print orthoPoly.dim
                            #print orthoPoly.shape
                            
                            #print sample.ravel().shape
                            #print len(data[0])
                            
                            
                            polynomialTime = pc.fitter_lr(orthoPoly, sample.T, data[0])    
                            polynomialAmp  = pc.fitter_lr(orthoPoly, sample.T, data[1])    
                            
                            extremaDict = {'Time':polynomialTime,'Amp':polynomialAmp}
                            polyDict[tag]= extremaDict
                            
                    endTime = time.clock()
                    polynomsT.append(polyDict)
                    
                print "    It took %1.2f seconds to create the polynoms with order %d\n" % (endTime-startTime,order)
                
                if not os.path.exists(polynomPath):
                    os.makedirs(polynomPath)    
                
                #save the polynom
                #create file with the network name in solution path directory
                FILE = open(polynomFile,"w")
                # store pickle and close file
                cPickle.dump(polynomsT,FILE,protocol=2)
                FILE.close()
                polynoms = [polynomsT]
      
        ################################################################################
        ################################################################################
        #--- 3.step post processing - sensitivity analysis
        
        ########################################################################################## 
        ### load polynoms 
        if ((calculateGPCE == False) or len(polynomsToPlotOrder) > 1.0) and postProcessing == True:
            
            print '   loading polynoms'
            startTime = time.clock()
            polynoms = []
            for polyToPlot in polynomsToPlotOrder:
                
                if polyToPlot != 0:
                    polynomPathTemp = ''.join([polynomPath,'polyChaos_solution_',dataNumber,'_order_',str(polyToPlot).zfill(2),'_method_',method,'_polynom','.pickle'])
                    #Load solution data 
                    loadFile = open(polynomPathTemp,"rb")
                else:
                    loadFile = open(polynomFile,"rb")
                # load pickle
                polynom = cPickle.load(loadFile)
                loadFile.close()
                polynoms.append(polynom)
            endTime = time.clock()
            print "    It took %1.2f seconds to load the polynoms\n" % (endTime-startTime)
    
        ##########################################################################################      
        
        if postProcessing == True:
            #load sim result without polychoas
            solutionDataPlain = [] # for the solution data
            if deterministicDataSetNumbers != []:
                endName1 = "_SolutionData_"
                endName2 = ".pickle"
                networkPath = str(cur+"/NetworkFiles/")
                startTime = time.clock()
                for setNumber in deterministicDataSetNumbers:
                    #standart solution path
                    filePath= ''.join([networkPath,networkName,'/',networkName,endName1,str(setNumber).zfill(3),endName2])
                    #create file with filename in soultionPath
                    FILE = open(filePath,"rb")
                    solutionDataLoad = cPickle.load(FILE)
                    solutionDataPlain.append(solutionDataLoad[1][0])
                    FILE.close()   
                endTime = time.clock()
                print "It took %1.2f seconds to load the orginal result\n" % (endTime-startTime) 
        
            # prepare plots of the peaks
            if plotPeaks == True: 
                peakDirs  = [plotDirectoryPeakP_bar,plotDirectoryPeakP_Conf,plotDirectoryPeakP_box,plotDirectoryPeakPLogs]
                for dirs in peakDirs:
                    if not os.path.exists(dirs):
                        os.makedirs(dirs) 
     
                #init log file
                logfile = open(''.join([plotDirectoryPeakPLogs,'Polychaos_peak_sensitivity','_','Pressure','_N.',str(order),'_method_',method,'_.txt']),'wb')                    
                logfile.write('')
                logfile.close()    
                logfile = open(''.join([plotDirectoryPeakPLogs,'Polychaos_peak_sensitivity','_','Flow','_N.',str(order),'_method_',method,'_.txt']),'wb')                    
                logfile.write('')
                logfile.close()    
                
            for pTC in polynomsToCalculate:  
                
                vesselId = pTC[0]
                girdPoint = pTC[1]
                
                count = polynomsToCalculate.index(pTC, )
                name = names[count]
                
                print 'Sensitivity analysis for point',name
                # load external data sets (no plotting methods implemented  (IKEA-prinzip)(do it your self) )
                deterministicSolution = None
                if solutionDataPlain != []:
                    timeValuesDet = []
                    timeValuesDetSplit = []
                    deterministicSolution = {'Pressure':None,'Flow':None,
                                             'Pressure_f':None,'Pressure_b':None,
                                             'Flow_f':None,'Flow_b':None}
                    detPres = []
                    detFlow = []
                    detPres_b = []
                    detFlow_b = []
                    detPres_f = []
                    detFlow_f = []
                    rho = vascularNetwork.globalFluid['rho']
                    
                    for solutionData in solutionDataPlain:
                        pCurr = solutionData['Pressure'][vesselId][:,[girdPoint]]
                        detPres.append(pCurr)
                        qCurr = solutionData['Flow'][vesselId][:,[girdPoint]]
                        detFlow.append(qCurr)
                        timeValuesDet.append(np.linspace(0,totalTime,len(pCurr)))
                        aCurr = solutionData['Area'][vesselId][:,[girdPoint]]
                        cCurr = solutionData['WaveSpeed'][vesselId][:,[girdPoint]]
                        pArray_f,pArray_b,qArray_f,qArray_b = linearWaveSplitting(pCurr,qCurr,aCurr,cCurr,rho)
                        detPres_b.append(pArray_b)
                        detFlow_b.append(qArray_b)
                        detPres_f.append(pArray_f)
                        detFlow_f.append(qArray_f)
                        timeValuesDetSplit.append(np.linspace(0,totalTime,len(pArray_f)))
                        
                    deterministicSolution['Pressure']   = detPres
                    deterministicSolution['Flow']       = detFlow
                    deterministicSolution['Pressure_b'] = detPres_b
                    deterministicSolution['Flow_b']     = detFlow_b
                    deterministicSolution['Pressure_f'] = detPres_f
                    deterministicSolution['Flow_f']     = detFlow_f
                    
                        
                ##########################################################################################          
                #### 3.1 mean and STD plots for the total P and Q signal
                if plotMeanSTD == True:
                    
                    print '   Plotting: mean and STD'
                    if not os.path.exists(plotDirectoryMeanVar):
                        os.makedirs(plotDirectoryMeanVar) 
                    
                    unitScale = {'Pressure': 1.0/133.32, 'Pressure_f':1.0/133.32, 'Pressure_b':1.0/133.32,
                                 'Flow':1e6,'Flow_f':1e6,'Flow_b':1e6}
                    
                    # load results from Reymond
                    DataReymondCode55      = cPickle.load(open(cur+'/Visualisation/PlotsForCodeValidations/outputDataReymondCode.v1dfExD','rb'))
                    p_Reymond55            = DataReymondCode55['SolutionData'][1]['Pressure']
                    
                    DataReymondCode55Visco = cPickle.load(open(cur+'/Visualisation/PlotsForCodeValidations/outputDataReymondCodeViscoElastic.v1dfExD','rb'))
                    p_Reymond55Visco       = DataReymondCode55Visco['SolutionData'][1]['Pressure']
                    
                    for quantity,unit in plotMeanStdQuantities.iteritems():
                        print '    Quantity: ',quantity
                        # polinomual solution 
                        
                        expected = pc.E(polynoms[0][count][quantity],distributions)
                        STD = np.sqrt(np.abs(pc.Var(polynoms[0][count][quantity],distributions)))
                        timeValuesPolynom = np.linspace(0,totalTime,len(expected))
                        
                        maxNumber = 2000
                        timeValuesInterpolated = np.linspace(0,totalTime,maxNumber)
                        
                        intFunctionExp = interp1d(timeValuesPolynom,expected)
                        expectedInterpolated = intFunctionExp(timeValuesInterpolated)
                        
                        intFunctionSD  = interp1d(timeValuesPolynom,STD)
                        SDinterpolated = intFunctionSD(timeValuesInterpolated)
                        
                        #for expi,sdi in zip(expectedInterpolated.tolist(),SDinterpolated.tolist()):
                        #    print expi*unitScale[quantity],sdi*unitScale[quantity]
                        
                                            
                        
                        # insert mean values of the distributions into the polynom (automatically)
                        distMean = pc.E(distributions)
                        #polynomWithMeans = polynoms[0][count][quantity](*distMean)
                        
                        if plotMeanConfidenceInterval == True:
                            print '     Calculate {} Confidence Intervals'.format(100-plotMeanConfidenceAlpha)
                            lenPolynoms = len(polynoms[0][count][quantity])
                            confIntervalLow = np.zeros(lenPolynoms)
                            confIntervalUp  = np.zeros(lenPolynoms)
                            index = 0
                            for polynom in polynoms[0][count][quantity]:
                                if polynom != 0:
                                    confInterval = pc.Perc(polynom, [plotMeanConfidenceAlpha/2., 100-plotMeanConfidenceAlpha/2.],distributions)
                                    if len(distributions) == 1: confInterval = confInterval[0]
                                    confIntervalLow[index] = confInterval[0]
                                    confIntervalUp[index]  = confInterval[1]
                                index += 1
                            print '     Done'            
                        partialPolyExpected = []
                        partialPolySTD =  []
                        #for indexD in number:
                            #calculate mean and STD with respect to a certain random varibale
                            #currDistMean = pc.E(distributions)
                            #currDistMean[indexD] = np.nan
                            # reduce polynomials
                            #currPolynomTime = polynoms[0][count][quantity](*currDistMean)
                            #partialPolyExpected.append(pc.E(currPolynomTime,distributions))
                            #partialPolySTD.append(np.sqrt(np.abs(pc.Var(currPolynomTime,distributions))))
                        from matplotlib import rc
                        from matplotlib import rcParams
                        from matplotlib.font_manager import FontProperties

                        fontLegend = FontProperties()
                        fontLegend.set_size('medium')
                        
                        rcParams['text.usetex']=True
                        rcParams['text.latex.unicode']=True
                        rcParams['font.family']='sans-serif'
                        #plt.rc('font', family='sans-serif')
                        
                        fig = plt.figure(1, edgecolor='k', dpi = 800)
                        fig.subplots_adjust(hspace=0.5)   
                        fig.subplots_adjust(bottom=0.2)   
                        fig.set_figwidth(8.27)
                        fig.set_figheight((11.69/3))
                        plt.clf()
                        
                        ax = plt.subplot(111)
                                                
                        #ax.plot(timeValuesPolynom, polynomWithMeans*unitScale[quantity],color='g' ,linestyle = ':', linewidth = 1.0, label = 'dist.means')
                        
                        if deterministicSolution:
                            try:
                                ax.plot(timeValuesDet[0], deterministicSolution[quantity][0]*unitScale[quantity],color='r' ,linestyle = '-', linewidth = 1.0, label = 'deterministic young case')
                            except:pass 
                            try:
                                ax.plot(timeValuesDet[1], deterministicSolution[quantity][1]*unitScale[quantity],color='g' ,linestyle = '-', linewidth = 1.0, label = 'deterministic old case')
                            except:pass
                            try:
                                ax.plot(timeValuesDetSplit[0], deterministicSolution[quantity][0]*unitScale[quantity],color='r' ,linestyle = '-', linewidth = 1.0, label = 'deterministic young case')
                            except:pass 
                            try:
                                ax.plot(timeValuesDetSplit[1], deterministicSolution[quantity][1]*unitScale[quantity],color='g' ,linestyle = '-', linewidth = 1.0, label = 'deterministic old case')
                            except:pass                      
                            
                        #plt.plot(timeValuesDet[0], detPres[0]*unitScale[quantity],color='r' ,linestyle = '--', linewidth = 1.5, label = 'det. young person')
                        #plt.plot(timeValuesDet[1], detPres[1]*unitScale[quantity],color='g' ,linestyle = ':', linewidth = 1.5, label = 'det. old person')
                        
                        #plt.plot(timeValuesPolynom, partialPolyExpected[0]*unitScale[quantity],color='r' ,linestyle = '-', linewidth = 1.5, label = 'aortic arch')
                        #plt.plot(timeValuesPolynom, partialPolyExpected[2]*unitScale[quantity],color='g' ,linestyle = '-', linewidth = 1.5, label = 'arm arteries')
                        #plt.plot(timeValuesPolynom, partialPolyExpected[1]*unitScale[quantity],color='m' ,linestyle = ':', linewidth = 2.5, label = 'head arteries')
                        #plt.plot(timeValuesPolynom, partialPolyExpected[3]*unitScale[quantity],color='b' ,linestyle = '--', linewidth = 1.5, label = 'leg arteries')
                        
                        
                        #plt.plot(timeValuesPolynom, confIntervalLow,color='g' ,linestyle = ':', linewidth = 0.5, label = 'lower bound')
                        #plt.plot(timeValuesPolynom, confIntervalUp,color='g' ,linestyle = ':', linewidth = 0.5, label = 'upper bound')
                        if plotMeanSigmaInterval:
                            factor = 1
                            ax.fill_between(timeValuesPolynom, expected*unitScale[quantity]-STD*factor*unitScale[quantity], expected*unitScale[quantity]+STD*factor*unitScale[quantity], facecolor='yellow', alpha=0.5, label='sigma')
                            
                            sigma = plt.Rectangle((0, 0), 0, 0, color='y', label=r'$E \pm SD$' )
                            
                            ax.add_patch(sigma)
                        if plotMeanConfidenceInterval:    
                            plt.fill_between(timeValuesPolynom, confIntervalLow*unitScale[quantity], confIntervalUp*unitScale[quantity], facecolor='#7DDAFF', alpha=0.5,label='plotPeaksConfidenceAlpha confidence')
                            
                            percent = str(100-plotMeanConfidenceAlpha)
                            #labelConf = ''.join([percent,'$\displaystyle\\% \\plotPeaksConfidenceAlpha$-interval'])
                            labelConf = ''.join([percent,'$\\%$ Confidence'])
                            
                            conf = plt.Rectangle((0, 0), 0, 0, color='#7DDAFF', label= labelConf)
                            ax.add_patch(conf)
                        
                        t_Reymond55New = np.linspace(0,totalTime,len(p_Reymond55))  
                        ax.plot(t_Reymond55New, p_Reymond55,color='r' ,linestyle = '--', linewidth = 2.0, label = '(R) - elastic') #all groups
                        
                        t_Reymond55New = np.linspace(0,totalTime,len(p_Reymond55Visco))  
                        ax.plot(t_Reymond55New, p_Reymond55Visco,color='k' ,linestyle = ':', linewidth = 2.0, label = '(R) - visco-elastic') #all groups
                        
                        ax.plot(timeValuesPolynom, expected*unitScale[quantity],color='b' ,linestyle = '-', linewidth = 1.0, label = '(S) expected') #all groups
                        
                                                
                        ax.set_xlabel(r'Time [s]')
                        plotNameQuantity = quantity
                        if "_" in quantity: plotNameQuantity = plotNameQuantity.replace('_','-')
                        ax.set_ylabel(' '.join([plotNameQuantity,unit]))
                        if limits[quantity] != [0,0]: plt.set_ylim(limits[quantity])
                        
                        ax.set_xlim(plotxTimeStart,plotxTimeEnd)
                        
                        # set x, y lables
                        plt.xticks(np.linspace(plotxTimeStart,plotxTimeEnd,len(plotXLablesTime)), plotXLablesTime)
                        plt.yticks(np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],plotYLabelNumber),np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],plotYLabelNumber)) 
                        
                        
                        ax.set_title(r'Biological variation of $\beta$ versus different numerical models')
                        #plt.grid(True,linestyle='-',color='0.1',linewidth = 0.4)    
                        ax.legend(loc='best', frameon=False, ncol=1, prop = fontLegend)
                        
                        ax.spines['top'].set_visible(False)
                        ax.tick_params(axis='x',top='off')
                        ax.spines['right'].set_visible(False)
                        ax.tick_params(axis='y',right='off')                       
                        
                        ax = fig.gca() 
                        #ax.yaxis.major.formatter.set_powerlimits((3,3))
                        additionalFilenameIdentifier = ''
                        if plotMeanSigmaInterval:         additionalFilenameIdentifier = 'SigmaInterval'
                        elif plotMeanConfidenceInterval:  additionalFilenameIdentifier = 'ConfidenceInterval'
                        
                        savePath = ''.join([plotDirectoryMeanVar,'Polychaos_point_',name,'_',quantity,'_N',str(order),'_method_',method,'_mean_',str(additionalFilenameIdentifier),plotFileType])
                        plt.savefig(savePath)
                        plt.clf()
                        ####################################################################################################################################
                        ### STD plots 
                        
                        ax = plt.subplot(111)
                        
                        ax.plot(timeValuesPolynom, STD*unitScale[quantity],color='k' ,linestyle = '-', linewidth = 1, label = 'STD')
                        
                        #ax.plot(timeValuesPolynom, STD, color='k' ,linestyle = '-', linewidth = 1, label = 'STD')
                        
                        #plt.plot(timeValuesPolynom, partialPolySTD[0]*unitScale[quantity],color='r' ,linestyle = '-', linewidth = 1.5, label = 'aortic arch')
                        #plt.plot(timeValuesPolynom, partialPolySTD[2]*unitScale[quantity],color='g' ,linestyle = '-', linewidth = 1.5, label = 'arm arteries')
                        
                        
                        #plt.plot(timeValuesPolynom, partialPolySTD[3]*unitScale[quantity],color='b' ,linestyle = '--', linewidth = 1.5, label = 'leg arteries')
                        #plt.plot(timeValuesPolynom, partialPolySTD[1]*unitScale[quantity],color='m' ,linestyle = ':', linewidth = 2.5, label = 'head arteries')
                        ax.set_xlabel(r'Time [s]')
                        plotNameQuantity = quantity
                        if "_" in quantity: plotNameQuantity = plotNameQuantity.replace('_','-')
                        ax.set_ylabel(' '.join([plotNameQuantity,unit]))
                        #plt.ylabel('Pressure '+unit)
                        if limits[quantity] != [0,0]: plt.set_ylim(limits[quantity])
                        ax.set_xlim(plotxTimeStart,plotxTimeEnd)
                        #ax.set_title ('STD')
                        #plt.grid(True,linestyle='-',color='0.1',linewidth = 0.4)    
                        ax.legend(loc='best', frameon=False, ncol=1, prop = fontLegend)
                        ax.spines['top'].set_visible(False)
                        ax.tick_params(axis='x',top='off')
                        ax.spines['right'].set_visible(False)
                        ax.tick_params(axis='y',right='off') 
                        
                        plt.xticks(np.linspace(plotxTimeStart,plotxTimeEnd,len(plotXLablesTime)), plotXLablesTime)
                        plt.yticks(np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],plotYLabelNumber),np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],plotYLabelNumber)) 
                        
                        # Shink current axis by 20%                                                
                        ax = fig.gca() 
                        #ax.yaxis.major.formatter.set_powerlimits((3,0))
                        savePath = ''.join([plotDirectoryMeanVar,'Polychaos_point_',name,'_',quantity,'_N',str(order),'_method_',method,'_STD',plotFileType])
                        plt.savefig(savePath)
                        plt.clf()
                        
                        fig.clear()
                 
                ##########################################################################################
                #### 3.2 sensitivity of peaks  
                if plotPeaks == True:
                    
                    print '   plotting: pressure peak deviation in time and amplitude'
                                        
                    #quantities to plot: analytic bool
                    quantitiesToPlotPeaks = {'extremaPressure':False,'extremaFlow':False,
                                             'extremaPressure_f':False,'extremaFlow_f':False,
                                             'extremaPressure_b':False,'extremaFlow_b':False}
                    quantitiesNames = { 'extremaPressure':'Pressure','extremaFlow':'Flow',
                                        'extremaPressure_f':'Pressure_f','extremaFlow_f':'Flow_f',
                                        'extremaPressure_b':'Pressure_b','extremaFlow_b':'Flow_b'}
                    quantitiesUnits = {'extremaPressure':'[mmHg]', 'extremaPressure_f':'[mmHg]', 'extremaPressure_b':'[mmHg]',
                                       'extremaFlow':'[ml/s]','extremaFlow_f':'[ml/s]','extremaFlow_b':'[ml/s]'}
                    unitScale = {'Pressure': 1.0/133.32, 'Pressure_f':1.0/133.32, 'Pressure_b':1.0/133.32
                                 ,'Flow': 1/1.e-6,'Flow_f':1/1.e-6,'Flow_b':1/1.e-6}
                    
                    approxError = 1.0e-4
                    
                    interSavedData = {}
                    
                    ### Find all Bifurcations parsing the tree
                    bifurcationList = []
                    viz = []
                    #find root
                    root = vascularNetwork.root    
                    # add root to the viz vessels if root has daughters:
                    if vascularNetwork.vessels[root].leftDaughter != None:
                        viz.append(root)
                    # loop through tree until all daughters are conected
                    while len(viz) != 0:
                        # get the mother vessel (already added) and add its daughters
                        motherVessel = viz.pop(0)
                        # find left daughter
                        leftDaughter  = vascularNetwork.vessels[motherVessel].leftDaughter
                        rightDaughter = vascularNetwork.vessels[motherVessel].rightDaughter
                        #append the mother to the calc list
                        curCalcList = [motherVessel]
                        
                        if leftDaughter != None:
                            #append the left daughter to the calc list
                            curCalcList.append(leftDaughter)
                           
                            if rightDaughter != None:
                                #append the left daughter to the calc list
                                curCalcList.append(rightDaughter)
                            else:
                                curCalcList.append(None)
                            # check if leftDaughter has also daughters 
                            if vascularNetwork.vessels[leftDaughter].leftDaughter != None:
                                viz.append(leftDaughter)
                                
                            if rightDaughter != None:
                                # check if rightDaughter has also daughters 
                                if vascularNetwork.vessels[rightDaughter].leftDaughter != None:
                                    viz.append(rightDaughter)
                        bifurcationList.append(curCalcList)
                    
                    quantitiesToPlotTemp = {}
                    for quantity,peaks in peaksToEvaluate[name].iteritems():
                        if len(peaks) != 0.0:
                            quantitiesToPlotTemp[quantity] = quantitiesToPlotPeaks[quantity]
                    
                    ###############################################################################################################
                    ###############################################################################################################
                    ### Peak Analysis: calculate E(x), STD and Confidence Intervals 
                    
                    if peakAnalysis is False: 
                        try:
                            FILE = open(peakAnalysisFile,"rb")
                            interSavedData = cPickle.load(FILE)
                            FILE.close()
                        except:
                            print 'ERROR try to load PeakAnlaysis File - Peak Anlaysis is not Performed'
                            input = raw_input('  Should Peak Anlaysis be performed now (y)? ')
                            if input == 'y': peakAnalysis = True
                            else: exit()
                    
                    
                    if peakAnalysis:
                        for quantity,analyticBool in quantitiesToPlotTemp.iteritems():
                            print '    run peak analysis for', quantitiesNames[quantity]
                            if len(polynoms[0][count][quantity]['Time']) == 0:
                                print 'no polynoms for ',quantity
                                break
                                
                            # sensitivity with respect to all random variables
                            peakTimeExpected = pc.E(polynoms[0][count][quantity]['Time'],distributions)
                            peakTimeSTD = np.abs(pc.Var(polynoms[0][count][quantity]['Time'],distributions))
                            peakTimeSTD= np.sqrt(peakTimeSTD) #[peakTimeSTD != 0.0] 
                            
                            timeValuesPeak = np.linspace(0,len(peakTimeExpected),len(peakTimeExpected))
                            
                            peakAmpExpected = pc.E(polynoms[0][count][quantity]['Amp'],distributions)
                            peakAmpSTD = np.abs(pc.Var(polynoms[0][count][quantity]['Amp'],distributions))
                            peakAmpSTD= np.sqrt(peakAmpSTD) #[peakAmpSTD != 0.0] 
                                
                            # 100-plotPeaksConfidenceAlpha% confidence intervals:
                            peakconfIntTime = (pc.Perc(polynoms[0][count][quantity]['Time'],  [plotPeaksConfidenceAlpha/2., 100-plotPeaksConfidenceAlpha/2.],distributions))
                            peakconfIntAmp = (pc.Perc(polynoms[0][count][quantity]['Amp'],  [plotPeaksConfidenceAlpha/2., 100-plotPeaksConfidenceAlpha/2.],distributions))
                            
                            # sensitivity with respect to a certain random variables
                            
                            #variables for E, STD and correlation 
                            meanTime = []
                            meanAmp = []
                            STDTime = []
                            STDAmp = []
                            correlationTime = []
                            correlationAmp = []
                            sobolIndexTime = []
                            sobolIndexAmp = []
                            confIntTime = []
                            confIntAmp = []
                            x = pc.variable(max(number)+1)
                            
                            timepeaks = np.linspace(0,len(polynoms[0][count][quantity]['Time'])-1,len(polynoms[0][count][quantity]['Time']))
                            amppeaks = np.linspace(0,len(polynoms[0][count][quantity]['Amp'])-1,len(polynoms[0][count][quantity]['Amp']))
                            
                            #calculate sensitivity sobol indices                                
                            sobolIndexTime = pc.Sens_m(polynoms[0][count][quantity]['Time'],distributions)
                            sobolIndexAmp  = pc.Sens_m(polynoms[0][count][quantity]['Amp'],distributions)
                        
                            for indexD in number:
                                #calculate mean and STD with respect to a certain random varibale
                                currDistMean = pc.E(distributions)
                                currDistMean[indexD] = np.nan
                                # reduce polynomials
                                currPolynomTime = polynoms[0][count][quantity]['Time'](*currDistMean)
                                currPolynomAmp = polynoms[0][count][quantity]['Amp'](*currDistMean)
                                
                                meanTime.append(pc.E(currPolynomTime,distributions))
                                STDTime.append(np.sqrt(np.abs(pc.Var(currPolynomTime, distributions))))
                                
                                meanAmp.append(pc.E(currPolynomAmp,distributions))
                                STDAmp.append(np.sqrt(np.abs(pc.Var(currPolynomAmp, distributions))))
                                                                
                                # calculate correlation
                                corrDistTime = []
                                corrDistAmp = []
                                for timepeakCount in timepeaks:
                                    timeCorr = pc.Corr([polynoms[0][count][quantity]['Time'][int(timepeakCount)],x[indexD]],distributions)[(0,1)]
                                    if abs(timeCorr) >= approxError: corrDistTime.append(timeCorr)
                                    else: corrDistTime.append(0.0) 
                                correlationTime.append(corrDistTime)
                                
                                for amppeakCount in amppeaks:
                                    ampCorr = pc.Corr([polynoms[0][count][quantity]['Amp'][int(amppeakCount)],x[indexD]],distributions)[(0,1)]
                                    if abs(ampCorr) >= approxError: corrDistAmp.append(ampCorr)
                                    else: corrDistAmp.append(0.0)
                                correlationAmp.append(corrDistAmp)
                                
                                # 100-plotPeaksConfidenceAlpha% confidence intervals:
                                confIntTime.append(pc.Perc(currPolynomTime,  [plotPeaksConfidenceAlpha/2., 100-plotPeaksConfidenceAlpha/2.],distributions))
                                confIntAmp.append(pc.Perc(currPolynomAmp,  [plotPeaksConfidenceAlpha/2., 100-plotPeaksConfidenceAlpha/2.], distributions))
                                
                            
#                             bifCount = max(number)                        
#                             for bif in bifurcationList:
#                                 
#                                 number.append(bifCount)
#                                 bifCount= bifCount+1
#                                 uncertainVarNames.append(str(bif))
#                                 #calculate mean and STD with respect to a certain random varibale
#                                 currDistMean = pc.E(distributions)
#                                 currDistMean[bif[0]] = np.nan
#                                 if bif[1] != None: currDistMean[bif[1]] = np.nan
#                                 if bif[2] != None: currDistMean[bif[2]] = np.nan
#                                 # reduce polynomials
#                                 currPolynomTime = polynoms[0][count][quantity]['Time'](*currDistMean)
#                                 currPolynomAmp = polynoms[0][count][quantity]['Amp'](*currDistMean)
#                                 
#                                 meanTime.append(pc.E(currPolynomTime,distributions))
#                                 STDTime.append(np.sqrt(np.abs(pc.Var(currPolynomTime, distributions))))
#                                 
#                                 meanAmp.append(pc.E(currPolynomAmp,distributions))
#                                 STDAmp.append(np.sqrt(np.abs(pc.Var(currPolynomAmp, distributions))))
#                                                            
#                                 # 100-plotPeaksConfidenceAlpha% confidence intervals:
#                                 confIntTime.append(pc.Perc(currPolynomTime,  [plotPeaksConfidenceAlpha/2., 100-plotPeaksConfidenceAlpha/2.],distributions))
#                                 confIntAmp.append(pc.Perc(currPolynomAmp,  [plotPeaksConfidenceAlpha/2., 100-plotPeaksConfidenceAlpha/2.], distributions))
                            # create peak sensitivity log file
                                               
                            logfile = open(''.join([plotDirectoryPeakPLogs,'Polychaos_peak_sensitivity','_',quantitiesNames[quantity],'_N.',str(order),'_method_',method,'_.txt']),'a')
                            logfile.write(''.join(['\n','point_',name,'\n']))
                            logfile.write(' Sensitivity of the peaks point %s, method %s, order %d, for %s \n'%(name,method,order,quantity))
                            logfile.write(' expected values: \n')
                            logfile.write(''.join([ '    time       ',' '.join(str(i) for i in peakTimeExpected),'\n','          STD: ',' '.join(str(i) for i in peakTimeSTD),'\n']))
                            logfile.write(''.join([ '            ConfInt',str(100-plotPeaksConfidenceAlpha),'%' ,str(peakconfIntTime),'\n']))
                            logfile.write(''.join([ '    amplitude  ',' '.join(str(i*unitScale[quantitiesNames[quantity]]) for i in peakAmpExpected),'\n','          STD: ',' '.join(str(i*unitScale[quantitiesNames[quantity]]) for i in peakAmpSTD), '\n']))
                            logfile.write(''.join([ '            ConfInt',str(100-plotPeaksConfidenceAlpha),'%' ,str(peakconfIntAmp*unitScale[quantitiesNames[quantity]]),'\n']))
                        
                            for indexD in number:
                                logfile.write(''.join(['\n    Sensitivity \n']))
                                logfile.write(''.join(['    ',uncertainVarNames[indexD],'\n']))
                                logfile.write(''.join(['    time       ',' \n']))
                                logfile.write(''.join([ '        Mean ',' '.join(str(i) for i in meanTime[indexD]),'\n']))
                                logfile.write(''.join([ '        STD    ',' '.join(str(i) for i in STDTime[indexD]),'\n']))
                                #logfile.write(''.join([ '        Corr   ',' '.join(str(i) for i in correlationTime[indexD]),'\n']))
                                logfile.write(''.join([ '        ConfInt',str(100-plotPeaksConfidenceAlpha),'%' ,' '.join(str(i) for i in confIntTime[indexD]),'\n']))
                                logfile.write(''.join(['    amplitude',' \n']))
                                logfile.write(''.join([ '        Mean ',' '.join(str(i*unitScale[quantitiesNames[quantity]]) for i in meanAmp[indexD]),'\n']))
                                logfile.write(''.join([ '        STD    ',' '.join(str(i*unitScale[quantitiesNames[quantity]]) for i in STDAmp[indexD]),'\n']))
                                #logfile.write(''.join([ '        Corr',' '.join(str(i) for i in correlationAmp[indexD]),'\n']))
                                logfile.write(''.join([ '        ConfInt',str(100-plotPeaksConfidenceAlpha),'%' ,' '.join(str(i*unitScale[quantitiesNames[quantity]]) for i in confIntAmp[indexD]),'\n']))
                            logfile.close() 
                        
                        
                            interSavedData[quantitiesNames[quantity]] = {}
                            #expected values
                            interSavedData[quantitiesNames[quantity]]['peakTimeExpected']=peakTimeExpected
                            interSavedData[quantitiesNames[quantity]]['peakTimeSTD']=peakTimeSTD
                            interSavedData[quantitiesNames[quantity]]['peakAmpExpected']=peakAmpExpected
                            interSavedData[quantitiesNames[quantity]]['peakAmpSTD']=peakAmpSTD
                            interSavedData[quantitiesNames[quantity]]['peakconfIntTime']=peakconfIntTime
                            interSavedData[quantitiesNames[quantity]]['peakconfIntAmp']=peakconfIntAmp
                           
                            interSavedData[quantitiesNames[quantity]]['meanTime']=meanTime
                            interSavedData[quantitiesNames[quantity]]['meanAmp']=meanAmp
                            interSavedData[quantitiesNames[quantity]]['STDTime']=STDTime
                            interSavedData[quantitiesNames[quantity]]['STDAmp']=STDAmp
                            interSavedData[quantitiesNames[quantity]]['sobolIndexTime'] = sobolIndexTime
                            interSavedData[quantitiesNames[quantity]]['sobolIndexAmp'] = sobolIndexAmp
                            interSavedData[quantitiesNames[quantity]]['correlationTime']=correlationTime
                            interSavedData[quantitiesNames[quantity]]['correlationAmp']=correlationAmp
                            interSavedData[quantitiesNames[quantity]]['confIntTime']=confIntTime
                            interSavedData[quantitiesNames[quantity]]['confIntAmp']=confIntAmp
                        
                        FILE = open(peakAnalysisFile,"w")
                        # store pickle and close file
                        cPickle.dump(interSavedData,FILE,protocol=2)
                        FILE.close()
                    
                    ###############################################################################################################
                    ###############################################################################################################
                    ### Creating Plots
                                         
                    print '    Creating Plots'                           
                    for quantity,analyticBool in quantitiesToPlotTemp.iteritems():     
                        
                        peakconfIntTime = interSavedData[quantitiesNames[quantity]]['peakconfIntTime']
                        peakconfIntAmp = interSavedData[quantitiesNames[quantity]]['peakconfIntAmp']
                        peakTimeExpected = interSavedData[quantitiesNames[quantity]]['peakTimeExpected']
                        peakTimeSTD = interSavedData[quantitiesNames[quantity]]['peakTimeSTD']
                        peakAmpExpected = interSavedData[quantitiesNames[quantity]]['peakAmpExpected']
                        peakAmpSTD = interSavedData[quantitiesNames[quantity]]['peakAmpSTD']
                        meanTime = interSavedData[quantitiesNames[quantity]]['meanTime']
                        meanAmp = interSavedData[quantitiesNames[quantity]]['meanAmp']
                        STDTime = interSavedData[quantitiesNames[quantity]]['STDTime']
                        STDAmp = interSavedData[quantitiesNames[quantity]]['STDAmp']
                        sobolIndexTime = interSavedData[quantitiesNames[quantity]]['sobolIndexTime']
                        sobolIndexAmp  = interSavedData[quantitiesNames[quantity]]['sobolIndexAmp']
                        correlationTime = interSavedData[quantitiesNames[quantity]]['correlationTime']
                        correlationAmp = interSavedData[quantitiesNames[quantity]]['correlationAmp']
                        confIntTime = interSavedData[quantitiesNames[quantity]]['confIntTime']
                        confIntAmp = interSavedData[quantitiesNames[quantity]]['confIntAmp']
                        
                        ###############################################################################################################
                        ### create peak sensitivity plots
                        if plotPeaksAnalyticSensitivity == True:
                            
                            N = 2*(max(number)+1)
                            ind = np.arange(N)    # the x locations for the groups
                            ind1 = np.arange(0,N-1,2)
                            ind14 = np.arange(0.4,N-0.6,2)
                            ind2 = np.arange(1,N,2)
                            ind24 = np.arange(1.4,N+0.4,2)
                            
                            width = 0.35       # the width of the bars: can also be len(x) sequence
                            
                            # fixing the error of the peak function
                            peakCounts = [0]
                            subNum = 1.
                            numPeak = ['1','2','2']
                            if len(correlationTime[0]) > 1:
                                if quantity == 'extremaPressure': peakCounts = [0,2]
                                if quantity == 'extremaFlow': peakCounts = [0,1]
                                subNum = 2.
                            
                            from matplotlib import rc
                            from matplotlib import rcParams
                            rcParams['text.usetex']=True
                            rcParams['text.latex.unicode']=True
                            fig = plt.figure(1, edgecolor='k', dpi = 800)
                            fig.subplots_adjust(hspace=0.5)   
                            fig.set_figwidth(8.27)
                            fig.set_figheight((11.69/3.)*subNum)
                            plt.clf()
                            
                            analyticTime = {'A':[[0.0,0.0,0.0],[],[1.0,0.0,0.0]],'B':[[1.0,0.0,0.0],[],[1.0,0.0,0.0]],'C':[[1.0,0.0,0.0]],'D':[[2.0/3.0,1.0/3.0,0.0]],'E':[[0.5,0.5,0.0]]}
                            analyticAmp = {'A':[[1.0,0.0,0.0],[],[0.25,0.375,0.375]],'B':[[1.0,0.0,0.0],[],[0.25,0.375,0.375]],'C':[[22.0/31.0,9.0/62.0,9.0/62.0]],'D':[[44.0/69.0,16.0/69.0,3.0/23.0]],'E':[[44.0/69.0,16.0/69.0,3.0/23.0]]}
                            
                            for peakCount in peakCounts:
                                plt.subplot(subNum,1,peakCounts.index(peakCount,)+1)
                                timeSensOfThePeak = np.empty((max(number)+1))
                                ampSensOfThePeak = np.empty((max(number)+1))
                                tickss = []
                                for indexD in number:  
                                    timeSensOfThePeak[indexD] = correlationTime[indexD][peakCount]
                                    ampSensOfThePeak[indexD] = correlationAmp[indexD][peakCount]
                                    tickss.append(uncertainVarNames[indexD])
                                
                                #calculate relative correlation
                                if sum(abs(timeSensOfThePeak)) == 0: sumTimeSensOfThePeak = 1.0
                                else: sumTimeSensOfThePeak = sum(abs(timeSensOfThePeak))
                                if sum(abs(ampSensOfThePeak)) == 0: sumAmpSensOfThePeak = 1.0
                                else: sumAmpSensOfThePeak = sum(abs(ampSensOfThePeak))
                                
                                p1 = plt.bar(ind1, timeSensOfThePeak,   width, color='r', label = 'occurence time')
                                p2 = plt.bar(ind2, ampSensOfThePeak,   width, color='b', label = 'amplitude')
                                if analyticBool: 
                                    p14 = plt.bar(ind14, [analyticTime[name][peakCount][0],analyticTime[name][peakCount][1],analyticTime[name][peakCount][2]],   width, color='m', label = 'analytic occ. time')
                                    p24 = plt.bar(ind24, [analyticAmp[name][peakCount][0],analyticAmp[name][peakCount][1],analyticAmp[name][peakCount][2]],   width, color='c', label = 'analytic amplitude')
                                else:
                                    p141 = plt.bar(ind14, abs(timeSensOfThePeak/sumTimeSensOfThePeak),   width, color='m', label = 'rel occ time')
                                    p242 = plt.bar(ind24, abs(ampSensOfThePeak/sumAmpSensOfThePeak),   width, color='c', label = 'rel amp')
                                
                                plt.plot([-0.5,6.5],[0,0], color='k')
                                plt.ylabel('Relative sensitivity in [-1,0]')
                                #plt.title('Relativ Sensitivity of '+quantitiesNames[quantity]+' - Peak '+numPeak[peakCount])
                                plt.xticks(0.9+ind1+width/2., tickss)
                                plt.axis([-0.5, 6.5,-1.0,2.0])
                                plt.yticks(np.arange(-1.0,1.2,0.2))
                                plt.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)
                                
                                savePath = ''.join([plotDirectoryPeakP_bar,'Polychaos_point_',name,'_',quantitiesNames[quantity],'_N',str(order),'_method_',method,'_sensitivity',plotFileType])
                                plt.savefig(savePath)
                                plt.clf()
                        
                        ############################################################################################################## 
                        # box plots of mean STD for peaks 
                        
                        a = True
                        if a == True:
                            peakCounts = peaksToEvaluate[name][quantity]
                            
                        for peakCount in peakCounts:
                            from matplotlib import rcParams
                            #rcParams['text.usetex']=True
                            #rcParams['text.latex.unicode']=True
                            fig = plt.figure(1, edgecolor='k',dpi=800)
                            fig.subplots_adjust(hspace=0.5)   
                            fig.set_figwidth(8.27)
                            if plotPeaksMeanSTDBoxPlotsSingle: fig.set_figheight((11.69/3.*1.1))
                            else: fig.set_figheight((11.69/3.*1.1)*len(number))
                            plt.clf()
                            
                            from matplotlib.font_manager import FontProperties

                            fontLegend = FontProperties()
                            fontLegend.set_size('small')
                                                        
                            ###################################################################################################
                            # calculate global y limits for time and amplitude
                            ymaxT = -1e64
                            yminT = 1e64
                            ymaxA = -1e64
                            yminA = 1e64
                                                        
                            for indexD in number:
                                ymaxT = max([ymaxT, max([meanTime[indexD][peakCount]+STDTime[indexD][peakCount],peakTimeExpected[peakCount]+peakTimeSTD[peakCount]])])
                                yminT = min([yminT, min([meanTime[indexD][peakCount]-STDTime[indexD][peakCount],peakTimeExpected[peakCount]-peakTimeSTD[peakCount]])])
                                ymaxA = max([ymaxA, max([(meanAmp[indexD][peakCount]+STDAmp[indexD][peakCount])*unitScale[quantitiesNames[quantity]],
                                                         (peakAmpExpected[peakCount]+peakAmpSTD[peakCount])*unitScale[quantitiesNames[quantity]]])])
                                yminA = min([yminA, min([(meanAmp[indexD][peakCount]-STDAmp[indexD][peakCount])*unitScale[quantitiesNames[quantity]],
                                                         (peakAmpExpected[peakCount]-peakAmpSTD[peakCount])*unitScale[quantitiesNames[quantity]]])])
                            
                            ycorrT = (ymaxT-yminT)*0.1
                            ymaxT  = ymaxT+ycorrT
                            yminT  = yminT-ycorrT
                            
                            ycorrA = (ymaxA-yminA)*0.1
                            ymaxA  = ymaxA+ycorrA
                            yminA  = yminA-2*ycorrA
                                                                              
                            subCount = 1
                            for indexD in number:
                                if plotPeaksMeanSTDBoxPlotsSingle: ax = plt.subplot(1,1,subCount)
                                else:                             
                                    ax = plt.subplot(len(number),1,subCount)
                                    subCount=subCount+1
                                
                                ###################################################################################################
                                ### Time
                                
                                mean = meanTime[indexD][peakCount]
                                STD = STDTime[indexD][peakCount]
                                expected = peakTimeExpected[peakCount]
                                expectedSTD = peakTimeSTD[peakCount]
                                                                
                                #expected 
                                ax.plot([0.55,0.95],[expected,expected],color='k' ,linestyle = '--',linewidth ='1.5', label='E(X) all')
                                ax.broken_barh([(0.6,0.3)],(expected-expectedSTD, 2.0*expectedSTD), linestyle='dashed',alpha=0.5,facecolors = 'yellow',edgecolor='grey',label="STD all" )    #this parameter
                                STDall = plt.Rectangle((0, 0), 0, 0, color='y', alpha=0.5, label='STD all' )
                                ax.add_patch(STDall)
                                # the current parameter
                                ax.broken_barh([(0.5,0.5)],(mean-STD, 2.0*STD), linestyle='solid',facecolors = 'None',edgecolor='b', label= 'STD parameter' )
                                
                                ax.plot([0.4,1.1],[mean,mean],color='r' ,linestyle = '-',linewidth ='2' ,label='E(X) parameter')
                                STDpara = plt.Rectangle((0, 0), 0, 0, color='b', label='STD parameter' )
                                ax.add_patch(STDpara)
                                                                
                                ax.set_ylabel('Time [s]',fontsize=13)
                                # Shink current axis's height by 10% on the bottom
                                box = ax.get_position()
                                ax.set_position([box.x0, box.y0 + box.height * 0.1,
                                                 box.width, box.height * 0.9])
                                
                                # Put a legend below current axis
                                ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),frameon=False, ncol=4,prop = fontLegend)
                                
                                ax.set_ylim(yminT,ymaxT)
                                
                                ###################################################################################################
                                ### Amplitude 
                                
                                ax2 = ax.twinx()
                                mean = meanAmp[indexD][peakCount]  *unitScale[quantitiesNames[quantity]]
                                STD = STDAmp[indexD][peakCount] *unitScale[quantitiesNames[quantity]]
                                expected = peakAmpExpected[peakCount] *unitScale[quantitiesNames[quantity]]
                                expectedSTD = peakAmpSTD[peakCount] *unitScale[quantitiesNames[quantity]]
                                
                                ## amplitude
                                #expected 
                                ax2.broken_barh([(1.6,0.3)],(expected-expectedSTD, 2.0*expectedSTD), linestyle='dashed',alpha=0.5,facecolors = 'yellow',edgecolor='grey',label='STD total')    #this parameter
                                #this parameter
                                ax2.broken_barh([(1.5,0.5)],(mean-STD, 2.0*STD), linestyle='solid',facecolors = 'None',edgecolor='b',label='E parameter' )
                                ax2.plot([1.4,2.1],[mean,mean],color='r' ,linestyle = '-',linewidth ='2' ,label='Mean parameter')
                                #expected 
                                ax2.plot([1.55,1.95],[expected,expected],color='k' ,linestyle = '--',linewidth ='1.5', label='E total')
                                                                
                                plotNameQuantity = quantitiesNames[quantity]
                                if "_" in plotNameQuantity: plotNameQuantity = plotNameQuantity.replace('_','-')
                                ax2.set_ylabel(' '.join([plotNameQuantity,quantitiesUnits[quantity]]),fontsize=13)

                                ax2.set_ylim(yminA,ymaxA)
                                
                                #ax2.yaxis.major.formatter.set_powerlimits((0,0)) 
                                plt.xlim(0,2.5)
                                plt.xticks([0.75,1.75], ['Occ. Time', 'Amplitude'])
                                #plt.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)
                                plt.title(' '.join(['Variable:',uncertainVarNames[indexD],'Peak:',str(peakCount)]),fontsize=13)  
                                
                                if plotPeaksMeanSTDBoxPlotsSingle == True:
                                    savePath = ''.join([plotDirectoryPeakP_box,'Polychaos_point_',name,'_',quantitiesNames[quantity],'_N',str(order),'_method_',method,'_mean_STD_peak',str(peakCount),'Variable:',uncertainVarNames[indexD],plotFileType])
                                    plt.savefig(savePath)
                                    plt.clf()
                               
                            if plotPeaksMeanSTDBoxPlotsSingle == False: 
                                savePath = ''.join([plotDirectoryPeakP_box,'Polychaos_point_',name,'_',quantitiesNames[quantity],'_N',str(order),'_method_',method,'_mean_STD_peak',str(peakCount),plotFileType])
                                #plt.savefig(savePath)
                                plt.clf()
                        
                            plt.clf()
                            
                            ############################################################################################
                            #-- plots with sobol Index STD-sensitivity
                            yValues = np.linspace(1,len(uncertainVarNames),len(uncertainVarNames))
                            rcParams['text.usetex']=True
                            rcParams['text.latex.unicode']=True
                            fig = plt.figure(1, edgecolor='k',dpi=800)
                            
                            fig.set_figheight((11.69/3.*1.3))
                            fig.set_figwidth(8.27*1.0)
                            
                            ax = plt.subplot(1,1,1)
                            
                            #fig.subplots_adjust(bottom=0.25)   
                            fig.subplots_adjust(top=0.90)
                            #fig.subplots_adjust(left=0.125)
                            fig.subplots_adjust(left=0.14)
                            fig.subplots_adjust(right=0.98) 
                            fig.subplots_adjust(bottom=0.12)   
                            fontSizeLabel = 18
                            
                            ax.axhline(y=0,color='k',linewidth = 1.5)
                                                                                   
                            #STDsensitivityAmp
                            plt.plot(yValues,sobolIndexTime.T[peakCount],color='k' ,linestyle = '',mfc = 'w',mec = 'k', marker =  '*',label=r'$Si^T$', linewidth = 10., markersize = 12) 
                            plt.plot(yValues,sobolIndexAmp.T[peakCount],color='k' ,linestyle = '', marker =  'o',label=r'$Si^A$', linewidth = 5., markersize = 7)                          
                            
                            ### for vessel data
                            uncertainVarVesselNames = []
                            for vesselId in vesselIds:
                                uncertainVarVesselNames.append(vascularNetwork.vessels[vesselId].name)
                            uncertainVarVesselNames = ['ascending aorta', 'aortic arch A', 'aortic arch B', 'thoracic A', 'thoracic B', 'abdominal A','abdominal B','abdominal C','abdominal D','abdominal E']
                            #uncertainVarVesselNames = ['asc.\n aorta','A','B','A','B','A','B','C','D','E']
                            plt.xticks(range(1,len(uncertainVarVesselNames)+1),uncertainVarVesselNames, fontsize = fontSizeLabel)
                            
                            ### for multiple
                            plt.xticks(range(1,len(uncertainVarNames)+1),uncertainVarNames, fontsize = fontSizeLabel)
                            
                            plt.yticks(fontsize = fontSizeLabel)
                            #sensitivity index 
                            ax.set_ylabel(r'$Si$', fontsize = fontSizeLabel)
                            ax.set_ylim(-0.03,1)
                            ax.set_xlim(0.,len(uncertainVarNames)+0.5)
                            ax.spines['top'].set_visible(False)
                            ax.tick_params(axis='x',top='off')
                            ax.spines['right'].set_visible(False)
                            ax.spines['bottom'].set_visible(False)
                            ax.tick_params(axis='y',right='off')
                                                      
                            box = ax.get_position()
                            ax.set_position([box.x0, box.y0 + box.height * 0.2,box.width, box.height * 0.8])
                            
                            quantitiesNameNow = { 'extremaPressure':'pressure wave','extremaPressure_b':'backward propagating pressure wave'}
                                                        
                            if peakCount == 0  : ax.set_title(quantitiesNameNow[quantity]+' - diastole ',fontsize = fontSizeLabel)
                            elif peakCount == 1 : ax.set_title(quantitiesNameNow[quantity]+' - systole ',fontsize = fontSizeLabel)
                            elif peakCount == 2 : ax.set_title(quantitiesNameNow[quantity]+' - shoulder ',fontsize = fontSizeLabel)
                            elif peakCount == 3 : ax.set_title(quantitiesNameNow[quantity]+'\\ - point of inflection',fontsize = fontSizeLabel)
                            
                            print quantitiesNameNow[quantity],peakCount
                            print 'time'
                            print sobolIndexTime.T[peakCount]
                            print 'amp'
                            print sobolIndexAmp.T[peakCount]
                            print 
                            
                            ax.legend(loc='upper right', frameon=False, ncol=1,fontsize = fontSizeLabel, numpoints=1)
                            
                            plt.setp(plt.xticks()[1], rotation=30, ha='right')
                                                        
                            savePath = ''.join([plotDirectoryPeakP_box,'Polychaos_point_',name,'_',quantitiesNames[quantity],'_N',str(order),'_method_',method,'_relative_STD_peak',str(peakCount),plotFileType])
                            plt.savefig(savePath,dpi=800)
                            savePath = ''.join([plotDirectoryPeakP_box,'Polychaos_point_',name,'_',quantitiesNames[quantity],'_N',str(order),'_method_',method,'_relative_STD_peak',str(peakCount),'.pdf'])
                            plt.savefig(savePath,dpi=800)
                            plt.clf()
                            
                            
                        ##############################################################################################################
                        # plots with mean and confidence intervall for peaks 
                                                
                        from matplotlib import rc
                        from matplotlib import rcParams
                        rcParams['text.usetex']=True
                        rcParams['text.latex.unicode']=True
                        fig = plt.figure(1, edgecolor='k')
                        fig.subplots_adjust(hspace=0.5)   
                        fig.set_figwidth(8.27)
                        fig.set_figheight((11.69/3.)*len(number))
                        plt.clf()
                        
                        # polynomial solution 
                        # signal with mean of current investigation point
                        distMean = pc.E(distributions)
                        polynomWithMeans = pc.E(polynoms[0][count][quantitiesNames[quantity]],distributions) #(*distMean)
                        timeValuesPolynom = np.linspace(0,totalTime,len(polynomWithMeans))
                        
                        print quantity
                        
                        linstyleO = '-'
                        # main signal if forward and backward
                        if '_f' in quantity or '_b' in quantity:
                            pureQuantity = quantitiesNames[quantity].split('_')[0]
                        
                        
                            polynomWithMeansTotal = pc.E(polynoms[0][count][pureQuantity],distributions) # (*distMean)
                            timeValuesPolynomTotal = np.linspace(0,totalTime,len(polynomWithMeansTotal))
                                                    
                            saveConfIntTimeCur = interSavedData[pureQuantity]['confIntTime']# interSavedData[pureQuantity]['peakconfIntTime'] # #interval due to parameter
                            saveConfIntAmpCur = interSavedData[pureQuantity]['confIntAmp'] #interSavedData[pureQuantity]['peakconfIntAmp'] # #interval due to parameter ## if changed back don forget [indexD] downstaris ;=
                            peakCountsPure = peaksToEvaluate[name][''.join(['extrema',pureQuantity])]
                            
                            peakCounts = peaksToEvaluate[name][quantity]
                            varColors =['r','m','g','b','r','m','g','b']
                            
                            initialPressureCorrection = polynomWithMeansTotal[0] *unitScale[quantitiesNames[quantity]] #initialPressureCorrection = vascularNetwork.globalFluid['initialPressure']*unitScale['Pressure']
                            
                            print number, peakCounts
                            for indexD in number:
                                ax = plt.subplot(max(number)+1,1,indexD+1)
                                
                                ax2 = ax.twinx()
                                ax.plot(timeValuesPolynomTotal, polynomWithMeansTotal *unitScale[quantitiesNames[quantity]],color='k' ,linestyle = '-', linewidth = 1.0, label = 'Expected values from simulation with uncertain parameters total signal')
                                ax2.plot(timeValuesPolynom, polynomWithMeans *unitScale[quantitiesNames[quantity]],color='k' ,linestyle = '--', linewidth = 1.0, label = 'Expected values from simulation with uncertain parameters')
                                    
                                
                                
                                
                                for peakCount in peakCounts:
                                    currAmpVal = []
                                    currTimeVal = []
                                    
                                    saveAmpVal  = []
                                    saveTimeVal = []
                                    
                                    #take the first ones
                                    for peakCountLarge in peakCounts:
                                        currAmpVal.append(confIntAmp[indexD][peakCountLarge][peakCount]*unitScale[quantitiesNames[quantity]])
                                        currTimeVal.append(confIntTime[indexD][peakCountLarge][peakCount])
                                        
                                        saveAmpVal.append( saveConfIntAmpCur[indexD][peakCountLarge][peakCount]*unitScale[quantitiesNames[quantity]])
                                        saveTimeVal.append(saveConfIntTimeCur[indexD][peakCountLarge][peakCount])
                                        
                                    #print currAmpVal
                                    #print 
                                    #peakExpTplot = [peakTimeExpected[peakCount],peakTimeExpected[peakCount]]
                                    #peakExpAplot = [peakAmpExpected[peakCount],peakAmpExpected[peakCount]]
                                    
                                    #currAmpVal = confIntAmp[indexD][peakCount]
                                    #currTimeVal = confIntTime[indexD][peakCount]
                                    
                                    ax2.plot(currTimeVal,currAmpVal, color=varColors[indexD], linewidth =1.5,marker = 'o',markersize=0.5)
                                    
                                    #for peakCount in peakCounts:                                        
                                    ax.plot(saveTimeVal,saveAmpVal, color=varColors[indexD], linewidth = 1)
                                    
                                    
                                ax.set_ylabel(' '.join([pureQuantity,quantitiesUnits[''.join(['extrema',pureQuantity])]]),fontsize=13)
                                ax2.set_ylabel(' '.join([plotNameQuantity,quantitiesUnits[quantity]]),fontsize=13)
                                                           
                                ax.set_xlabel('Time [s]')
                                
                                ax.set_title (uncertainVarNames[indexD])
                                
                                ax.set_xlim(plotxTimeStart,plotxTimeEnd)
                                ax2.set_xlim(plotxTimeStart,plotxTimeEnd)
                                
                                ## Align y1 and y2 axis for Flow
                                if pureQuantity == "Flow":
                                    y1,y2 = ax.get_ylim()
                                    y3,y4 = ax2.get_ylim()
                                    yt1 = min([y1,y3])
                                    yt2 = max([y2,y4])
                                    ax.set_ylim(yt1,yt2)
                                    ax2.set_ylim(yt1,yt2)
                                    if limits[quantitiesNames[quantity]] != [0,0]: ax.set_ylim(limits[quantitiesNames[quantity]])
                                
                                ## Align y1 and y2 axis for pressure
                                else:
                                    y1,y2 = ax.get_ylim()
                                    y3,y4 = ax2.get_ylim()
                                    yt1 = min([y1,initialPressureCorrection+y3])
                                    yt2 = max([y2,initialPressureCorrection+y4])
                                    ax.set_ylim(yt1,yt2)
                                    ax2.set_ylim(yt1-initialPressureCorrection,yt2-initialPressureCorrection)
                                    if limits[quantitiesNames[quantity]] != [0,0]: ax.set_ylim(limits[quantitiesNames[quantity]])
                                
                                ax.spines['top'].set_visible(False)
                                ax.tick_params(axis='x',top='off')
                                    
                                #ax2.set_ylim(y1,y2)
                                    
                            savePath = ''.join([plotDirectoryPeakP_Conf,'Polychaos_point_',name,'_',quantitiesNames[quantity],'_N.',str(order),'_method_',method,'_',str(100-plotPeaksConfidenceAlpha),'%_confidence_randomVariable_',str(indexD),plotFileType])
                            plt.savefig(savePath)
                            plt.clf()
        
#            ##########################################################################################  
#            if plotSensitiviy == True:
#                if not os.path.exists(plotDirectorySensitivity):
#                    os.makedirs(plotDirectorySensitivity) 
#                
#                print '   plotting: sensitivity'
#            ##                
#            ##                polyBeta1 = polynoms[0][count]['Pressure'](y=796020.0/orderOfMagnitudeFactor1,z=796020.0/orderOfMagnitudeFactor1)
#            ##                SensBeta1Exp = pc.E(polyBeta1,distributions)
#            ##                SensBeta1STD = np.sqrt(np.abs(pc.Var(polyBeta1,distributions)))
#            ##                
#            ##                polyBeta2 = polynoms[0][count]['Pressure'](x=324970.0/orderOfMagnitudeFactor1,z=796020.0/orderOfMagnitudeFactor1)
#            ##                SensBeta2Exp = pc.E(polyBeta2,distributions)
#            ##                SensBeta2STD = np.sqrt(np.abs(pc.Var(polyBeta2,distributions)))
#            ##                
#            ##                polyBeta3 = polynoms[0][count]['Pressure'](x=324970.0/orderOfMagnitudeFactor1,y=796020.0/orderOfMagnitudeFactor1)
#            ##                SensBeta3Exp = pc.E(polyBeta3,distributions)
#            ##                SensBeta3STD = np.sqrt(np.abs(pc.Var(polyBeta3,distributions)))
#            ##                
#            ##                timeValuesPolynom = np.linspace(0,vascularNetwork.simulationContext['totalTime'],len(SensBeta1STD))
#            ##                
#            ##                fig = plt.figure(1, edgecolor='k')
#            ##                fig.subplots_adjust(hspace=0.5)   
#            ##                fig.set_figwidth(8.27)
#            ##                fig.set_figheight((11.69/3))
#            ##                
#            ##                plt.plot(timeValuesPolynom, SensBeta1STD ,color='b' ,linestyle = '-', linewidth = 1, label = 'Sensitivity beta 1 (STD)')
#            ##                plt.plot(timeValuesPolynom, SensBeta2STD ,color='r' ,linestyle = '-', linewidth = 1, label = 'Sensitivity beta 2 (STD)')
#            ##                plt.plot(timeValuesPolynom, SensBeta3STD ,color='r' ,linestyle = '-', linewidth = 1, label = 'Sensitivity beta 3 (STD)')
#            ##         
#            ##                
#            ##                plt.xlabel('Time (s)')
#            ##                plt.ylabel('Pressure (Pa)')
#            ##                plt.axis([ 0,0.5,-5.0,14])
#            ##                #plt.title ('Pressure ')
#            ##                plt.grid(True,linestyle='-',color='0.1',linewidth = 0.4)    
#            ##                plt.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)
#            ##                ax = fig.gca() 
#            ##                ax.yaxis.major.formatter.set_powerlimits((0,0))
#            ##                savePath = ''.join([plotDirectorySensitivity,'Polychaos_point_',name,'_N.',str(order),'_method_',method,'_sensitivity.png'])
#            ##                plt.savefig(savePath)
#            ##                fig.clear()
            plt.close()
            # free memory
            vascularNetwork.__del__()
            del(vascularNetwork)

def printIntervalData( printList, listName, printArguments, printTrue = True):
    '''
    This Function creats nice console output of data
    '''
    number,nameLength,numberVessel,vesselIds,varNames,numberBC,BCvesselID,BCtype,numberFluid,FluidType,numberNetworkSolver,uncertainParameterNetworkSolver,numberMulti,multiVesselIDS,uncertainParameter,uncertainVarNames = printArguments
    
    if len(printList) != len(number):
        print " WARNING: print List has wrong length!"
        return
    string1 = ''.join(['Nr.'.center(3),'|','Id'.center(3),'|','Variable'.ljust(nameLength),'|',listName,'\n'])
    if printTrue: print 'Nr.'.center(3),'|','Id'.center(3),'|','Variable'.ljust(nameLength),'|',listName
    string2 = ''

    # Vessel data
    for num in numberVessel:
        stringCur = ''.join([str(number[num]).center(3),'|',str(vesselIds[num]).center(3),'|',str(varNames[num]).ljust(nameLength),'|',str(printList[num]),'\n'])
        string2 = ''.join([string2,stringCur])
        if printTrue: print str(number[num]).center(3),'|',str(vesselIds[num]).center(3),'|',str(varNames[num]).ljust(nameLength),'|',str(printList[num])
    
    # BC conditions
    splitNumber = len(numberVessel)
    for num in numberBC:
        stringCur = ''.join([str(number[num+splitNumber]).center(3) ,'|',str(BCvesselID[num]).center(3),'|',str(BCtype[num]).ljust(nameLength),  '|',str(printList[num+splitNumber]),'\n'])
        string2 = ''.join([string2,stringCur])
        if printTrue: print str(number[num+splitNumber]).center(3),'|',str(BCvesselID[num]).center(3),'|',str(BCtype[num]).ljust(nameLength),'|',str(printList[num+splitNumber])
    splitNumber = splitNumber+ len(numberBC)
    
    # global FLUID properties
    for num in numberFluid:
        stringCur = ''.join([str(number[num+splitNumber]).center(3) ,'|','gF'.center(3),'|',str(FluidType[num]).ljust(nameLength),  '|',str(printList[num+splitNumber]),'\n'])
        string2 = ''.join([string2,stringCur])
        if printTrue: print str(number[num+splitNumber]).center(3),'|','gF'.center(3),'|',str(FluidType[num]).ljust(nameLength),'|',str(printList[num+splitNumber])
    splitNumber = splitNumber+ len(numberFluid)   
        
    #networkSolverCount
    # global NETWORK SOLVER properties
    for num in numberNetworkSolver:
        stringCur = ''.join([str(number[num+splitNumber]).center(3) ,'|','nS'.center(3),'|',str(uncertainParameterNetworkSolver[num]).ljust(nameLength),  '|',str(printList[num+splitNumber]),'\n'])
        string2 = ''.join([string2,stringCur])
        if printTrue: print str(number[num+splitNumber]).center(3),'|','nS'.center(3),'|',str(uncertainParameterNetworkSolver[num]).ljust(nameLength),'|',str(printList[num+splitNumber])
    splitNumber = splitNumber+ len(numberNetworkSolver)
    
    # Multi random variables
    for num in numberMulti:
        for multiVesselID in multiVesselIDS[num]:
            stringCur = ''.join([str(number[num+splitNumber]).center(3),'|',''.join(['Multiple ',str(numberMulti[num])]).center(3),'|',str(uncertainVarNames[num+splitNumber]).ljust(nameLength),'|',str(printList[num+splitNumber]),'\n'])
            string2 = ''.join([string2,stringCur])
        if printTrue: print str(number[num+splitNumber]).center(3),'|',''.join(['Multiple ',str(numberMulti[num])]).center(3),'|',str(uncertainVarNames[num+splitNumber]).ljust(nameLength),'|',str(printList[num+splitNumber])        
        
    if printTrue: print "------------------------------------------------- \n"
    return ''.join([string1,string2,"------------------------------------------------- \n\n"])


if __name__ == '__main__':
    main()