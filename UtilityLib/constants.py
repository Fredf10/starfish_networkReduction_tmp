
## Units
import numpy as np
##########################################################################################
#----------------------------------------------------------------------------------------# 
# Conversion
def unitConversion(unitsDict,value,unit):
    try: 
        value = float(value)
    except:
        if value is None: return value
        else: print 'ERROR: unitConversion not feasible: value {} not convertable to float'.format(value)
    if ' ' in unit:
        unit = unit.split(' ')
        for un in range (0,len(unit),1): value = value*unitsDict[unit[un]]
    else: 
        value = value*unitsDict[unit]
    return value

#----------------------------------------------------------------------------------------# 
#convert unit to SI system
unitsDictSI = {'m':1.0,
            'cm':1.0E-2,
            'mm':1.0E-3,
            'm-1':1.0,
            'cm-1':1.0E2,
            'mm-1':1.0E3,
            'm2':1.0,
            'cm2':1.E-4,
            'mm2':1.E-6,
            'm-2':1.0,
            'cm-2':1.E4,
            'mm-2':1.E6,
            'm3':1.0,
            'cm3':1E-6,
            'mm3':1E-9,
            'm-3':1.0,
            'cm-3':1.0E6,
            'mm-3':1.0E9,
            'Pa':1.0,
            'Pa-1':1.0,
            'mmHg':133.32,
            'mmHg-1':1.0/133.32,
            'ml': 1E-6,
            'ml-1':1E6,
            'N':  1.0,
            'N-1':1.0,
            's':  1.0,
            's2': 1.0,
            's-2':1.0,
            's-1':1.0,
            'min': 60.0,
            'min-1': 1.0/60.0,
            'min2': 3600.0,
            'min-2': 1.0/3600.0,
            'e6':1.E6,
            'e5':1.E5,
            'e4':1.E4,
            'e3':1.0E3, 
            'e2':1.0E2,
            'e-6':1.E-6,
            'e-5':1.E-5,
            'e-4':1.E-4,
            'e-3':1.0E-3, 
            'e-2':1.0E-2,
            '' : 1.0,
            ' ':1.0,
            'kg': 1.0,
            'g': 1.E-3,
            'mg': 1.E-6,
            'rad':1.0,
            'deg':np.pi/180.
            }

#----------------------------------------------------------------------------------------# 
#convert unit to medical system  Medical
unitsDictMed = {'m':1.0E2,
                'cm':1.0,
                'mm':1.0E-1,
                'm-1':1.0E-2,
                'cm-1':1.0,
                'mm-1':1.0E1,
                'm2':1.0E4,
                'cm2':1.0,
                'mm2':1.E-2,
                'm-2':1.E-4,
                'cm-2':1.0,
                'mm-2':1.E2,
                'm3':1.E6,
                'cm3':1.0,
                'mm3':1.E-3,
                'm-3':1.E-6,
                'cm-3':1.0,
                'mm-3':1.E3,
                'Pa':1.0/133.32,
                'Pa-1':133.32,
                'mmHg':1.0,
                'mmHg-1':1.0,
                'ml': 1.0,
                'ml-1':1.0,
                'N':  1.0,
                'N-1':1.0,
                's':  1.0,
                's2': 1.0,
                's-2':1.0,
                's-1':1.0,
                'min': 60.0,
                'min-1': 1.0/60.0,
                'min2': 3600.0,
                'min-2': 1.0/3600.0,
                'e6':1.E6,
                'e5':1.E5,
                'e4':1.E4,
                'e3':1.0E3, 
                'e2':1.0E2,
                'e-6':1.E-6,
                'e-5':1.E-5,
                'e-4':1.E-4,
                'e-3':1.0E-3, 
                'e-2':1.0E-2,
                '' : 1.0,
                ' ':1.0,
                'kg': 1.0E3,
                'g': 1.0,
                'mg': 1.E-3,
                'rad':1.0,
                'deg':np.pi/180.0
                }

#----------------------------------------------------------------------------------------# 
# units dict for all variables in the system
# based in SI units
variableUnitsSI = { # SimulationContext
                    'CFL': '',
                    'totalTime': 's',
                    'NumScheme': '',
                    'EqSystem': '',
                    # initialisationControls
                    'initMeanFlow':'m3 s-1',
                    'initMeanPressure':'Pa',
                    'compTotalSys':'m3 Pa-1',
                    'centralVenousPressure': 'Pa',
                    'minimumVenousPressure': '',
                    'gravityConstant': 'm s-2',
                    'angleXSystem': '',
                    'gravitationalField': '',
                   # BoundaryConditions
                    'Rt': '',
                    'RtOpen': '',
                    'Topen1': 's',
                    'Topen2': 's',
                    'RtClosed': '',
                    'Tclosed1': 's',
                    'Tclosed2':'s',
                    'Rc':'Pa s m-3',
                    'C':'m3 Pa-1',
                    'R':'Pa s m-3',
                    'RT':'Pa s m-3',
                    'Z':'Pa s m-3',
                    'Rtotal':'Pa s m-3',
                    'amp':'m3 s-1',
                    'ampConst':'m3 s-1',
                    'Tpulse': 's',
                    'Tspace': 's',
                    'Npulse': '',
                    'scale': '',
                    'freq': 's-1',
                    'runtimeEvaluation': '',
                    'T' : 's',
                    'Emax' : 'Pa m-3',
                    'Emin' : 'Pa m-3',
                    'Tpeak' : 's',
                    'V0' : 'm3',
                    'K'  : 's m-3',
                    'filePathName' :'',
                   # Vessel
                    'id': '',
                    'name': '',
                    'start': '',
                    'end': '',
                    'leftDaughter': '',
                    'rightDaughter': '',
                    'angleYMother': '',
                    'geom': '', 
                    'length': 'm',
                    'radiusA': 'm',
                    'radiusB': 'm',
                    'N': '',
                    'complianceType': '', 
                    'Pfunc': '',
                    'Ps': 'Pa',
                    'As': 'm2',
                    'wallThickness': 'm', 
                    'youngModulus': 'Pa',
                    #'beta': 'Pa m-1',
                    'beta': '',
                    # Fluid properties
                    'my': 'Pa s',
                    'rho': 'kg m-3',
                    'gamma': '',
                    # Plot Variables,
                    'Flow': 'm3 s-1',
                    'Pressure':'Pa'}

#----------------------------------------------------------------------------------------# 
# based on medical units
variableUnitsMed  = {
                    # SimulationContext
                    'CFL': '',
                    'totalTime': 's',
                    'NumScheme': '',
                    'EqSystem': '',
                    # initialisationControls
                    'initMeanFlow':'ml s-1',
                    'initMeanPressure':'mmHg',
                    'compTotalSys':'cm3 mmHg-1',
                    'centralVenousPressure': 'Pa',
                    'minimumVenousPressure': '',
                    'gravityConstant': 'm s-2',
                    'angleXSystem': '',
                    'gravitationalField': '',
                    # BoundaryConditions
                    'Rt': '',
                    'RtOpen': '',
                    'Topen1': 's',
                    'Topen2': 's',
                    'RtClosed': '',
                    'Tclosed1': 's',
                    'Tclosed2':'s',
                    'Rc':'mmHg s cm-3',
                    'C':'cm3 mmHg-1',
                    'R':'mmHg s cm-3',
                    'RT':'mmHg s cm-3',
                    'Z':'mmHg s cm-3',
                    'Rtotal':'mmHg s cm-3',
                    'amp':'ml s',
                    'ampConst':'ml s',
                    'Tpulse': 's',
                    'Tspace': 's',
                    'Npulse': '',
                    'scale': '',
                    'freq': 's-1',
                    'runtimeEvaluation': '',
                    'T' : 's',
                    'Emax' : 'mmHg cm-3',
                    'Emin' : 'mmHg cm-3',
                    'Tpeak' : 's',
                    'V0' : 'cm3',
                    'K'  : 's cm-3',
                    'filePathName' :'',
                   # Vessel
                    'id': '',
                    'name': '',
                    'start': '',
                    'end': '',
                    'leftDaughter': '',
                    'rightDaughter': '',
                    'angleYMother': '',
                    'geom': '', 
                    'length': 'cm',
                    'radiusA': 'cm',
                    'radiusB': 'cm',
                    'N': '',
                    'complianceType': '', 
                    'Pfunc': '',
                    'Ps': 'mmHg',
                    'As': 'cm2',
                    'wallThickness': 'cm', 
                    'youngModulus': 'mmHg',
                    'beta': 'mmHg cm-1',
                    # Fluid properties
                    'my': 'mmHg s',
                    'rho': 'g cm-3',
                    'gamma': '',
                    'centralVenousPressure': 'mmHg',
                    'minimumVenousPressure': 'mmHg',
                    # Plot Variables,
                    'Flow': 'ml s-1',
                    'Pressure':'mmHg'}

##########################################################################################
#### definition of all input/outut variables


variablesDict = {## class Vascular Network
                 # simulation context
                 'totalTime'                : {'type':'float', 'unitSI': 's',  'strCases': None, 'multiVar': False},
                 'CFL'                      : {'type':'float', 'unitSI': None, 'strCases': None, 'multiVar': False},
                 'cycleMode'                : {'type':'bool',  'unitSI': None, 'strCases': None, 'multiVar': False},
                 'numberCycles'             : {'type':'int',   'unitSI': None, 'strCases': None, 'multiVar': False},
                 'numberSavedCycles'        : {'type':'int',   'unitSI': None, 'strCases': None, 'multiVar': False},
                 'gravitationalField'       : {'type':'bool',  'unitSI': None,      'strCases': None, 'multiVar': False},
                 'gravityConstant'          : {'type':'float', 'unitSI': 'm s-2',   'strCases': None, 'multiVar': False},
                 'centralVenousPressure'    : {'type':'float',      'unitSI': 'Pa',      'strCases': None, 'multiVar': False},
                 'minimumVenousPressure'    : {'type':'float None', 'unitSI': 'Pa',      'strCases': None, 'multiVar': False},
                 # solver calibration
                 'solvingSchemeField'       : {'type':'str',   'unitSI': None, 'strCases': ['MacCormack'], 'multiVar': False},
                 'solvingSchemeConnections' : {'type':'str',   'unitSI': None, 'strCases': ['Linear'], 'multiVar': False},
                 'rigidAreas'               : {'type':'bool',  'unitSI': None, 'strCases': None, 'multiVar': False},
                 'simplifyEigenvalues'      : {'type':'bool',  'unitSI': None, 'strCases': None, 'multiVar': False},
                 'riemannInvariantUnitBase' : {'type':'str',   'unitSI': None, 'strCases': ['Flow', 'Pressure'], 'multiVar': False},
                 'automaticGridAdaptation'  : {'type':'bool',  'unitSI': None, 'strCases': None, 'multiVar': False},
                 # initialisationControl    
                 'initialsationMethod'          : {'type':'str',        'unitSI': None,      'strCases': ['Auto','MeanFlow','MeanPressure','ConstantPressure'], 'multiVar': False},
                 'initMeanFlow'                 : {'type':'float',      'unitSI': 'm3 s-1',  'strCases': None, 'multiVar': False},
                 'initMeanPressure'             : {'type':'float',      'unitSI': 'Pa',      'strCases': None, 'multiVar': False},
                 'estimateWindkesselCompliance' : {'type':'str',        'unitSI': None,      'strCases': ['Tree','Wk3','System','No'], 'multiVar': False},
                 'compPercentageWK3'            : {'type':'float',      'unitSI': None,      'strCases': None, 'multiVar': False},
                 'compPercentageTree'           : {'type':'float',      'unitSI': None,      'strCases': None, 'multiVar': False},
                 'compTotalSys'                 : {'type':'float',      'unitSI': 'm3 Pa-1', 'strCases': None, 'multiVar': False},
                 # global fluid
                 'my'                : {'type':'float',      'unitSI': 'Pa s',   'strCases': None, 'multiVar': False},
                 'rho'               : {'type':'float',      'unitSI': 'kg m-3', 'strCases': None, 'multiVar': False},
                 'gamma'             : {'type':'float None', 'unitSI':  None,    'strCases': None, 'multiVar': False},
                 ## class Communicators
                 'comType'           : {'type':'str', 'unitSI': None, 'strCases': ['CommunicatorRealTimeViz','CommunicatorBaroreceptor'], 'multiVar': False},
                 'comId'             : {'type':'int', 'unitSI': None, 'strCases': None, 'multiVar': False},
                 'vesselId'          : {'type':'int', 'unitSI': None, 'strCases': None, 'multiVar': False},
                 'node'              : {'type':'int', 'unitSI': None, 'strCases': None, 'multiVar': False},
                 'dn'                : {'type':'int', 'unitSI': None, 'strCases': None, 'multiVar': False},
                 'quantitiesToPlot'  : {'type':'str', 'unitSI': None, 'strCases': ['Pressure','Area','Flow'], 'multiVar': True},
                 ## boundary conditions
                 # flux bc
                 'amp'               : {'type':'float',  'unitSI': 'm3 s-1',    'strCases': None, 'multiVar': False},
                 'ampConst'          : {'type':'float',  'unitSI': 'm3 s-1',    'strCases': None, 'multiVar': False},
                 'Npulse'            : {'type':'float',  'unitSI': None,        'strCases': None, 'multiVar': False},
                 'Tpulse'            : {'type':'float',  'unitSI': 's',         'strCases': None, 'multiVar': False},
                 'freq'              : {'type':'float',  'unitSI': 's-1',       'strCases': None, 'multiVar': False},
                 'Tspace'            : {'type':'float',  'unitSI': 's',         'strCases': None, 'multiVar': False},
                 'runtimeEvaluation' : {'type':'bool',   'unitSI': None,        'strCases': None, 'multiVar': False},
                 'filePathName'      : {'type':'str',    'unitSI': None,        'strCases': ['anything'], 'multiVar': False},
                 'prescribe'         : {'type':'str',    'unitSI': None,        'strCases': ['influx','total'], 'multiVar': False},
                 'gaussC'            : {'type':'float',  'unitSI': None,        'strCases': None, 'multiVar': False},
                 # reflection coefficients
                 'Rt'                : {'type':'float',  'unitSI': None,        'strCases': None, 'multiVar': False},
                 'RtOpen'            : {'type':'float',  'unitSI': None,        'strCases': None, 'multiVar': False},
                 'Topen1'            : {'type':'float',  'unitSI': 's',         'strCases': None, 'multiVar': False},
                 'Topen2'            : {'type':'float',  'unitSI': 's',         'strCases': None, 'multiVar': False},
                 'RtClosed'          : {'type':'float',  'unitSI': None,        'strCases': None, 'multiVar': False},
                 'Tclosed1'          : {'type':'float',  'unitSI': 's',         'strCases': None, 'multiVar': False},
                 'Tclosed2'          : {'type':'float',  'unitSI': 's',         'strCases': None, 'multiVar': False},
                 # windkessel 
                 'Rtotal'            : {'type':'float None',      'unitSI': 'Pa s m-3', 'strCases': None, 'multiVar': False},
                 'Rc'                : {'type':'float str None',  'unitSI': 'Pa s m-3', 'strCases': ['VesselImpedance'],  'multiVar': False},
                 'Z'                 : {'type':'float str',       'unitSI': 'Pa s m-3', 'strCases': ['VesselImpedance'], 'multiVar': False},
                 'C'                 : {'type':'float',           'unitSI': 'm3 Pa-1',  'strCases': None, 'multiVar': False},
                 # heart model
                 'T'                 : {'type':'float',  'unitSI': 's',         'strCases': None, 'multiVar': False},
                 'Emax'              : {'type':'float',  'unitSI': 'Pa m-3',    'strCases': None, 'multiVar': False},
                 'Emin'              : {'type':'float',  'unitSI': 'Pa m-3',    'strCases': None, 'multiVar': False},
                 'Tpeak'             : {'type':'float',  'unitSI': 's',         'strCases': None, 'multiVar': False},
                 'V0'                : {'type':'float',  'unitSI': 'm3',        'strCases': None, 'multiVar': False},
                 'K'                 : {'type':'float',  'unitSI': 's m-3',     'strCases': None, 'multiVar': False},
                 ## class vessel
                 # attributes
                 'name'              : {'type':'str',  'unitSI': None,     'strCases': ['anything'], 'multiVar': False},
                 'Id'                : {'type':'int',  'unitSI': None,     'strCases': None, 'multiVar': False},
                 # topology
                 'leftDaughter'      : {'type':'int None',  'unitSI': None,     'strCases': None, 'multiVar': False},
                 'rightDaughter'     : {'type':'int None',  'unitSI': None,     'strCases': None, 'multiVar': False},
                 'angleYMother'      : {'type':'float None', 'unitSI': 'rad',    'strCases': None, 'multiVar': False},
                 # geometry
                 'geometryType'      : {'type':'str',    'unitSI': None,    'strCases': ['uniform','cone','constriction'], 'multiVar': False},
                 'length'            : {'type':'float',  'unitSI': 'm',     'strCases': None, 'multiVar': False},
                 'radiusProximal'    : {'type':'float',  'unitSI': 'm',     'strCases': None, 'multiVar': False},
                 'radiusDistal'      : {'type':'float',  'unitSI': 'm',     'strCases': None, 'multiVar': False},
                 'N'                 : {'type':'float',  'unitSI': None,    'strCases': None, 'multiVar': False},
                 # compliances
                 'complianceType'    : {'type':'str',    'unitSI': None,   'strCases': ['Laplace', 'Laplace2', 'Exponential', 'Hayashi', 'Reymond'], 'multiVar': False},
                 'constantCompliance': {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},
                 'externalPressure'  : {'type':'float',  'unitSI': 'Pa',   'strCases': None, 'multiVar': False},
                 'Ps'                : {'type':'float',  'unitSI': 'Pa',   'strCases': None, 'multiVar': False},
                 'As'                : {'type':'float None','unitSI': 'm2','strCases': None, 'multiVar': False},
                 'wallThickness'     : {'type':'float',  'unitSI': 'm',    'strCases': None, 'multiVar': False},
                 'youngModulus'      : {'type':'float',  'unitSI': 'm2',   'strCases': None, 'multiVar': False},
                 'beta'              : {'type':'float',  'unitSI': 'Pa',   'strCases': None, 'multiVar': False},
                 'betaLaplace'       : {'type':'float',  'unitSI': 'Pa',   'strCases': None, 'multiVar': False},
                 'betaExponential'   : {'type':'float',  'unitSI': None,     'strCases': None, 'multiVar': False},
                 'betaHayashi'       : {'type':'float',  'unitSI': 'Pa',     'strCases': None, 'multiVar': False},
                 'distensibility'    : {'type':'float',  'unitSI': 'Pa-1', 'strCases': None, 'multiVar': False},
                 # fluid 
                 'applyGlobalFluid'  : {'type':'bool',   'unitSI': None,   'strCases': None, 'multiVar': False},
                 }

###### XML description version 4.0
##########################################################################################
## class VascularNetwork

# simulationContext
simulationContextElements = ['totalTime',
                             'CFL',
#                             'cycleMode'
#                             'numberCycles',
#                             'numberSavedCycles'
                            'gravitationalField',
                            'gravityConstant',
                            'centralVenousPressure',
                            'minimumVenousPressure']

# initialisation controls
solverCalibrationElements =  [
#                              'solvingSchemeField',
#                              'solvingSchemeConnections',
                              'rigidAreas',
                              'simplifyEigenvalues',
                              'riemannInvariantUnitBase',
                              'automaticGridAdaptation'  ]

# initialisation controls
initialisationControlsElements = ['initialsationMethod',
                                  'initMeanFlow',
                                  'initMeanPressure',
                                  'estimateWindkesselCompliance',
                                  'compPercentageWK3',
                                  'compPercentageTree',
                                  'compTotalSys']

#----------------------------------------------------------------------------------------# 
# global fluid
globalFluidElements = ['my',
                       'rho',
                       'gamma']


vascularNetworkElements = ['simulationContext',
                           'solverCalibration',
                           'initialisationControls']

##########################################################################################
## Communicator objects

communicatorRealTimeViz  = ['comType', 'comId', 'vesselId', 'node', 'dn', 'quantitiesToPlot' ]
CommunicatorBaroreceptor = ['comType', 'comId', 'vesselId', 'node']

communicatorReference = {'CommunicatorRealTimeViz' : communicatorRealTimeViz,
                         'CommunicatorBaroreceptor': CommunicatorBaroreceptor
                         }

##########################################################################################
## class BoundaryConditions

#----------------------------------------------------------------------------------------# 
# Dictionary with boundary class references! Name_in_xml : class in classBoundaryConditions.py
bcTagsClassReferences = { # BoundaryConditions normal
                          'Flow-PhysiologicalFunction'       :'PhysiologicalFunction',
                          'Flow-Sinus'                       :'Sinus',
                          'Flow-Sinus2'                      :'Sinus2',
                          'Flow-RampMean'                    :'RampMean',
                          'Flow-Fourier'                     :'Fourier',
                          'Flow-PhysiologicalData'           :'PhysiologicalData',
                          'Flow-FromFile'                    :'FlowFromFile',
                          'Velocity-Gaussian'                :'expVelocity',
                          'Pressure-Sinus'                   :'Sinus',
                          'Pressure-Sinus2'                  :'Sinus2',
                          'Pressure-RampMean'                :'RampMean',
                          'ReflectionCoefficient'            :'ReflectionCoefficient',
                          'ReflectionCoefficientTimeVarying' :'ReflectionCoefficientTimeVarying',
                          'Resistance'                       :'Resistance',
                          'Windkessel-2Elements'             :'Windkessel2',
                          'Windkessel-3Elements'             :'Windkessel3',
                          'L-network'                        :'L_network',
                          'VaryingElastanceHeart'            :'VaryingElastance',
                          # BoundaryCondition names if 1 Vessel '_' == end-positionreference 
                          '_Flow-PhysiologicalFunction'       :'PhysiologicalFunction',
                          '_Flow-Sinus'                       :'Sinus',
                          '_Flow-Sinus2'                      :'Sinus2',
                          '_Flow-RampMean'                    :'RampMean',
                          '_Flow-Fourier'                     :'Fourier',
                          '_Flow-PhysiologicalData'           :'PhysiologicalData',
                          '_Flow-FromFile'                    :'FlowFromFile',
                          '_Velocity-Gaussian'                :'expVelocity',
                          '_Pressure-Sinus'                   :'Sinus',
                          '_Pressure-Sinus2'                  :'Sinus2',
                          '_Pressure-RampMean'                :'RampMean',
                          '_ReflectionCoefficient'            :'ReflectionCoefficient',
                          '_ReflectionCoefficientTimeVarying' :'ReflectionCoefficientTimeVarying',
                          '_Resistance'                       :'Resistance',
                          '_Windkessel-2Elements'             :'Windkessel2',
                          '_Windkessel-3Elements'             :'Windkessel3',
                          '_L-network'                        :'L_network',
                          '_VaryingElastanceHeart'            :'VaryingElastance'}

#----------------------------------------------------------------------------------------# 
# Tag dictionary defines the variables a boundaryCondition needs as input via xml!
boundaryConditionElements = {
          # Type 1 attribute
          # BoundaryConditions normal
          'ReflectionCoefficient'       :['Rt'],
          'ReflectionCoefficientTimeVarying' : ['RtOpen','Topen1','Topen2','RtClosed','Tclosed1','Tclosed2'],
          'Resistance'                  :['Rc'],
          'Windkessel-2Elements'        :['Rc','C'],
          'Windkessel-3Elements'        :['Rc','Rtotal','C','Z'],
          'Flow-PhysiologicalFunction'  :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-Sinus'                  :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-Sinus2'                 :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-RampMean'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Velocity-Gaussian'           :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','gaussC','prescribe'],
          'Pressure-Sinus'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Pressure-Sinus2'             :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Pressure-RampMean'           :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-Fourier'                :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','scale','prescribe'],
          'Flow-PhysiologicalData'      :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          'Flow-FromFile'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','filePathName','prescribe'],
          'L-network'                   :['Z','C'],
          'VaryingElastanceHeart'       :['T', 'Emax', 'Emin', 'Tpeak', 'V0', 'K'],
          # BoundaryCondition names if 1 Vessel '_' == end-positionreference 
          '_Flow-PhysiologicalFunction' :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-Sinus'                 :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-Sinus2'                :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-RampMean'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-Fourier'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','scale','prescribe'],
          '_Flow-PhysiologicalData'     :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Flow-FromFile'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','filePathName','prescribe'],
          '_Velocity-Gaussian'          :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','C','prescribe'],
          '_Pressure-Sinus'             :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Pressure-Sinus2'            :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_Pressure-RampMean'          :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
          '_ReflectionCoefficient'      :['Rt'],
          '_ReflectionCoefficientTimeVarying' : ['RtOpen','Topen1','Topen2','RtClosed','Tclosed1','Tclosed2'],
          '_Resistance'                 :['Rc'],
          '_Windkessel-2Elements'       :['Rc','C'],
          '_Windkessel-3Elements'       :['Rc','Rtotal','C','Z'],
          '_L-network'                  :['Z','C'],
          'None'                        :['']} 
##########################################################################################
## class Vessel

vesselAttributes = ['Id',
                    'name']

vesselTopologyElements = ['leftDaughter',
                          'rightDaughter',
                          'angleYMother']
   
vesselGeometryElements = ['geometryType',
                          'length',
                          'radiusProximal',
                          'radiusDistal',
                          'N']

vesselComplianceElements = {'Laplace'     :['complianceType','constantCompliance','externalPressure','Ps','As','betaLaplace'],
                            'Laplace2'    :['complianceType','constantCompliance','externalPressure','Ps','As','wallThickness','youngModulus'],
                            'Exponential' :['complianceType','constantCompliance','externalPressure','Ps','As','betaExponential'],
                            'Hayashi'     :['complianceType','constantCompliance','externalPressure','Ps','As','betaHayashi'],
                            'Reymond'     :['complianceType','constantCompliance','externalPressure','Ps','As','distensibility']
                            }

vesselFluidElements = ['applyGlobalFluid',
                       'my',
                       'rho',
                       'gamma']

vesselElementReference = {'topology'     :vesselTopologyElements,
                          'geometry'     :vesselGeometryElements,
                          'compliance'   :vesselComplianceElements,
                          'fluid'        :vesselFluidElements}

vesselElements = ['topology',
                  'geometry',
                  'compliance',
                  'fluid']


##########################################################################################
## XML file elements

xmlElements  = [ 'simulationContext',
                 'solverCalibration',
                 'initialisationControls',
                 'globalFluid',
                 'communicators',
                 'boundaryConditions',
                 'vessels' ]

xmlElementsReference = {'simulationContext'     : simulationContextElements,
                        'solverCalibration'     : solverCalibrationElements,
                        'initialisationControls': initialisationControlsElements,
                        'boundaryConditions'    : boundaryConditionElements,
                        'globalFluid'           : globalFluidElements,
                        'communicators'         : communicatorReference,
                        'vessels'               : vesselElements }


# ######## old xml version 3.0 to be deleted soon
# #----------------------------------------------------------------------------------------# 
# # dictionary to set the output-du matrix for common flow/pressure bc-function as Sinus,Sinus^2 and Mean
# import numpy as np
# bcDuMatrix = { 'Pressure'   : np.array([1,0]),
#                'Flow'       : np.array([0,1])}
# 
# bcTags = {# Type 1 attribute
#           # BoundaryConditions normal
#           'ReflectionCoefficient'       :['Rt'],
#           'ReflectionCoefficientTimeVarying' : ['RtOpen','Topen1','Topen2','RtClosed','Tclosed1','Tclosed2'],
#           'Resistance'                  :['Rc'],
#           'Windkessel-2Elements'        :['Rc','C'],
#           'Windkessel-3Elements'        :['Rc','Rtotal','C','Z'],
#           'Flow-PhysiologicalFunction'  :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           'Flow-Sinus'                  :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           'Flow-Sinus2'                 :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           'Flow-RampMean'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           'Velocity-Gaussian'           :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','gaussC','prescribe'],
#           'Pressure-Sinus'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           'Pressure-Sinus2'             :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           'Pressure-RampMean'           :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           'Flow-Fourier'                :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','scale','prescribe'],
#           'Flow-PhysiologicalData'      :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           'Flow-FromFile'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','filePathName','prescribe'],
#           'L-network'                   :['Z','C'],
#           'VaryingElastanceHeart'       :['T', 'Emax', 'Emin', 'Tpeak', 'V0', 'K'],
#           # BoundaryCondition names if 1 Vessel '_' == end-positionreference 
#           '_Flow-PhysiologicalFunction' :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           '_Flow-Sinus'                 :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           '_Flow-Sinus2'                :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           '_Flow-RampMean'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           '_Flow-Fourier'               :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','scale','prescribe'],
#           '_Flow-PhysiologicalData'     :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           '_Flow-FromFile'              :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','filePathName','prescribe'],
#           '_Velocity-Gaussian'          :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','C','prescribe'],
#           '_Pressure-Sinus'             :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           '_Pressure-Sinus2'            :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           '_Pressure-RampMean'          :['amp','ampConst','Npulse','Tpulse','freq','Tspace','runtimeEvaluation','prescribe'],
#           '_ReflectionCoefficient'      :['Rt'],
#           '_ReflectionCoefficientTimeVarying' : ['RtOpen','Topen1','Topen2','RtClosed','Tclosed1','Tclosed2'],
#           '_Resistance'                 :['Rc'],
#           '_Windkessel-2Elements'       :['Rc','C'],
#           '_Windkessel-3Elements'       :['Rc','Rtotal','C','Z'],
#           '_L_network'                  :['Z','C'],
#           'None'                        :['']} 


##########################################################################################
## Vascular Polynomial Chaos

polynomialChaosDistributions = ['Uniform','Normal']

vPCconfigurationTemplate ={ 'createDistributions' : True,
                            'createOrthoPoly'  : True,
                            'createSample'     : True,
                            'runSimulations'   : True,
                            'calculateGPCE'    : True,
                            'preProcessData'   : True,
                            'plotMinMaxPoints' : True,
                            'postProcessing'   : True,
                            'plotMeanSTD'      : True,
                            'plotPeaks'        : True,
                            'polynomialOrders' : [3],
                            'sampleMethod'     : 'M',
                            'linearWaveSplit'  : True,
                            'velocityProfileCoefficient' : 2.0,
                            'polynomsToCalculate' : [[0,0], [0,2]],
                            'names' : ['A' ,'B'],
                            'delta'  :{'A': 
                                        {'Pressure' :1.0,'Flow':1.0e-8,
                                        'Pressure_f':1.0,'Flow_f':1.0e-8,
                                        'Pressure_b':1.0,'Flow_b':1.0e-8} ,
                                        'B':
                                        {'Pressure' :2.7,'Flow':0.5e-7,
                                        'Pressure_f':2.7,'Flow_f':0.5e-7,
                                        'Pressure_b':2.7,'Flow_b':0.5e-7} },
                            'peaksToEvaluate' : {'A': 
                                                {'extremaPressure':[0],'extremaFlow':[0],
                                                'extremaPressure_f':[0],'extremaFlow_f':[0],
                                                'extremaPressure_b':[0],'extremaFlow_b':[0]} ,
                                                'B':
                                                {'extremaPressure':[0],'extremaFlow':[0],
                                                'extremaPressure_f':[0],'extremaFlow_f':[0],
                                                'extremaPressure_b':[0],'extremaFlow_b':[0]} },
                            'plotDirectory' : 'PolynomialChaos',
                            'deterministicDataSetNumbers' : [],
                            'polynomsToPlotOrder' : [0],
                            'plotMeanConfidenceInterval' : False,
                            'plotMeanConfidenceAlpha' : 1.0,
                            'plotMeanSigmaInterval' : False,
                            'peakAnalysis' : True,
                            'plotPeaksConfidenceAlpha' : 1.0,
                            'plotPeaksAnalyticSensitivity' : False,
                            'plotPeaksMeanSTDBoxPlotsSingle' : True                            }

