
# coding: utf-8

# # [Starfish](http://www.ntnu.no/starfish)
# A python application for the simulation and vizualization of 1d arterial networks. Capabilities include coupling with a wide variety of boundary conditions, usage of a number of material models for the arterial compliance, and stochastic simulation accounting for uncertainty in any model parameters.
# 
# ## Demo
# 1. Example simulation of a 55 vessel arterial network
# 2. Analysis of calibration of the same network to noninvasive clinical imaging
# 
# ## Outline of this tutorial
# 
# 1. We will construct a simple 1 segment arterial network with a windkessel boundary condition
# 2. We will simulate this network and look at the pressure, flow and area of the vessel
# 3. We will vary geometrical and material parameters and plot the effects
# 4. We will uncertain parameters and evaluate the uncertainty and sensitivity of pressure, flow and area
# 
# ## Basic results for waves in compliant tubes
# 
# Fluid flow in compliant tubes exhibits wave properties, thus it is of interest to note the fundamental relationship between the material parameters and the resulting wave speed
# 
# $$
# c_{pw}^2 = \frac{A}{\rho C} 
# $$
# where $A$ is the cross sectional area and $C=\frac{\partial A}{\partial P}$ is called the compliance.
# 
# We will employ a Hookean type material model linearized about a reference pressure $P_s$
# $$
# A = (\frac{P-P_s}{\beta} + \sqrt{A_s})^2
# $$
# and
# $$
# \beta = \frac{4}{3} \sqrt{\pi} E h
# $$
# where $E$ is the elastic (Young's) modulus and $h$ is the thickness of the arterial wall.
# From the function $A(P)$ we see that
# $$
# C = \frac{2}{\beta} (\frac{P-P_s}{\beta} + \sqrt{A_s})
# $$

# In[1]:

import matplotlib.pyplot as plt


# In[3]:

# In the real world
# 1. run singleVesselUQSA.xml with standard simulator?
# 2. run VascularPolynomialChaos.py  on the singleVesselUQSA.xml to produce vpcCase
# 3. run VascularPolynomialChaos.py on the resulting vpc.xml to evaluate all
import sys, os
import shutil
import numpy as np
import math
import starfish.UtilityLib.moduleFilePathHandler as mFPH
import starfish.UtilityLib.moduleXML as mXML
import starfish.NetworkLib.classVascularNetwork as cVN
import starfish.NetworkLib.classBoundaryConditions as cBC
import starfish.SolverLib.class1DflowSolver as c1dFS


# In[4]:

unit_mmhg_pa = 133.3
unit_pa_mmhg = 1./unit_mmhg_pa
unit_cm2_m2 = 1. / 100. / 100.
unit_cm_m= 1. / 100.
unit_m2_cm2 = 1. / unit_cm2_m2


# In[5]:

vascularNetwork = cVN.VascularNetwork()


# In[6]:

vesselId = 1
vascularNetwork.addVessel(vesselId)
vessel = vascularNetwork.vessels[vesselId]


# In[7]:

vessel.geometryType = 'uniform'
vessel.length = 0.1
vessel.radiusDistal = (2.5/2*1e-2) # This is not used by uniform
vessel.radiusProximal =(2.5/2*1e-2)


# In[8]:

vessel.complianceType = "Hayashi"
#vessel.constantCompliance = False
#vessel.externalPressure = 0.0
vessel.betaHayashi = 3.67 # Approximate value for a human aorta
vessel.As = None # Use Reference area 
vessel.Ps = 100*133.32# Reference Pressure at area As
vessel.N = 10


# In[9]:

inlet = cBC.Sinus()
inlet.name = "Flow-Sinus2" # Necessary to set as an inflow boundary condition
inlet.amp = 50*1e-6
inlet.ampConst = 50*1e-6
inlet.freq = 1
inlet.Npulse =3


# In[10]:

outlet = cBC.ReflectionCoefficient()
outlet.name = "_ReflectionCoefficient"
outlet.Rt  = 0.2


# In[11]:

vascularNetwork.boundaryConditions[vessel.Id]=[inlet, outlet]


# In[12]:

# Save network
networkName = "tutorial"
vascularNetwork.name = networkName
vascularNetwork.description=''
mXML.writeNetworkToXML(vascularNetwork)


# In[18]:

dataNumber = "000"


# In[13]:


import starfish.VascularPolynomialChaosLib.classDistributionManager as cDistMng
import starfish.VascularPolynomialChaosLib.moduleFilePathHandlerVPC as mFPH_VPC
import starfish.VascularPolynomialChaosLib.moduleBatchSimulationManager as mBSM
import starfish.VascularPolynomialChaosLib.classUqsaCase as cUqsaCase
import starfish.VascularPolynomialChaosLib.classSampleManager as cSM
import starfish.VascularPolynomialChaosLib.classRandomInput as cRI
import starfish.VascularPolynomialChaosLib.classRandomInputManager as cRIM

import starfish.VascularPolynomialChaosLib.classUqsaMethods as cUQM
import starfish.VascularPolynomialChaosLib.classLocationOfInterestManager as cLOIM


# In[14]:

polychaos_method = cUQM.UqsaMethodPolynomialChaos()
polychaos_method.polynomialOrder = 3
polychaos_method.sampleFactor = 2


# In[15]:

loi_manager = cLOIM.LocationOfInterestManager()
locationId = "vessel_1"
locationName = "root"
quantitiesOfInterestToProcess = ["Pressure", "Flow"]
xVal = 8*unit_cm_m
confidenceAlpha = 5
loi_manager.addLocationOfInterest(locationId, locationName,
                                  quantitiesOfInterestToProcess, xVal, confidenceAlpha)


# In[16]:

rim = vascularNetwork.randomInputManager
rim = vascularNetwork.randomInputManager = cRIM.RandomInputManager()


# In[26]:

uqsaCaseFile = mFPH_VPC.getFilePath('uqsaCaseXmlFile', networkName, dataNumber, 'write') # Creates path
uqsaCase = cUqsaCase.UqsaCase() 
if False and os.path.exists(uqsaCaseFile):
    uqsaCase.loadXMLFile(uqsaCaseFile)
    print(uqsaCase.uqsaMethods)
else:
    pass
# 1. specify random variables
rim = vascularNetwork.randomInputManager = cRIM.RandomInputManager()
a = vessel.radiusProximal
pct_unc = 0.05
b = pct_unc*vessel.radiusProximal #TODO is this symmetric for uniform?
rvType = "Uniform"
parameter = "vessel_1_radiusProximal"
rim.addRandomInput("Z1", a, b, rvType, parameter)

a = vessel.betaHayashi
pct_unc = 0.05
b = pct_unc*vessel.betaHayashi #TODO is this symmetric for uniform?
rvType = "Uniform"
parameter = "vessel_1_betaHayashi"
rim.addRandomInput("Z2", a, b, rvType, parameter)

sample_manager = uqsaCase.sampleManager = cSM.SampleManager()
sample_manager.dependentCase = False
sample_manager.samplingMethod = "S"

# 2. specify quantities of interest
uqsaCase.locationOfInterestManager = loi_manager = cLOIM.LocationOfInterestManager()
locationId = "vessel_1"
locationName = "root"
quantitiesOfInterestToProcess = ["Pressure", "Flow"]
xVal = 8*unit_cm_m
confidenceAlpha = 5
loi_manager.addLocationOfInterest(locationId, locationName,
                                  quantitiesOfInterestToProcess, xVal, confidenceAlpha)

# 3. specify uqsa method
#     3.a gpc order, number of samples, sampling method

polychaos_method = cUQM.UqsaMethodPolynomialChaos()
polychaos_method.polynomialOrder = 3
polychaos_method.sampleFactor = 2
uqsaCase.uqsaMethods ={"PC-3":polychaos_method}

# 4.
uqsaCase.initialize(networkName,dataNumber)



import starfish.VascularPolynomialChaos as vpc

uqsaCase.multiprocessing = True
uqsaCase.simulateEvaluations = True
uqsaCase.locationOfInterestManager.evaluateSimulationTime = True
vpc.run_uqsa_case(uqsaCase, vascularNetwork=vascularNetwork)

