import numpy as np 
import csv
from math import exp
import ODESolver as OD

import os,sys
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )

class BoundaryCondition(object):
	'''
	Base-class for all boundary conditions	
	'''
	position = None
	name = None
	
	##
	# if False prescribe influx of Flow and Pressure 
	# if True  prescribe total Values of Pessure and Flow
	prescribeTotalValues = False
	##
	# conditionQuantity <string> should be either 'Flow' or 'Pressure'
	conditionQuantity = None  
	
	
	def update(self,bcDict):
		'''
		updates the updateBoundaryDict data using a dictionary in form of 
		bcDict = {'variableName': value}
		'''
		for key,value in bcDict.iteritems():
			try:
				self.__getattribute__(key)
				self.__setattr__(key,value)
			except: 
				print 'Error boundaryConditions.update(): wrong key: %s, could not set up boundaryCondition' %key
				
	def setPosition(self, position):
		'''Set the position of the boundaryCondition '''
		self.position = position
	
	def getVariableValue(self,variableName):
		'''
		Returns value of variable with name : variableName
		States Error if not such variable
		'''
		try:
			return self.__getattribute__(variableName)
		except: 
			print "ERROR BoundaryCondition.getVariable() : BoundaryCondition has no variable {}".format(variableName)
		
		
class BoundaryConditionType2(BoundaryCondition):
	def setPosition(self, position):
		'''
		Set the position of the boundaryCondition
		and determines the return function
		based on the position of the bC (0,-1)
		'''
		self.position = position
		if self.position == 0:
			self.returnFunction = self.funcPos0
		elif self.position ==-1:
			self.returnFunction = self.funcPos1
		else:
			self.returnFunction = None
		
########################################################################################
# Type 1 boundary Conditions
########################################################################################

class BoundaryConditionType1(BoundaryCondition):
	"""	
	Boundary profile - type 1
	
		creates frame for a periodic signal 
		
		signal type can be quantified with function1 / function2 
		which are a function of (t)
		
		call function input (n,dt)
		gives amplitude value for pressure/flow back
		for the function delta g
	"""
	def __init__(self):
		
		BoundaryCondition.__init__(self)
		
		self.type = 1
		self.prescribe = 'influx'
		
		## variables from xml
		self.amp	  = 0
		self.ampConst = 0 
		self.Npulse   = 1.
		self.freq	  = 1.
		self.Tpulse   = 0.
		self.Tspace   = 0.
		self.runtimeEvaluation = False
				
		## evaluated values
		self.lastU       = 0.0
		self.TmeanFlow   = 0.0
		self.MeanFlow    = 0.0
		self.Tperiod     = 0.0
		self.pulseTime   = []
		self.initialStep = False
		self.duMatrix    =  np.ones(2)
		self.duVector    = np.empty(2)
		
		# values to manipulate period during runtime
		self.updateTime = -10.0
		self.freqNew    = self.freq
		self.TspaceNew  = self.Tspace
	
	def initialize(self,bcDict):
		'''
		updates - the updateBoundaryDict data using a dictionary in from of 
				  bcDict = {'variableName': value}
				- time period of one pulse
				- pulseTimes
		'''
		
		self.update(bcDict)
		
		self.Tperiod = self.Tspace+1.0/self.freq
				
		self.pulseTime = []
		for puls in np.arange(self.Npulse):
			#self.pulseTime.append([self.Tpulse + (self.Tperiod)*puls,
			#						   self.Tpulse + 1.0/self.freq+ (self.Tperiod)*puls - self.TmeanFlow,
			#						   self.Tpulse + self.Tperiod +(self.Tperiod)*puls  + self.TmeanFlow,
			#						   self.Tpulse + self.Tperiod +(self.Tperiod)*puls  + self.TmeanFlow])
			
			
			if puls == 0:
				self.pulseTime.append([self.Tpulse + (self.Tperiod)*puls,
									   self.Tpulse + 1.0/self.freq+ (self.Tperiod)*puls - self.TmeanFlow,
									   self.Tpulse + self.Tperiod +(self.Tperiod)*puls ,
									   self.Tpulse + self.Tperiod +(self.Tperiod)*puls ])
				
			else:
				self.pulseTime.append([self.Tpulse + (self.Tperiod)*puls,
									   self.Tpulse + 1.0/self.freq+ (self.Tperiod)*puls,
									   self.Tpulse + self.Tperiod +(self.Tperiod)*puls,
									   self.Tpulse + self.Tperiod +(self.Tperiod)*puls ])
		
		self.initialStep = False
		
		self.lastU = self.calculateOneStep(0,0)
		
		if self.prescribe == 'total':
			self.prescribeTotalValues = True
		elif self.prescribe == 'influx':
			self.prescribeTotalValues = False
			
		if 'Flow' in self.name:
			self.duMatrix = np.array([0,1])
			self.conditionQuantity = 'Flow'
		elif 'Pressure' in self.name:
			self.duMatrix = np.array([1,0])
			self.conditionQuantity = 'Pressure'
		
		## load file for Flow-From File
		if self.name == 'Flow-FromFile':
			if self.loadedFile == False:
				self.loadFile()
			
	def updatePeriodRuntime(self, TperiodNew, updateTime):
		'''
		This function updates the freq, Tspace and pulsTime for  new given Tperiod
		'''
				
		freq    = self.freq
		Tspace  = self.Tspace
		Tperiod = self.Tperiod
		## calculate new values
		self.freqNew   = Tperiod*freq/TperiodNew
		self.TspaceNew = TperiodNew*Tspace/Tperiod
		
		manipulateOtherPulses = False
		findFirstPulse = True
		
		for pulse in self.pulseTime:
			
			if  updateTime >= pulse[0] and updateTime <= pulse[2] and findFirstPulse == True:
				pStart = pulse[3]
				manipulateOtherPulses = True
				
				self.updateTime = pulse[2]				
								
			elif updateTime >= pulse[2] and updateTime <= pulse[3]:
				updateTime = pulse[3] 
					
			if manipulateOtherPulses == True:
				if findFirstPulse == False:
					pulse[0] = pStart
					pulse[1] = pStart + 1.0/self.freqNew  #- self.TmeanFlow
					pulse[2] = pStart + TperiodNew        #- self.TmeanFlow
					pulse[3] = pStart + TperiodNew 				
					pStart = pulse[3]
				else: findFirstPulse = False
		
	def calculateDuVector(self,Tsteps,dt):
		'''
		Function calculates the duVector for a given number of time steps Tsteps and dt
		return: self.duVector
		'''
		self.duVector = np.zeros((Tsteps,2))
		
		lastStep = self.calculateOneStep(0,dt)
		for nt in range(len(self.duVector)):
			nextStep = self.calculateOneStep(nt+1,dt)
			self.duVector[nt] = nextStep - lastStep
			lastStep = nextStep
			
		return self.duVector
		
	def calculateDu(self,n,dt):
		'''
		Function calculates the du for a given time step and dt
		return: du
		'''
		nextStep = self.calculateOneStep(n+1,dt)	
		du = nextStep - self.lastU
		self.lastU = nextStep
		
		
		if n*dt > self.updateTime-dt and n*dt < self.updateTime + dt :
			
			self.freq    = self.freqNew
			self.Tspace  = self.TspaceNew		
			self.Tperiod = self.TspaceNew + 1.0/self.freqNew
			
			self.updateTime = -10
			#print "bcType1 updatedTime", n*dt, self.updateTime + dt, self.updateTime - dt, self.freq , self.Tperiod
			
		return du
		
	def calculateOneStep(self,n,dt):
		'''
		calculates the amplitude value for one time step n and dt
		
		return: array([ampP,ampQ])
		'''
		t = n*dt
		ampT = self.ampConst
				
		timeInPulseTime = False
		# get time slots
		for pTA in self.pulseTime:
			if pTA[0] <= t and t <= pTA[3]:
				pulsNum = self.pulseTime.index(pTA)
				pTA1 = pTA[0]
				pTA2 = pTA[1]
				pTA3 = pTA[2]
				pTA4 = pTA[3]
				timeInPulseTime = True
				break
		if timeInPulseTime == False:
			#print "WARNING BoundaryCondtitionType1.calculateOneStep(): t = {} not in pulseTime, return constant amp!".format(t)
			return ampT
		
		if pulsNum == 0: tmeanFlowShift = self.TmeanFlow
		else : tmeanFlowShift =  0
				
		try:
			if pTA1 <= t and t <= pTA2:# and self.initialStep == True:
				t = t + tmeanFlowShift#self.TmeanFlow
				t0 = pTA1 #self.Tpulse+pulsNum*self.Tperiod
								
				## function1 according to the signal type
				ampT = ampT + self.function1(t,t0,pulsNum)
				
			elif pTA2 <= t and t < pTA3:
				t = t+ tmeanFlowShift#self.TmeanFlow
				t0 = pTA1 #self.Tpulse+pulsNum*self.Tperiod
				## function2 according to the signal type in most cases return self.ampConst
				ampT = ampT + self.function2(t,t0,pulsNum)
			
			elif pTA3 <= t and t <= pTA4: 
				t = t-self.Tperiod++tmeanFlowShift #+self.TmeanFlow#
				t0 = pTA1 #self.Tpulse+pulsNum*self.Tperiod
				## function3 according to the signal type in most cases same as function 1
				ampT = ampT + self.function3(t,t0,pulsNum)
			
			t = n*dt
			if self.Tpulse <= t and t < self.Tpulse+1.5*dt and self.initialStep == False:
				self.initialStep = True
				
			return ampT*self.duMatrix
		
		except:
			return ampT*self.duMatrix	
		
	def findMeanFlowAndMeanTime(self,givenMeanFlow = None, quiet = False):
		'''
		This function calculates the mean flow of the signal self.MeanFlow
		and the first occurence time of the mean flow self.TmeanFlow
		'''
		#find meanFlow
		self.initialize({})
		
		period = self.Tperiod
		totalTime = self.Tperiod+self.Tpulse
		dt = 0.001
		nTsteps = int(np.round(totalTime/dt, 0))
		nTstepsStart = int(np.round(self.Tpulse/dt, 0))
		integral = 0.0
		bcFunction = np.zeros(nTsteps-nTstepsStart-1)
		time = np.zeros(nTsteps-nTstepsStart-1)
		
		for n in xrange(nTstepsStart,nTsteps-1):
			flow = self.calculateOneStep(n,dt)
			integral = integral + (flow+self.calculateOneStep(n+1,dt))/2.0
			bcFunction[n-nTstepsStart] = flow[1]
			time[n-nTstepsStart]= (n-nTstepsStart)*dt
			
		if givenMeanFlow == None:
			self.MeanPressure = (integral*dt/period)[0]
			self.MeanFlow = (integral*dt/period)[1]
		else:
			self.MeanFlow = givenMeanFlow
			difference = abs(givenMeanFlow - (integral*dt/period)[1])
			if quiet == False:
				print '''\n  WARNING:  given meanFlow given {} differs from meanFlow of boundaryCondition {}
	            (evaluated with integral over one period). 
	            The difference is {} ml s-1 \n'''.format(givenMeanFlow*1.e6,(integral*dt/period)[1]*1.e6,difference*1.e6) 
						
		## finde meanFlow time
		from scipy import interpolate as interpolate
		numberIntPoints = 1000000
		inflowFunction = interpolate.interp1d(time,bcFunction)
		newTime = np.linspace(0,time[-1], numberIntPoints)
		newValues = inflowFunction(newTime)
		self.TmeanFlow = 0
		for n in np.linspace(0,numberIntPoints-1,numberIntPoints ): 
			ti = newTime[n]
			qi = newValues[n]
			if abs(self.MeanFlow-qi) <= 1.e-5:
				notFound = True
				while notFound == True:
					self.TmeanFlow = ti
					n = n+1
					if n == numberIntPoints: break
					if abs(self.MeanFlow-newValues[n+1]) > abs(self.MeanFlow-newValues[n]):
						self.TmeanFlow = newTime[n]
						break
				break
		
		if quiet == False:
			print '====================================='
			print '___BC _Type1: mean flow evaluation___'
			print 'meanFlow evaluated (ml/s)  {:.6}'.format(str((integral*dt/period)[1]).ljust(5))
			print 'meanFlowTime (s)           {:.6}'.format(str(self.TmeanFlow).ljust(5))
			print 'total volume/period (ml)   {:.6}'.format(str(integral[1]*dt*1e6).ljust(5))
		
		if self.TmeanFlow != 0.0: self.Tpulse = 0.0
		self.initialize({})
		
	def function1(self,t,t0,pulsNum):
		'''
		amplitude function caracterized by signal type
		'''
		pass
	def function2(self,t,t0,pulsNum):
		'''
		amplitude function caracterized by signal type
		'''
		return 0
	
	def function3(self,t,t0,pulsNum):
		'''
		amplitude function caracterized by signal type
		'''
		return self.function1(t,t0,pulsNum)

class RampMean(BoundaryConditionType1):
	"""	
	Boundary profile - type 1
	
		ramps to a mean amplitude self.amp starting from self.ampConst
	"""
		
	def function1(self,t,t0,pulsNum):
		return self.amp*pow(np.sin(np.pi*(t - self.Tpulse)/(1./self.freq*2.0)),2.0)
	
	def function2(self,t,t0,pulsNum):
		return self.amp
	
	def function3(self,t,t0,pulsNum):
		return self.amp
	
class Sinus(BoundaryConditionType1):
	"""	
	Boundary profile - type 1
	
		creates in a periodic sinus signal
	"""
	def function1(self,t,t0,pulsNum):
		return self.amp*np.sin(2*np.pi*self.freq*(t - t0))
	

class Sinus2(BoundaryConditionType1):
	"""	
	Boundary profile - type 1
	
		creates in a periodic sinus-squared signal
	"""
	def function1(self,t,t0,pulsNum):
		return self.amp*pow(np.sin(np.pi*self.freq*(t - t0 ) ),2.0)

class PhysiologicalFunction(BoundaryConditionType1):
	"""	
	Boundary profile - type 1
	
		creates a similar heart-outflow signal as found in
		Stergiopulos et al. 1992
		
		set together from 4 continous functions sin2,sin2,cos,sin2
		to lead to a continous function
	"""
	def __init__(self):
		
		BoundaryConditionType1.__init__(self)
				
		## additional variables for this function
		self.lowPoint = 0.1739
		self.fracSin2 = 0.36
		self.fracCos = 0.43
		self.fracRes = 1.0 - (self.fracSin2+self.fracCos)
			
	def function1(self,t,t0,pulsNum):
		if t < self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq)+self.fracSin2*1/self.freq:
			ampT = self.amp*pow(np.sin(np.pi*self.freq/(2*self.fracSin2)*(t - (self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq)))),2.0)
		elif t < self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq)+(self.fracSin2+self.fracCos)*1/self.freq:
			ampT = self.amp*(np.cos((np.pi*self.freq/(2.0*self.fracCos)*(t-self.fracSin2/self.freq- (self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq))))) )
		elif t < self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq)+(self.fracSin2+self.fracCos+self.fracRes*3.0/8.0)*1/self.freq:
			ampT =self.amp*(pow(np.sin((t-(self.fracSin2+self.fracCos)/self.freq-(self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq))-(1.0/(self.freq)*self.fracRes)*3./8.)*np.pi*self.freq/(1.5*self.fracRes)),2.0)*2*self.lowPoint-self.lowPoint)
		else:
			ampT =self.amp*(pow(np.sin((t-(self.fracSin2+self.fracCos)/self.freq-(self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq))-(1.0/(self.freq)*self.fracRes)*3./8.0)*np.pi*self.freq/(self.fracRes*1.25)),2.0)*self.lowPoint-self.lowPoint)		
		return ampT

class expVelocity(BoundaryConditionType1):
	"""	
	Boundary profile - type 1
	
		creates a single gaussian peak signal as found
		in D.Xiu and S.P.Sherwin 2007
	"""
	def __init__(self):
		
		BoundaryConditionType1.__init__(self)
				
		## additional variables for this function
		self.gaussC = 5000
		self.area   = 1
		
		## set to duMatrix
		self.duMatrix[1] = 1
		
		
	def function1(self,t,t0,pulsNum):
		return self.amp * exp(-self.gaussC*(t-t0+0.5/self.freq)**2)* self.area

class Fourier(BoundaryConditionType1):
	"""	
	Boundary profile - type 1
	
		creates fourier signal similar to the flow of the heart
	"""
	def __init__(self):
		
		BoundaryConditionType1.__init__(self)
		
		## additional variables for this function
		self.scale  = 1
		self.harm = [[0.2465,0.1975,0.2624,0.2132,0.0424,0.0722,0.0411,0.0093,0.0141,0.0044],
					 [0.0,2.2861,4.5723,6.8584,9.1445,11.4307,13.7168,16.0029,18.2891,20.5752],
					 [0.0,2.5010,-2.9986,-0.2689,1.4904,-2.9854,-0.0317,1.5292,-2.5394,0.5771]]
		
		## set variables
		self.freq =  self.harm[1][1]
				
	def function1(self,t,t0,pulsNum):
		amp = 0
		for i in range(len(self.harm[0])):
			amp += self.harm[0][i]*np.sin(2.0*np.pi*self.harm[1][i]*(t - t0) + self.harm[2][i] + np.pi/2.0)
		#amp = amp*self.scale
		return amp
		

class PhysiologicalData(BoundaryConditionType1):
	"""	
	Boundary profile - type 1
	
		creates signal based on measured physiological values
		(source unknown: values in NetworkLib/physiologicalData.py)
	"""
	def __init__(self):
		
		BoundaryConditionType1.__init__(self)
		
		## additional variables for this function
		from physiologicalData import aorticFlowPressure
		self.data = aorticFlowPressure()
		
		tDataMax = max(self.data.t)
		tDataMin = min(self.data.t)
		self.timeData = np.linspace(tDataMin,tDataMax,2000)
		
		
	def function1(self,t,t0,pulsNum):
		self.timeSim = np.linspace(0.0,1./self.freq,2000)
		tReset = t - t0
		tInter = np.interp(tReset,self.timeSim,self.timeData) 
		Q = np.interp(tInter,self.data.t,self.data.Q)*1.e6
		#P = np.interp(t-pulseNumber*self.tmax,self.data.t,self.data.P)
		return Q*self.duMatrix
			
class FlowFromFile(BoundaryConditionType1):
	"""	
	Boundary profile - type 1
	
		creates signal based on data values stored in *.csv file
		saved in network directory 
		
		csv delimiter ;
		first line of the colums must be t and Q respectively: 
		t;Q
		..; .. 
		
		expected units: t [s]; Q[m3 / s]
		
	"""
	def __init__(self):
		
		BoundaryConditionType1.__init__(self)
		
		## additional variables fill in with data in the readXML file
		self.filePathName = ''
		self.dataTime = []
		self.dataFlow = []
		self.loadedFile = False
		
	def loadFile(self):
		try:
			# set the path relative to THIS file not the executing file!
			if '.csv' not in self.filePathName: self.filePathName = self.filePathName.join(['','.csv'])
			pathAndFilename = ''.join([cur,'/../NetworkFiles/',self.filePathName])
			reader = csv.DictReader(open(pathAndFilename,'rb'),delimiter=';')
		except:
			print 'ERROR: boundaryConditions.FlowFromFile could not open file <<{}>> with boundary values, system exit'.format(self.filePathName.split('/')[-1])
			exit()
		try:
			dataTime = []
			dataFlow = []
			for row in reader:
				dataTime.append(float(row['t']))
				dataFlow.append(float(row['Q']))
			self.dataTime = np.asarray(dataTime)
			self.dataFlow = np.asarray(dataFlow) 
			self.loadedFile = True
		except: pass
		
		
	def function1(self,t,t0,pulsNum):
		# set up simulation time (needed here if freq changes after init known to slow down)
		self.timeData = np.linspace(min(self.dataTime),max(self.dataTime),2000)
		self.timeSim = np.linspace(0.0,1./self.freq,2000)
		# reset time to fit in array [0 , 1/freq]
		tReset = t - t0
		# interpolate time to dataTime
		tInter = np.interp(tReset,self.timeSim,self.timeData) 
		# interpolate the Q value from data and dataTime
		Q = np.interp(tInter,self.dataTime,self.dataFlow)
		#P = np.interp(t-pulseNumber*self.tmax,self.data.t,self.data.P)
		return Q*self.duMatrix
		
		
	
		
########################################################################################
# Type 2 boundary Conditions
########################################################################################
class PrescribedInflux(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	calculates omega if influx of pressure or flow is prescribed with a type1 condition
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	
	def __init__(self):
		self.type = 2
		self.name = "standardType2"

		self.returnFunction = None
		self.omegaNew = np.empty((2))
		self.dQInOut = np.empty((2))
		
	def __call__(self,_domegaField_,duPrescribed,R,L,n,dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_,duPrescribed,L,R)

	def funcPos0(self,_domegaField,duPrescribed,L,R):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		domegaPrescribed_ = np.dot(L[0],duPrescribed)

		self.omegaNew[0] = domegaPrescribed_
		self.omegaNew[1] =  _domegaField
		
		self.dQInOut = R[:][1]*self.omegaNew 	
						
		return np.dot(R,self.omegaNew),self.dQInOut	
		
	
	def funcPos1(self,domegaField_,duPrescribed,L,R):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		_domegaPrescribed = np.dot(L[1],duPrescribed)
		
		self.omegaNew[0] = domegaField_
		self.omegaNew[1] =  _domegaPrescribed
		
		self.dQInOut = (R[:][1]*self.omegaNew)[::-1].copy()	
		
		return np.dot(R,self.omegaNew) ,self.dQInOut
		
	
class PrescribedTotalFlow(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	calculates omega if total flow is prescribed with a type1 condition
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	
	def __init__(self):
		self.type = 2
		self.name = "standardType2"

		self.R = np.empty((2,2))
		self.returnFunction = None
		self.omegaNew = np.empty((2))
		self.duNew = np.empty((2)) 
		self.dQInOut = np.zeros((2))
	
	def __call__(self,_domegaField_,duPrescribed,R,L,n,dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_,duPrescribed,L)

	
	def funcPos0(self,_domegaField,duPrescribed,L):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		'''
			# 1. of L set values l11 = 0 l12 = 1.
			# 2. and make inverse
			# det = l11*l22 - l12*l21 = -l21 
			self.R[0][0] =  L[1][1] / -L[1][0] 
			self.R[0][1] =  1./L[1][0] 	#-1 / -L[1][0] 
			self.R[1][0] =  1.			#-L[1][0] / -L[1][0]
			self.R[1][1] =  0           	
			# 3. set omega vector
			self.omegaNew[0] = duPrescribed[1]
			self.omegaNew[1] = _domegaField
			# 4. du = self.R * self.omegaNew
			du = np.dot(self.R,self.omegaNew)
		How it is done:
			inserted R in dot product
		'''
		self.duNew[0] = - L[1][1] / L[1][0] * duPrescribed[1] + 1./L[1][0] * _domegaField
		self.duNew[1] = duPrescribed[1]
		self.dQInOut[0] = duPrescribed[1]
		return self.duNew, self.dQInOut
	
	def funcPos1(self,domegaField_,duPrescribed,L):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		'''
		How to do it:
			# 1. of L set values l21 = 0 l22 = 1.
			# 2. and make inverse
			# det = l11*l22 - l12*l21 = l11 
			self.R[0][0] =  1. / L[0][0] # 1/l11
			self.R[0][1] =  -L[0][1]/L[0][0] # -l12/l11
			self.R[1][0] =  0			
			self.R[1][1] =  1.      
			# 3. set omega vector
			self.omegaNew[0] = domegaField_
			self.omegaNew[1] = duPrescribed[1] 
			# 4. du = self.R * self.omegaNew
			du = np.dot(self.R,self.omegaNew)
		How it is done:
			inserted R in dot product
		'''
		self.duNew[0] = domegaField_ / L[0][0]  - L[0][1] / L[0][0] * duPrescribed[1]
		self.duNew[1] = duPrescribed[1]
		
		self.dQInOut[1] = duPrescribed[1]		
		return self.duNew, self.dQInOut

class PrescribedTotalPressure(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	calculates omega if total pressure is prescribed with a type1 condition
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	
	def __init__(self):
		self.type = 2
		self.name = "standardType2"

		self.returnFunction = None
		self.omegaNew = np.empty((2))
		self.R = np.empty((2,2))
		self.duNew = np.empty((2)) 
		self.dQInOut = np.zeros(2)
		
	def __call__(self,_domegaField_,duPrescribed,R,L,n,dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_,duPrescribed,L)

	
	def funcPos0(self,_domegaField,duPrescribed,L):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		'''
		How to do it:
			# 1. of L set values l11 = 0 l12 = 1.
			# 2. and make inverse
			# det = l11*l22 - l12*l21 = l22 
			self.R[0][0] =  1.
			self.R[0][1] =  0 	
			self.R[1][0] =  -L[1][0]/L[1][1]			
			self.R[1][1] =  1./L[1][1]	           
			# 3. set omega vector
			self.omegaNew[0] = duPrescribed[0]
			self.omegaNew[1] = _domegaField
			# 4. du = self.R * self.omegaNew
			du = np.dot(self.R,self.omegaNew)
		How it is done:
			inserted R in dot product
		'''
		self.duNew[0] = duPrescribed[0]
		self.duNew[1] = -L[1][0]/L[1][1] * duPrescribed[0] + _domegaField / L[1][1] 
		self.dQInOut[0] = self.duNew[1]
		return self.duNew, self.dQInOut
	
	def funcPos1(self,domegaField_,duPrescribed,L):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		'''
		How to do it:
			# 1. of L set values l11 = 0 l12 = 1.
			# 2. and make inverse
			# det = l11*l22 - l12*l21 = -l12 
			self.R[0][0] =  0. 
			self.R[0][1] =  1.	
			self.R[1][0] =  1./L[0][1]
			self.R[1][1] =  L[0][0] / -L[0][1]           
					
			# 3. set omega vector
			self.omegaNew[0] = domegaField_
			self.omegaNew[1] = duPrescribed[0]
			# 4. du = self.R * self.omegaNew
			du = np.dot(self.R,self.omegaNew)
		How it is done:
			inserted R in dot product
		'''
		self.duNew[0] = duPrescribed[0]
		self.duNew[1] = domegaField_ / L[0][1] - L[0][0] / L[0][1]  * duPrescribed[0] 
		self.dQInOut[1] = self.duNew[1]
		return self.duNew, self.dQInOut

class ReflectionCoefficient(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	Terminal reflection 
	only in combination with prescribed influx condition (Type1) or alone
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	def __init__(self):
		self.type = 2
		self.Rt = 0
	
		self.returnFunction = None
		
		self.omegaNew = np.empty((2))
	
	def __call__(self, _domegaField_, duPrescribed, R, L, n, dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_,duPrescribed,L,R)
		
	def funcPos0(self,_domegaField,duPrescribed,L,R):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		domegaPrescribed_ = np.dot(L[0],duPrescribed)
		domegaReflected_ = _domegaField*self.Rt + domegaPrescribed_
		
		self.omegaNew[0] = domegaReflected_
		self.omegaNew[1] = _domegaField
		
		self.dQInOut = R[:][1]*self.omegaNew 
		
		return np.dot(R,self.omegaNew),self.dQInOut
		

	def funcPos1(self,domegaField_,duPrescribed,L,R):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		_domegaPrescribed = np.dot(L[1],duPrescribed)
		_domegaReflected = domegaField_*self.Rt + _domegaPrescribed
		
		self.omegaNew[0] = domegaField_
		self.omegaNew[1] = _domegaReflected
		
		self.dQInOut = (R[:][1]*self.omegaNew)[::-1].copy()
		
		return np.dot(R,self.omegaNew),self.dQInOut
	

class ReflectionCoefficientTimeVarying(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	Terminal reflection 
	only in combination with prescribed influx condition (Type1) or alone
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	def __init__(self):
		self.type = 2
		
		## reflection Values
		self.RtOpen   = 0.2
		self.RtClosed = 0.8
		
		## Timing
		self.Topen1 =  0.0 #verschiebung nach rechts von 0.5/freq
		self.Topen2 =  0.0 #verschiebung nach links von 0.5/freq
		
		self.Tclosed1 = -0.0 #verschiebung nach rechts von Tspace anfang
		self.Tclosed2 =  0.0 # verschiebung nach links von Tspace ende
		
		### need this from boundaryConditionType1 ####
		## is set with update method in vascularNetwork in calculateInitialValues
		self.Tpulse = 0.0
		self.Tspace = 0.0
		self.freq   = 1.
		self.Npulse	= 1.
		self.TmeanFlow = 0.0
		###
		
		self.pulseTimeR = []
		
		self.returnFunction = None
		
		self.omegaNew = np.empty((2))
		self.dQInOut = np.empty((2))
	
	def __call__(self, _domegaField_, duPrescribed, R, L, n, dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_,duPrescribed,L,R,n,dt)
		
	def funcPos0(self,_domegaField, duPrescribed, L, R, n, dt):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		Rt = self.calcRt(n,dt)
		domegaPrescribed_ = np.dot(L[0],duPrescribed)
		domegaReflected_ = _domegaField*Rt + domegaPrescribed_
		
		self.omegaNew[0] = domegaReflected_
		self.omegaNew[1] = _domegaField
				
		self.dQInOut = R[:][1]*self.omegaNew 
				
		return np.dot(R,self.omegaNew), self.dQInOut
		

	def funcPos1(self,domegaField_, duPrescribed, L, R, n, dt):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		Rt = self.calcRt(n,dt)
		_domegaPrescribed = np.dot(L[1],duPrescribed)
		_domegaReflected = domegaField_*Rt + _domegaPrescribed
		
		self.omegaNew[0] = domegaField_
		self.omegaNew[1] = _domegaReflected
		
		self.dQInOut = (R[:][1]*self.omegaNew)[::-1].copy()	
		
		return np.dot(R,self.omegaNew), self.dQInOut
		


	def update(self,bcDict):
		'''
		updates the updateBoundaryDict data using a dictionary in from of 
		bcDict = {'variableName': value}
		'''
		for key,value in bcDict.iteritems():
			try:
				self.__getattribute__(key)
				self.__setattr__(key,value)
			except: pass
				#print 'ValueError: wrong key: %s, could not set up boundaryCondition' %key	

		self.type = 2
		self.name = 'ReflectionCoefficientTimeVarying'

		self.pulseTimeR = []
		
		#let the valve be open 20 + of the pulse
		#self.Topen1 = - (0.5/self.freq)*0.1 
		#self.Topen2 =   (0.5/self.freq)*0.1 
		
		self.Tperiod = self.Tspace+1.0/self.freq
		
		for puls in np.arange(self.Npulse):
			self.pulseTimeR.append([self.Tpulse + (self.Tperiod)*puls + self.Tclosed2 - self.TmeanFlow,
									self.Tpulse + 0.5/self.freq+ (self.Tperiod)*puls - self.TmeanFlow + self.Topen1,
									self.Tpulse + 0.5/self.freq+ (self.Tperiod)*puls - self.TmeanFlow + self.Topen2,
									self.Tpulse + 1.0/self.freq+ (self.Tperiod)*puls - self.TmeanFlow + self.Tclosed1,
									self.Tpulse + self.Tperiod +(self.Tperiod)*puls  - self.TmeanFlow + self.Tclosed2,
									self.Tpulse + self.Tperiod +(self.Tperiod)*puls + self.Tclosed2 ])
	

	def calcRt(self,n,dt):
		t = n*dt
		RC = self.RtClosed
		
		# get time slots
		for pTA in self.pulseTimeR:
			if pTA[0]< t and t < pTA[5]:
				#pulsNum = self.pulseTimeR.index(pTA)
				pTA1 = pTA[0]
				pTA2 = pTA[1]
				pTA3 = pTA[2]
				pTA4 = pTA[3]
				pTA5 = pTA[4]
				pTA6 = pTA[5]
				break
		try:
			if (pTA1 < t and t < pTA2):
				RC = (self.RtOpen-self.RtClosed)/(pTA2-pTA1)*(t-pTA1) +self.RtClosed
			elif pTA2 <= t and t <= pTA3: 
				RC = self.RtOpen
			elif pTA3 < t and t < pTA4:
				RC = -(self.RtOpen-self.RtClosed)/(pTA4-pTA3)*(t-pTA3) + self.RtOpen
			elif (pTA5 < t and t < pTA6):
				#RC = self.RtOpen
				RC = (self.RtOpen-self.RtClosed)/(pTA2-pTA1)*(t-pTA5) +self.RtClosed
			else: pass
			return RC
		except:
			return RC
		

class Resistance(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	signle Resistance element
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	
	def __init__(self):
		self.type = 2
		self.Rc = 1
		
		self.venousPressure = 20.*133# not needed
		
		self.returnFunction = None
		self.omegaNew = np.empty((2))
		self.dQInOut = np.empty((2))
		
	def __call__(self, _domegaField_, duPrescribed, R, L, n, dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_,R,Z1)
	
	def funcPos0(self,_domegaField, R, Z1):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		r11,r12,r21,r22 = R[0][0],R[0][1],R[1][0],R[1][1]
		
		Rc = self.Rc
		if Rc == 'VesselImpedance':
			Rc = Z1
		
		domega_ = -(r22*Rc + r12)/(r21*Rc + r11)*_domegaField

		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domegaField
		
		self.dQInOut = R[:][1]*self.omegaNew 
				
		return np.dot(R,self.omegaNew), self.dQInOut
			
	
	def funcPos1(self,domegaField_, R, Z1):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		r11,r12,r21,r22 = R[0][0],R[0][1],R[1][0],R[1][1]
		
		Rc = self.Rc
		if Rc == 'VesselImpedance':
			Rc = Z1
		
		_domega = (r11 - r21*Rc)/(r22*Rc - r12)*domegaField_
			
		self.omegaNew[0] = domegaField_
		self.omegaNew[1] = _domega
		
		self.dQInOut = (R[:][1]*self.omegaNew)[::-1].copy()	
		
		return np.dot(R,self.omegaNew), self.dQInOut
			

class Windkessel2(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	2 Element Windkessel 
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	
	def __init__(self):
		self.type = 2
		
		self.Rc = 1
		self.C  = 0
		
		self.venousPressure = 7.*133
		
		self.returnFunction = None
		self.omegaNew = np.empty((2))
		self.dQInOut = np.empty((2))

	def __call__(self, _domegaField_, duPrescribed, R, L, n, dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_, R, dt, P, Q, n)
	
	def funcPos0(self,_domegaField, R, dt, P, Q, n):
		"""return function for position 0 at the start
		of the vessel
		"""
		
		"""Knut Petter: The old version of the WK2 scheme is shown in the comments so you can compare it to the new version,
		 the new scheme is based on the "half step central difference"-scheme shown in my master thesis, which was tested and gave good results.
		
		Old version:
		P0 = du[0] #/2.take this out is old ##??
		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]
		
		if self.Rc == None: self.Rc = -1./R[1][1]	
		tau = 2.*self.Rc*self.C
		taudt = tau/dt
		denom = 1.+taudt+r21*self.Rc
		domega_ = ((-1.-taudt-r22*self.Rc)*_domega_ + taudt*(dw_+_dw) + P0)/denom
		return np.array([domega_ , _domega_])

		"""
		''' Version does not take hole R matrix into account
		r21,r22 = R[1][0],R[1][1]
		taudt = 2*self.Rc*self.C/dt
		a = 1 + taudt + self.Rc*r21
		b = 1 + taudt + self.Rc*r22
		c = 2*(self.Rc*Q + P - self.venousPressure)
		
		domega_ = (-b*_domegaField_ - c)/a
		#return np.array([domega_ , _domega])
		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domegaField_
		return self.omegaNew	
		'''
		r11,r12,r21,r22 = R[0][0],R[0][1],R[1][0],R[1][1]
		
		taudt = self.Rc*self.C/dt
		a = -self.Rc*r21 - (1. + 2.*taudt)*r11
		b = (2.*taudt+1.)*r12 + self.Rc*r22
		
		domega_ = (2*(self.Rc*Q + (P - self.venousPressure[n])) + b*_domegaField)/a
		
		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domegaField
		
		self.dQInOut = R[:][1]*self.omegaNew 
		
		return np.dot(R,self.omegaNew), self.dQInOut
		
		
	def funcPos1(self,domegaField_, R, dt, P, Q, n):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		
		"""Old version:
		P0 = du[0] #/2.take this out is old ##??
		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]
		
		if self.Rc == None: self.Rc = 1./R[1][0]
		tau = 2.*self.Rc*self.C
		taudt = tau/dt
		denom = 1.+taudt-r22*self.Rc
		_domega = ((-1.-taudt+r21*self.Rc)*_domega_ + taudt*(dw_+_dw) + P0)/denom
		return np.array([_domega_ , _domega])
		"""
		''' Version does not take hole R matrix into account
		r21,r22 = R[1][0],R[1][1]
		taudt = 2*self.Rc*self.C/dt
		a = taudt - self.Rc*r21 - 1
		b = taudt - self.Rc*r22 - 1
		c = 2*(self.Rc*Q - (P - self.venousPressure))
		
		_domega = (-a*domega_ + c)/b
		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domega
		return self.omegaNew	
		#return np.array([domega_ , _domega])
		'''
		r11,r12,r21,r22 = R[0][0],R[0][1],R[1][0],R[1][1]
		
		taudt = self.Rc*self.C/dt
		
		a = self.Rc*r21 - (1. + 2.*taudt)*r11
		b = (2.*taudt+1.)*r12 - self.Rc*r22
		
		_domega = (2*(self.Rc*Q - (P - self.venousPressure[n])) + a*domegaField_)/b
		
		self.omegaNew[0] = domegaField_
		self.omegaNew[1] = _domega
		
		self.dQInOut = (R[:][1]*self.omegaNew)[::-1].copy()
		
		return np.dot(R,self.omegaNew), self.dQInOut
			
		
class Windkessel3(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	3 Element Windkessel 
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	
	def __init__(self):
		self.type = 2
		
		# parameters in xml
		self.Rc  = None
		self.Z   = 'VesselImpedance'
		self.C   = 1
		self.Rtotal = None
		# parameters for calculation
		self.Z0    = None # vesselimpedance at first step
		self.RcNum = None
			
		self.venousPressure = 0 #7.*133.
		self.returnFunction = None
		self.omegaNew = np.empty((2))
		self.dQInOut = np.empty((2))
		
		self.firstRun = False
		
	def __call__(self, _domegaField_, duPrescribed, R, L, n, dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_, R, dt, P, Q, Z1, Z2, n)
	
	def funcPos0(self,_domegaField, R, dt, P, Q, Z1, Z2, n):
		'''
		return function for position 0 at the start
		of the vessel
		
		Old version
		venousPressure = P

		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]
		
		if self.R1 == None: self.R1 = -1./r22
		if self.Rc == None: self.Rc = self.RT - self.R1 
		tau = 2.*self.Rc*self.C
		taudt = tau/dt
		denom = taudt + 1. + r21*(self.R1*taudt + self.Rc + self.R1)
		domega_ = ((1.+r21*self.R1)*taudt*dw_ + (1.+r22*self.R1)*taudt*_dw + (-1.-taudt-r22*(self.R1*taudt + self.Rc + self.R1))*_domega_ + P0)/denom
		
		Old Version: does not take hole R matrix into account
		r21,r22 = R[1][0],R[1][1]
		
		Z = self.Z 
		if self.Z == 'VesselImpedance':# and self.firstRun == False: 
			Z = Z1 #+Z2
			#self.Z = Z
			#self.firstRun = True
		Rc = self.Rc
		if self.Rc == None:
			Rc = self.Rtotal - Z
			## Rc not time-varying activate this line:
			#self.Rc = Rc
		C = self.C
			
		taudt = 2*Rc*C/dt		
		
		a = taudt + 1
		b = (taudt +1)*Z + Rc
		c = 2*Q*(Rc + Z) - 2*(P-self.venousPressure)
		
		domega_ = -((a + b*r21)*_domega_ - c)/(a + b*r22)
		
		
		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domega_
		return self.omegaNew	
		'''
		r11,r12,r21,r22 = R[0][0],R[0][1],R[1][0],R[1][1]
		
		
		
		
		if self.Z == 'VesselImpedance' and self.firstRun == False: 
			Z = Z1 
			## Z not time-varying activate this lines:
			self.Z0 = Z1
			self.firstRun = True
		elif self.firstRun == True:
			Z = self.Z0
		else:
			Z = self.Z 
				
		Rc = self.Rc
		if self.Rc == None:
			Rc = self.Rtotal - Z
			## Rc not time-varying activate this line:
			#self.Rc = Rc
			##
		
		C = self.C
		
		a = Z*Rc*C+0.5*dt*(Z+Rc)
		b = 0.5*dt+Rc*C 
		
		domega_ = (dt*(P-self.venousPressure[n])+dt*(Z+Rc)*Q+_domegaField*(a*r22+b*r12))/(-a*r21-b*r11)
		
		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domegaField
		
		self.dQInOut = R[:][1]*self.omegaNew 
				
		return np.dot(R,self.omegaNew), self.dQInOut
			
		
	
	def funcPos1(self, domegaField_, R, dt, P, Q, Z1, Z2, n):
		'''
		return function for position -1 at the end
		of the vessel  (newer version from Knut Petter ??? )
		
		Knut Petter: Yes this is a newer version based on the "half-step central difference"-scheme (shown for WK2 in my master thesis, but not WK3)
		I am also putting the old version back into the comments so the two can be compared.
		
	
		Old Version:
		P0 = du[0] #/2.take this out is old ##??

 		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]

		if self.R1 == None: self.R1 = 1./r21
		if self.Rc == None: self.Rc = self.RT - self.R1 
		tau = 2*self.Rc*self.C
		taudt = tau/dt	
		denom =  taudt + 1. - r22*(self.R1*taudt + self.Rc + self.R1)
		_domega = ((1.-r21*self.R1)*taudt*dw_ + (1.-r22*self.R1)*taudt*_dw + (-1.-taudt+r21*(self.R1*taudt + self.Rc + self.R1))*_domega_ + P0)/denom
		return np.array([_domega_ , _domega])
		
		Old Version Not taking total R-matrix into account:
		taudt = 2*Rc*C/dt
		
		a = taudt + 1
		b = (taudt +1)*Z + Rc
		c = 2*Q*(Rc + Z) - 2*(P-self.venousPressure)
		
		_domega = -((a + b*r21)*domega_ - c)/(a + b*r22)
		
		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domega
		return self.omegaNew	
		'''
		r11,r12,r21,r22 = R[0][0],R[0][1],R[1][0],R[1][1]
		
		
		if self.Z == 'VesselImpedance' and self.firstRun == False: 
			Z = Z1 
			## Z not time-varying activate this lines:
			self.Z0 = Z1
			self.firstRun = True
		elif self.firstRun == True:
			Z = self.Z0
		else:
			Z = self.Z 
		
		Rc = self.Rc
		if self.Rc == None:
			Rc = self.Rtotal - Z
			## Rc not time-varying activate this line:
			#self.Rc = Rc
		
		C = self.C
		
		a = Z*Rc*C+0.5*dt*(Z+Rc)
		b = 0.5*dt+Rc*C 
		
		_domega = (dt*(P-self.venousPressure[n])-dt*(Z+Rc)*Q+domegaField_*(b*r11-a*r21))/(a*r22-b*r12)
		
		self.omegaNew[0] = domegaField_
		self.omegaNew[1] = _domega
		
		self.dQInOut = (R[:][1]*self.omegaNew)[::-1].copy()	
		
		return np.dot(R,self.omegaNew), self.dQInOut
		


class L_network(BoundaryConditionType2):
	"""	
	Boundary profile - type 2
	
	L-network
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	return:
	the domega-vector with (domega_ , _domega) based on the input values and its returnFunction
	"""
	
	def __init__(self):
		self.type = 2
		
		self.C  = 0
		self.R1 = 1

		self.returnFunction = None
		self.omegaNew = np.empty((2))

	def __call__(self, _domegaField_, duPrescribed, R, L, n, dt, P, Q, A, Z1, Z2):
		return self.returnFunction(_domegaField_, duPrescribed, R, L, n, dt)
	
	def funcPos0(self,_domega_, du, R, dt):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		print "ERROR: boundaryCondition Lnet is not implemented correct!"
		exit()
		
		dQ0 = du[1]/2.0

		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = 23,23 ## dO[self.position][0],dO[self.position][1]
		if self.R1 == None: self.R1 = -1./R[1][1]	
		tau = 2.*self.R1*self.C
		taudt = tau/dt	
		denom = taudt + taudt*r21*self.R1 + r21*self.R1
		domega_ = ((1.+r21*self.R1)*taudt*dw_ + (1.+r22*self.R1)*taudt*_dw + (-taudt - taudt*r22*self.R1 - r22*self.R1)*_domega_ + self.R1*dQ0)/denom
		#return np.array([domega_ , _domega_])
		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domega_
		
		return np.dot(R,self.omegaNew)
	
	def funcPos1(self,_domega_, du, R, dt):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		
		print "ERROR: boundaryCondition Lnet is not implemented correct!"
		exit()
		
		dQ0 = du[1]/2.0
		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = 23,23 #dO[self.position][0],dO[self.position][1]
		if self.R1 == None: self.R1 = 1./R[1][0]	
		tau = 2.*self.R1*self.C
		taudt = tau/dt
		denom = -taudt + taudt*r22*self.R1 + r22*self.R1
		_domega = ((-1.+r21*self.R1)*taudt*dw_ + (-1.+r22*self.R1)*taudt*_dw + (taudt - taudt*r21*self.R1 - r21*self.R1)*_domega_ + self.R1*dQ0)/denom
		
		self.omegaNew[0] = _domega_
		self.omegaNew[1] = _domega
		
		return np.dot(R,self.omegaNew)
	
class VaryingElastanceOld(BoundaryConditionType2):
	"""
	An implementation of a time-varying elastance model of the left ventricle,ing based on the modfied varying elastance equation including a source resistance K
	(as proposed by Shroff), and a parametrized time varying elastance function (as given by Stergiopulos).
	
	The general shape of the elastance function is given by three shape parameters. Various conditions of heart contractility are then created by scaling this 
	function using the parameters:
	T - Heart period
	Emax - Maximum elastance
	Emin - Minimum elastance
	Tpeak - Time to peak elastance
	
	The equation also requires the two constants:
	V0 - Volume axis intercept
	K - Source resistance
	
	NB! The source resistance (K) was introduced to the implementation as an experiment. It modifies the elastance curve depending on ventricular outflow,
	so that it becomes dependent on the afterload of the heart. Introducing the source resistance did produce a load dependence 
	(shown by curved isochrones in p-v loops), however the results are in no way to be trusted since the modified varying elastance equation was not intended
	to be used together with the specific varying elastance curve shape used here (Stergiopulos). A proper implementation of the source resistance 
	requires a different curve shape. The parameter K is therefore set to zero by default, but it should perhaps be removed from the code altogether??
	
	Currently only the return method "def funcPos0" has been implemented so that the boundary condition can only be put at the proximal end of a blood vessel.
	It is fairly straightforward to implement funcPos1 if necessary, this does however require a lot of duplicated code.   """
	def __init__(self):
		self.type = 2
		
		self.subiterations = 0
		
		self.omegaNew = np.empty((2))
		
		#Default parameters
		self.T     = 1
		self.Emax  = 2.31 * 133.3e6 
		self.Emin  = 0.06 * 133.3e6 
		self.Tpeak = 0.4
		
		self.V0 = 20e-6 
		
		self.K  =  0.0
		
		"""Shape parameters"""
		self.alpha = 1.672
		self.n1 = 1.32
		self.n2 = 21.9

		#n-1 values
		self.aorticFlowPreviousTimestep = None
		
		self.system = {'both open':np.array([0,1,2]), 'mitral open': np.array([0,1]), 'aortic open':np.array([1,2])} 
		
				
		self.cycleNumber = 0
		self.atriumPressure = 7.5 * 133.32 #Pressure in the atrium ## venouse pressure?!
		
		self.x0 = np.array([0.0, 0.0, 0.0]) #Initial values for the iterative solver
		
		self.mitral = None # mitral valve
		self.aortic = None # aortic valve
		self.initializeValves() # intialize valves
		
		self.dQInOut = np.empty((2))
		
		
	def initializeSolutionVectors(self, Tsteps):
		"""Initializes some solution vectors storing pressure, flow and volume of the ventricle, as well as opening and closing state 
		
		NB! This method is not called from the class constructor, but is called externally by the initializeSolutionMatrices method in the solver,
		this is a bit messy, but was the easiest way to do it since the BC is initiated before the number of time steps is known.
		"""
		
		""" Initialize Solution Vectors """
		
		print """ Initialize Solution Vectors """
		
		self.mitral.initializeSolutions(Tsteps)
		self.aortic.initializeSolutions(Tsteps)
	
		self.pressure = np.zeros(Tsteps)
		self.volume = np.zeros(Tsteps)
		self.mitralQ = np.zeros(Tsteps)
		self.Elastance = np.zeros(Tsteps) #New
		self.Flow = np.zeros(Tsteps)
		self.DtFlow = np.zeros(Tsteps)
		self.Turb=np.zeros(Tsteps)
		self.Inert=np.zeros(Tsteps)
		self.InbyTurb=np.zeros(Tsteps)
		self.deltaP=np.zeros(Tsteps)
		self.aortaP=np.zeros(Tsteps)
		
		""" Initial conditions in the ventricle"""
		self.pressure[0] = self.atriumPressure
		self.volume[0] = self.atriumPressure/self.E(0) + self.V0
				
		
	def initializeValves(self):
		"""Mitral valve parameters"""
		
		self.mitral_annulus_area = 0.0007
		
		mitral_M_st          = 0.4
		mitral_M_rg = 0.0
		mitral_delta_p_open = 0
		mitral_delta_p_close = 0 
		mitral_K_v_open = 0.3
		mitral_K_v_close = 0.4

		
		"""Aortic valve parameters"""
		aortic_M_st          = 1
		aortic_M_rg          = 0.00
		aortic_delta_p_open  = 0*133.32
		aortic_delta_p_close = 2*133.32 # 2mmHg
		aortic_K_v_open      = 0.12
		aortic_K_v_close     = 0.12

		
		#Create valves
		self.mitral = Valve(mitral_M_st, mitral_M_rg, mitral_delta_p_open, \
						mitral_delta_p_close, mitral_K_v_open, mitral_K_v_close)
		self.aortic = Valve(aortic_M_st, aortic_M_rg, aortic_delta_p_open, \
						aortic_delta_p_close, aortic_K_v_open, aortic_K_v_close)
	
	def __call__(self, _domegaField_, duPrescribed, R, L, n, dt, P, Q, A, Z1, Z2):
	
		self.updateValves(P, n, dt)                     # Update the state of the mitral and aortic valve at timestep n + 1
		self.startNewCycleIfCriteriaIsMet(n, dt)
		self.funcPos0(_domegaField_, R, n, dt, P, Q, A)      # Compute the riemann variant going into the vessel save in omegaNew
		
		self.dQInOut = R[:][1]*self.omegaNew 
		# calculate du and return this!
		return np.dot(R,self.omegaNew),self.dQInOut


	def updateValves(self, P, n, dt):
		mitralPressureDifference = self.atriumPressure - self.pressure[n]
		self.mitral.updateValveState(mitralPressureDifference, n, dt)
		
		aorticPressureDifference  = self.pressure[n] - P
		self.aortic.updateValveState(aorticPressureDifference, n, dt)


	def getCycleTime(self, n, dt):
		return n*dt - self.T*self.cycleNumber
	
	def startNewCycleIfCriteriaIsMet(self, n, dt):
		if self.getCycleTime(n+1, dt) > self.T:
			self.cycleNumber += 1

	def funcPos0(self, _domega, R, n, dt, Pn, Qn, A):
		
		# Qn1 == value at old time step
		# change to self.aorticPressurePreviousTimestep ...
		
		Qn1 = self.aorticFlowPreviousTimestep
		
		
		r11,r12,r21,r22 =  R[0][0],R[0][1],R[1][0],R[1][1]
		L = self.aortic.computeL(A, n+1)
		LdivB = self.aortic.LdivideB(A, n+1)
		B = self.aortic.computeB(A, n+1)
		mitrL = self.mitral.computeL(self.mitral_annulus_area, n+1)
		mitrLdivB = self.mitral.LdivideB(A, n+1)
		mitrB = self.mitral.computeB(self.mitral_annulus_area, n+1) #
		mitrQn = self.mitralQ[n]
		mitrQn1 = self.mitralQ[n-1]
		venoP = self.atriumPressure
		t = self.getCycleTime(n+1, dt)
# 		ttemp = t-dt
		E = self.E(t)
		Vn = self.volume[n]
		self.Elastance[n+1]=E/133.3e6
		self.Flow[n]=Qn*1e6
		self.aortaP[n]=Pn
#		self.DtFlow[n]=(Qn-Qnold)/dt
		ventrPn = self.pressure[n]
		print "B is: ", B
		print "pn is: ", Pn
		print "ventrPn is: ", ventrPn
		print "ventricle pressure should be:", E*(Vn-self.V0)
		print "ventrPn-Pn is:", ventrPn-Pn
# 		print "volume is", Vn
# 		print "E is,", E
		print "t is:",t
		if B:
			self.DtFlow[n]=1e6*(Qn-Qn1)/dt
			self.Turb[n]=abs(Qn)*Qn*B/133
			self.Inert[n]=L*(Qn-Qn1)/(dt*133)
			self.deltaP[n]=abs(Qn)*Qn*B/133+L*(Qn-Qn1)/(dt*133)
#			self.aortaP[n]=self.pressure[n]-self.deltaP[n]
#			self.aortaP[n]=Pn
			try:
				self.InbyTurb[n]=(abs(Qn)*Qn*B)/(L*(Qn-Qn1)/dt)
				
			except: pass
			print "Qn is: ", Qn
			print "dp1 (B part) is:", abs(Qn)*Qn*B
			print "dp2 (L part)is:", L*(Qn-Qn1)/dt
			print "dP total is:", abs(Qn)*Qn*B+L*(Qn-Qn1)/dt
		
		""" Because of the large differences in magnitude between pressure and flow in the currently used dimensions some attemts were made to scale the
		variables using for example these scaling factors for B,p, q,and omega"""
		B_ref = 1060/(2*A**2)
		n_p = 1.#self.Emax*self.V0
		n_q = 1.#(n_p/B_ref)**0.5
		n_o = 1.#n_q/r21
		if Qn<-20e-6:
			B=None
			
		args = dt, mitrLdivB, mitrB,LdivB, L, mitrL, B, mitrQn1, mitrQn, ventrPn, venoP, E, Vn, Qn, Qn1, r11, r12, r21, r22, Pn, _domega, n_q, n_p, B_ref
		
		
		"""The following section computes the increment domega_ which goes into the vessel from the ventricle, """
		
		if not B and not mitrB:
			"""both valves are closed: """
			
			
			self.mitralQ[n+1] = 0 
			domega_ = _domega
			self.volume[n+1] = self.volume[n] - 0.5*(Qn - mitrQn)*dt
			self.pressure[n+1] = E*(self.volume[n+1] - self.V0)
			print "both closed"
			if Qn == 0:
				domega_ = _domega
			else:
				domega_ = (-0.5*Qn - r22*_domega)/r21 #Correction for non-zero flow at full aortic valve closure
				
			self.x0 = np.array([0,0,domega_])
		else:
			if not mitrB:
				"""only the aortic valve is open"""
				x = self.newtonSolver(self.x0,args, partialSystem='aortic open')
				self.mitralQ[n+1] = 0
				print "aortic open"
				self.pressure[n+1] = self.pressure[n] + x[0]*n_p
				domega_  = x[1]*n_o
				self.x0 = np.concatenate((np.array([0]), x))
				
			elif not B:
				"""only  the mitral valve is open"""
				x = self.newtonSolver(self.x0,args, partialSystem='mitral open')
				self.mitralQ[n+1] = self.mitralQ[n] + x[0]*n_q
				self.pressure[n+1] = self.pressure[n] + x[1]*n_p
				
			
				if Qn == 0:
					domega_ = _domega
				else:
					domega_ = (-0.5*Qn - r22*_domega)/r21 #Correction for non-zero flow at full aortic valve closure
					
				self.x0 = np.concatenate((x, np.array([0])))	
			else:
				"""both valves are open"""
				x = self.newtonSolver(self.x0,args)
				print "both open"
			
				self.mitralQ[n+1] = self.mitralQ[n] + x[0]*n_q
				self.pressure[n+1] = self.pressure[n] + x[1]*n_p
				domega_ = x[2]*n_o
				self.x0 = x
				
			dQ = r21*domega_ + r22*_domega
			self.volume[n+1] = Vn - (Qn + 0.5*(-self.mitralQ[n] - self.mitralQ[n+1] + dQ))*dt
			
		self.aorticFlowPreviousTimestep = Qn
#		Qnold=Qn
		
		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domega
		
		


	def funcPos1(self, _domega, R, L, n, dt, P, Q, A):
		pass
	
	def E(self, t):
		"""Computes the value of the elastance at time t, according to the shape parameters given by Stergiopolus and scaled
		   according to Tpeak, T, Emax and Emin. """
		a1 = 0.708*self.Tpeak
		a2 = 1.677*a1
		
		n1, n2 = self.n1, self.n2
		shapeFunction1 = (t/(a1))**n1/(1+(t/(a1))**n1)
		shapeFunction2 = (1 + (t/(a2))**n2)**(-1)	
		return (self.Emax-self.Emin)*self.alpha*shapeFunction1*shapeFunction2 + self.Emin

	def newtonSolver(self, x0, args,  partialSystem = 'both open'):
		"""Solves the partial or full equation system of the varying elastance model"""
		maxIterations = 20 
		iterations = 0

		xn = x0[self.system[partialSystem]]
		res = self.solverResiduals(xn,*args, partialSystem = partialSystem)
		#error = np.linalg.norm(res, 2)
								
		while True:
			#print x0
			iterations +=1
			J_inv = self.solverInverseJacobian(xn, *args, partialSystem = partialSystem)
			
 			
			
			x = xn - np.dot(J_inv, res).T
			
						
			error = np.linalg.norm(x - xn, 2)/np.linalg.norm(xn, 2)
			if error < 0.0001:
				break
			xn = x
			res = self.solverResiduals(x, *args, partialSystem = partialSystem)

			if iterations > maxIterations: 
				x *= 0
				break

		return x
		
	def solverResiduals(self, x_partial, dt, mitrLdivB, mitrB, LdivB, L, mitrL, B, mitrQn1, mitrQn,ventrPn, atrP, E, Vn, Qn, Qn1, r11, r12, r21, r22, Pn, _domega, n_q, n_p, B_ref, partialSystem = np.array([0,1,2])):#dt, mitrL, mitrB,L,B, mitrQn1, mitrQn, venoP, E, Vn, Qn, Qn1, r21, r22, Pn, _domega):
		"""Computes  are the resisduals of the functions f1,f2 and f3, they are defined as functions that are only called when they are needed. The
		argument partialSystem determines which of the residuals are computed and returned."""
		
		x = np.array([0.0, 0.0, 0.0])
		x[self.system[partialSystem]] += x_partial
		dQm, dPv, domega_ = x 
		
# 		Knut Petters version:
		"""		
		def f1(): 
			a = mitrQn/n_q + dQm
			return a*abs(a) + mitrLdivB/(2*n_q*dt)*(3*dQm + (mitrQn1 - mitrQn)/n_q) + (n_p*dPv + ventrPn - atrP)/(mitrB*n_q**2)
		def f2(): 
			return E/n_p*(Vn - (Qn - mitrQn + 0.5*(n_q*domega_ + r22*_domega - n_q*dQm))*dt - self.V0)*(1-self.K*(Qn + n_q*domega_ + r22*_domega)) - ventrPn/n_p - dPv
		def f3():
			a = (Qn + r22*_domega)/n_q + domega_
			return a*abs(a) + LdivB/(2*n_q*dt)*(3*domega_ + (3*r22 - Qn + Qn1)/n_q) + (n_q/r21*domega_ + _domega + Pn - ventrPn - n_p*dPv)/(B*n_q**2)
		"""
		
		def f1():
			a = mitrQn/n_q + dQm
			return mitrB*a*abs(a)+(0.5/dt)*mitrL*(3*dQm + (mitrQn1 - mitrQn)/n_q)+(n_p*dPv + ventrPn - atrP)
		def f2():
			return E/n_p*(Vn - (Qn - mitrQn + 0.5*(n_q*domega_*r21 + r22*_domega - n_q*dQm))*dt - self.V0)*(1-self.K*(Qn + n_q*domega_*r21 + r22*_domega)) - ventrPn/n_p - dPv
		
		def f3():
			a = (Qn + r22*_domega)/n_q + domega_*r21
			return B*a*abs(a)+(0.5/dt)*L*(3*domega_*r21 + (3*r22*_domega - Qn + Qn1)/n_q)+ (r11*domega_ + _domega*r12 + Pn - ventrPn - n_p*dPv)
		
		def f1simple():
			a= mitrQn +dQm
			return mitrB*a*abs(a)+0.5*mitrL*(3*dQm+mitrQn1-mitrQn)/dt-atrP+ventrPn+dPv
		
		def f2simple():
			return E*(Vn-0.5*dt*(r21*domega_+r22*_domega+2*Qn-2*mitrQn-dQm)-self.V0)-ventrPn-dPv
		
		def f3simple():
			a=r21*domega_+r22*_domega+Qn
			return B*a*abs(a)+(L/(2*dt))*(3*r21*domega_+3*r22*_domega-Qn+Qn1)-ventrPn -dPv + Pn + r11*domega_ + _domega*r12
			
		functions = np.array([f1simple, f2simple, f3simple])
		return np.array([f() for f in functions[self.system[partialSystem]]])
	
	def solverInverseJacobian(self, x_partial, dt, mitrLdivB, mitrB,LdivB, L, mitrL, B, mitrQn1, mitrQn,ventrPn, venoP, E, Vn, Qn, Qn1, r11, r12, r21, r22, Pn, _domega, n_q, n_p,B_ref, partialSystem = np.array([0,1,2])):
		"""Computes the inverse Jacobian of the system. The components a1, a2, a3, a4, a5 and a6 are declared using strings, and evaluated using eval()
		 only when needed. (Not sure how smart this is). Reduces lines of code, but is probably slower."""
		
		
		x = np.array([0, 0, 0])
		x[self.system[partialSystem]] += x_partial 
		dmQ, dvP, domega_ = x
		
		if partialSystem == 'mitral open':
			
			a1 = mitrB*2*(mitrQn/n_q + dmQ)*np.sign(mitrQn/n_q + dmQ)+1.5*mitrL/dt
			a2 = 0.5*E*dt*(1-self.K*(r21*domega_ + r22*_domega + Qn))
			
			a1simple = mitrB*2*abs(mitrQn+dmQ)+1.5*mitrL/dt
			a2simple = 0.5*E*dt
			
			J_inv = np.array([[1,  1], 
							  [a2simple, -a1simple]])/(a1simple+a2simple)
			return J_inv
		elif partialSystem == 'aortic open':
			a3 = -0.5*E*dt*(1-self.K*(r21*domega_ + r22*_domega + Qn))-r21*E*self.K*(Vn -(Qn - mitrQn + 0.5*(n_q*domega_*r21 + r22*_domega - n_q*dmQ))*dt - self.V0)
			a4 = B*r21*2*(Qn + r21*domega_+ r22*_domega)*np.sign(Qn + r21*domega_+ r22*_domega)+1.5*L*r21/dt+r11
			
			a3simple = -0.5*E*dt*r21
			a4simple = 2*r21*B*abs(r21*domega_+r22*_domega+Qn)+1.5*L*r21/(dt)+r11
			
			J_inv =np.array([[a4simple, -a3simple], 
							 [1, -1]])/(-a4simple +a3simple)
			return J_inv
		else:
			a1 = mitrB*2*(mitrQn/n_q + dmQ)*np.sign(mitrQn/n_q + dmQ)+1.5*mitrL/dt
			a2 = 0.5*E*dt*(1-self.K*(r21*domega_ + r22*_domega + Qn))
			a3 = -0.5*E*dt*(1-self.K*(r21*domega_ + r22*_domega + Qn))-r21*E*self.K*(Vn -(Qn - mitrQn + 0.5*(n_q*domega_*r21 + r22*_domega - n_q*dmQ))*dt - self.V0)
			a4 = B*r21*2*(Qn + r21*domega_+ r22*_domega)*np.sign(Qn + r21*domega_+ r22*_domega)+1.5*L*r21/dt+r11

			J_inv = np.array([[ -a4+a3,  -a4,  a3  ],
						      [ -a2*a4, a1*a4, -a1*a3 ],
						      [  -a2, a1, a1+a2 ]])/(a1*a3 - a1*a4 - a2*a4)
			

			return J_inv 
		
		
		
		"""
		expressions = np.array([
		'2*(mitrQn/n_q + dmQ)*np.sign(mitrQn/n_q + dmQ) + 1.5*mitrLdivB/(n_q*dt)',
		'B_ref/mitrB',
		'0.5*n_q/n_p*E*dt*(1-self.K*(Qn + n_q*domega_ + r22*_domega))',
		'-0.5*n_q/n_p*E*dt*(1-self.K*(Qn + n_q*domega_ + r22*_domega)) - E*self.K*n_q/n_p*(Vn -(Qn - mitrQn + 0.5*(n_q*domega_ + r22*_domega - n_q*dmQ))*dt - self.V0)',
		'-B_ref/B',
		'2*((Qn + r22*_domega)/n_q + domega_)*np.sign((Qn + r22*_domega)/n_q + domega_) + 1.5*LdivB/(n_q*dt) + 1/(r21*B*n_q)'])




		

		if partialSystem == 'mitral open':
			a1, a2, a3 = [eval(e) for e in expressions[[0,1,2]]]
			
			J_inv = np.array([[-1,  -a2], 
							  [-a3, a1]])/(-a1-a2*a3)
			return J_inv
		elif partialSystem == 'aortic open':
			a4, a5,a6 = [eval(e) for e in expressions[[3,4,5]]]
			J_inv =np.array([[a6, -a4], 
							 [-a5, -1]])/(-a6 - a4*a5)
			return J_inv
		else:
			a1, a2, a3, a4, a5, a6 = [eval(e) for e in expressions]

			J_inv = np.array([[ a6+a4*a5,  a2*a6,  -a2*a4   ],
						      [    a3*a6, -a1*a6,   a1*a4   ],
						      [   -a3*a5,  a1*a5,   a1+a2*a3]])/(a3*a2*a6 + a1*a6 + a1*a4*a5)
			

			return J_inv 
		
		"""
		

class Valve:
	""" A valve model meant to be used with the varying elastance boundary condition. Two instances of the class are then made representing the mitral and 
	aortic valves. The valve instance has a variable state, which is updated every timestep using the updateValveState method and is then able to return
	updated values of the coefficients B and L, used by the varying elastance solver. 
	"""
	
	def __init__(self, M_st, M_rg, delta_p_open, delta_p_close, K_v_open, K_v_close):
		
		"""
		M_st - Determines the effective orifice area (EOA) at maximum opening. M_st ~= 1 => Healthy valve
																			   M_st  < 1 => Stenosed valve
																			   (If set to 1 there is no pressure difference, meaning that the valve can not close
																			   once it has fully opened, 0.99 has been used)
               																			   
		M_rg - Determines the EOA at minimum opening. M_rg = 0 => Healthy valve
													  M_rg > 0 => Regurgitant valve
															  														  
		"""
		self.M_st = M_st
		self.M_rg = M_rg
		self.delta_p_open = delta_p_open
		self.delta_p_close = delta_p_close
		self.K_v_open = K_v_open
		self.K_v_close = K_v_close
		
		self.rho = 1060 ## default value is update while start up in classVascularNetwork.initialize()
	
	def initializeSolutions(self, Tsteps):
		"""This method is called by the initializeSolutionVectors method in the VaryingElastance-instance to initialize the vector containing the state variable,
		the state variable is stored so that it can be used in visualization to show the opening and closing of the valve """
		self.state = np.zeros(Tsteps)
	
	def computeB(self, A, n):
		""" Returns the turbulent resistance coefficient B, used in computing the pressure difference across the valve"""
		A_eff = self.effectiveOrificeArea(A,n)
		

		
		if A_eff == 0:
			B = None
		elif A/A_eff > 1e4:
			B = None
		else:
			B = 5*0.5*self.rho*(1/A_eff - 1/A)**2
#			B = 0.5*self.rho*(1/A_eff - 1/A)  #test
		return B

	def computeL(self, A, n):
		""" Returns the inertance coefficient L,  used in computing the pressure difference across the valve"""
		A_s = self.effectiveOrificeArea(A,n)	
		if A_s == 0:
			L = None
		elif A/A_s > 1e4:
			L = None
		else:
			L = 4*np.pi*self.rho*(1/A_s - 1/A)**0.5
		return L

	
	def LdivideB(self, A, n):
		A_s = self.effectiveOrificeArea(A,n)
		if A_s == 0:
			return None
		elif A/A_s > 1e4:
			return None
		else:
			return 4*np.pi*(1/A_s - 1/A)**(-1.5)

	def effectiveOrificeArea(self, A, n):
		""" Computes the effective orifice area (A_eff) by interpolating between maximum and minimum area."""
		A_max = self.M_st * A
		A_min = self.M_rg * A
		
		return (A_max - A_min) * self.state[n]  + A_min
	
	def A_max(self, A):
		return self.M_st * A
	
	def A_min(self, A):
		return self.M_rg * A
	
	def updateValveState(self, delta_p, n, dt):
		"""The opening state at timestep n+1 is computed by explicit numerical solution of the opening and closing rate equations (Mynard) """
		
		if delta_p > self.delta_p_open:
			"""The valve starts to open if the pressure difference is positive and above the opening threshold pressure"""
			if self.state[n] == 1.0:
				self.state[n+1] = 1.0
			else:
				self.state[n+1] = self.state[n] + (1 - self.state[n])*self.K_v_open*(delta_p - self.delta_p_open)*dt
				if self.state[n+1] > 1.0:
					self.state[n+1] = 1.0
					
		elif delta_p < -self.delta_p_close:
			""" The valve starts to close """
			if self.state[n] == 0.0:
				self.state[n+1] = 0.0
			else:
				self.state[n+1] = self.state[n] + self.state[n]*self.K_v_close*(delta_p - self.delta_p_close)*dt
				if self.state[n+1] < 0:
					self.state[n+1] = 0.0
		else:
			""" The case where -delta_p_close < delta_p < delta_p_open, the valve state stays unchanged """
			self.state[n+1] = self.state[n]
			
class VaryingElastance(BoundaryConditionType2):
	"""
	An implementation of a time-varying elastance model of the left ventricle,ing based on the modfied varying elastance equation including a source resistance K
	(as proposed by Shroff), and a parametrized time varying elastance function (as given by Stergiopulos).
	
	The general shape of the elastance function is given by three shape parameters. Various conditions of heart contractility are then created by scaling this 
	function using the parameters:
	T - Heart period
	Emax - Maximum elastance
	Emin - Minimum elastance
	Tpeak - Time to peak elastance
	
	The equation also requires the two constants:
	V0 - Volume axis intercept
	K - Source resistance
	
	NB! The source resistance (K) was introduced to the implementation as an experiment. It modifies the elastance curve depending on ventricular outflow,
	so that it bself.mitral = None # mitral valve
		self.aortic = None # aortic valve
		self.initializeValves() # intialize valvesecomes dependent on the afterload of the heart. Introducing the source resistance did produce a load dependence 
	(shown by curved isochrones in p-v loops), however the results are in no way to be trusted since the modified varying elastance equation was not intended
	to be used together with the specific varying elastance curve shape used here (Stergiopulos). A proper implementation of the source resistance 
	requires a different curve shape. The parameter K is therefore set to zero by default, but it should perhaps be removed from the code altogether??
	
	Currently only the return method "def funcPos0" has been implemented so that the boundary condition can only be put at the proximal end of a blood vessel.
	It is fairly straightforward to implement funcPos1 if necessary, this does however require a lot of duplicated code.   """
	def __init__(self):
		self.type = 2
		
		self.subiterations = 0
		
		self.omegaNew = np.empty((2))
		
		#Default parameters
		self.T = 1
		self.Emax  = 2.31 * 133.3e6 
		self.Emin = 0.06 * 133.3e6 
		self.Tpeak = 0.4
		
		self.V0 = 20e-6 
		
		self.K  =  0.0
		self.Rv = 0.02*133/(10**-6)
		
		"""Shape parameters"""
		self.alpha = 1.672
		self.n1 = 1.32
		self.n2 = 21.9
		self.R11 = None
		self.R12 = None
		self.R21 = None
		self.R22 = None
		self.DtW2 = None

		

		#n-1 values
#		self.aorticFlowPreviousTimestep = None
		
#		self.system = {'both open':np.array([0,1,2]), 'mitral open': np.array([0,1]), 'aortic open':np.array([1,2])} 
		
				
		self.cycleNumber = 0
		self.atriumPressure = 7.5 * 133.32 #Pressure in the atrium ## venouse pressure?!
		
#		self.x0 = np.array([0.0, 0.0, 0.0]) #Initial values for the iterative solver
		
		
		self.dQInOut = np.empty((2))
		
		
	def initializeSolutionVectors(self, Tsteps):
		"""Initializes some solution vectors storing pressure, flow and volume of the ventricle, as well as opening and closing state 
		
		NB! This method is not called from the class constructor, but is called externally by the initializeSolutionMatrices method in the solver,
		this is a bit messy, but was the easiest way to do it since the BC is initiated before the number of time steps is known.
		"""
		
		""" Initialize Solution Vectors """
		
		print """ Initialize Solution Vectors """
		
	
		self.pressure = np.zeros(Tsteps)
		self.volume = np.zeros(Tsteps)
		self.mitralQ = np.zeros(Tsteps)
		self.Elastance = np.zeros(Tsteps) #New
		self.Flow = np.zeros(Tsteps)
		self.DtFlow = np.zeros(Tsteps)
		self.deltaP=np.zeros(Tsteps)
		self.aortaP=np.zeros(Tsteps)
		
		""" Initial conditions in the ventricle"""
		self.pressure[0] = self.atriumPressure
		self.volume[0] = self.atriumPressure/self.E(0) + self.V0
				
		
	
	def __call__(self, _domegaField_, duPrescribed, R, L, n, dt, P, Q, A, Z1, Z2):
	
#		self.updateValves(P, n, dt)                     # Update the state of the mitral and aortic valve at timestep n + 1
		self.startNewCycleIfCriteriaIsMet(n, dt)
		self.funcPos0(_domegaField_, R, n, dt, P, Q, A)      # Compute the riemann variant going into the vessel save in omegaNew
		
		self.dQInOut = R[:][1]*self.omegaNew 
		# calculate du and return this!
		return np.dot(R,self.omegaNew),self.dQInOut





	def getCycleTime(self, n, dt):
		return n*dt - self.T*self.cycleNumber
	
	def startNewCycleIfCriteriaIsMet(self, n, dt):
		if self.getCycleTime(n, dt) > self.T:
			self.cycleNumber += 1

	def funcPos0(self, _domega, R, n, dt, Pn, Qn, A):
		
		# Qn1 == value at old time step
		# change to self.aorticPressurePreviousTimestep ...
		
#		Qn1 = self.aorticFlowPreviousTimestep
		
		L=np.linalg.inv(R)
		L11, L12, L21, L22 = L[0][0],L[0][1],L[1][0],L[1][1]
		omegaprevious_ = L11*Pn + L12*Qn
		r11,r12,r21,r22 =  R[0][0],R[0][1],R[1][0],R[1][1]
	#	deltatdiff = 0.00001
#		mitrQn = self.mitralQ[n]
#		mitrQn1 = self.mitralQ[n-1]
		venoP = self.atriumPressure
		t = self.getCycleTime(n+1, dt)
		t2 = self.getCycleTime(n, dt)
# 		ttemp = t-dt
		E = self.E(t)
		e2 = self.E(t)
	#	dE= (self.E(t+deltatdiff) -E)/deltatdiff
		Vn = self.volume[n]
		self.R11 = r11
		self.R12 = r12
		self.R21 = r21
		self.R22 = r22
		self.DtW2 = _domega/dt

		
		self.Elastance[n+1]=E/133.3e6
		self.Flow[n]=Qn*1e6
		self.aortaP[n]=Pn
#		self.DtFlow[n]=(Qn-Qnold)/dt
		ventrPn = self.pressure[n]
		
		print "pn is: ", Pn/133
		print "ventrPn is: ", ventrPn/133
		print "ventricle pressure should be:", E*(Vn-self.V0)/133
		print "ventrPn-Pn is:", ventrPn-Pn
		print "volume is", Vn*10**6
		print "Qn is", Qn*10**6
		print "start volume was", self.volume[0]*10**6
		print "t is:",t
		print "E is,", E/133.3e6
		print "e is", e2/133.3e6
		
		
		def diastole(u,t):
		
			"""Differential equations during diastole u[0]=V, u[1]=Pv, return dV/dt and dP/dt"""
			
			deltatdiff = 0.00001
			E= self.E(t)
			dE= (self.E(t+deltatdiff) -E)/deltatdiff
			DV = (venoP-u[1])/self.Rv
			DP = (dE*u[0]-self.V0)+E*(venoP-u[1])/self.Rv
			return [DV,DP]
	
	
		def systole(u,t):
		
			deltatdiff = 0.00001
			Etemp= self.E(t)
			dE= (self.E(t+deltatdiff) -Etemp)/deltatdiff
			
			DV = -u[2]
			DP = (dE*(u[0]-self.V0) - Etemp*u[2] )
	#		#DW1 = ((dE*(u[0]-self.V0) - Etemp*u[3] - self.R12*self.DtW2)/self.R11)
			DQ = (self.R21/self.R11)*(dE*(u[0]-self.V0) - Etemp*u[2] ) +self.DtW2*(self.R22-(self.R12*self.R21)/self.R11)
			return[DV,DP,DQ]
			
		solverSys = OD.RungeKutta4(systole) #Runga Kutta solver for ejection phase, Systole

		solverdias = OD.RungeKutta4(diastole) #Runga Kutta solver for diastole
			
			
		if (venoP>ventrPn):
			solverdias.set_initial_condition([Vn,ventrPn])
			if (t2)<0:
				t_pointsd = np.linspace(0,dt,2)
				
			else:
				t_pointsd = np.linspace(t2,t2+dt,2)
			print t_pointsd
			print self.cycleNumber

			udiastole,td = solverdias.solve(t_pointsd)
			udiastole = udiastole[1]
			V = udiastole[0]
			Pv = udiastole[1]
			self.volume[n+1]=V
			self.pressure[n+1]=Pv
			if Qn ==0:
				domega_=_domega
			else:
				domega_ = (-0.5*Qn - r22*_domega)/r21

			print "diastole"
			
		elif (Qn>=0 and ventrPn-Pn>(-0.5*133.32)):
			solverSys.set_initial_condition([Vn,ventrPn,Qn])
			t_pointss = np.linspace(t2,t2+dt,2)
			uSystole,ts = solverSys.solve(t_pointss)
			uSystole = uSystole[1]
			print t_pointss
			print t_pointss.shape
			V = uSystole[0]
			P = uSystole[1]
			Q =uSystole[2]
			omeganew_ = L11*P + L12*Q
			print "systole"
			
			self.volume[n+1]=V
			self.pressure[n+1]=P
			domega_= omeganew_ - omegaprevious_
			
		else:
			self.volume[n+1]=Vn
			self.pressure[n+1]=E*(Vn-self.V0)
			
			domega_=_domega
			print "iso"
			

		self.omegaNew[0] = domega_
		self.omegaNew[1] = _domega


	def funcPos1(self, _domega, R, L, n, dt, P, Q, A):
		pass



	def E(self, t):
		"""Computes the value of the elastance at time t, according to the shape parameters given by Stergiopolus and scaled
		   according to Tpeak, T, Emax and Emin. """
		a1 = 0.708*self.Tpeak
		a2 = 1.677*a1
		
		n1, n2 = self.n1, self.n2
		shapeFunction1 = (t/(a1))**n1/(1+(t/(a1))**n1)
		shapeFunction2 = (1 + (t/(a2))**n2)**(-1)	
		return (self.Emax-self.Emin)*self.alpha*shapeFunction1*shapeFunction2 + self.Emin
# 	def diastole(self,u,t):
# 		
# 		"""Differential equations during diastole u[0]=V, u[1]=Pv, return dV/dt and dP/dt"""
# 		
# 		deltatdiff = 0.00001
# 		E= self.E(t)
# 		dE= (self.E(t+deltatdiff) -E)/deltatdiff
# 		DV = (self.atriumPressure-u[1])/self.Rv
# 		DP = (dE*u[0]-self.V0)+E*(self.atriumPressure-u[1])/self.Rv
# 		return [DV,DP]
# 	
# 	
# 	def systole(self,u,t):
# 		
# 		deltatdiff = 0.00001
# 		E= self.E(t)
# 		dE= (self.E(t+deltatdiff) -E)/deltatdiff
# 		
# 		DV = -u[3]
# 		DP = (dE*(u[0]-self.V0) - E*u[3] )
# 		DW1 = ((dE*(u[0]-self.V0) - E*u[3] - self.R12*self.DtW2)/self.R11)
# 		DQ = (self.R21/self.R11)*(dE*(u[0]-self.V0) - E*u[3] ) +self.DtW2*(self.R22-(self.R12*self.R21)/self.R11)
# 		return[DV,DP,DW1,DQ]


