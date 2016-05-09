from VascularPolynomialChaosLib.testBaseClass import TestBaseClass 
import numpy as np
from scipy import optimize

# class MeasurmentRoutine(TestBaseClass):
#        
#     def __init__(self):
#         '''
#         Class for adapting the vascularn etwork to a patient specific condition with measurements
#         of cardiovascular parameter defined in this class.
#         Includes also algorithm to adapt the defined network parameter to patient specific condtions
#         given by the measurements
#         '''
#     
#     def adaptationToPatientSpecificCondition(self, vascularNetwork):
#         '''
#         This method starts the adapation of the given vessels to the patient specific conditions
#         defined by the measurments
#         
#         Args: vascularNetwork (class): instance of the vascularNetwork to adapt.
#         '''
#         pass
#     
    
class MeasurmentRoutine(TestBaseClass):
   
    externVariables      = {'vesselId': TestBaseClass.ExtValue([int],    unit = ''),
                            'diastolicPressure'      : TestBaseClass.ExtValue([float],  unit = 'Pa'),
                            'diastolicArea'      : TestBaseClass.ExtValue([float], unit = 'm2' ),
                            'waveSpeed_ff_carotidFemoral'    : TestBaseClass.ExtValue([float],  unit = 'm s-1' ),
                            #'carotidArteryId': TestBaseClass.ExtValue([int],    unit = ''),
                            #'gemoralArteryId': TestBaseClass.ExtValue([int],    unit = ''),
                            'vesselIdsAorticToCarotid' : TestBaseClass.ExtValue([int],   multiVar=True, unit = ''),
                            'vesselIdsAorticToFemoral' : TestBaseClass.ExtValue([int],   multiVar=True, unit = '')
                            }
    
    externXmlAttributes  = []
    externXmlElements    = ['vesselId',
                            'diastolicPressure',
                            'diastolicArea',
                            'waveSpeed_ff_carotidFemoral',
                            #'carotidArteryId',
                            #'gemoralArteryId',
                            'vesselIdsAorticToCarotid',
                            'vesselIdsAorticToFemoral']
    
    def __init__(self):
        '''
        This class has three measurments defined: diastolic pressure, diastolic area at a given location,
        and the diastolic pulse wave velocity from carotid to femoral artery.
        
        Based on these three measurements, the wall model, nameley As,Ps, and the "stiffness coefficients"
        are adapted such that the vascular netowrk reflects the measured values.
        
        The wall models which are currently supported are Laplace, Hayashi, Reymond.
        
        Attributes:
            vesselId (int): vessel id of the measurement side 
            diastolicPressure (float): diastolic pressure measured at the measurement location
            diastolicArea (float): diastolic area measured at the measurement location
            waveSpeed_ff_carotidFemoral (float): foot-to-foot wave speed measured from caritid to femoral artery
            carotidArteryId (int): vesselId of the carotid artery
            gemoralArteryId (int): vesselId of the femoral artery
            aorticToCarotidVesselIds (list) of (int): vesselIds from aortic arch to carotid artery used to calculate the distance for c_ff
            aorticToFemoralVesselIds (list) of (int):  vesselIds from aortic arch to femoral artery used to calculate the distance for c_ff
            
        '''
        self.vesselId                    = None
        self.diastolicPressure           = None
        self.diastolicArea               = None
        self.waveSpeed_ff_carotidFemoral = None
        #self.carotidArteryId             = None
        #self.gemoralArteryId             = None
        self.vesselIdsAorticToCarotid    = None
        self.vesselIdsAorticToFemoral    = None
        
        
    def adaptationToPatientSpecificCondition(self, vascularNetwork):
        '''
        This method starts the adapation of the given vessels to the patient specific conditions
        defined by the measurments.
        Algorithm:
        
        1. adapt Ps in all vessels relative to vessel where diastolicPressure is measurement 
        and adapt As in all vessels relative to vessel where diastolicArea is measurement 
        2. adapt "wall stiffness coefficients" such that measured waveSpeed_ff_carotidFemoral matches with
        the intial c_ff given by c(Ps,As) in the carotid to femoral artery.
        
        
        Args: vascularNetwork (class): instance of the vascularNetwork to adapt.
        '''
        # 1. adapt Ps  and 2. adapt As
        # check if vesselId is valid Id
        if self.vesselId not in vascularNetwork.vessels.keys():
            raise ValueError('MeasurmentRoutine_Pd_Ad_pvwff.adaptationToPatientSpecificCondition(): the measurement side defined vesselId {} is not represented as a vessel in the netowrk'.format(self.vesselId))
        # 1.1 find out ratios alpha_p and alpha_a Ps and As is changed in the measurement vessel
        measurmentVessel = vascularNetwork.vessels[self.vesselId]
        alpha_p  = self.diastolicPressure / measurmentVessel.compliance.Ps
        alpha_a  = self.diastolicArea / np.mean(measurmentVessel.compliance.As)        
        # 1.2 loop through vessels and adjust Ps, and As with alpha_1 and alpha_2
        for vessel in vascularNetwork.vessels.itervalues():
            vessel.compliance.Ps = vessel.compliance.Ps*alpha_p
            vessel.compliance.As = vessel.compliance.As*alpha_a
            
        # 2. adapt wall stiffness coefficients
        # check if vesselIds are valid Ids
        existingVessels = vascularNetwork.vessels.keys()
        for vesselId in self.vesselIdsAorticToCarotid+self.vesselIdsAorticToFemoral:
            if vesselId not in existingVessels:
                raise ValueError('MeasurmentRoutine_Pd_Ad_pvwff.adaptationToPatientSpecificCondition(): vesselId {} connecting carotid and femoral artery is not represented as a vessel in the netowrk'.format(vesselId))
        
        # 2.1 evaluate alpha_c iteratively such that c_ff from network matches with measurment waveSpeed_ff_carotidFemoral
        alpha_c = 1.0
        args = [self.waveSpeed_ff_carotidFemoral,vascularNetwork.vessels]
        #args = [c_ff_analytic,wallModelType,dataDict]
        sol = optimize.root(self.estimateCff, alpha_c, args)
        alpha_c = sol.x[0]
        
        print 
        print "MeasurmentRoutine"
        print 
        print "estimated alpha_p:", alpha_p
        print "estimated alpha_a:", alpha_a
        print "estimated alpha_c:", alpha_c, type(alpha_c)
        print         
        
        # 2.2 loop through vessels and adjust all coefficients with alpha_c
        for vessel in vascularNetwork.vessels.itervalues():
            vessel.compliance.adaptMaterialCoefficient(alpha_c)
            
    def estimateCff(self, alpha_c, args):
        '''
        function to estimate Carotid Femoral wave speed for a given alpha_c and wallmodel
        method of wave speed evaluation as the analytic way reported in Vardoulis 2013
        
        Note: 
        to use a wall model it need an update function .adaptMaterialCoefficient(alpha_c)
        which updates the wall model coefficient "beta" such that the equation
        wave speed at reference pressure:
        c1 = alpha_c c_wallModel(Ps) 
        and 
        
        wallmodel.adaptMaterialCoefficient(alpha_c)
        c2 = c_wallModel(Ps)
        are the same!
        '''
        c_ff,vessels = args
            
        distanceAorticBifFemoralRight = 0
        distanceCarotidAorticBif = 0
        
        ## evaluate wave speed analytic as reported in Vardoulis 2013
        timeCarotid_analytic = 0
        for i in self.vesselIdsAorticToCarotid:
            # use the mean c
            vessel = vessels[i]
            pressure = np.ones(vessel.N) * vessel.compliance.Ps
            c_i = alpha_c*np.mean(vessel.waveSpeed(vessel.A(pressure),vessel.C(pressure)))
            x_i = vessel.length
            
            timeCarotid_analytic = timeCarotid_analytic + x_i/c_i
            distanceCarotidAorticBif = distanceCarotidAorticBif + x_i
            
        timeFemoral_analytic = 0
        for i in self.vesselIdsAorticToFemoral:
            # use the mean c
            vessel = vessels[i]
            pressure = np.ones(vessel.N) * vessel.compliance.Ps
            c_i = alpha_c*np.mean(vessel.waveSpeed(vessel.A(pressure),vessel.C(pressure)))
            x_i = vessel.length
            
            timeFemoral_analytic = timeFemoral_analytic + x_i/c_i
            distanceAorticBifFemoralRight = distanceAorticBifFemoralRight+x_i
            
        c_ff_model = (distanceAorticBifFemoralRight-distanceCarotidAorticBif)/(timeFemoral_analytic-timeCarotid_analytic)

        return abs(c_ff_model-c_ff)
                
                
    