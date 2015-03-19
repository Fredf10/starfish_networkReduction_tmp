########################################################################################
#                            Vascular Polynomial Chaos 0.3
########################################################################################
## 
# created by Vinzenz Eck vinzenz.eck@mytum.de
# uses polynomial Chaos toolbox from Jonathan Feinberg, Simula Center Oslo
##

import sys,os
cur = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append('/'.join([cur,'VascularPolynomialChaosLib']))
from classVpcConfiguration import VpcConfiguration
from classVpcDistribution import VpcDistribution

sys.path.append('/'.join([cur,'UtilityLib']))
import moduleStartUp as mSU
import moduleXML

import chaospy as cp
import pprint



def vascularPolyChaos():
    '''
    Perform vascular polynomial chaos 
    or MonteCarlo analysis for STARFiSh simulation case
    # steps
    # 1. load vpc case and configuration 
    
    # 2. create distributions
    
    # 3. add correlation if existent
    
    # 4. create samples
    
    # 5. create Orthogonal polynomials
    
    # 6. evaluate model / on local machine or on server
    
    # 7. postprocess evaluated data, peak finding etc
    
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
    vascularNetwork = moduleXML.loadNetworkFromXML(''.join([networkName,'_vpc']), dataNumber)
    
    # 2. create distributions
    vpcDistributions = VpcDistribution()
    # create distributions defined in vascularNetwork or try to load
    if vpcConfiguration.createDistributions == True:
        vpcDistributions.createDistributions(vascularNetwork)
    else:
        vpcDistributions.loadDistributions()
        
    # 3. add correlation if existent
    
    # 4. create samples
        
    # 5. create Orthogonal polynomials
    
    # 6. evaluate model / on local machine or on server
    
    # 7. postprocess evaluated data, peak finding etc
    
    # 8. calculate polynomial chaos expansion
    
    # 9. uncertainty quantfication, sensitivity analysis


if __name__ == '__main__':
    vascularPolyChaos()