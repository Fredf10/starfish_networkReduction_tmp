#!/usr/bin/env python
# -*- coding: utf-8 -*-

from testBaseClass import TestBaseClass 

from classLocationOfInterestManager import LocationOfInterestManager
from classVpcConfiguration import VpcConfiguration

class ConfigurationUQSA(TestBaseClass):
    
    externVariables      = {'vpcConfiguration' : TestBaseClass.ExtObject({'VpcConfiguration':VpcConfiguration}),
                            'locationOfInterestManager' : TestBaseClass.ExtObject({'LocationOfInterestManager':LocationOfInterestManager}),
                           } 
    externXmlAttributes  = []
    externXmlElements    = ['vpcConfiguration',
                            'locationOfInterestManager']
    
    
    def __init__(self):
        
        self.vpcConfiguration          = None
        self.locationOfInterestManager = None
        
    
        
        
    
        
    
        
        
        
        
        
        
        
        