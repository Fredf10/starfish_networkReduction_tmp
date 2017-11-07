'''
Created on Sep 28, 2016

@author: fredrik
'''


import sys

import starfish.NetworkTests.classConvergenceCase as cConvergenceCase

convergenceXmlCaseFile = sys.argv[1]


cCCase = cConvergenceCase.ConvergenceCase()

cCCase.loadXMLFile(convergenceXmlCaseFile)
cCCase.createEvaluationCaseFiles(convergenceXmlCaseFile.replace('.xml', '.p'))