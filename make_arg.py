'''
Created on Sep 28, 2016

@author: fredrik
'''

### Script thant launches reductionCase/batch
import sys,os
cur = os.path.dirname(os.path.realpath('__file__'))


import starfish.VascularNetworkReductionLib.classReductionCase as cRedCase

convergenceXmlCaseFile = sys.argv[1]

CPUTimeFile = convergenceXmlCaseFile[:-4] + "_CPU.txt"

redCase = cRedCase.ReductionCase(CPUTimeFile=CPUTimeFile)

redCase.loadXMLFile(convergenceXmlCaseFile)

redCase.loadBatchDatalist()
redCase.loadoptimizeParamsDataFile()


redCase.createEvaluationCaseFiles()
