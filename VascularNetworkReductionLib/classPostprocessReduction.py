'''
Created on Sep 28, 2016

@author: fredrik
'''


import sys, os
import h5py
import numpy as np
import multiprocessing
cur = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur + '/../')



import moduleFilePathHandlerVNR as mFPHVNR
import UtilityLib.moduleXML as mXML


class PostprocessReduction:
    
    def __init__(self, batchDataList=None):
        
        """batchData = {'networkName': newNetworkName, 
                        'dataNumber': dataNumber, 
                        'networkXmlFileLoad': newNetworkXmlFile,
                        'networkXmlFileSave': solutionFileXML,
                        'pathSolutionDataFilename' : solutionFilehd5}"""
        
        self.batchData = batchDataList
        
    def postprocessSingle(self, batchData):

        #simulationIndex          = batchData['simulationIndex']
        networkName              = batchData['networkName']
        dataNumber               = batchData['dataNumber']
        networkXmlFileLoad       = batchData['networkXmlFileLoad']
        networkXmlFileSave       = batchData['networkXmlFileSave']
        pathSolutionDataFilename = batchData['pathSolutionDataFilename']
        
        self.processVascularNetwork(networkName, dataNumber, networkXmlFileSave, pathSolutionDataFilename)
        

    def postprocessMulti(self, batchDataList, numberWorkers = None, quiet = False):
        '''
        Run a set of simulations on one core without multiprocessing
        
        Input:
            batchDataList <list> := with data for each batch job [batchData1, batchData .. ]
                batchData <dict> := dict with {simulationIndex: , networkName: , dataNumber: , networkXmlFile: , pathSolutionDataFilename: }
            
        '''
        if numberWorkers == None: numberWorkers = multiprocessing.cpu_count()

        print '------Multiprocessing Batch Job Postprocess------'
        print 'numberWorkers:   {}'.format(numberWorkers)
        print 'numberOfEval.:   {}'.format(len(batchDataList))
        pool = multiprocessing.Pool(numberWorkers)
        results = pool.imap(self.postprocessSingle, batchDataList)
        pool.close() 
        while (True):
            completed = results._index
            if (completed == len(batchDataList)): break
        pool.join()

        
    def processVascularNetwork(self, networkName, dataNumber, networkXmlFile, pathSolutionDataFilename):
        
        vascularNetwork = mXML.loadNetworkFromXML(networkName, dataNumber = dataNumber, 
                                                  networkXmlFile = networkXmlFile, pathSolutionDataFilename = pathSolutionDataFilename)
        
        vascularNetwork.initialize(initializeForSimulation = True)

        if pathSolutionDataFilename == None:
            pathSolutionDataFilename = mFPHVNR.getFilePath('solutionFile', self.name, self.dataNumber, 'read')
        # TODO, what if this fails? do we know?
        solutionDataFile = h5py.File(pathSolutionDataFilename, "r")
        pathSolutionDataFilenameNew = pathSolutionDataFilename + "tmp"
        newSolutionDataFile = h5py.File(pathSolutionDataFilenameNew, "w")

        for groupName, group in solutionDataFile.iteritems():
            if groupName == 'VascularNetwork':
                dsetGroup = newSolutionDataFile.create_group('VascularNetwork')
                simulationTime = group['simulationTime'][:]
                dsetGroup.create_dataset('simulationTime', data=simulationTime)
            elif groupName == 'vessels': # or '-' in groupName: # '-' is loads older hdf5 data files
                newVesselDataGroup = newSolutionDataFile.create_group('vessels')
                for vesselh5pyName, vesselh5pyGroup in group.iteritems():
                    vesselId = int(vesselh5pyName.split(' - ')[-1])
                    vessel = vascularNetwork.vessels[vesselId]
                    
                    self.processVesselData(vessel, newVesselDataGroup, vesselh5pyName, vesselh5pyGroup, simulationTime)
        
        solutionDataFile.close()
        newSolutionDataFile.close()
        os.remove(pathSolutionDataFilename)
        os.rename(pathSolutionDataFilenameNew, pathSolutionDataFilename)
        
        del vascularNetwork
        
    def processVesselData(self, vessel, newVesselDataGroup, vesselh5pyName, vesselh5pyGroup, simulationTime, Npoints=3):
        """ method that reduces number of nodes to Npoints, and also create datasets of wavesplit for each vessel in hdf file"""
        
        newVesselh5pyGroup = newVesselDataGroup.create_group(vesselh5pyName)
        
        Psol = vesselh5pyGroup["Psol"]
        Asol = vesselh5pyGroup["Asol"]
        Qsol = vesselh5pyGroup["Qsol"]
         
        l = vessel.length
        N = vessel.N
        
        zOld = np.linspace(0, l, N)
        zNew = np.linspace(0, l, Npoints)
        
        PsolNew = self.myInterpolater(zNew, zOld, Psol)
        QsolNew = self.myInterpolater(zNew, zOld, Qsol)
        AsolNew = self.myInterpolater(zNew, zOld, Asol)
        CsolNew = self.myInterpolaterC(zNew, zOld, Psol, vessel)
        #print CsolNew
        gamma = vessel.gamma
        
        PsolForward = np.zeros((len(simulationTime), Npoints))
        PsolBackward = np.zeros((len(simulationTime), Npoints))
        QsolForward = np.zeros((len(simulationTime), Npoints))
        QsolBackward = np.zeros((len(simulationTime), Npoints))
        
        for n in range(Npoints):
            P = PsolNew[:, n]
            Q = QsolNew[:, n]
            A = AsolNew[:, n]
            C = CsolNew[:, n]
            U = Q/A
            c = vessel.waveSpeed(A, C)
            
            P_forward, P_backward, Q_forward, Q_backward = self.waveSplit(P, Q, A, U, c, C, gamma)
            
            PsolForward[:, n] = P_forward
            PsolBackward[:, n] = P_backward

            QsolForward[:, n] = Q_forward
            QsolBackward[:, n] = Q_backward
            
        newVesselh5pyGroup.create_dataset('Psol', data=PsolNew)
        newVesselh5pyGroup.create_dataset('Qsol', data=QsolNew)
        newVesselh5pyGroup.create_dataset('Asol', data=AsolNew)
        
        newVesselh5pyGroup.create_dataset('Psol_f', data=PsolForward)
        newVesselh5pyGroup.create_dataset('Psol_b', data=PsolBackward)
        
        newVesselh5pyGroup.create_dataset('Qsol_f', data=QsolForward)
        newVesselh5pyGroup.create_dataset('Qsol_b', data=QsolBackward)
        
        
    def waveSplit(self, P, Q, A, U, c, C, gamma):
        """ nonlinear wavesplit using omega, R and L matrices."""
        alpha = (gamma + 2.)/(gamma + 1.)
        
        lambda1 = alpha*U + np.sqrt(c**2 + alpha*(alpha - 1)*U**2) 
        lambda2 = alpha*U - np.sqrt(c**2 + alpha*(alpha - 1)*U**2)
        
        Z1 = 1./(lambda1*C)
        Z2 = -1./(lambda2*C)
        
        deltaP = P[1:]- P[:-1]
        deltaQ = Q[1:]- Q[:-1]
        
        deltaW1 = deltaP + Z2[:-1]*deltaQ
        deltaW2 = deltaP - Z1[:-1]*deltaQ
        
        deltaP_forward = Z1[:-1]*deltaW1/(Z1[:-1] + Z2[:-1])
        deltaP_backward = Z2[:-1]*deltaW2/(Z1[:-1] + Z2[:-1])
        

        deltaQ_forward = deltaW1/(Z1[:-1] + Z2[:-1])
        deltaQ_backward = -deltaW2/(Z1[:-1] + Z2[:-1])
        
        P_forward = np.cumsum(np.append(P[0], deltaP_forward))
        P_backward = np.cumsum(np.append(P[0], deltaP_backward))
        
        Q_forward = np.cumsum(np.append(Q[0], deltaQ_forward))
        Q_backward = np.cumsum(np.append(Q[0], deltaQ_backward))
        
        return P_forward, P_backward, Q_forward, Q_backward
        
        
    def myInterpolater(self, zNew, zOld, Xsol):
        tpoint, zpoint = np.shape(Xsol)
        
        XsolNew = np.zeros((tpoint, len(zNew)))
        for i in range(tpoint):
            
            f = Xsol[i, :].copy()
            
            fnew = np.interp(zNew, zOld, f)
            XsolNew[i, :]= fnew
    
        return XsolNew
    
    def myInterpolaterC(self, zNew, zOld, Psol, vessel):
        """ methods that """
        
        tpoint, zpoint = np.shape(Psol)
        
        CsolNew = np.zeros((tpoint, len(zNew)))
        
        for i in range(tpoint):
            P = Psol[i, :]
            C = vessel.C(P)
             
            fnew = np.interp(zNew, zOld, C)
            CsolNew[i, :]= fnew
        
        return CsolNew
    

        
    
