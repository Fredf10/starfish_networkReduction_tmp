import numpy as np 


class System(object):
    
    def __init__(self,vessel,simplifyEigenvalues,riemannInvariantUnitBase,n,dt):
        '''
        Constructor of System with the SystemEqations
        
        Input:
            vessel <classVessel>              : vessel for which the system Equations should be build
            simplifyEigenvalues <bool>        : TRUE  == simplified eigenvalues lambda1/2 = +- c
                                                FALSE == lambda1/2 = dlt*v +- sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            riemannInvariantUnitBase <String> : 'Pressure' == calculate R,L based on riemann invariants with unit
                                                              of Pressure
                                                'Flow'     == calculate R,L based on riemann invariants with unit
                                                              of Flow
        Variables (most important):
           
            m12 <np.array> : 12 value of system Matrix M
            m21 <np.array> : 21 value of system Matrix M
            m22 <np.array> : 22 value of system Matrix M
            b2  <np.array> : 2 value of righthandside vector B
            LAMBDA   <nested list> = [<np.array>,<np.array>] : Eigenvalues at each boundary node 
            R        <nested list> = [<np.array>,<np.array>] : Right-Eigenvector matrix at each boundary node 
            L        <nested list> = [<np.array>,<np.array>] : Left-Eigenvector matrix at each boundary node 
            Z        <nested list> = [<np.array>,<np.array>] : Impedances (Z1, Z2) (forward/backward) at each boundary node 
        
        '''        
        self.name = ' '.join(['system equations', str(vessel.Id)])
        
        # compliance properties
        self.C    = vessel.C
        self.c    = vessel.c
        
        self.C_nID = vessel.C_nID
        
        #Fluid properties
        self.my     = vessel.my
        self.rho    = vessel.rho
        self.gamma  = vessel.gamma
        
        # vessel net gravity
        self.netGravity = vessel.netGravity
        
        # properties for Riemann invariant calculation
        self.z = vessel.z
        self.dt = dt 
        self.n = n
                
        # system matrices
        self.m12    = None
        self.m21    = None
        self.m22    = None
        self.M      = None
        self.b2     = None
        self.B      = None
        self.LAMBDA = [np.ones(2),np.ones(2)]  # list with 2 objects one for each boundary
        self.R      = [np.ones((2,2)),np.ones((2,2))]  # list with 2 objects one for each boundary
        self.L      = [np.ones((2,2)),np.ones((2,2))]  # list with 2 objects one for each boundary
        self.Z      = [np.ones(2),np.ones(2)] # [forward Z, backward Z ] list with 2 objects one for each boundary
        self.domega = [0.0,0.0] # [omega2(position0),omega1(position-1)] list with 1 domega at each boundary
        self.du     = np.empty(2)
        
        #calculate dlt
        self.dlt = (self.gamma+2.0)/(self.gamma+1.0)
        
        # eigen values
        if simplifyEigenvalues   == True:
            if riemannInvariantUnitBase == 'Flow':
                self.updateLARL = self.updateLARLSys0InvariantFlow
            elif riemannInvariantUnitBase == 'Pressure':
                self.updateLARL = self.updateLARLSys0InvariantPressure
            
        elif simplifyEigenvalues == False:
            if riemannInvariantUnitBase == 'Flow':
                self.updateLARL = self.updateLARLSys1InvariantFlow
            elif riemannInvariantUnitBase == 'Pressure':
                self.updateLARL = self.updateLARLSys1InvariantPressure
    
        self.Re = 3000.0
                    
    def updateSystem(self,P,Q,A,pi=np.pi):
        '''
        Update all the system-equation and matrices
        
        Input:
            P <np.array> : current pressure values of the vessel
            Q <np.array> : current flow values of the vessel
            A <np.array> : current area values of the vessel
        
        '''
        n = self.n[0]
        # calculate needed values
        C = self.C(P) 
        c = self.c(A,C)
        v = Q/A
        
# #         Re = np.max(np.sqrt(A/np.pi)*2.0*v/self.my*self.rho)    
# #         if Re > self.Re: # > 3000.0
# #             print "WARNING approaching high reynoldsnumber for vessel", self.name, Re
# #             self.Re = Re
#             
        dlt = self.dlt
        
        m12 = (1. / C)[1:-1]
        m21 = (C * (c**2 - dlt * v**2))[1:-1] # vinz expression #||# (A/self.rho)[1:-1] # paul expression  # (C*(c**2-self.dlt*v**2))[1:-1] # vinz expression
        m22 = (2 * dlt * v)[1:-1]  #vinz expression #||# self.dlt #paul expression #(2*self.dlt*v)[1:-1]
        # set up B matrix
        b2 = (-2.0 / self.rho * pi * v * (self.gamma+2) * self.my + A * self.netGravity[n])[1:-1] #Correction with net gravity force
        
        ## set up LAMBDA, L and R matrices
        #self.updateLARL(P,Q,A,Ct=C,ct=c)
        
        return m12,m21,m22,b2
            
    def updateLARLSys0InvariantFlow(self,P,Q,A,idArray=[0,-1],Ct=None,ct=None):
        '''
        Update LAMBDA,R,L,Z of the system equations
        
        Special terms:
            sys0 :          uses simplified eigenvalues lambda = c (wave speed) 
                            as result Z1=Z2=Zc
            invariantFlow:  calculates R,L using riemannInvariants based on unit FLow
                            as in Paper of Paul Roger Leinan but with R2*(-1) correction to
                            have additive invariants for pressure
        
        Input:
            P                 <np.array> : current pressure values of the vessel
            Q                 <np.array> : current flow values of the vessel
            A                 <np.array> : current area values of the vessel
            idArray = [0,-1]  <list>     : define which boundary should be updated (default = both)
            update  = 'all'   <string>   : set to 'L' if only L and not R should be updated
            Ct = None         <np.array> : Compliance C if avialiable otherwise it will be calculated
            ct = None         <np.array> : waveSpeed  c if avialiable otherwise it will be calculated
        '''
        n = self.n[0]
        for Id in idArray:
            
            
            
            Aid = A[Id]       
            
            if Ct == None: C = self.C_nID(P,Id)
            else:  C = Ct[Id]
            if ct == None: c = self.c(A[Id],C)
            else:  c = ct[Id]
                        
            # Z1 = Z2 = Zc
            Zc =  1.0/(C*c)
            
            self.Z[Id][0] = Zc
            self.Z[Id][1] = Zc
            
            self.LAMBDA[Id][0] = c
            self.LAMBDA[Id][1] = -c
                  
            ## left eigenvalue matrix
            # riemannInvariants with unit flow    
            self.L[Id][0][0] =  0.5/(Zc)
            self.L[Id][0][1] =  0.5
            self.L[Id][1][0] =  0.5/(Zc)
            self.L[Id][1][1] = -0.5

            ## right eigenvalue matrix
            # riemannInvariants with unit flow 
            self.R[Id][0][0] =  Zc
            self.R[Id][0][1] =  Zc
            self.R[Id][1][1] = -1.0
            
            L = self.L[Id]
            R = self.R[Id]
            ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
                            
            ### check consistency: calculate sum(R*L) == sum(I) == 2.0 
            errorIdentity = abs((L[0][0]+L[1][0])*(R[0][0]+R[0][1])+(L[0][1]+L[1][1])*(R[1][0]+R[1][1])-2.0)
            if errorIdentity > 5.e-16:
                print "WARNING: SystemEquations, inverse of L and R differ, error {} > 5.e-16".format(errorIdentity)
            
            # calculate omegas
            du = self.du
            
            zi = self.z[Id] - [l2,l1][Id] * self.dt
                        
            du[0]  = np.interp(zi,self.z,P) - P[Id]
            du[-1] = np.interp(zi,self.z,Q) - Q[Id]
            
            self.domega[Id] = np.dot( self.L[Id][1+Id], du) + self.dt * self.L[Id][1+Id][1] * Aid * self.netGravity[n] 
                
             
    def updateLARLSys0InvariantPressure(self,P,Q,A,idArray=[0,-1],Ct=None,ct=None):
        '''
        Update LAMBDA,R,L,Z of the system equations
        
        Special terms:
            sys0 :          uses simplified eigenvalues lambda = c (wave speed)
                            as result Z1=Z2=Zc
            invariantFlow:  calculates R,L using riemannInvariants based on unit Pressure
                            as derieved by Leif Rune Hellevik
        
        Input:
            P                 <np.array> : current pressure values of the vessel
            Q                 <np.array> : current flow values of the vessel
            A                 <np.array> : current area values of the vessel
            idArray = [0,-1]  <list>     : define which boundary should be updated (default = both)
            update  = 'all'   <string>   : set to 'L' if only L and not R should be updated
            Ct = None         <np.array> : Compliance C if avialiable otherwise it will be calculated
            ct = None         <np.array> : waveSpeed  c if avialiable otherwise it will be calculated
        '''
        n = self.n[0]
        for Id in idArray: 
            
            Aid = A[Id]      
            
            if Ct == None: C = self.C_nID(P,Id)
            else:  C = Ct[Id]
            if ct == None: c = self.c(A[Id],C)
            else:  c = ct[Id]
            
            self.LAMBDA[Id][0] = c
            self.LAMBDA[Id][1] = -c
            
            # Z1 = Z2 = Zc
            Zc =  1.0/(C*c)
                        
            self.Z[Id][0] = Zc
            self.Z[Id][1] = Zc
            
            ## left eigenvalue matrix
            # riemannInvariants with unit pressure
            self.L[Id][0][1] =  Zc
            self.L[Id][1][1] = -Zc

            ## right eigenvalue matrix
            # riemannInvariants with unit pressure
            self.R[Id][0][0] =  0.5
            self.R[Id][0][1] =  0.5
            self.R[Id][1][0] =  0.5/(Zc)
            self.R[Id][1][1] = -0.5/(Zc)

            L = self.L[Id]
            R = self.R[Id]
            ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
                            
            ### check consistency: calculate sum(R*L) == sum(I) == 2.0 
            errorIdentity = abs((L[0][0]+L[1][0])*(R[0][0]+R[0][1])+(L[0][1]+L[1][1])*(R[1][0]+R[1][1])-2.0)
            if errorIdentity > 5.e-16:
                print "WARNING: SystemEquations, inverse of L and R differ, error {} > 5.e-16".format(errorIdentity)
            
            # calculate omegas
            du = self.du
            
            zi = self.z[Id] - [l2,l1][Id] * self.dt
                        
            du[0]  = np.interp(zi,self.z,P) - P[Id]
            du[-1] = np.interp(zi,self.z,Q) - Q[Id]
            
            self.domega[Id] = np.dot( self.L[Id][1+Id], du) + self.dt * self.L[Id][1+Id][1] * Aid * self.netGravity[n]    
            
    def updateLARLSys1InvariantFlow(self,P,Q,A,idArray=[0,-1],Ct=None,ct=None,sqrt=np.sqrt):
        '''
        Update LAMBDA,R,L,Z of the system equations
        
        Special terms:
            sys1 :          uses correct eigenvalues lambda1/2 = dlt*v +- sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            invariantFlow:  calculates R,L using riemannInvariants based on unit FLow
                            as in Paper of Paul Roger Leinan but with R2*(-1) correction to
                            have additive invariants for pressure
        
        Input:
            P                 <np.array> : current pressure values of the vessel
            Q                 <np.array> : current flow values of the vessel
            A                 <np.array> : current area values of the vessel
            idArray = [0,-1]  <list>     : define which boundary should be updated (default = both)
            update  = 'all'   <string>   : set to 'L' if only L and not R should be updated
            Ct = None         <np.array> : Compliance C if avialiable otherwise it will be calculated
            ct = None         <np.array> : waveSpeed  c if avialiable otherwise it will be calculated
        '''
        dlt = self.dlt  
        n = self.n[0]
        
        for Id in idArray:       
               
            Aid = A[Id]
               
            if Ct == None: C = self.C_nID(P,Id)
            else:  C = Ct[Id]
            if ct == None: c = self.c(Aid,C)
            else:  c = ct[Id]
            
            v = Q[Id]/Aid
            
            l1 = dlt*v+sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            l2 = dlt*v-sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            
            self.LAMBDA[Id][0] = l1
            self.LAMBDA[Id][1] = l2
            
            Z1 =  1.0/(C*l1)
            Z2 = -1.0/(C*l2)
                        
            self.Z[Id][0] = Z1
            self.Z[Id][1] = Z2
            
            ## left eigenvalue matrix
            # riemannInvariants with unit flow    
            self.L[Id][0][0] =  1.0/(Z1+Z2)
            self.L[Id][0][1] =  Z2/(Z1+Z2)
            self.L[Id][1][0] =  1.0/(Z1+Z2)
            self.L[Id][1][1] = -Z1/(Z1+Z2)

            ## right eigenvalue matrix
            # riemannInvariants with unit flow 
            self.R[Id][0][0] =  Z1
            self.R[Id][0][1] =  Z2
            self.R[Id][1][1] = -1.0

            L = self.L[Id]
            R = self.R[Id]
            ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
                            
            ### check consistency: calculate sum(R*L) == sum(I) == 2.0 
            errorIdentity = abs((L[0][0]+L[1][0])*(R[0][0]+R[0][1])+(L[0][1]+L[1][1])*(R[1][0]+R[1][1])-2.0)
            if errorIdentity > 5.e-16:
                print "WARNING: SystemEquations, inverse of L and R differ, error {} > 5.e-16".format(errorIdentity)
            
            # calculate omegas
            du = self.du
            
            zi = self.z[Id] - [l2,l1][Id] * self.dt
                        
            du[0]  = np.interp(zi,self.z,P) - P[Id]
            du[-1] = np.interp(zi,self.z,Q) - Q[Id]
            
            self.domega[Id] = np.dot( self.L[Id][1+Id], du) + self.dt * self.L[Id][1+Id][1] * Aid * self.netGravity[n]   
              
              
    def updateLARLSys1InvariantPressure(self,P,Q,A,idArray=[0,-1],Ct=None,ct=None,sqrt=np.sqrt):
        '''
        Update LAMBDA,R,L,Z of the system equations
        
        Special terms:
            sys1 :          uses correct eigenvalues lambda1/2 = dlt*v +- sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            invariantFlow:  calculates R,L using riemannInvariants based on unit Pressure
                            as derieved by Leif Rune Hellevik
        
        Input:
            P                 <np.array> : current pressure values of the vessel
            Q                 <np.array> : current flow values of the vessel
            A                 <np.array> : current area values of the vessel
            idArray = [0,-1]  <list>     : define which boundary should be updated (default = both)
            update  = 'all'   <string>   : set to 'L' if only L and not R should be updated
            Ct = None         <np.array> : Compliance C if avialiable otherwise it will be calculated
            ct = None         <np.array> : waveSpeed  c if avialiable otherwise it will be calculated
        '''
        dlt = self.dlt  
        n = self.n[0]
        for Id in idArray:
            
            Aid = A[Id]       
               
            if Ct == None: C = self.C_nID(P,Id)
            else:  C = Ct[Id]
            if ct == None: c = self.c(A[Id],C)
            else:  c = ct[Id]
            
            v = Q[Id]/A[Id] 
            
            l1 = dlt*v+sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            l2 = dlt*v-sqrt(c**2.+dlt*(dlt-1.)*(v)**2.0)
            
            self.LAMBDA[Id][0] = l1
            self.LAMBDA[Id][1] = l2
            
            Z1 =  1.0/(C*l1)
            Z2 = -1.0/(C*l2)
                                    
            self.Z[Id][0] = Z1
            self.Z[Id][1] = Z2
            
            ## left eigenvalue matrix
            # riemannInvariants with unit poressure
            self.L[Id][0][1] =  Z2
            self.L[Id][1][1] = -Z1

            ## right eigenvalue matrix
            # riemannInvariants with unit pressure
            self.R[Id][0][0] =  Z1/(Z1+Z2)
            self.R[Id][0][1] =  Z2/(Z1+Z2)
            self.R[Id][1][0] =  1.0/(Z1+Z2)
            self.R[Id][1][1] = -1.0/(Z1+Z2)
            
            L = self.L[Id]
            R = self.R[Id]
            ## calculate riemann invariants w1 at pos = -1 and w2 at pos = 0
                            
            ### check consistency: calculate sum(R*L) == sum(I) == 2.0 
            errorIdentity = abs((L[0][0]+L[1][0])*(R[0][0]+R[0][1])+(L[0][1]+L[1][1])*(R[1][0]+R[1][1])-2.0)
            if errorIdentity > 5.e-16:
                print "WARNING: SystemEquations, inverse of L and R differ, error {} > 5.e-16".format(errorIdentity)
            
            # calculate omegas
            du = self.du
            
            zi = self.z[Id] - [l2,l1][Id] * self.dt
                        
            du[0]  = np.interp(zi,self.z,P) - P[Id]
            du[-1] = np.interp(zi,self.z,Q) - Q[Id]
            
            self.domega[Id] = np.dot( self.L[Id][1+Id], du) + self.dt * self.L[Id][1+Id][1] * Aid * self.netGravity[n]  
    