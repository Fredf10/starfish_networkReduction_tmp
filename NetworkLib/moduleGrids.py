import numpy as np 
    
def uniform(length, radiusA, radiusB, N):
    """
    Function to calculate a Uniform grid 
    with a single radius over the hole vessel length
    
    Input:  vessel-length, radiusA, radiusB, number of Gridpoints
    Output:  z  = gridpoint-array 
             dz = spacing array
             A0 = area-array with given radius
    """         
    A0 = np.ones(N)*np.pi*radiusA**2.
    z  = np.linspace(0.,length,N)
    dz = np.diff(z)
    
    return z,dz,A0

def cone(length, radiusA, radiusB, N):
    """
    Function to calculate a Cone grid 
    with a radiusA at the proximal end and radiusB at the distal end
    
    Input:  vessel-length, radiusA, radiusB, number of Gridpoints
    Output:  z  = gridpoint-array 
             dz = spacing array
             A0 = area-array with given radius
    """         
    z  = np.linspace(0.,length,N)
    dz = np.diff(z)
    A0 = np.pi*(radiusA + ((radiusB-radiusA)*np.linspace(0.0,1.0,N)))**2
    
    return z,dz,A0

def constriction(length, radiusA, radiusB, N):
    """
    Function to calculate a constricted grid 
    with a radiusA as proximal and distal radius and radiusB
    as constriction radius; the constriction is calculated via a cos-function
    
    Input:  vessel-length, radiusA, radiusB, number of Gridpoints
    Output:  z  = gridpoint-array 
             dz = spacing array
             A0 = area-array with given radius
    """             
    z  = np.linspace(0.,length,N)
    dz = np.diff(z)
    R0 = radiusA + (radiusB-radiusA)*(1.0 + np.cos(np.linspace(-1.,1.,N)*np.pi))/2.
    A0 = np.pi*R0**2
    print "WARNING: constricted vessel geometry was never verified and vaildated!"
    return z,dz,A0

def stenosed(length, radiusA, radiusB, N, position, constriction, stenosisLength):
    """
    Function to calculate a Cone grid 
    with a radiusA at the proximal end and radiusB at the distal end
    
    In addition a stenoses is inserted with a defined stenoses level, length and position
    
    Input:  vessel-length, radiusA, radiusB, number of Gridpoints
    Output:  z  = gridpoint-array 
             dz = spacing array
             A0 = area-array with given radius
    """         
    z  = np.linspace(0.,length,N)
    dz = np.diff(z)
    A0 = np.pi*(radiusA + ((radiusB-radiusA)*np.linspace(0.0,1.0,N)))**2
    
    # find indices where the stenoses is
    target = [position,position+stenosisLength]
    idx = z.searchsorted(target)
    idx = np.clip(idx, 1, len(z)-1)
    left = z[idx-1]
    right = z[idx]
    idx -= target - left < right - target
    
    
    A0[idx[0]:idx[1]+1] = A0[idx[0]:idx[1]+1] * constriction
    
    return z,dz,A0