"""
In order to rotate the sample to perform VPSC simulation for
the 2011 Numisheet benchmark problem.
"""
import euler
import numpy as np
from vpsc_param import __makenp__ 
euler=euler.euler
import upf

def main(gr=None, A=None, x=[1,0,0], y=[0,1,0], z=[0,0,1], ifig=1):
    """
    gr : grains,
    A : rotation matrix, which transform the given coordinate system
    to the desired system.
    """
    if A==None:
        A = np.zeros((3, 3))
        x, y, z = __makenp__(x, y, z)
        A[0] = x
        A[1] = y
        A[2] = z
        A = A.transpose()
        pass
    
    gr = gr.copy()
    ngr = gr.shape[0]
    
    for i in range(ngr):
        #phi1 , phi, phi2
        a = euler(ph=gr[i][0], th=gr[i][1], tm=gr[i][2], echo=False)
        Ainv=A.transpose()
        newa = np.dot(a,Ainv) #newa = a A^-1
        phi1, phi, phi2 = euler(a=newa, echo=False)
        gr[i][0] = phi1
        gr[i][1] = phi
        gr[i][2] = phi2
        pass

    ## pole figure plotting
    temp = upf.polefigure(grains=gr, csym='cubic')
    temp.pf(ifig=ifig)
    
    
    return gr
