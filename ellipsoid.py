"""
Plotting an ellipsoid
and rotate and project it onto a proper plane to see
how a morphology of a grain is.

One has in total, 6 variables.
Three principal axes of the ellipsoid, and its Euler angles representing
its orientation with respect to a reference coordinate system.

Try to have def in the outer of the class to facilitate
the debugging and accessibility from other modules.

"""
import numpy as np
import euler
from pylab import *
from enthought.mayavi import mlab
p3d = mlab.points3d
pi = np.pi
eul = euler.euler

def rotation(theta = 30.):
    sth = np.sin(theta*pi/180.)
    cth = np.cos(theta*pi/180.)
    return np.array([[cth,sth],[-sth,cth]])

def rotmat(theta=30.):
    sth = np.sin(theta*pi/180.)
    cth = np.cos(theta*pi/180.)
    return np.matrix([[cth,sth],[-sth,cth]])

def xy(x,y,theta):
    rst = []
    for i in range(len(x)):
        rst.append(np.dot(rotation(theta=theta),[x[i],y[i]]))
    xnew = []
    ynew = []
    for i in range(len(rst)):
        xnew.append(rst[i][0])
        ynew.append(rst[i][1])
    return xnew, ynew

def plot(x,y,x0=0,y0=0,ifig=1):
    pylab.figure(ifig)
    plot(x+x0,y+y0)


def el_3d(npoint=10000, a=1., b=1.2, c=0.9,
          phi1=90., phi=90., phi2=90., x0=0, y0=0, z0=0, ireturn =False):
    x,y,z=[],[],[]
    beta = []
    lam = []
    for i in range(npoint):
        beta.append(np.random.rand()*2*pi -pi)
        lam.append(np.random.rand()*2*pi - pi)
    for i in range(npoint):
        x.append(a*np.cos(beta[i])*np.cos(lam[i]))
        y.append(b*np.cos(beta[i])*np.sin(lam[i]))
        z.append(c*np.sin(beta[i]))
    x = np.array(x); y = np.array(y); z= np.array(z)
    x = x + x0; y = y + y0; z= z + z0
    beta = np.array(beta); lam = np.array(lam)
    rot = eul(ph=phi1, th = phi, tm=phi2,echo=False)
    rot = np.array(rot)
    xnew, ynew, znew = [], [], []
    
    for i in range(len(x)):
        rst = np.dot(rot,[x[i],y[i],z[i]])
        xnew.append(rst[0]);ynew.append(rst[1]);znew.append(rst[2])
        
    xnew = np.array(xnew); ynew=np.array(ynew); znew=np.array(znew)
    #p3d(x,y,z,mode='point')
    p3d(xnew,ynew,znew,mode='point')
    if ireturn == True: return xnew,ynew,znew

def el_3d_(npoint = 1000, a=1., b=1.2, c=0.9,
           phi1=90., phi=90., phi2=90., 
           x0=0.,y0=0.,z0=0., ireturn = False):
    x, y, z = [], [], []
    psi, ph = [],[]
    for i in range(npoint):
        psi.append(np.random.rand()*pi)
        ph.append(np.random.rand()*2*pi)
    psi = np.array(psi); ph = np.array(ph)
    for i in range(npoint):
        x.append(a*np.sin(psi[i])*np.cos(ph[i]))
        y.append(b*np.sin(psi[i])*np.sin(ph[i]))
        z.append(c*np.cos(psi[i]))
    x = np.array(x); y = np.array(y); z = np.array(z)
    x = x + x0; y = y + y0; z = z + z0
    rot = eul(ph=phi1, th = phi, tm = phi2, echo=False)
    rot = np.array(rot)
    xnew, ynew, znew = [],[],[]
    for i in range(len(x)):
        rst = np.dot(rot,[x[i],y[i],z[i]])
        xnew.append(rst[0]); ynew.append(rst[1]); znew.append(rst[2])
    xnew = np.array(xnew); ynew=np.array(ynew); znew=np.array(znew)
    p3d(xnew,ynew,znew, mode='point')
    if ireturn == True: return xnew,ynew,znew

def el_(a=1.,b=1.,c=1.,phi1=90.,phi=90.,phi2=90.,x0=0.,y0=0.,z0=0.):
    """
    Plots a rotated and elongated principal axes of the ellipoid.
    """
    
    pass
    


def proj(x,y,z):
    
    pass
    


class ellips:
    def __init__(self, axes=[1.,1.,1.], angles =[90.,90.,90.]):
        
        pass

