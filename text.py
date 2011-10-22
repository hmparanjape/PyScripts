"""
Gaussian texture components.

With given components (final aim, will start from given Euler angle set)
"""
import random
import math
import os
import pf
gauss = random.gauss

sin = math.sin
cos = math.cos
atan = math.atan
atan2 = math.atan2
asin = math.asin
acos = math.acos
pi = math.pi
sqrt = math.sqrt
#coef= math.atan(1.0)/45. #*math.pi/180.


class textur:
    """
    Design an aggregate for the given distribution parameters
    Arguments :

    Euler angles
    p1         : phi1
    p2         : phi2
    p          : PHI

    w0

    ngrain     : Number of to be created grains
    filename = 'temp.txt'
    dist = 'g','n','e','l'

    iplot=True
    """
    def __init__(self, p1, p2, p, w0, ngrain,
                 filename='temp.txt', dist='g',
                 iplot=True, idot=False):
        path =os.getcwd()
        f=open(path+'\\'+filename, 'w')
        f.writelines('Designed texture using probability distributions \n')
        f.writelines('Given Euler angles are as below \n')
        f.writelines('ph1,   phi2,   PHI = ' + str(p1) + str(p2)+ str(p) )
        
        f.write('\nB   '+ str(ngrain))
        for i in range(ngrain):
            txt = text(p1=p1, p2=p2, p=p, w0=w0, dist=dist)
            angle = txt.angle

            f.write('\n %9.2f %9.2f %9.2f %9.3f'%(angle[0],angle[1],angle[2], 0.1))
        f.close()
        if iplot==True:
            pf.pf(ftex='temp.txt',idot=idot)

class text:
    """
    Gaussian texture for the given angle
    Arguments:

    p1
    p2
    p

    w0
    dist = 'g','n','e','l'
    """
    def __init__(self, p1=45., p2=45., p=45., w0=15., dist='g'):
        """
        Arguments :
        p1
        p2
        p
        w0
        """
        ## generates a random rotation axis, about which
        ## the given grain rotates
        delta, phi = self.rot_axis()
        
        if dist=='g': 
            B = self.transform_matrix(delta=delta,
                                      phi=phi,
                                      w=self.gaussian(w0))
        if dist=='e':
            B = self.transform_matrix(delta=delta,
                                      phi=phi,
                                      w=self.expo(w0))
        if dist=='l': 
            B = self.transform_matrix(delta=delta, 
                                      phi=phi,
                                      w=self.lognorm(w0))
        if dist=='n': 
            B = self.transform_matrix(delta = delta,
                                    phi = phi,
                                    w=self.normal(w0))

        A = eul(iopt = 2, tm = p1, ph = p2, th = p) #rotation matrix of the grain to the sample
        
        c = [[0,0,0],[0,0,0],[0,0,0]]
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    c[i][j] = c[i][j]+ B[m][i] * A[m][j] # rot_{ij} = B_{mi}A_{mj}
                    
        p1, p, p2 = eul(iopt=1, a=c)
        self.angle=[p1,p,p2]
        
    def transform_matrix(self, delta, phi, w):
        """
        delta, phi, w = delta*math.pi/180. , phi*math.pi/180. , w*math.pi/180.
        Generates the rotation matrix
        """
        w = w * math.pi/180.
        d1 = sin(delta) * cos(phi)
        d2 = sin(delta) * sin(phi)
        d3 = cos(delta)
        cw = cos(w)
        sw = sin(w)
        p=[[None,None,None],[None,None,None],[None,None,None]]
        p[0][0] = ( 1 - d1**2 ) * cw + d1**2
        p[0][1] = d1 * d2 * ( 1 - cw ) + d3 * sw
        p[0][2] = d1 * d3 * ( 1 - cw ) - d2 * sw
        p[1][0] = d1 * d2 * ( 1 - cw ) - d3 * sw
        p[1][1] = ( 1 - d2**2 ) * cw + d2**2
        p[1][2] = d2 * d3 * ( 1 - cw ) + d1 * sw
        p[2][0] = d1 * d3 * ( 1 - cw ) + d2 * sw
        p[2][1] = d2 * d3 * ( 1 - cw ) - d1 * sw
        p[2][2] = ( 1 - d3**2 ) * cw + d3**2
        return p
    
    def rot_axis(self):
        """
        Random rotation axis generator
        """
        delta = random.uniform(0.,1.)
        delta = acos(delta)
        phi = random.uniform(0.,2.*math.pi)
        return delta, phi   #radians...
    
    def gaussian(self, w0):
        dp = gauss(mu=0., sigma=w0)
        return dp
    
    def expo(self, w0):
        dp = random.expovariate(w0)
        return dp
    
    def lognorm(self, w0):
        dp = random.lognormvariate(mu=0, sigma=w0)
        return dp
    
    def normal(self,w0):
        dp = random.normalvariate(mu=0, sigma=w0)
        return dp

def eul(iopt=1, ph=None, th=None, tm=None, a=None):
    """
    Calculate the euler angle associated with the tansformation
    maxtrix a(i,j) if iopt = 1 and viceversa if iopt=2
    a(i,j) transforms from system sa to system ca.
    ph,th,om are the euler angles(in degrees) of ca referred to sa.
    """
    pi = math.pi
    if iopt==1:
        if a == None:
            print 'ERR: input a is missed!!'
            return None
        th = acos(a[2][2])
        if (abs(a[2][2]) > 0.9999):
            tm = 0.
            ph = atan2(a[0][1],a[0][0])
        else:
            sth = sin(th)
            tm = atan2(a[0][2]/sth,a[1][2]/sth)
            ph = atan2(a[2][0]/sth,-a[2][1]/sth)
        th = th * 180./pi
        ph = ph * 180./pi
        tm = tm * 180./pi
        return tm, th, ph
    elif iopt==2:
        a = [[ None, None,None ],[ None, None,None ], [None,None ,None ] ]
        sph=sin(ph*pi/180.)
        cph=cos(ph*pi/180.)
        sth=sin(th*pi/180.)
        cth=cos(th*pi/180.)
        stm=sin(tm*pi/180.)
        ctm=cos(tm*pi/180.)
        a[0][0]=ctm*cph-sph*stm*cth
        a[1][0]=-stm*cph-sph*ctm*cth
        a[2][0]=sph*sth
        a[0][1]=ctm*sph+cph*stm*cth
        a[1][1]=-sph*stm+cph*ctm*cth
        a[2][1]=-sth*cph
        a[0][2]=sth*stm
        a[1][2]=ctm*sth
        a[2][2]=cth
        return a

def miller2euler(hkl=None,uvw=None):
    """
    Provided (hkl)// ND and [uvw]//RD,
    calculates and returns the Euler angles.

    * Euler angles are using Bunge nomenclature
    """
    import numpy as np
    if hkl==None or uvw==None:
        print "both hkl and uvw must be given"
        raise IOError
    hkl = np.array(hkl)
    uvw = np.array(uvw)
    #hkl = hkl/sqrt(((hkl.copy())**2).sum())
    #uvw = uvw/sqrt(((uvw.copy())**2).sum())
    print 'uvw = ', uvw
    print 'hkl = ' , hkl
    
    h = hkl.copy()[0]; k = hkl.copy()[1]; l = hkl.copy()[2]
    u = uvw.copy()[0]; v = uvw.copy()[1]; w = uvw.copy()[2]
    
    phi= acos(l/math.sqrt(h**2. + k**2. + l**2.))
    p2 = asin(h/math.sqrt(h**2. + k**2.))
    p2_= acos(k/math.sqrt(h**2. + k**2.))

    temp = w    / sqrt( u**2. + v**2. + w**2.)
    temp = temp * sqrt( h**2. + k**2. + l**2. )
    temp = temp / sqrt( h**2. + k**2.)
    p1 = asin(temp)

    return np.array([p1,phi,p2])*180./pi

    
    
