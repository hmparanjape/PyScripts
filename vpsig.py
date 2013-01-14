"""
Stress implicit increment for VPSC
"""
from vp_f2py import vpsc
import numpy as np
import os
from vpsc_in import histmaker
from math import *

### guess initial (strain incremental)
### derivative ? (if possible for accelerating optimization)
def main(epsincr=None, incrsig=None, eqincr=0.005,
         jobold=0, jobnew=1, save=False):
    """
    save will be flag for saving the postmort file
    for starting from the last run by the implicit loop
    """
    if None in [epsincr, incrsig]:
        raise IOError, 'You should input epsincr and incrsig'

    fnpostold = 'POSTMORT_%s.OUT'%str(jobold).zfill(4)
    fnpostnew = 'POSTMORT_%s.IN'%str(jobnew).zfill(4)
    os.rename(fnpostold, fnpostnew) #rename

    ## optimization loop
    #tmp = vpsc(prcs=prcs, 
    ##
    return #return suitable strain rate increment...

def optimization(x0=None, #parameter
                 ftex=None, fsx=None,
                 xtol=0.000000001,
                 maxiter=2000,jobid=200,eqincr=0.005,thet0=45., #target
                 ):
    """
    thet0=target

    x0    :[degree] initial value
    thet0 :[degree] target
    """
    from scipy.optimize import fmin_bfgs, fmin
    xopt = fmin(
        ipbiaxial, #defunc to be minimized
        x0=[x0],
        args=(eqincr, ftex, fsx, jobid,thet0), #thet0 is the target
        maxiter=maxiter,xtol=xtol,
        maxfun=maxiter, retall=True
        )
    return xopt

def ipbiaxial(x0, eqincr, ftex, fsx, jobid, theta0):
    """
    in-plane biaxial

    theta in radian
    x0 in radian
    """
    filename = 'ipbiaxial_attemp_temp_%s'%str(jobid).zfill(4)
    r = 1
    x0 = x0[0]
    x = cos(x0)*r
    y = sin(x0)*r
    x = round(x,9)
    y = round(y,9)
    z = -round(( x + y ),9)

    histmaker(
        filename=filename,
        nstep = 1, ictrl = 7, eqincr= eqincr,
        iudot=[[1,1,1],
               [0,1,1],
               [0,0,0]],
        udot=[[x,0,0],
              [0,y,0],
              [0,0,z]],
        iscau=[0,0,1,1,1,1],
        scauchy=[0,0,0,0,0,0]
        )
    prcs = ['0',filename]
    job = vpsc(prcs=prcs, texture=ftex, fsx=fsx,stp=1,jobid=jobid)
    tmp = job.run()
    dsx0=tmp['sbar'][0][0][0,0] - tmp['sbar'][0][0][2,2]
    dsy0=tmp['sbar'][0][0][1,1] - tmp['sbar'][0][0][2,2]
    dsx1=tmp['sbar'][0][1][0,0] - tmp['sbar'][0][1][2,2]
    dsy1=tmp['sbar'][0][1][1,1] - tmp['sbar'][0][1][2,2]

    dsx = dsx1-dsx0
    dsy = dsy1-dsy0
    
    os.system('clear')
    print 'edotx=',x
    print 'edoty=',y    
    print 'dsig_xx =', dsx
    print 'dsig_yy =', dsy

    theta1 = atan2(dsy, dsx) #radian
    print 'theta1 (degree) = ', theta1*180./np.pi
    return theta1-theta0*np.pi/180.
