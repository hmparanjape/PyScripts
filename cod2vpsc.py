"""
Sensitivity of grain number in RVE upon VPSC modeling is performed
such that the RVE is prepared for individually for given COD grid file


Planned: r-values at 0., 45., 90.
         yield stresses at 0., 45., 90.
"""
###
# Over the preparation for the MSMSE
###
import numpy as np
import randomEuler; reload(randomEuler)
ss = randomEuler.sampling_simplex
#ss = randomEuler.sampling_simplex(odf=None, iang=10., ngrain, maxiter, xtol=10)...
import vp_f2py
import vpsc_examples
import matplotlib.pyplot as plt
from vpsc_param import __makenp__


def multir(niter=100, ngrain=500, interaction=3,
           odf=None, header='tmp', fsx=None):
    """
    """
    r0=[]; r45=[]; r90=[]; r135=[]; r180=[]
    ys0=[]; ys45=[]; ys90=[]; ys135=[]; ys180=[]

    for i in range(niter):
        rvalues, ysvalues = r(ngrain=ngrain, odf=odf, header=header,
                              fsx=fsx, interaction=interaction)
        #print 'r, ys as below'
        #print rvalues, ysvalues; raw_input()
        r0.append(rvalues[0])
        r45.append(rvalues[1])
        r90.append(rvalues[2])
        r135.append(rvalues[3])
        r180.append(rvalues[4])
        
        ys0.append(ysvalues[0])
        ys45.append(ysvalues[1])
        ys90.append(ysvalues[2])
        ys135.append(ysvalues[3])
        ys180.append(ysvalues[4])        
        pass

    r0, r45, r90, r135, r180 = __makenp__(r0, r45, r90, r135, r180)
    ys0, ys45, ys90, ys135, ys180 = __makenp__(ys0, ys45, ys90, ys135, ys180)

    #print 'r0', r0;raw_input()
    print 'r0 =',   r0.mean(),   r0.std()
    print 'r45 =',  r45.mean(),  r45.std()
    print 'r90 =',  r90.mean(),  r90.std()
    print 'r135 =', r135.mean(), r135.std()
    print 'r180 =', r180.mean(), r180.std()
    
    rmeans = np.array([r0.mean(), r45.mean(), r90.mean(), r135.mean(), r180.mean()])
    rstds = np.array([r0.std(), r45.std(), r90.std(), r135.std(), r180.std()])

    ysmeans = np.array([ys0.mean(), ys45.mean(), ys90.mean(), ys135.mean(), ys180.mean()])
    ysstds  = np.array([ys0.std(), ys45.std(), ys90.std(), ys135.std(), ys180.std()])

    ang = np.arange(0., 180.01, 45)
    
    return ang, rmeans, rstds, ysmeans, ysstds


def r(ngrain=500, odf=None, header='tmp', fsx=None, interaction=None):
    """
    R-values
    """
    import os
    ss(odf=odf, iang=10, ngrain=ngrain, header=header)
    
    ftex = '%s_%s.cmb'%(header,str(ngrain).zfill(5))
    a = vpsc_examples.prob(dang=45, texture=ftex, fsx=fsx,
                           interaction=interaction,
                           ifig=3)
    a.__run__()
    a.ang; a.r; a.ys 
    os.remove(ftex)   
    return a.r, a.ys


def plot(niter=2, ngrain=[100,200], interaction=3,
         odf=None, header='tmp', fsx=None, ifig=3):
    
    fig = plt.figure(ifig)
    ax0 = fig.add_subplot(111)
    fig = plt.figure(ifig+1)
    ax1 = fig.add_subplot(111)


    masterx = np.array(ngrain)

    margin = 0.01 * max(ngrain)
    r0mean =[]
    r45mean=[]
    r90mean=[]
    ys0mean=[]
    ys45mean=[]
    ys90mean=[]
    r0err=[]
    r45err=[]
    r90err=[]
    ys0err=[]
    ys45err=[]
    ys90err=[]
    
    for i in range(len(ngrain)):
        ang, rmean, rstds, ysmean, ysstds = multir(
            niter=niter, ngrain=ngrain[i], interaction=interaction,
            odf=odf, header=header, fsx=fsx)
        currentx = masterx[i]
        
        r0mean.append(rmean[0])
        r45mean.append(rmean[1])
        r90mean.append(rmean[2])
        ys0mean.append(ysmean[0])
        ys45mean.append(ysmean[1])
        ys90mean.append(ysmean[2])
        
        r0err.append(rstds[0])
        r45err.append(rstds[1])
        r90err.append(rstds[2])
        ys0err.append(ysstds[0])
        ys45err.append(ysstds[1])
        ys90err.append(ysstds[2])        
        pass
        
        #RD
    ax0.errorbar(x=masterx-margin, y=r0mean, yerr=r0err,
                 marker='o', color='black', ls=' ', mfc='None',
                 label='RD',
                 ms = 10.,
                 )
        #DD
    ax0.errorbar(x=masterx, y=r45mean, yerr=r45err,
                 marker='<', color='black', ls=' ', mfc='None',
                 ms=10.,
                 label='DD',
                 )
        #TD
    ax0.errorbar(x=masterx+margin, y=r90mean, yerr=r90err,
                 marker='d', color='black', ls=' ', mfc='None',
                 ms=10.,
                 label='TD',
                 )
    
        #RD
    ax1.errorbar(x=masterx-margin, y=ys0mean, yerr=ys0err,
                 marker='o', color='black', ls=' ', mfc='None',
                 label='RD',                 
                 ms=10.,
                 )
        #DD
    ax1.errorbar(x=masterx, y=ys45mean, yerr=ys45err,
                 marker='<', color='black', ls=' ', mfc='None',
                 ms=10.,
                 label='DD',
                 )
        #TD
    ax1.errorbar(x=masterx+margin, y=ys90mean, yerr=ys90err,
                 marker='d', color='black', ls=' ', mfc='None',
                 label='TD',
                 ms=10.,
                 )
    


    ax0.set_xlabel('number of grains')
    ax0.set_xticks(masterx)
    ax0.set_ylabel('R', dict(fontsize=20))
    ax1.set_xlabel('number of grains', dict(fontsize=20))
    ax1.set_xticks(masterx)
    ax1.set_ylabel(r'$\sigma^{YS}$ [MPa]',dict(fontsize=20))
    ax0.legend(loc='best')
    ax1.legend(loc='best')
    pass
