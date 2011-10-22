"""
yield surface plotting predcited by a crystal plasticity code
"""
from vpsc_in import histmaker
import os
import numpy as np
import vp_f2py
import matplotlib.pyplot as plt
from math import atan2, sqrt,tan
cos = np.cos; sin = np.sin; pi=np.pi
tiny = 10**-8
### history file maker and returns the filename
def xxyy(jobid=0, n=10, ftex=None, fsx=None, interaction=3, ifig=1, init=0, fin=360):
    """
    Yield locus in the x-y plane
    """
    theta = np.linspace(init,fin,n) * np.pi / 180.
    filename = []
    tmp = 0
    for i in range(n):
        while True:
            fn = 'hist/xxyy_%s_%s.hist'%(
                str(jobid).zfill(4),
                str(tmp).zfill(4))
            if not(os.path.exists(fn)): break
            tmp = tmp + 100
            pass
        tmp = tmp + 1
        filename.append(fn)
        pass
    r = 1.0
    x = cos(theta) * r
    y = sin(theta) * r
    z = -(x+y)
    prcs = [ ]
    for i in range(len(filename)):
        prcs.append('0')
        prcs.append(filename[i])
        x[i] = round(x[i],9)
        y[i] = round(y[i],9)
        z[i] = -x[i]-y[i]
        histmaker(
            filename = filename[i],
            nstep = 1, ictrl = 7, eqincr=0.005,
            iudot=[[1,0,0],
                   [1,1,0],
                   [1,1,0]],
            udot = [[x[i],    0,  0],
                    [   0, y[i],  0],
                    [   0,    0,  z[i]]],
            iscau   = [0, 0, 1, 1, 1, 1], 
            scauchy = [0, 0, 0, 0, 0, 0],
            )
        pass
    x = [ ] #sig_{xx}
    y = [ ] #sig_{yy}
    job=vp_f2py.vpsc(
        texture=ftex, fsx=fsx, iupdate=[0,0,0,0], stp=0,
        interaction=interaction,prcs=prcs, jobid=jobid)
    tmpmaster = job.run()
    sbar = tmpmaster['sbar']
    dbar = tmpmaster['dbar']
    sxx = [] ; dxx=[]
    syy = [] ; dyy=[]
    for i in range(len(filename)):
        sxx.append(sbar[i][0][0,0] - sbar[i][0][2,2])
        syy.append(sbar[i][0][1,1] - sbar[i][0][2,2])
        dxx.append(dbar[i][0][0,0])
        dyy.append(dbar[i][0][1,1])
        pass
    sxx.append(sxx[0])
    syy.append(syy[0])
    d = [] #strain rate angle (normal to yield surface)
    for i in range(len(filename)):
        d.append(atan2(dyy[i],dxx[i]))
        pass
    ## plotting
    fig = plt.figure(ifig); ax=fig.add_subplot(111)
    ax.plot(sxx,syy,'.',label=r'$\sigma_{xx}-\sigma_{yy}$')
    ax.set_xlabel(r'$\sigma_{xx}$', dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{yy}$', dict(fontsize=20))
    r = 30
    for i in range(len(d)):
        ax.plot(
            [sxx[i],sxx[i]+r*cos(d[i])], #[x0,x1]
            [syy[i],syy[i]+r*sin(d[i])], #[y0,y1]
            color='gray', alpha=0.5, ls='-')
        pass
    ## plotting
    for fn in filename:
        os.remove(fn)
        pass
    return sxx,syy

def xxxy(jobid=0, n=10,ftex=None,fsx=None,interaction=3, ifig=1):
    """
    plane consisting of sig_{xx} and tau_{xy}
    
    Yield locus in the xx-xy plane
    """
    theta = np.linspace(0,360,n) * np.pi / 180.
    filename=[]
    tmp = 0
    for i in range(n):
        while True:
            fn = 'hist/xxyy_%s_%s.hist'%(
                str(jobid).zfill(4),
                str(tmp).zfill(4))
            if not(os.path.exists(fn)): break
            tmp = tmp + 100
            pass
        tmp = tmp + 1
        filename.append(fn)
        pass
    r = 1.0
    x = cos(theta) * r #xx
    y = sin(theta) * r #xy
    z = -(x+y)
    prcs = [ ]
    for i in range(len(filename)):
        prcs.append('0')
        prcs.append(filename[i])
        x[i] = round(x[i],4)
        y[i] = round(y[i],4)
        z[i] = round(-x[i],4)
        if x[i]+z[i]!=0: raise IOError,'Non zero'
        histmaker(
            filename = filename[i],
            nstep = 1, ictrl = 7, eqincr=0.005,
            iudot=[[1,1,0],
                   [1,0,0],
                   [1,1,0]],
            udot = [[x[i], y[i],  0],
                    [   0,    0,  0],
                    [   0,    0,  z[i]]],
            iscau   = [0, 1, 1, 1, 1, 0], 
            scauchy = [0, 0, 0, 0, 0, 0],
            )
        pass

    x = [ ] #sig_{xx}
    y = [ ] #sig_{yy}

    job=vp_f2py.vpsc(
        texture=ftex, fsx=fsx, iupdate=[0,0,0,0], stp=0,
        interaction=interaction,prcs=prcs, jobid=jobid)
    tmpmaster = job.run()
    sbar = tmpmaster['sbar']
    dbar = tmpmaster['dbar']

    sxx = [] ; dxx=[]
    syy = [] ; dyy=[]
    
    for i in range(len(filename)):
        sxx.append(sbar[i][0][0,0] - sbar[i][0][2,2])
        syy.append(sbar[i][0][0,1])
        dxx.append(dbar[i][0][0,0])
        dyy.append(dbar[i][0][0,1])        
    sxx.append(sxx[0])
    syy.append(syy[0])
    d=[]
    for i in range(len(filename)):
        d.append(atan2(dyy[i],dxx[i]))
        pass
    ## plotting
    fig = plt.figure(ifig); ax=fig.add_subplot(111)
    ax.plot(sxx,syy,'.',label=r'$\sigma_{xx}-\sigma_{xy}$')
    ax.set_xlabel(r'$\sigma_{xx}$', dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{xy}$', dict(fontsize=20))
    r = 30
    for i in range(len(d)):
        ax.plot(
            [sxx[i],sxx[i]+r*cos(d[i])], #[x0,x1]
            [syy[i],syy[i]+r*sin(d[i])], #[y0,y1]
            color='gray', alpha=0.5, ls='-')
        pass
    ## plotting
    
    #return sbar, dbar
    return sxx,syy

def yyxy(jobid=0, n=10,ftex=None,fsx=None,interaction=3,ifig=1):
    """
    plane consisting of sig_{yy} and tau_{xy}
    
    Yield locus in the yy-xy plane
    """
    theta = np.linspace(0,360,n) * np.pi / 180.
    filename=[]
    tmp = 0
    for i in range(n):
        while True:
            fn = 'hist/xxyy_%s_%s.hist'%(
                str(jobid).zfill(4),
                str(tmp).zfill(4))
            if not(os.path.exists(fn)): break
            tmp = tmp + 100
            pass
        tmp = tmp + 1
        filename.append(fn)
        pass
    r = 1.0
    x = cos(theta) * r #yy
    y = sin(theta) * r #xy
    z = -(x+y)
    prcs = [ ]
    for i in range(len(filename)):
        prcs.append('0')
        prcs.append(filename[i])

        x[i] = round(x[i],4)
        y[i] = round(y[i],4)
        z[i] = round(-x[i],4)
        if x[i]+z[i]!=0: raise IOError,'Non zero'
        histmaker(
            filename = filename[i],
            nstep = 1, ictrl = 7, eqincr=0.005,
            iudot=[[0,1,0],
                   [1,1,0],
                   [1,1,0]],
            udot = [[0, y[i],  0],
                    [   0,   x[i],  0],
                    [   0,    0,  z[i]]],
            iscau   = [1, 0, 1, 1, 1, 0], 
            scauchy = [0, 0, 0, 0, 0, 0],
            )
        pass

    x = [ ] #sig_{xx}
    y = [ ] #sig_{yy}

    job=vp_f2py.vpsc(
        texture=ftex, fsx=fsx, iupdate=[0,0,0,0], stp=0,
        interaction=interaction,prcs=prcs, jobid=jobid)
    tmpmaster = job.run()
    sbar = tmpmaster['sbar']
    dbar = tmpmaster['dbar']

    sxx = [] ; dxx=[]
    syy = [] ; dyy=[]
    for i in range(len(filename)):
        sxx.append(sbar[i][0][1,1] - sbar[i][0][2,2])
        syy.append(sbar[i][0][0,1])
        dxx.append(dbar[i][0][1,1])
        dyy.append(dbar[i][0][0,1])
        pass        
    sxx.append(sxx[0])
    syy.append(syy[0])
    d = [] #strain rate angle (normal to yield surface)
    for i in range(len(filename)):
        d.append(atan2(dyy[i],dxx[i]))
        pass    
    ## plotting
    fig = plt.figure(ifig); ax=fig.add_subplot(111)
    ax.plot(sxx,syy,'.',label=r'$\sigma_{yy}-\sigma_{xy}$')
    ax.set_xlabel(r'$\sigma_{yy}$', dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{xy}$', dict(fontsize=20))
    r = 30
    for i in range(len(d)):
        ax.plot(
            [sxx[i],sxx[i]+r*cos(d[i])], #[x0,x1]
            [syy[i],syy[i]+r*sin(d[i])], #[y0,y1]
            color='gray', alpha=0.5, ls='-')
        pass
    ## plotting
    #return sbar, dbar
    return sxx,syy

def xxxz(jobid=0, n=10,ftex=None,fsx=None,interaction=3,ifig=1):
    """
    plane consisting of sig_{xx} and tau_{xz}
    
    Yield locus in the xx-xz plane
    """
    theta = np.linspace(0,360,n) * np.pi / 180.
    filename=[]
    tmp = 0
    for i in range(n):
        while True:
            fn = 'hist/xxyy_%s_%s.hist'%(
                str(jobid).zfill(4),
                str(tmp).zfill(4))
            if not(os.path.exists(fn)): break
            tmp = tmp + 100
            pass
        tmp = tmp + 1
        filename.append(fn)
        pass
    r = 1.0
    x = cos(theta) * r #xx
    y = sin(theta) * r #xz
    z = -(x+y)
    prcs = [ ]
    for i in range(len(filename)):
        prcs.append('0')
        prcs.append(filename[i])

        x[i] = round(x[i],4)
        y[i] = round(y[i],4)
        z[i] = round(-x[i],4)
        if x[i]+z[i]!=0: raise IOError,'Non zero'
        histmaker(
            filename = filename[i],
            nstep = 1, ictrl = 7, eqincr=0.005,
            iudot=[[1,0,1],
                   [1,0,0],
                   [1,1,0]],
            udot = [[x[i], 0,  y[i]],
                    [   0,    0,  0],
                    [   0,    0,  z[i]]],
            iscau   = [0, 1, 1, 1, 0, 1], 
            scauchy = [0, 0, 0, 0, 0, 0],
            )
        pass

    x = [ ] #sig_{xx}
    y = [ ] #sig_{yy}

    job=vp_f2py.vpsc(
        texture=ftex, fsx=fsx, iupdate=[0,0,0,0], stp=0,
        interaction=interaction,prcs=prcs, jobid=jobid)
    tmpmaster = job.run()
    sbar = tmpmaster['sbar']
    dbar = tmpmaster['dbar']

    sxx = [] ; dxx=[]
    syy = [] ; dyy=[]
    for i in range(len(filename)):
        sxx.append(sbar[i][0][0,0] - sbar[i][0][2,2])
        syy.append(sbar[i][0][0,2])
        dxx.append(dbar[i][0][0,0])
        dyy.append(dbar[i][0][0,2])
    sxx.append(sxx[0])
    syy.append(syy[0])
    d = [] #strain rate angle (normal to yield surface)
    for i in range(len(filename)):
        d.append(atan2(dyy[i],dxx[i]))
        pass

    ## plotting
    fig = plt.figure(ifig); ax=fig.add_subplot(111)
    ax.plot(sxx,syy,'.',label=r'$\sigma_{xx}-\sigma_{xz}$')
    ax.set_xlabel(r'$\sigma_{xx}$', dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{xz}$', dict(fontsize=20))
    r = 30
    for i in range(len(d)):
        ax.plot(
            [sxx[i],sxx[i]+r*cos(d[i])], #[x0,x1]
            [syy[i],syy[i]+r*sin(d[i])], #[y0,y1]
            color='gray', alpha=0.5, ls='-')
        pass
    ## plotting    

    #return sbar, dbar
    return sxx,syy

def yyxz(jobid=0, n=10,ftex=None,fsx=None,interaction=3,ifig=1):
    """
    plane consisting of sig_{yy} and tau_{xz}
    
    Yield locus in the yy-xz plane
    """
    theta = np.linspace(0,360,n) * np.pi / 180.
    filename=[]
    tmp = 0
    for i in range(n):
        while True:
            fn = 'hist/xxyy_%s_%s.hist'%(
                str(jobid).zfill(4),
                str(tmp).zfill(4))
            if not(os.path.exists(fn)): break
            tmp = tmp + 100
            pass
        tmp = tmp + 1
        filename.append(fn)
        pass
    r = 1.0
    x = cos(theta) * r #xx
    y = sin(theta) * r #xz
    z = -(x+y)
    prcs = [ ]
    for i in range(len(filename)):
        prcs.append('0')
        prcs.append(filename[i])

        x[i] = round(x[i],4)
        y[i] = round(y[i],4)
        z[i] = round(-x[i],4)
        if x[i]+z[i]!=0: raise IOError,'Non zero'
        histmaker(
            filename = filename[i],
            nstep = 1, ictrl = 7, eqincr=0.005,
            iudot=[[1,0,1],
                   [1,0,0],
                   [1,1,0]],
            udot = [[x[i], 0,  y[i]],
                    [   0,    0,  0],
                    [   0,    0,  z[i]]],
            iscau   = [0, 1, 1, 1, 0, 1], 
            scauchy = [0, 0, 0, 0, 0, 0],
            )
        pass

    x = [ ] #sig_{xx}
    y = [ ] #sig_{yy}

    job=vp_f2py.vpsc(
        texture=ftex, fsx=fsx, iupdate=[0,0,0,0], stp=0,
        interaction=interaction,prcs=prcs, jobid=jobid)

    tmpmaster = job.run()
    sbar = tmpmaster['sbar']
    dbar = tmpmaster['dbar']

    sxx = []; dxx=[]
    syy = []; dyy=[]
    for i in range(len(filename)):
        sxx.append(sbar[i][0][0,0] - sbar[i][0][2,2])
        syy.append(sbar[i][0][0,2])
        dxx.append(dbar[i][0][0,0])
        dyy.append(dbar[i][0][0,2])
        pass
    
    sxx.append(sxx[0])
    syy.append(syy[0])
    d = [] #strain rate angle (normal to yield surface)
    for i in range(len(filename)):
        d.append(atan2(dyy[i],dxx[i]))
        pass

    ## plotting
    fig = plt.figure(ifig); ax=fig.add_subplot(111)
    ax.plot(sxx,syy,'.',label=r'$\sigma_{yy}-\sigma_{xz}$')
    ax.set_xlabel(r'$\sigma_{yy}$', dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{xz}$', dict(fontsize=20))
    r = 30
    for i in range(len(d)):
        ax.plot(
            [sxx[i],sxx[i]+r*cos(d[i])], #[x0,x1]
            [syy[i],syy[i]+r*sin(d[i])], #[y0,y1]
            color='gray', alpha=0.5, ls='-')
        pass
    ## plotting    
    return sxx,syy

def xxyz(jobid=0, n=10,ftex=None,fsx=None,interaction=3,ifig=3):
    """
    plane consisting of sig_{xx} and tau_{yz}
    
    Yield locus in the xx-yz plane
    """
    theta = np.linspace(0,360,n) * np.pi / 180.
    filename=[]
    tmp = 0
    for i in range(n):
        while True:
            fn = 'hist/xxyy_%s_%s.hist'%(
                str(jobid).zfill(4),
                str(tmp).zfill(4))
            if not(os.path.exists(fn)): break
            tmp = tmp + 100
            pass
        tmp = tmp + 1
        filename.append(fn)
        pass
    r = 1.0
    x = cos(theta) * r #xx
    y = sin(theta) * r #xz
    z = -(x+y)
    prcs = [ ]
    for i in range(len(filename)):
        prcs.append('0')
        prcs.append(filename[i])

        x[i] = round(x[i],4)
        y[i] = round(y[i],4)
        z[i] = round(-x[i],4)
        if x[i]+z[i]!=0: raise IOError,'Non zero'
        histmaker(
            filename = filename[i],
            nstep = 1, ictrl = 7, eqincr=0.005,
            iudot=[[1,0,0],
                   [1,0,1],
                   [1,1,0]],
            udot = [[x[i],  0,  0],
                    [   0,    0,  y[i]],
                    [   0,  0,  z[i]]],
            iscau   = [0, 1, 1, 0, 1, 1], 
            scauchy = [0, 0, 0, 0, 0, 0],
            )
        pass

    x = [ ] #sig_{xx}
    y = [ ] #sig_{yy}

    job=vp_f2py.vpsc(
        texture=ftex, fsx=fsx, iupdate=[0,0,0,0], stp=0,
        interaction=interaction,prcs=prcs, jobid=jobid)
    tmpmaster = job.run()
    sbar = tmpmaster['sbar']
    dbar = tmpmaster['dbar']

    sxx = [];dxx=[]
    syy = [];dyy=[]
    for i in range(len(filename)):
        sxx.append(sbar[i][0][0,0] - sbar[i][0][2,2])
        syy.append(sbar[i][0][1,2])
        dxx.append(dbar[i][0][0,0])
        dyy.append(dbar[i][0][1,2])
        pass
    
    sxx.append(sxx[0])
    syy.append(syy[0])

    d=[]
    for i in range(len(filename)):
        d.append(atan2(dyy[i],dxx[i]))
        pass

    ## plotting
    fig = plt.figure(ifig); ax=fig.add_subplot(111)
    ax.plot(sxx,syy,'.',label=r'$\sigma_{xx}-\sigma_{yz}$')
    ax.set_xlabel(r'$\sigma_{xx}$', dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{yz}$', dict(fontsize=20))
    r = 30
    for i in range(len(d)):
        ax.plot(
            [sxx[i],sxx[i]+r*cos(d[i])], #[x0,x1]
            [syy[i],syy[i]+r*sin(d[i])], #[y0,y1]
            color='gray', alpha=0.5, ls='-')
        pass
    ## plotting
    
    #return sbar, dbar
    return sxx,syy


def yyyz(jobid=0, n=10,ftex=None,fsx=None,interaction=3,ifig=1):
    """
    plane consisting of sig_{yy} and tau_{yz}
    
    Yield locus in the xx-xz plane
    """
    theta = np.linspace(0,360,n) * np.pi / 180.
    filename=[]
    tmp = 0
    for i in range(n):
        while True:
            fn = 'hist/xxyy_%s_%s.hist'%(
                str(jobid).zfill(4),
                str(tmp).zfill(4))
            if not(os.path.exists(fn)): break
            tmp = tmp + 100
            pass
        tmp = tmp + 1
        filename.append(fn)
        pass
    r = 1.0
    x = cos(theta) * r #xx
    y = sin(theta) * r #xz
    z = -(x+y)
    prcs = [ ]
    for i in range(len(filename)):
        prcs.append('0')
        prcs.append(filename[i])

        x[i] = round(x[i],4)
        y[i] = round(y[i],4)
        z[i] = round(-x[i],4)
        if x[i]+z[i]!=0: raise IOError,'Non zero'
        histmaker(
            filename = filename[i],
            nstep = 1, ictrl = 7, eqincr=0.005,
            iudot=[[0,0,0],
                   [1,1,1],
                   [1,1,0]],
            udot = [[0, 0,  0],
                    [   0, x[i],  y[i]],
                    [   0,    0,  z[i]]],
            iscau   = [1, 0, 1, 0, 1, 1], 
            scauchy = [0, 0, 0, 0, 0, 0],
            )
        pass

    x = [ ] #sig_{xx}
    y = [ ] #sig_{yy}

    job=vp_f2py.vpsc(
        texture=ftex, fsx=fsx, iupdate=[0,0,0,0], stp=0,
        interaction=interaction,prcs=prcs, jobid=jobid)
    tmpmaster = job.run()
    sbar = tmpmaster['sbar']
    dbar = tmpmaster['dbar']

    sxx = []; dxx=[]
    syy = []; dyy=[]
    for i in range(len(filename)):
        sxx.append(sbar[i][0][1,1] - sbar[i][0][2,2])
        syy.append(sbar[i][0][1,2])
        dxx.append(dbar[i][0][1,1])
        dyy.append(dbar[i][0][1,2])
        pass

    sxx.append(sxx[0])
    syy.append(syy[0])
    d = [] #strain rate angle (normal to yield surface)
    for i in range(len(filename)):
        d.append(atan2(dyy[i],dxx[i]))
        pass    
    ## plotting
    fig = plt.figure(ifig); ax=fig.add_subplot(111)
    ax.plot(sxx,syy,'.',label=r'$\sigma_{yy}-\sigma_{yz}$')
    ax.set_xlabel(r'$\sigma_{yy}$', dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{yz}$', dict(fontsize=20))
    r = 30
    for i in range(len(d)):
        ax.plot(
            [sxx[i],sxx[i]+r*cos(d[i])], #[x0,x1]
            [syy[i],syy[i]+r*sin(d[i])], #[y0,y1]
            color='gray', alpha=0.5, ls='-')
        pass
    ## plotting
    #return sbar, dbar
    return sxx,syy




### -----------------------------------------------------
# yield surface probing (xx-yy-xy) stress component space
def xxyyxy(
    jobid=0, n=10, m=3, ftex=None,
    fsx=None, interaction=3, ifig=1
    ):
    """
    Yield locus in the xx-yy-xy space
    """
    from mpl_toolkits.mplot3d import axes3d, Axes3D

    r=1.0
    theta = np.linspace(0,360,n) * np.pi / 180.
    phi = np.linspace(-89.99,89.99,m) * np.pi / 180.
    filename=[]
    tmp = 0
    for i in range(n-1):
        for j in range(m-1):
            while True:
                fn = 'hist/xxyy_%s_%s.hist'%(
                    str(jobid).zfill(4),
                    str(tmp).zfill(4))
                if not(os.path.exists(fn)): break
                tmp = tmp + 100
                pass
            tmp = tmp + 1
            filename.append(fn)
            pass
        pass
    for j in range(2):
        while True:
            fn = 'hist/xxyy_%s_%s.hist'%(
                str(jobid).zfill(4),
                str(tmp).zfill(4))
            if not(os.path.exists(fn)): break
            tmp = tmp + 100
            pass
        tmp = tmp + 1
        filename.append(fn)
        pass
    
    ## Velocity gradient imposition ----------------------------
    xy=[];x=[];y=[]
    for i in range(m+1): #Loop over xy axis (phi)

        if i==0 or i==m:
            if i==0: xy.append(1)
            elif i==m: xy.append(-1)
            x.append(0.)
            y.append(0.)
            pass
        else:
            temp = sin(phi[i]) * r
            dd = cos(phi[i]) * r
            for j in range(n-1):
                xy.append(temp)
                x.append(cos(theta[j])*dd)
                y.append(sin(theta[j])*dd)
                pass
            pass
        pass
    x = np.array(x); y = np.array(y); xy = np.array(xy)
    prcs = [ ]
    for i in range(len(filename)):
        prcs.append('0')
        prcs.append(filename[i])
        x[i] = round(x[i],9)
        y[i] = round(y[i],9)
        z = -x[i]-y[i]
        histmaker(
            filename = filename[i],
            nstep = 1, ictrl = 7, eqincr=0.005,
            iudot=[[1,1,0],
                   [1,1,0],
                   [1,1,0]],
            udot = [[x[i], xy[i],  0],
                    [   0,  y[i],  0],
                    [   0,     0,  z]],
            iscau   = [0, 0, 1, 1, 1, 0], 
            scauchy = [0, 0, 0, 0, 0, 0],
            )
        pass
    ## -------------------------------------------------------
    x = [ ] #sig_{xx}
    y = [ ] #sig_{yy}
    z = [ ] #sig_{xy}
    job=vp_f2py.vpsc(
        texture=ftex, fsx=fsx, iupdate=[0,0,0,0], stp=0,
        interaction=interaction,prcs=prcs, jobid=jobid)
    tmpmaster = job.run()
    sbar = tmpmaster['sbar']; dbar = tmpmaster['dbar']
    sxx = [] ; dxx=[]
    syy = [] ; dyy=[]
    sxy = [] ; dxy=[]
    for i in range(len(filename)):
        sxx.append(sbar[i][0][0,0] - sbar[i][0][2,2])
        syy.append(sbar[i][0][1,1] - sbar[i][0][2,2])
        sxy.append(sbar[i][0][0,1])
        dxx.append(dbar[i][0][0,0])
        dyy.append(dbar[i][0][1,1])
        dxy.append(dbar[i][0][0,1])
        pass

    # sxx.append(sxx[0])
    # syy.append(syy[0])
    # sxy.append(sxy[0])
    d = [] #strain rate angle (normal to yield surface)
    t = []
    for i in range(len(filename)):
        norm = sqrt(dxx[i]**2+dyy[i]**2)
        d.append(atan2(dyy[i],dxx[i]))
        t.append(atan2(dxy[i],norm))
        pass

    ## plotting
    fig = plt.figure(ifig); #ax=fig.add_subplot(111)
    ax = Axes3D(fig)
    ax.plot(sxx, syy, sxy, '.', label=r'$\sigma_{xx}-\sigma_{yy}-\sigma_{xy}$')
    ax.set_xlabel(r'$\sigma_{xx}$', dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{yy}$', dict(fontsize=20))
    ax.set_zlabel(r'$\sigma_{xy}$', dict(fontsize=20))

    ax.set_axis_off()
    
    r = 30
    # for i in range(len(d)):
    #     ax.plot(
    #         [sxx[i],sxx[i]+cos(d[i])], #[x0,x1]
    #         [syy[i],syy[i]+sin(d[i])], #[y0,y1]
    #         [sxy[i],sxy[i]+tan(t[i])], #[z0,z1]
    #         color='gray', alpha=0.5, ls='-')
    #     pass
    ## plotting
    
    return sxx,syy,sxy


def uni(ifig=1, ftex=None, fsx=None, dang=None, jobid=20):
    """
    in-plane variation of uniaxial yield stress mapped into
    xx-yy-xy stress space
    """
    from vpsc_examples import prob
    R = np.zeros((3,3)) #rotation matrix
    
    myprob = prob(
        dang=dang, texture=ftex,
        fsx=fsx, jobid=jobid)
    
    ang, ys, r , nys = myprob.__run__()

    xx=[];yy=[];xy=[]
    for i in range(len(ang)):
        # rotation matrix construction
        rang = ang[i]*np.pi/180.
        R[0,0] = cos(rang)
        R[0,1] = sin(rang)
        R[1,0] = -R[0,1]
        R[1,1] = R[0,0]
        R[2,2] = 1.
        ## --
        ## - uniaxial tension stress-state
        sig = np.zeros((3,3))
        sig[0,0] = ys[i]
        ## --
        ## rotation of the stress tensor
        #sig = np.dot(np.dot(R,sig),R.transpose())
        sig = np.dot(np.dot(R.transpose(),sig),R)
        xx.append(sig[0,0])
        yy.append(sig[1,1])
        xy.append(sig[0,1])
        pass
    return xx,yy,xy

def coder(ifig=1, ftex=None, fsx=None,m=None,n=None):
    """
    Rotate the axis and code it as a movie file
    """
    
    
    import matplotlib.pyplot as plt
    import subprocess
    # xxyyxy yield loci with different level of xy
    sxx, syy, sxy = xxyyxy(jobid=0, n=n,m=m,ifig=ifig, ftex=ftex,fsx=fsx)
    # in-plane biaxial
    bxx, byy = xxyy(jobid=22, n=n, ftex=ftex, fsx=fsx,
                   ifig=ifig+2, init=-45,
                   fin = 135)
    bzz = np.zeros((len(byy),))
    # in-plane uniaxial yield locus fit into the same xxyyxy stress space
    ux, uy, uxy = uni(ifig=ifig+1, ftex=ftex, fsx=fsx, dang=5)
    ax = plt.figure(ifig).axes[0]

    ## in-plane uniaxial and biaxial tests --------------
    ax.plot(ux, uy, uxy, label='in-plane uniaxial')
    upto = len(bxx)-1
    ax.plot(bxx[0:upto], byy[0:upto], bzz[0:upto], label='in-plane biaxial')
    ## --------------------------------------------------

    data = np.array([sxx,syy,sxy]).transpose()
    np.savetxt('xyz',data) #save the data for latter use.
    fig = plt.figure(ifig)
    ax = fig.axes[0]

    ## rotation of the view
    azimuth = np.linspace(-180,180,1000)
    elevation=np.array([45.])
    os.system('rm -f *.png') #remove the existing png files..
    k = 0
    for i in elevation:
        for j in azimuth:
            filename=str('%s'%str(k).zfill(5))+'.png'
            ax.view_init(i,j)
            fig.savefig(filename,dpi=100)
            k = k + 1
            pass
        pass
    ## -------------------
    
    command= ('mencoder',
              'mf://*.png',
              '-mf',
              'type=png:w=1200:h=900:fps=40',
              '-ovc',
              'lavc',
              '-lavcopts',
              'vcodec=mpeg4',
              '-oac',
              'copy',
              '-o',
              'output.avi')
    subprocess.check_call(command)
    print 'movie is written'
    command=(
    'rm',
    '-f',
    '*.png'
    )
    subprocess.check_call(command)
    print 'All png files are removed'
    command=(
        'mplayer',
        'output.avi')
    subprocess.check_call(command)
    print 'now the movie is playing'
    
                    
    
    
    
    

