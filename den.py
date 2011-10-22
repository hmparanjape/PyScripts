"""
Density
post-processing module
"""
import glob
import numpy as np
class Densities:

    def __init__(self):
        """
        The wild card of the files is
        fixed to be 'denf_*', 'denr_*', and 'dent_*'
        """
        
        ffilenames = glob.glob('denf_*.out')
        rfilenames = glob.glob('denr_*.out')
        tfilenames = glob.glob('dent_*.out')
        ffilenames.sort()
        rfilenames.sort()
        tfilenames.sort()
        self.ngr, self.nsm = self.__info__(filename=ffilenames[0])
        self.nstep = len(ffilenames)
        self.fmd = np.zeros( # forward density master data
            (self.ngr, self.nsm, self.nstep))
        self.rmd = np.zeros( # forward density master data
            (self.ngr, self.nsm, self.nstep))
        self.tmd = np.zeros( # forward density master data
            (self.ngr, self.nsm, self.nstep))



        for istp in range(self.nstep):
            fdata = np.loadtxt(ffilenames[istp])
            rdata = np.loadtxt(rfilenames[istp])
            tdata = np.loadtxt(tfilenames[istp])
            for ism in range(self.nsm):
                for igr in range(self.ngr):
                    self.fmd[igr,ism,istp] = fdata[igr,ism]
                    self.rmd[igr,ism,istp] = rdata[igr,ism]
                    self.tmd[igr,ism,istp] = tdata[igr,ism]
                    pass
                pass
            pass
        pass

    def __info__(self, filename):
        block = np.loadtxt(filename)
        ngr = block.shape[0]
        nsm = block.shape[1]
        return ngr, nsm

    def grain(self,igr,sm=None):
        """
        Return grain igr's denf, denr, dent
        """
        if sm==None: fdata = np.zeros(
            (self.nsm, self.nstep))
        else: fdata = np.zeros((self.nstep))
        rdata = fdata.copy()
        tdata = fdata.copy()
        
        for istp in range(self.nstep):
            if sm==None:
                for ism in range(self.nsm):
                    fdata[ism,istp] = self.fmd[igr,ism,istp]
                    rdata[ism,istp] = self.rmd[igr,ism,istp]
                    tdata[ism,istp] = self.tmd[igr,ism,istp]
                    pass
                pass
            else:
                fdata[istp] = self.fmd[igr,sm,istp]
                rdata[istp] = self.rmd[igr,sm,istp]
                tdata[istp] = self.tmd[igr,sm,istp]
                pass
            pass
        return fdata, rdata, tdata

    
    def __info__(self, filename):
        block = np.loadtxt(filename)
        try:
            ngr = block.shape[0]
            nsm = block.shape[1]
        except:
            raise IOError, 'Need more than 1 grain'
        return ngr, nsm
    pass # end of the class Densities

def den_gr(igr=0, den=None):
    fd, rd, td = den.grain(igr)
    return fd, rd, td
def gr_plot(igr=None, ngr=1, ifig=1,
            figname='gr_densities.pdf',
            ism=None, nsm=None):
    """
    Plot the evolution of dislocatin densities
    """
    import matplotlib.pyplot as plt
    try: plt.ioff()
    except TclError:
        print 'Tcl Error Happened'
        pass
    fig = plt.figure(ifig)
    fig2 = plt.figure(ifig+1)
    fig.clf(); fig2.clf()
    ax = fig2.add_subplot(111)
    ax00 = fig.add_subplot(231)
    ax01 = fig.add_subplot(232)
    ax02 = fig.add_subplot(233)
    ax10 = fig.add_subplot(234)
    ax11 = fig.add_subplot(235)
    ax12 = fig.add_subplot(236)
    
    myDen = Densities()
    y = []

    if igr==None:
        igr = []
        if ngr==1 and myDen.ngr==1: igr=[0]
        else:
            for i in range(ngr):
                igr.append(
                    np.random.randint(
                        low=0, high=myDen.ngr)
                    )
                pass
            pass
        pass
    else: ngr = len(igr)
    
    if ism==None:
        ism = []
        if nsm==None: nsm = myDen.nsm
        if nsm==1 and myDen.nsm==1: ism=[0]
        else:
            for i in range(nsm):
                ism.append(
                    np.random.randint(
                        low=0, high=nsm)
                    )
                pass
            pass
        pass
    else:
        if any(ism[i] not in range(myDen.nsm)):
            print "Requested slip mode %i is not in the list"%ism[i]
            print "Possible slip mode is as below"
            print range(myDen.nsm)
            pass
        else: nsm = len(ism)
        pass
    maxy, miny = 0, 0
    yf = []
    yr = []
    yt = []
    for i in range(len(igr)):
        alinef, aliner, alinet = den_gr(
            igr=igr[i], den = myDen)
        yf.append(alinef)
        yr.append(aliner)
        yt.append(alinet)
        alinef = np.array(alinef)
        aliner = np.array(aliner)
        alinet = np.array(alinet)

        for j in ism:
            x = np.arange(len(alinef[j]))
            ax00.semilogy(x,alinef[j], '-')
            try:
                ax01.semilogy(x,aliner[j], 'o')
                pass
            except:pass
            ax02.semilogy(x,alinet[j], '+')
            ax10.plot(x,alinef[j], '-')
            ax11.plot(x,aliner[j], 'o')
            ax12.plot(x,alinet[j], '+')
            ax.semilogy(x,alinef[j],'-',color='r')
            try:
                ax.semilogy(x,aliner[j],'o',color='b')
                pass
            except: pass
            ax.semilogy(x,alinet[j],'+',color='k')
        pass

    fig.savefig(figname)
    fig2.savefig('%s_%i.pdf'%(figname.split('.pdf')[0],2))
    print 'grains',igr
    print 'slimodes',
    print '%s saved'%figname
    return myDen, igr, ism

def gr_plotavg(figname='gr_densities_avg.pdf', ifig=3):
    """
    Plotting for average dislocation densities through all grains
    and their slip modes.
    """
    import matplotlib.pyplot as plt
    try: plt.ioff()
    except TclError:
        print 'Tcl Error, try to enable X11'
        return -1
    fig = plt.figure(ifig)
    ax = fig.add_subplot(111)
    myDen = Densities()
    yf, yr, yt = [], [], []
    df = np.zeros((myDen.nstep,))
    dr = df.copy()
    dt = df.copy()
    for i in range(myDen.ngr):
        f, r, t = den_gr(igr=i, den=myDen)
        ayf = np.array(f)
        ayr = np.array(r)
        ayt = np.array(t)
        for j in range(myDen.nsm):
            for k in range(myDen.nstep):
                df[k] = df[k] + ayf[j][k]
                dr[k] = dr[k] + ayr[j][k]
                dt[k] = dt[k] + ayt[j][k]
                pass
            pass
        pass
    df = df/myDen.nsm/myDen.ngr
    dr = dr/myDen.nsm/myDen.ngr
    dt = dt/myDen.nsm/myDen.ngr
    ax.semilogy(df,'-')
    ax.semilogy(dr,'o')
    ax.semilogy(dt,'+')
    fig.savefig(figname)
    print '%s saved'%figname
    return df,dr,dt

    
    
