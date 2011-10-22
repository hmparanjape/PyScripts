"""
VPSC post-process script on analysis on gamdot(is,igr)
Created on July 2011
Author: Youngung Jeong
"""

import glob
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os

#gammadot class
class GammaDot:
    def __init__(self):
        """
        The wildcard of the files is
        fixed to be 'gam_*.out'
        """
        filenames = 'gam_*.out'
        filenames = glob.glob(filenames)
        ## sorting the list of filenames!
        filenames.sort()
        self.ngr, self.nsm = self.__info__(filename=filenames[0])
        self.nstep = len(filenames)
        print 'Current gam_*.out files summary'
        print 'Number of grain: %i'%self.ngr
        print 'Number of slip mode; %i'%self.nsm
        print 'Number of deformation steps: %i'%self.nstep
        self.masterdata = np.zeros(
            (self.ngr, self.nsm, self.nstep))
        for istp in range(self.nstep):
            print 'istp: %i'%(istp+1)
            data = np.loadtxt(filenames[istp], skiprows=2)
            for ism in range(self.nsm):
                for igr in range(self.ngr):
                    self.masterdata[
                        igr,ism,istp] = data[igr,ism]
                    pass
                pass
            pass
        pass
    
    def __info__(self,filename):
        """
        Returns information of the resulting gam_00000.out file
        """
        block = np.loadtxt(
            filename,skiprows=2)
        ngr = block.shape[0] #number of grains
        nsm = block.shape[1] #number of slip modes
        return ngr, nsm
    
    def grain(self,igr,sm=None):
        """
        Return grain igr's gamdot

        if sm is not given, the data will be 2D (sm, stp)
        else it will be 1D (stp)
        """
        if sm==None: data = np.zeros(
            (self.nsm, self.nstep))
        else: data = np.zeros((self.nstep))
        for istp in range(self.nstep):
            if sm==None:
                for ism in range(self.nsm):
                    data[ism,istp] = self.masterdata[
                        igr,ism,istp]
                    pass
                pass
            else:
                data[istp] = self.masterdata[igr,sm,istp]
                pass
            pass
        return data
    pass #end of class GamDot

def gam_gr(igr=0, relative=False, gamdot=None):
    """
    Returns gam
    the variation of the slip system of a grain
    with respect to deformation step

    argument: igr, relative=False, gamdot: a GamDot class

    'data[i] = data[i]/sum(np.abs(data[i]))' is
    returned if relative is True
    """
    data = gamdot.grain(igr) # (ism, istp)
    if relative: #if relative value is required
        data = data.transpose() #(istp, ism)
        for i in range(len(data)): # loop over step
            data[i] = data[i]/sum(np.abs(data[i]))
            pass
        data = data.transpose() #(ism, istp)
        pass
    return data

def gr_plot(igr=None, ngr=1, relative=False,
            ifig=1, figname='temp.pdf', ism=None, nsm=None):
    """
    Plot the gammadot changes and saves it.
    Arguments:
    igr=None, ngr=1, relative=False,
    ifig=1, figname='temp.pdf', ism=None, nsm=None
    """
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm

    try: plt.ioff() # turn of the interactive
    except TclError:
        print 'Tcl Error happened'
        
       
    fig = plt.figure(ifig)
    ax = fig.add_subplot(111)
    
    cmap = ListedColormap(['r','g'])
    norm = BoundaryNorm([-1,0,1], cmap.N)
    
    myGamDot = GammaDot()
    y = []

    ## grain pickup
    if igr==None:
        ## randomly picked up ngr number of grains
        ## and plot them in the same figure
        igr = []
        for i in range(ngr):
            igr.append(
                np.random.randint(
                    low=0, high=myGamDot.ngr)
                )
            pass
        pass
    else: ngr = len(igr)
    ## grain pickup ends

    ## slip mode pickup
    if ism==None:
        ism = []
        if nsm==None: nsm = myGamDot.nsm
        for i in range(nsm):
            ism.append(
                np.random.randint(
                    low=0, high=nsm)
                )
            pass
        pass
    else:
        if any(ism[i] not in range(myGamDot.nsm)):
            print "Requested slip mode %i is not in the list"%ism[i]
            print "Possible slip mode is as below"
            print range(myGamDot.nsm)
            pass
        else: nsm = len(ism)
        pass
    ## slip mode pickup ends
    
    maxy,miny = 0,0
    maxx,minx = 0,0
    for i in range(len(igr)):
        aline = gam_gr(
            igr=igr[i], relative=relative,
            gamdot=myGamDot)
        y.append(aline) #y[igr, ism, istp]
        aline = np.array(aline)
        for j in ism: #each slip mode
            x = np.arange(len(aline[j]))
            points = np.array([x,aline[j]]).T.reshape(-1,1,2)
            segments = np.concatenate([points[:-1],
                                       points[1:]], axis=1)
            lc = LineCollection(segments, cmap=cmap, norm=norm)
            lc.set_array(aline[j])
            lc.set_linewidth(1.)
            ax.add_collection(lc)
            if max(aline[j])>maxy: maxy = max(aline[j])* 1.1
            if min(aline[j])<miny: miny = min(aline[j])* 1.1
            if max(x)>maxx: maxx = max(x)*1.1
            if min(x)<minx: minx = min(x)*1.1
            pass
        pass
    ax.set_ylim(miny, maxy)
    ax.set_xlim(minx, maxx)
    ax.set_ylabel(r'$\dot{\gamma}$'); ax.set_xlabel('step')
    fig.savefig(figname)
    print '%s saved'%figname
    return myGamDot, igr, ism



