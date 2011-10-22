import os
os.sys.path.append('/home/youngung/Documents/py') ## path
import vp_f2py
reload(vp_f2py)
from vp_f2py import *
import numpy as np; import numpy
import scipy as sp; import scipy
import matplotlib.pyplot as plt
import vpsc_in, sx_maker, chg_basis
import glob, shutil
import math, crss, euler, pf
from vpsc_param import __texf__
from vpsc_param import __makenp__
from vpsc_param import fout
from vpsc_param import fplot
from vpsc_param import __interpolate__
import scipy.integrate as integrate

eul = euler.euler
pf = pf.pf
v_in = vpsc_in.vpsc_in
chg_basis=chg_basis.chg_basis
norm = numpy.linalg.norm

### examples ##
# 1. in-plane variation of the yield stress and R-value
# 2. multipath : RD-TD-DD-SS-BU  stress-strain curves and R-value vs strain
# 2-2 in_plane_var: In-plane variation of the stress-strain curves and R-values..
# 3. What FB requested on the inverse pole figure observation
# 4. Example 4 : ICOTOM paper (simulative work on hardening tendency
#                              influenced by slip system couplings)
# 5. in-plane stress work-contour construction
#    (counterpart to in-plane biaxial testing results)
# 6. simple-shear cyclic test (ICOTOM 16 - 2)
# 7. YS surface using parallel python package (not completed)
# 8. Optimization of the VOCE hardening parameters using SciPy opimization package
# 9. Checking the dislocation density based hardening implementation (2011-July~)
# 10. Run with RVES having different number of grains to fit the experiments
#   to find the best linerization scheme. (response to the MSMSE2011 reviews)

###
### Example 1
###   In-plane variation of the yield stress and R-value
###

"""
>>> myprob = vp_f2py.prob(texture='texture/00500.cmb',
                          fsx='sx/fcc.sx',
                          dang=45.
                          th_0=0., th_1=180.0001,
                          eqc=0.005, stp=1, ict=7,
                          interaction=3)
>>> myprob.__pp__(); myprob.__write__();
    --> __write__ results in YSprob, Rprob, prob files
       in prob the order of columns is ang, r, and ys

    if the above were performed in which an plt display is
    enabled, one can plot the results directly with using
    myprob.__plot__(ifig=3) with ifig being the id of the
    figure.
"""

def __inplanerotation__(ang=0., gr=None):
    """
    Rotates the polycrystalline aggregate sample 
    """
    th = ang * np.pi / 180.
    sth = np.sin(th) ; cth = np.cos(th)
    rmat = np.array([[cth, sth, 0.],[-sth, cth, 0.],[0., 0., 1.]])
    newgr = gr.copy()
    for igr in range(len(gr)):
        crot = eul(ph=gr[igr][0], th=gr[igr][1], tm=gr[igr][2], echo=False) #convert it into rotation matrix
        newrot = np.dot(crot, rmat)
        ph, th, tm = eul(a=newrot, echo=False)
        newgr[igr][0] = ph
        newgr[igr][1] = th
        newgr[igr][2] = tm
        newgr[igr][3] = gr[igr][3] #volume fraction
        pass
    return newgr

def __writegrtofile__(gr=None, fn=None):
    """
    write polycrystal aggregate into a file
    """
    f = open(fn, 'w')
    f.writelines('dummy\ndummy\ndummy\nB %i\n'%len(gr))
    for i in range(len(gr)):
        f.writelines('%13.9f %13.9f %13.9f %13.8e\n'%(
                gr[i][0], gr[i][1], gr[i][2], gr[i][3]))

class prob:
    """
    class probing for in-plane variation of YS and R
    
    The hardening flag is set to be [0,0,0,0] so that
    none of the state is updated after any loading process.

    Arguments:
        dang = 45
        texture
        fsx
        th_0=0.0
        th_1= 180.001
        eqc=0.05,
        stp=1,
        ict=7,
        interaction=3,
        ifig=3
    """
    def __init__(self, dang=45, texture=None, fsx=None,
                 th_0=0., th_1=180.0001,
                 eqc=0.005, stp = 1, ict=7, interaction=3,
                 ifig=3,jobid=1):
        if interaction==0: 
            print " Cannot Impose FC"; raise IOError
        self.prcs=[]
        angs = np.arange(th_0, th_1, dang)
        self.ang = angs
        vpsc_in.histmaker(filename='hist/temp0', nstep=stp,
                          ictrl=7, eqincr=eqc,)

        vpsc_in.rotfile(filename='temp0',ang = dang)
        for i in range(len(angs)):
            # """
            # print "*******************"
            # print "R_prob's job #%i"%i
            # print "PROCESS PREPARATION"
            # print "*******************"
            # """
            self.prcs.append('0')
            self.prcs.append('hist/temp0')
            if i==len(angs)-1: pass #No last rotation
            else:
                self.prcs.append('4')
                self.prcs.append('rot/temp0')
                
        self.job = vpsc(texture=texture, fsx=fsx, 
                        iupdate=[0,0,0,0], prcs = self.prcs,
                        interaction = interaction,
                        jobid=jobid)
        #self.__run__()
        #if ifig==None: pass
        #else: self.__plot__(ifig)
        #print "vpsc.so is completed"
    def __run__(self,):
        """ run! """
        self.job.run()
        self.__pp__()
        return self.ang, self.ys, self.r, self.ys/self.ys[0]

    def __pp__(self,):
        self.job.pp()
        self.ys = []
        self.r = []
        for i in range(len(self.prcs)/2):
            if self.prcs[i*2]!='0': pass
            else:
                YS=self.job.datamaster['scau'][i][0][0,0]
                ew=self.job.datamaster['dbar'][i][0][1,1]
                et=self.job.datamaster['dbar'][i][0][2,2]
                R = ew/et
                self.ys.append(YS)
                self.r.append(R)
        return self.job.datamaster

    def __write__(self,r_filename='Rprob',
                  ys_filename='YSprob',
                  filename='prob'):
        """
        Writes the ang-R, ang-YS, and ang-R-YS, respectively.
        """
        fout(r_filename,self.ang,self.r)
        fout(ys_filename,self.ang,self.ys)
        fout(filename, self.ang, self.r, self.ys) ## R and YS altogether
        pass
    
    def __plot__(self,ifig=None):
        """
        plot
        """
        fig = plt.figure(ifig)
        ax = fig.add_subplot(111)
        self.nys = self.ys/self.ys[0]
        ax.set_title('R and YS profile') #R
        ##
        ax.set_xlim(-5,185)
        ax.set_xticks(np.arange(0.,180.001,45.))
        ax.grid('on')
        ax_ = ax.twinx() #YS
        ax_.plot(self.ang,self.nys)
        ax_.plot(self.ang,self.r)

        ### deco
        ax_.set_xlabel(r'$\theta$'+'from +RD', 
                       dict(fontsize=20))
        ax_.set_ylabel('Normalized YS',
                       dict(fontsize=20))
        ax.set_ylabel('R value',
                      dict(fontsize=20))

def probplot(ifig):
    """
    """
    fig = plt.figure(ifig)
    ax = fig.add_subplot(111)
    ang, y = fin('prob',1,0,1,2)
    #filename, nhead, ix, *args for y
    ang = np.array(map(float,ang))
    r = np.array( map(float,y[0]))
    ys = np.array(map(float,y[1]))
    nys = ys/ys[0]
    #nys = self.ys/self.ys[0]
    ax.set_title('R and YS profile') #R
    ax.grid('on')
    ax_ = ax.twinx() #YS
    ax.plot(ang, r, '--')
    ax_.plot(ang, nys,'-')
    ### deco
    ax.set_xlabel(r'$\theta$', dict(fontsize=20))
    ax_.set_ylabel('Normalized YS', dict(fontsize=20))
    ax.set_ylabel('R value', dict(fontsize=20))
    ax.set_xlim(-5,185)
    ax.set_xticks(np.arange(0.,180.001,45))

def prob_files(wildcard='texture/*.cmb',
               fsx=None, dang=45.):
    """
    """
    files = glob.glob(wildcard)
    jobs = []
    for i in range(len(files)):
        jobs.append(
            prob(texture=files[i],
                 fsx=fsx, dang=dang)
            )
        pass
        
    for i in range(len(jobs)):
        jobs[i].__pp__()
        jobs[i].__write__(r_filename='%s%s'%
                          ('Rprob',str(i).zfill(3)),
                          ys_filename='%s%s'%
                          ('YSprob',str(i).zfill(3)),
                          filename='%s%s'%
                          ('prob',str(i).zfill(3))
                          )
        pass
    pass


def examples1(files=None, dang=None, fsx=None, interaction=None,
              tex=None, ifig=3, labels=None, iint=False,jobid=1):
    """
    In-plane variation of YS and R with changing the number of grains...
    in prepartion of MSMSE submission.
    
    --Arguments-------------------------------------
    files, dang, fsx, interaction, tex, ifig, labels
    ------------------------------------------------


    if iint==True:
      in-plane variation over different linearization schems...
    else:
      over the files given 

    example:
    >>> reload(vpsc_examples);
    d=vpsc_examples.examples1(
      dang=5,files=['304_surf_00500.cmb',
        '304_0.25t_00500.cmb',
        '304_0.5t_00500.cmb'], fsx='304.sx', interaction=3,
        ifig=1, labels=['surface',r'$\frac{1}{4}$t', r'$\frac{1}{2}$t'])

    >>> reload(vpsc_examples);d=vpsc_examples.examples1(
    dang=5,files=['304_surf_00500.cmb'], fsx='304.sx',
    ifig=1, labels=['surface',r'$\frac{1}{4}$t', r'$\frac{1}{2}$t'])
    """
    l = 0.15  #left blank
    b = 0.15  #bottom blank
    w = 0.75  #width
    h = 0.80  #height
    nax = 1.
    
    fig_width=7
    fig_height=5
    fig_rel_size= 0.9 
        
    fig_width, fig_height = np.array(
        [fig_width,fig_height]) * fig_rel_size

    fig  = plt.figure(ifig,   [fig_width*nax, fig_height])
    fig1 = plt.figure(ifig+1, [fig_width*nax, fig_height])
    el = l/nax
    ew = w/nax
    fullscale = 1.0
    incr = fullscale / nax
    
    ax0 = fig.add_axes((el,      b,ew,h))
    ax1 = fig1.add_axes((el     ,b,ew,h))

    markers = ['o','^','d','p','h','*','+',
               'o','^','d','p','h','*','+']
    linestyle = ['-','--','-.',':','-','--',
                 '-.',':','-','--','-.',':',]
    if iint==False:
        for i in range(len(files)):
            a = prob(
                texture = files[i],
                dang = dang, stp = 10,
                interaction = interaction,
                ifig=None, jobid=jobid)
            jobid = jobid + 1
            ang, ys, r, nys = a.__run__()
            
            ax0.plot(ang, r, color='black', mfc='None',
                     ls=linestyle[i],
                     marker = markers[i],
                     label = labels[i])
            ax1.plot(ang, nys, color = 'black', mfc = 'None',
                     ls=linestyle[i],
                     marker = markers[i],
                     label=labels[i])
            
            ax0.set_xlabel(r'$\theta$ from RD', dict(fontsize=20))
            ax0.set_ylabel(r'$R_{\theta}$', dict(fontsize=17))        
            ax1.set_xlabel(r'$\theta$ from RD', dict(fontsize=20))
            ax1.set_ylabel(r'$\sigma^{YS}_{\theta}/\sigma^{YS}_{\theta=0^{\circ}}$',
                           dict(fontsize=20))
            ax0.set_xticks(np.arange(0.,180.01,15))
            ax1.set_xticks(np.arange(0.,180.01,15))
            ax0.grid(); ax1.grid()

            pass
        pass
    elif iint==True:
        interaction=[1,2,3,4]
        labels=['affine','secant',
                r'$n^{eff}=10$','tangent']
        for i in range(len(interaction)):
            a = prob(
                texture = files[0],
                dang = dang, stp = 10,
                interaction = interaction[i],
                ifig=None, jobid=jobid)
            jobid = jobid + 1
            ang, ys, r, nys = a.__run__()
            ax0.plot(ang, r, color='black', mfc='None',
                     ls=linestyle[i],
                     marker = markers[i],
                     label = labels[i])
            ax1.plot(ang, nys, color = 'black', mfc = 'None',
                     ls=linestyle[i],
                     marker = markers[i],
                     label=labels[i])
            
            ax0.set_xlabel(r'$\theta$ from RD', dict(fontsize=20))
            ax0.set_ylabel('R-value', dict(fontsize=17))        
            ax1.set_xlabel(r'$\theta$ from RD', dict(fontsize=20))
            ax1.set_ylabel(r'$\sigma^{YS}/\sigma^{YS}_{RD}$ from RD',
                           dict(fontsize=20))
            
            ax0.set_xticks(np.arange(0.,180.01,15))
            ax1.set_xticks(np.arange(0.,180.01,15))
            ax0.grid(); ax1.grid()

            pass        
        pass
    ax0.legend(loc='best'); ax1.legend(loc='best')
        
###
### using pp  # still not completed (2011-02-21)
###
try: import pp
except: print "pp was not found"
else:
    print 'parellel python has been successfully imported'
    pass

def ppex():
    #job1 = vpsc(texture='texture/08000.cmb', mode='lankf')
    jobs = []
    jobs.append(prob(texture='texture/00500.cmb', dang=5.))
    jobs.append(prob(texture='texture/01000.cmb', dang=5.))
    jobs.append(prob(texture='texture/02000.cmb', dang=5.))    
    jobs.append(prob(texture='texture/04000.cmb', dang=5.))
    jobs.append(prob(texture='texture/08000.cmb', dang=5.))

    # jobs.append(vpsc(texture='texture/00500.cmb', mode='lankf'))
    # jobs.append(vpsc(texture='texture/00500.cmb', mode='ten_ang',ang=5.))
    # jobs.append(vpsc(texture='texture/00500.cmb', mode='ten_ang',ang=15.))
    # jobs.append(vpsc(texture='texture/00500.cmb', mode='ten_ang',ang=10.))

    job_server = pp.Server()
    results =[]
    for i in range(len(jobs)):
        #results.append(job_server.submit(jobs[i].run)())
        results.append(  job_server.submit(jobs[i].__run__)()  )
    job_server.wait()
    #f1 = job_server.submit(job2.run)
    #f2 = job_server.submit(job1.pp)
    #f3 = job_server.submit(job2.pp)
    job_server.print_stats()

    for i in range(len(jobs)):
        crs = results[i]
        fout('prob_pp_%s'%str(i).zfill(3),crs[0],crs[1],crs[2])
    return results
    
###
### example 2
### RD-TD-DD  str-str  r-str
###           stress-work R-work
class multipath:
#     def __init__(self, texture ='texture/00500.cmb',
#                  stp=3, eqc=0.005, ict=7):
#         self.modes = ['RD','DD','TD','SS','BU']
    def __init__(self,texture='texture/00500.cmb',
                 stp=3, eqc=0.005, ict=7,fsx=None,
                 modes=['RD','DD','TD','SS','BU'] ):
        """
        Assigns multiple jobs into myjobs object.
        End-user can add more jobs to the down below block.

        Arugments
        texture ='texture/00500.cmb'
        * args  : modes for class vpsc  'RD','DD','TD','SS'...
        modes = ['RD','DD','TD','SS','BU']
        stp = 3 # of deformation step
        eqc = 0.005 # equivalent incremental strain
        ict = 7 #strain control
        """
        #flag indicating whether the pp was performed
        self.waspp = False 
        
        ## hardwired modes are used.
        #self.modes = ['RD','DD','TD','SS','BU'] 
        self.modes = modes
        self.myjobs = []

        ## Shoots Information on to the screen
        print "\n***************************************************"
        print "You have total %i number"%len(self.modes),
        print " of different loading paths"
        print "%3s %15s"%('id','mode name')
        for i in range(len(self.modes)):
            print "%3i %15s"%(i,self.modes[i])
        print "***************************************************"

        for i in range(len(self.modes)):
            """
            Unlike vp_param.vpsc class
            vp_f2py.vpsc class needs the post-processing right after
            each run is finished since it, somehow, has its local
            variables initialized and does not save it independently
            as of 2011-02-21. 
            """
            self.myjobs.append(
                vpsc(mode=self.modes[i],
                     texture=texture, fsx=fsx,
                     stp=stp, eqc=eqc, ict=ict,jobid=i)
                )
            self.myjobs[i].run()
        self.waspp = True
        pass

    def __pp__(self):
        """
        Defines work and calculates work using cumtrapz
        """
        self.work = []
        self.sig = []
        self.eps = []
        self.edot = []
        for i, mode in zip(range(len(self.modes)), self.modes):
            #epstot from vpsc7 f2py module
            eps = self.myjobs[i].datamaster['eps'][0]  
            sig = self.myjobs[i].datamaster['sbar'][0]
            edot = self.myjobs[i].datamaster['dbar'][0]
            if mode=='DD' or mode=='dd' or mode=='TD' or mode=='td':
                eps = self.myjobs[i].datamaster['eps'][1]
                sig = self.myjobs[i].datamaster['sbar'][1]
                edot = self.myjobs[i].datamaster['dbar'][1]
            eqincr = self.myjobs[i].eqincr
            cdat = [eps, sig, edot]
            #return eps, sig

            if any(cdat[i]==None for i in range(len(cdat))):
                pass #no wonder. it is supposed to be 
                     #a rigid rotation process
            else:
                wrk = np.zeros((3,3))
                E, S = [], []
                for stp in range(len(eps)):
                    ceps = eps[stp].flatten()
                    csig = sig[stp].flatten()
                    cbar = edot[stp].flatten()
                    E.append(ceps)
                    S.append(csig)
                    pass

                E, S = np.array(E), np.array(S) #each in dimension of (stp,9)
                cwrk = []
                E = E.transpose()
                S = S.transpose()  # each in dimension of (9,stp)

                for icmp in range(len(E)): #for each component
                    temp_work = integrate.cumtrapz(x=E[icmp], y=S[icmp])
                    temp_work = np.insert(temp_work,[0],[0])
                    cwrk.append(temp_work)
                    pass
                cwrk = np.array(cwrk)
                cwork = []
                cwrk = cwrk.transpose() # now cwrk in (stp,9)
                #for each step of deformation
                for istp in range(len(E.transpose())): 
                    cwork.append(cwrk[istp].sum())
                    pass
                self.work.append(cwork)
                self.sig.append(sig)
                self.eps.append(eps)
                self.edot.append(edot)
                pass
            pass

        self.work = np.array(self.work)
        self.sig = np.array(self.sig)
        self.eps = np.array(self.eps)
        self.edot = np.array(self.edot)

    def __ppwrite__(self,iplot=True, ifig=94):
        """
        Writes down the important results in a relevant format
        for each loading case

        + and plotts!
        if iplot==True
        """

        l = 0.15  #left blank
        b = 0.15  #bottom blank
        w = 0.75  #width
        h = 0.80  #height 

        fig_width=7
        fig_height=5
        fig_rel_size= 0.9 
        
        fig_width, fig_height = np.array([fig_width,fig_height]) * fig_rel_size


        markers = ['o','^','d','p','h','*','+']
        
        if iplot==True:
            #fig = plt.figure(ifig, [fig_width *len(self.modes),fig_height])
            fig = plt.figure(ifig,   [fig_width,fig_height])
            fig1 = plt.figure(ifig+1, [fig_width,fig_height])
            # each figure's length and width
            el = l
            ew = w
            fullscale = 1.0

            ax0 = fig.add_axes((el,b,ew,h))
            ax1 = fig1.add_axes((el,b,ew,h))
            
            # for i in range(len(self.modes)):
            #     current_left = el + incr * i
            #     current_b    = b
            #     current_wid  = ew #+ incr * i                
            #     current_h    = h
            #     axes.append(
            #         fig.add_axes(
            #             (current_left, current_b,
            #              current_wid, current_h),label='ax #%i'%i
            #             )
            #         )
            #     pass
            pass
        
        ## writes into files. Work, sig, eps
        for i, mode in zip(range(len(self.modes)), self.modes):
            ## finds the effective stress and strain
            sig__ = self.sig[i].copy()
            eps__ = self.eps[i].copy()
            wrk__ = self.work[i].copy()
            edt__ = self.edot[i].copy()

            sig11 = []; sig12 = []; sig13 = []
            sig21 = []; sig22 = []; sig23 = []
            sig31 = []; sig32 = []; sig33 = []
            eps11 = []; eps12 = []; eps13 = []
            eps21 = []; eps22 = []; eps23 = []
            eps31 = []; eps32 = []; eps33 = []
            edt11 = []; edt12 = []; edt13 = []
            edt21 = []; edt22 = []; edt23 = []
            edt31 = []; edt32 = []; edt33 = []                        

            for j in range(len(sig__)):
                sig11.append(sig__[j][0,0])
                sig12.append(sig__[j][0,1])
                sig13.append(sig__[j][0,2])
                sig21.append(sig__[j][1,0])
                sig22.append(sig__[j][1,1])
                sig23.append(sig__[j][1,2])
                sig31.append(sig__[j][2,0])
                sig32.append(sig__[j][2,1])
                sig33.append(sig__[j][2,2])

            sig11 = np.array(sig11)
            sig12 = np.array(sig12)
            sig13 = np.array(sig13)
            sig21 = np.array(sig21)
            sig31 = np.array(sig31)
            sig32 = np.array(sig32)
            sig33 = np.array(sig33)

            for j in range(len(eps__)):
                eps11.append(eps__[j][0,0])
                eps12.append(eps__[j][0,1])
                eps13.append(eps__[j][0,2])
                eps21.append(eps__[j][1,0])
                eps22.append(eps__[j][1,1])
                eps23.append(eps__[j][1,2])
                eps31.append(eps__[j][2,0])
                eps32.append(eps__[j][2,1])
                eps33.append(eps__[j][2,2])

            eps11 = np.array(eps11)
            eps12 = np.array(eps12)
            eps13 = np.array(eps13)
            eps21 = np.array(eps21)
            eps31 = np.array(eps31)
            eps32 = np.array(eps32)
            eps33 = np.array(eps33)


            for j in range(len(edt__)):
                edt11.append(edt__[j][0,0])
                edt12.append(edt__[j][0,1])
                edt13.append(edt__[j][0,2])
                edt21.append(edt__[j][1,0])
                edt22.append(edt__[j][1,1])
                edt23.append(edt__[j][1,2])
                edt31.append(edt__[j][2,0])
                edt32.append(edt__[j][2,1])
                edt33.append(edt__[j][2,2])

            edt11 = np.array(edt11)
            edt12 = np.array(edt12)
            edt13 = np.array(edt13)
            edt21 = np.array(edt21)
            edt31 = np.array(edt31)
            edt32 = np.array(edt32)
            edt33 = np.array(edt33)            


            if mode=='rd' or mode=='RD':
                sig__ = sig11 - sig33
                eps__ = eps11
                r = edt22/edt33
                fout('%s%s'%('wrk_str_eps',mode),
                     wrk__, sig__,eps__, r)

            elif mode=='td' or mode=='TD':
                sig__ = sig11 - sig33
                eps__ = eps11
                r = edt22/edt33
                fout('%s%s'%('wrk_str_eps',mode),
                     wrk__, sig__,eps__,r)
            elif mode=='dd' or mode=='DD':
                sig__ = sig11 - sig33
                eps__ = eps11
                r = edt22/edt33                
                fout('%s%s'%('wrk_str_eps',mode),
                     wrk__, sig__,eps__,r)
            elif mode=='bu' or mode=='BU':
                sig__ = sig33 - sig11
                sig__ = abs(sig__)
                fout('%s%s'%('wrk_str_eps',mode),
                     wrk__, sig__)
            elif mode=='ss' or mode=='SS':
                sig__ = sig12
                fout('%s%s'%('wrk_str_eps',mode),
                     wrk__, sig__)
            else: 
                print "Unexpected mode is passed to __ppwrite__"
                raise IOError

            if iplot==True:
                ax0.plot(eps__, sig__,markers[i],
                         color='black',
                         label='%s'%mode, mfc='None')
                ax1.plot(eps__, r, markers[i],
                         color='black',
                         label='%s'%mode, mfc='None')
                pass

            ax0.set_xlabel(r'$\varepsilon^{pl}$', dict(fontsize=20))
            ax0.set_ylabel(r'$\sigma$ [MPa]', dict(fontsize=20))
            #ax0.set_ylim(0.,)
            ax1.set_xlabel(r'$\varepsilon^{pl}$', dict(fontsize=20))
            ax1.set_ylabel(r'R-value', dict(fontsize=15))
            #ax1.set_ylim(0.,)
            ax1.legend(loc='best')
            ax0.legend(loc='best')            
            #fout('%s%s'%('wrk_str_eps',mode), wrk__, sig__)

## plotting defunc for the above class
def wrk_str_plot():
    files = glob.glob('wrk_str_eps*')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in files:
        x,y = fplot(filename = i, ix=0, iy=1, nhead=1)
        ax.plot(x[0],y[0], label = i[11::1])
        ax.set_ylim(0.,)
        pass
    pass


##
## Example 2-2 another multipath (in-plane variation)
##
class in_plane_var:
    """
    >>> import vpsc_examples
    >>> a = vpsc_examples.in_plane_var(texture='texture.cmb', ang=45.,
                                       stp=20, eqc=0.01, smode='indv',
                                       hl=1.4, hb=-0.1, fsx=None,
                                       interaction=3)
    """
    def __init__(self,texture=None, fsx=None, ang=45., eqc=0.01, stp=20,
                 hs=1.0, hl=1.0, hb=1.0, smode='indv', interaction=3):
        """
        ang: incremental angle
        eqc
        stp : deformation step number
        """
        self.dang = ang
        if self.dang >0: self.angs = np.linspace(0., 180.001, 1800./self.dang+1)
        else: self.angs = np.linspace(0.,-90., (-90./self.dang)+1)
        self.jobs = []

        if fsx==None:
            ## bcc case
            tau0=1.; tau1=0.5; thet0=1.0; thet1=0.5            
            self.fsx = 'temp_bcc.sx'
            sx_maker.cubic(filename = self.fsx,  hii=hs, hij=hl, hb=hb,
                           tau0=tau0,  tau1=tau1, thet0=thet0, thet1=thet1,
                           iopsysx=0,
                           dependency=smode, #'iso' or 'indv'
                           b=[[1,1,1],[1,1,1],[1,1,1]], n=[[1,1,0],[1,1,2],[1,2,3]])
            ## fcc case
            # self.fsx = 'temp_fcc.sx'
            # sx_maker.cubic(filename = 'temp.sx',  hii=1.0, hij=1.0, hb=1.0,
            #                tau0=tau0,  tau1=tau1, thet0=thet0, thet1= thet1,
            #                iopsysx=0,
            #                dependency='indv',
            #                b=[[1,1,0]], n=[[1,1,1]])
            pass
        else: self.fsx=fsx
        
        ## rotate the texture file, and save it, and use it.
        self.tout = 'temp0.tex'
        for i in range(len(self.angs)):

            # rotation of the given file and saves it in self.tout
            self.__rot__(infile=texture, outfile=self.tout, ang=self.angs[i])

            # submit the job ---------------------------------------------------------------
            self.jobs.append(
                vpsc(texture=self.tout, stp=stp, eqc=eqc, fsx=self.fsx,
                     interaction=interaction, mode='RD')
                )
            self.jobs[i].run()
            # ------------------------------------------------------------------------------
            pass

        ###
        """post-process"""
        self.__pp__()

    def __pp__(self):
        """
        post-processing
        """
        self.datamaster = {}
        fig = plt.figure(1)
        # ax0 = fig.add_subplot(211)
        # ax1 = fig.add_subplot(212)
        fig = plt.figure(1)        
        ax0 = fig.add_axes((0.15,0.15,0.57,0.8), label=r'$\sigma$ $\varepsilon^{pl}_{axial}$')
        fig = plt.figure(2)
        ax1 = fig.add_axes((0.15,0.15,0.57,0.8), label=r'R vs $\varepsilon_{axial}$')

        for i in range(len(self.jobs)):
            eps = self.jobs[i].datamaster['eps'][0]
            dbar = self.jobs[i].datamaster['dbar'][0]
            sbar = self.jobs[i].datamaster['sbar'][0]
            e11 = []
            s11 = []
            R = []
            for j in range(len(dbar)):
                if j==0: e11.append(0)
                else: e11.append(eps[j-1,0,0])
                r = dbar[j,1,1] / dbar[j,2,2]
                R.append(r)
                s11.append(sbar[j,0,0] - sbar[j,2,2])
            e11, s11, R = __makenp__(e11, s11, R)
            ax0.plot(e11, s11, label=r'$\theta$: %4.1f'%self.angs[i])
            ax1.plot(e11, R, label=r'$\theta$: %4.1f'%self.angs[i])
            pass

        ax0.set_ylim(0.,); ax0.set_xlim(0.,)
        ax1.set_ylim(0.,); ax1.set_xlim(0.,)

        ax0.set_xlabel(r'$\varepsilon^{pl}_{axial}$', dict(fontsize=25))
        ax0.set_ylabel(r'$\sigma_{axial}$',dict(fontsize=25))
        ax1.set_xlabel(r'$\varepsilon^{pl}_{axial}$', dict(fontsize=25))
        ax1.set_ylabel(r'R')
        
        ax0.legend(bbox_to_anchor=(1.05, 1), loc=2)
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2)
        pass

    def __rot__(self, infile=None, outfile='temp.tex', ang=0.):
        """
        Imposes an in-plane rotation(as of an angle denoted as ang)
        to the texture file and saves it elsewhere (outfile)
        """
        if infile==None:raise IOError


        ## save header of the in texture file
        FILE = open(infile, 'r')
        header = FILE.readlines()[0:4]
        FILE.close()
        
        ## rotation        
        aa = np.loadtxt(infile, skiprows=4)
        aa = aa.transpose()
        aa[0] = aa[0] + ang
        for i in range(len(aa[0])):
            if aa[0,i]>360.: aa[0,i] = aa[0,i] - 360.
            else: pass
            pass

        ## writes the rotated texture!
        fout = open(outfile, 'w')
        for i in range(len(header)):
            fout.write(header[i]) #writes the header 
            pass

        aa = aa.transpose()  #transpose the matrix again
        for i in range(len(aa)):
            fout.writelines('%8.3e  %8.3e  %8.3e  %8.3e \n'%
                            (aa[i][0], aa[i][1], aa[i][2], aa[i][3]))
            pass
        fout.close()
        
###
### example 3
### What FB requested#3
###
 # task1 : inv pole figure rendering of the random texture --> became an independent module
 # task2 : and see how they moves in that configuration.
 # task3 : Assign a particular orientation and see if
 #         we have a bifurcation point in between (100) and (111) line.
 
##########################################################################
# Because of this, an object-oriented pole figure                        #
# plotting software has been developed (2011-1-Mar to 2011-5-Mar)        #
# That is saved in the local development folder of my GNU/Linux          #
# Ubuntu system, "/the_ultimate_polefigure". What a childish name!       #
#   --> Long time ago, it was abbreviated as upf (2011-Aug-16)           #
#                                                                        #
# The strategy is simple. Provided the random texture, elongates it      #
# along a certain sample axis, e.g. (001) z-axis. Then tracks of certain #
# inverse poles and measures its rate. The definition of the             #
# crystallographic orientation rate will be like 'theta/strain'.         #
# Here theta is the minimum rotation angle between two points. The angle #
# can be calculated by finding the axis around which the pole is rotate  #
# in the minimum angle to its next destination.                          #
# Until a better name comes up, I would name the class fb3...            #
##########################################################################

import upf  #upf is also a lousy name due to its childish attribute...
import randomEuler # for making a random file
re = randomEuler.randomEuler
from text import miller2euler

class fb3:
    def __init__(self,):
        #initial texture, note that it can increase if
        #additional grain is assigned.
        self.rt = re(ngrain=1000) 
        self.ftx = 'random.tex'
        self.rt.write(self.ftx)

        npf=20
        for i in np.linspace(20.,30., npf):
            self.__addgr__(filename=self.ftx, agr=[i, 44.98,0.])
        
        # Running vpsc ----------------------------------------------
        self.nstp = 30
        self.vp = vpsc(mode='RD', texture=self.ftx, stp=self.nstp, eqc=0.005)
        self.vp.run()
        #------------------------------------------------------------
        
        self.gr = self.vp.gr   #nph/nprc/ step/ ngr/ 4

        mypf = upf.polefigure(grains=self.gr[0][0][0])

        for i in range(len(np.linspace(20.,40., npf))):
            for step in range(self.nstp):
                if step ==0: color='r'
                else : color='k'
                mypf.dotplot(pole=[1,0,0,], ifig=1, proj='ipf', csym='cubic',
                             agrain=self.gr[0][0][step][-1-i], color=color)
            
        # mypf = upf.polefigure(grains=self.gr[0][0][0])
        # mypf.ipf(pole=[1,0,0], color='r')
        # mypf = upf.polefigure(grains=self.gr[0][0][1])
        # mypf.ipf(pole=[1,0,0]) 
        # mypf = upf.polefigure(grains=self.gr[0][0][2])
        # mypf.ipf(pole=[1,0,0])         
        # mypf = upf.polefigure(grains=self.gr[0][0][3])
        # mypf.ipf(pole=[1,0,0]) 
        # mypf = upf.polefigure(grains=self.gr[0][0][4])
        # mypf.ipf(pole=[1,0,0]) 
        
    def __addgr__(self, filename=None, agr=None, wgt=None):
        """
        Adds a specific grain to the existing e texture file.
        Note that, if two many grains, compared to the existing grains,
        are added, the representativeness that it used to have can be
        affected!

        The volume fraction of the grain will be average, unless given.
        vol(gr) = total_vol/total_ngr
        --------
        Argument
        agr = None
        """
        ## find the Euler angle nomenclature
        FILE = open(filename, 'r')
        lines = FILE.readlines()
        nomen = lines[3].split()[0]
        FILE.close()
        
        #add one,
        agg = np.genfromtxt(filename, skiprows=4) #current aggregate
        aggt = agg.copy().transpose() #agg[phi1, phi2, phi3, inten]

        #if wgt is given through wgt as well as agr[3]
        #conflict may occur. Thus, here the convention
        #is using the wgt case
        
        if wgt!=None and len(agr)==4:
            print "weight of grain is given both ways"
            print "As a convention, wgt is being used"
            print "in this case."

        #calculates the avg weight in case wgt isn't given
        if wgt==None and len(agr)!=4:
            wgt_tot = aggt[3].sum()
            wgt_avg = wgt_tot/len(agg)
            wgt = wgt_avg
        elif wgt==None and len(agr)==4: wgt = agr[3]
        elif wgt!=None:                 wgt = wgt
        else:
            print "Something is wrong"
            raise IOError
        
        agr = [agr[0], agr[1], agr[2], wgt] #assign the wgt
        
        grains = []
        for i in range(len(agg)):
            grains.append(agg[i])
        grains.append(agr)
        grains = np.array(grains) #make it numpy-arrayed.
        FILE = open(filename,'w')
        FILE.writelines('dummy\ndummy\ndummy\n')
        FILE.writelines('%s  %6i\n'%(nomen,len(grains)))
        for i in range(len(grains)):
            FILE.writelines('%9.3f   %9.3f   %9.3f  %12.5e\n'%
                            (grains[i][0],grains[i][1],
                             grains[i][2],grains[i][3]))
        pass
###
### Example 4
### For microstructural and crystallographic aspects of yield surface evolution
###  In prepartion for the ICOTOM paper #1 - Crystallographic aspects of yield surface evolution

def cr(bins=20):
    """
    Reads exisiting crss files and process the data
    """
    cr = crss.crss()               #path is None means PWD
    ngr = len(cr.master_source[0]) #number of grains
    nmode = cr.nmode               #number of slip systems
    nstep = cr.nstep               #number of strain increment
    
    masterdd = []
    for g in range(ngr):   #grain
        for s in range(nmode):  #slip system
            masterdd.append(cr.data(igrain=g, imode=s, istep=-1))
    masterdd = np.array(masterdd)
    #print masterdd.shape;raw_input()
    return np.histogram(masterdd, bins=bins)[0]  #normed=True)[0]
    
import randomEuler; reload(randomEuler)
re = randomEuler.randomEuler
import upf; reload(upf)

class YsParam:
    """
    Yield locus parameterization study
    - in-plane variation probing on YS and Lankford-coefficient (R-value)
    - Yield locus evolution with normalized with its direction
    
    Makes use of POSTMORT.IN!

    isave = nstep  :    writes grain states in POSTMOR.OUT at step 'isave'
    irecover = 1 (read) or 0 (not) from POSTMORT.IN

    # Communication through postmorph must envolve the proper
    renaming procedure on POSTMORT.out and POSTMORT.in

    2011-Aug-15~ attempts to add more complex case for prestrain path
    in connection with Rauch et al.'s dislocation hardening.
    (Look at the icase)
    2011-Oct-20~ ilat is added to test latent-hardening(in the sx_maker.py as well)
    
    --> the distortion of yield surface in case of reverse loading?
    """
    def __init__(self, eqc=0.005, nprob=32, stp=5, ncycle=2, interaction=3,
                 ngrain=100, ifig=1, pfmode='contourf',
                 iupdate=[1,1,1,0],
                 hii=1.0, hij=1.0, hb=1.0, hp=1.0,
                 bins=20, dang=15.,
                 mode='RD', u11=None, u22=None,

                 icase=1, #flag on the parametric case study
                 
                 ftex=None, fsx=None, # texture file and single crystal file
                                      # they are complementary to h-parameters
                                      # and ngrain
                 
                 ysproj = '2D', #yield surface projection
                 ihardlaw=0,  # 0: Voce -2:Rauch

                 ## dislocation hardening 'Rauch' case
                 idislmode=1,  #idislmode for all is active(1) or active-inactive mode(0)
                 ilat=0,       #ilat(1): Latent hardening, ilat(0): no interplanar interactions
                 ibau=0,       #ibau(1):sudden drop, ibau(0): no such drop
                 pp=0.8,       #used to drop the sudden density level.
                 iopsysx=0,     #iopsysx(1): a slip mode includes two glide direction (+ and -)
                               #       (0): one slip mode for one glide direction [110] != [-1-10]
                 sxinteraction ='Bauschinger Latent',
                 dependency='indv',  #or iso
                 copl=False,  #Tell coplanar slip mode from other slip modes
                                #to impose a different latent hardening rates.

                 ## associated dislocation hardening parameters
                 grsize=3.0e-05,    #grain size
                 burgers=2.46e-10,  #burgers vector
                 fdisln=180,        #
                 shearm=8.5e+04,    #maximum shear modulous
                 ftau0=38.,         #lattice friction (initial crss)
                 portion=0.8,       #dislocatin exchange portion between denf&denr
                 ftherm=2.8         #thermal coefficient in the recovery term
                 ):
        """
        ---------
        Arguments
        ---------
        
        eqc : equivalent strain increment
        nprob : number of YS probing
        stp : number of strain increment for each cycle
        ncycle : number of (YS/Prestrain) cycles
        interaction=3 #linearization for the power law constitutive model
        ngrain: number of the random texture grains
        ifig=1
        pfmode='contourf' :pole figure plotting mode
        iupdate=[1,1,1,0] # orient, shape, hardening, itran
        hii=1.0 : self-hardening coefficient
        hij=1.0 : latent-hardening coefficient
        hb =1.0 : Bauschinger hardening coefficient
        hp =1.0 : latent-hardening coefficient on coplanar slip systems.
        copl = False: coplanar latent hardening
        bins=20 :Number of bins used in the CRSS distribution curve
        dang=15 :increment of angle in probing YS and R in-plane variation
        mode='RD' : mode of vpsc loadings (RD, TD, SS, Bulge, inplanebiaxialfc...)
        u11, u22 : Only valid in the case 'inplanebiaxial' or 'inplanebiaxialfc'
        ftex=None, fsx=None : texture and single crystal file all None
        ysproj = '2D' or '3D'

        ihardlaw=0,  # 0: vode -2:rauch

        ## dislocation hardening 'rauch' case
        idislmode=1,  #idislmode for all is active(1) or active-inactive mode(0)
        ilat=0,       #ilat(1): Latent hardening, ilat(0): no interplanar interactions
        ibau=0,       #ibau(1):sudden drop, ibau(0): no such drop
        pp=0.8,       #used to drop the sudden density level.
        iopsysx=0,     #iopsysx(1): a slip mode includes two glide direction (+ and -)
                      #       (0): one slip mode for one glide direction [110] != [-1-10]

        ## associated dislocation hardening parameters
        grsize=3.0e-05,    #grain size
        burgers=2.46e-10,  #burgers vector
        fdisln=180,        #
        shearm=8.5e+04,    #maximum shear modulous
        ftau0=38.,         #lattice friction (initial crss)
        portion=0.8,       #dislocatin exchange portion between denf&denr
        ftherm=2.8         #thermal coefficient in the recovery term

        icase = 1, 2, 3, 4  #case study as described below:

        case #1:
           Prestrain monotonically. (completed)
        case #2:
           Prestrain and grandually reverse-loads (add flow stress curve!)
        case #3:
           Prestrain along RD and change loading path
           calculates theta (orthogonality of loading path change)
        """
        import os

        self.ifig = ifig; self.ncycle = ncycle; self.pfmode = pfmode; self.fsx = fsx
        self.interaction = interaction; self.ihardlaw = ihardlaw;
        self.nprob = nprob #  #of ys-tensile cycle
        self.dang = dang; self.mode = mode; self.eqc = eqc
        self.stp = stp; self.bins = bins
        self.u11 = u11; self.u22 = u22; self.iupdate = iupdate

        ## general single crystal file parameters
        self.hii = hii
        self.hij = hij
        self.hb = hb
        self.hp = hp
        self.iopsysx = iopsysx
                 ## dislocation hardening 'Rauch' case
        self.idislmode = idislmode  #idislmode for all is active(1) or active-inactive mode(0)
        self.ilat = ilat
        self.ibau = ibau       #ibau(1):sudden drop ibau(0): no such drop
        self.pp = pp

        self.sxinteraction = sxinteraction, #Bauschinger or Bauschinger latent or latent ...
        self.dependency=dependency #'indv' or 'iso'
        self.copl=copl
        
                 ## associated dislocation hardening parameters
        self.grsize = grsize    # grain size
        self.burgers = burgers  # burgers vector
        self.fdisln = fdisln    # 
        self.shearm = shearm    # maximum shear modulous
        self.ftau0 = ftau0      # lattice friction (initial crss)
        self.portion = portion  # dislocatin exchange portion between denf&denr
        self.ftherm = ftherm    # thermal coefficient in the recovery term
        
        del ifig, ncycle, pfmode, fsx, interaction, ihardlaw
        del nprob, dang, mode, eqc, stp, bins, u11, u22
        del idislmode, ilat, ibau, pp
        del sxinteraction, dependency, copl, iopsysx,
        del grsize, burgers, fdisln, shearm, ftau0, portion, ftherm
        
        ## Interactiveness switch off --------------------
        ## detects if connected through ssh
        log = os.popen('lastlog -u youngung')
        port = log.read(); port = port.split('\n')[1]
        port = port.split()[1][0:3]
        if port=='pts':
            # me ssh connected. turn off interactiveness
            plt.ioff() 
        else: plt.ion() 
        ## ----------------------------------------------
        
        self.dm = 7.5; self.dn = 7.5 # default resolution of bin
                                     # In the spherical space to which
                                     # grain population is interpolated.
                                     # during pole figure calc.
        
        # texture file or random texture generator ------
        if ftex!=None: self.ftex=ftex
        elif ftex==None:
            # random texture generation
            tex = re(ngrain=ngrain)
            self.ftex = '%i.tex'%ngrain
            tex.write(self.ftex)
            pass
        # -----------------------------------------------
        self.__figprep__()   ## all matplotlib canvases are prepared here!
        mypf, ini = self.__pf__(fsx=self.fsx) # pole figure preparation
        ## --------------------------------------------------

        #single crystal file making
        self.__sxmaker__()
        
        #####  ------------------------------------------------------ #######
        #                     Core of the script                            #
        self.Rmax=0.; self.Smax=0. #Max of R-value and flows stress

        print "\n***************************************"
        print "* Parametric yield surface prediction *"
        print "* Case :#  %i                         *"%icase
        print "***************************************\n"
        
        for i in range(self.ncycle):
            """
            Tensile and YS probing
            (2011-Aug-15, Some other prestrain path to be added)
            """
            ## Initial Yield locus probing and initial in-plane variations. 
            if i==0:
                self.__initys__() #probing of initial yield surface (undeformed)
                self.irecover=0 # no recover from postmort
            else:
                self.irecover=1 # recovery 
                pass
            
            ##-------------------------------------------------------------
            ## prestrain using def tensile

            if icase==1:
                lgr = self.__pretensile__(i=i) # rst contains grains after deform

            elif icase==2:
                ## reverse loading
                if i==0:
                    # figures
                    figr1 = plt.figure(); figr2 = plt.figure();figr1.clf();figr2.clf()
                    axra1 = figr1.add_subplot(111); axra2 = figr2.add_subplot(111)                    
                    e11, s11, r, lgr = self.__pretensile__(i=i,icase=icase) #prestrain 0
                    axra1.plot(e11, s11, label='Initial')
                    axra2.plot(e11,r, label='Initial')                    
                elif i==1:                              #reverse loading
                                                        # with reducted steps
                    # gradually reverse loading
                    # self.eqc = self.eqc/5.
                    self.stp = self.stp / 5
                    self.mode = '-%s'%self.mode
                    print 'Mode is %s'%self.mode
                    """
                    RD: uniaxial tension, -RD:unicompression along RD'
                    """
                    
                else:
                    if i==2:
                        eps=[]; sig=[]; rv=[]
                        eps0 = 0
                        pass
                    e11, s11, r, lgr = self.__pretensile__(i=i, icase=2)
                    e11 = abs(np.array(e11)) + eps0
                    axra1.plot(e11,abs(s11))
                    axra2.plot(e11,r)
                    eps0 = abs(e11[-1]) #
                    if i==self.ncycle - 1: # at the last cycle.
                        figr1.savefig('rauch_case2_flowstress.pdf')
                        figr2.savefig('rauch_case2_rvalue.pdf')
                        pass
                    pass
                        
            elif icase==3:
                ## Prestrain along RD and in-plane rotation then reload
                # check if ncycle==2
                if self.ncycle!=2:
                    print 'currently only self.ncycle==2 case is accepted'
                    raise IOError,'Unexpected ncycle'
                if i==0: #pre tensile for first step
                    # figures
                    figr1 = plt.figure(); figr2 = plt.figure();figr1.clf();figr2.clf()
                    axra1 = figr1.add_subplot(111); axra2 = figr2.add_subplot(111)
                    # Prestrain
                    e11,s11,r = self.__pretensile__(i=i, icase=icase)
                    axra1.plot(e11, s11, label='Initial')
                    axra2.plot(e11,r, label='Initial')
                    pass
                else: #
                    # now self.irecover equals to 1
                    self.angs = np.arange(0.,180., self.dang)
                    self.stp = self.stp * 2
                    for iang in range(len(self.angs)):
                        # explicit rotation of the texture: self.ftex
                        # newgr = __inplanerotation__(ang=self.angs[iang], gr = lgr)
                        # os.remove(self.ftex)
                        # __writegrtofile__(gr=newgr, fn=self.ftex)
                        self.mode='ten_ang'
                        e11, s11, r = self.__pretensile__(
                            i=i,
                            isave=False, #no save.
                            icase=3, ang=self.angs[iang],
                            )
                        axra1.plot(e11, s11, label=r'$\theta$=%4.1f'%self.angs[iang])
                        axra2.plot(e11, r, label=r'$\theta$=%4.1f'%self.angs[iang])
                        pass

                    pass
                pass
            elif icase==4:
                ##
                if i==0:
                    # figures
                    figr1 = plt.figure(); figr2 = plt.figure();figr1.clf();figr2.clf()
                    axra1 = figr1.add_subplot(111); axra2 = figr2.add_subplot(111)                    
                    e11, s11, r, lgr = self.__pretensile__(
                        i=i, icase=icase, isave=True, norm_component='s11',
                        ) #prestrain 0
                    eps0=e11[-1]
                    axra1.plot(e11, s11, label='Initial')
                    axra2.plot(e11,r, label='Initial')
                    self.__ysinvar__()
                    pass
                else:
                    self.irecover=True
                    self.mode='td_without_rotation'
                    e11, s11, r, lgr = self.__pretensile__(
                        i=1, isave=True, icase=icase,
                        norm_component='s22', # --> YS/sig_22
                        )
                    ## Accummulated strain
                    #e11 = e11 + eps0
                    eps0 = e11[-1]
                    ## --------------------
                    
                    axra1.plot(e11, s11, 'o', mec='black', mfc='None')
                    axra2.plot(e11,   r, 'o', mec='black', mfc='None')
                    self.__ysinvar__() #--> exists out below.
                    pass
                pass
            
            ## ------------------------------------------------------------
            ## pole figure plotting after prestrain
            if icase==3: pass
            else:
                mypf = upf.polefigure(
                    grains = lgr, csym=self.csym,
                    cdim=self.cdim, cang=self.cang)
                
                pf_ = mypf.pf(
                    pole=self.poles,  #poles
                    mode=self.pfmode,
                    ifig=self.ifig+1, dm=self.dm, dn=self.dn,
                    axes=[self.figpf.add_subplot(1,self.ncycle+1,2+i)]
                    )
                
                ####################################################
                ## ys - in-plane-variation calculation and plotting
                ## Note that this process does not update
                ## the state variables.
                if icase==2 and i==1 :pass
                elif icase==4: pass
                else:
                    self.__ysinvar__()
                ####################################################
                pass

            pass
        #                                                                   #
        #####  ------------------------------------------------------ #######

        ## detailing the plots --------------------------------------
        self.__detailingplots__()
        if icase in[3,4]:
            axra1.legend(loc='best'); axra2.legend(loc='best')
            axra1.set_ylim(0.,); axra2.set_ylim(0.,)
            figr1.savefig('rauch_00.pdf'); figr2.savefig('rauch_01.pdf')        
            pass
        self.fig.savefig('ys.pdf')       # 2D yield locus
        self.fig3ys.savefig('3ys.pdf')   # 3D yield surface (s11,s22,s12)
        self.fig3ysn.savefig('3ysn.pdf') # Normalized 3D yield surface
        self.figpf.savefig('pf.pdf')     # pole figure evolution
        self.figinp.savefig('invar.pdf') # in-plane variation
        #self.fcm.savefig('fcm.pdf')      # Sig and Eps all components
        self.fid.savefig('fid.pdf')      # CRSS distribution
        self.fig.clf()
        self.fig3ys.clf()
        self.fig3ysn.clf()
        self.figpf.clf()
        self.figinp.clf()
        self.fcm.clf()
        self.fid.clf()


        ## delete out files to ease the syncronizing pdf files via Dropbox
        # import glob
        # fouts = glob.glob('*.out')
        # for f in fouts: os.remove(f)
        ## 

        pass

    def __pretensile__(
        self, i=None, jobid=0, isave=True,
        icase=None, ang=None,
        norm_component = 's11'
        ):
        """
        prestrain as a uniaxial tensile

        1. Uniaxially tensiles
        2. Calculates accummulative strains (11, and VM)
        3. Initiates the figures when i==0
        4. Returns the texture after the last deformation step

        i: i in range(self.ncycle)
        joid: jobid for vp_f2py.vpsc for tensile 'mode'
        isave=True or False: flag whether or not to save the postmort
        icase: the global icase
        ang: angular term if icase is given such that it needs the 'ten_ang' mode
        norm_component: to decide which stress component to norm the yield locus.
        """
        rst = self.tensile(
            irecover=self.irecover, stp=self.stp, eqc=self.eqc,
            isave=isave, interaction=self.interaction,
            iupdate=self.iupdate, mode=self.mode, ihardlaw=self.ihardlaw,
            ftex=self.ftex, fsx=self.fsx,
            bins=self.bins, u11=self.u11, u22=self.u22,
            jobid=jobid, ang=ang)
        ## again normalized component is hardwired as s11
        self.norm = rst[norm_component][-1] - rst['s33'][-1]
        # the last element [-1] after pretensil is saved.
        
        dist = rst['cr'] # crss
        evm = rst['evm'] # von Mises strain
        svm = rst['svm'] # von Mises Stress
        
        ## Axial strain ------
        if i==0: #at the very first step (initial step)
            self.acc = 0.
            self.acc_vm = 0.
            self.acc_svm = 0.
            #E11 = np.array([0.,])
            self.gcf = plt.figure(self.ifig+2, figsize=(15, 4.8))
            self.gcf.clf()
            self.ax3 = self.gcf.add_subplot(131)  #R-value
            self.ax4 = self.gcf.add_subplot(132)  #Sig-eps
            self.fid = plt.figure(self.ifig+4)
            self.fid.clf()
            
            #Sig and eps all components                
            self.fcm = plt.figure(self.ifig+5)
            self.fcm.clf()
            self.ax_sg = self.fcm.add_subplot(131)
            self.ax_sg.set_ylabel(r'$\sigma$',dict(fontsize=18))
            self.ax_sg.set_xlabel(r'$w^{pl}$',dict(fontsize=18))
            
            self.ax_ep = self.fcm.add_subplot(132)
            self.ax_ep.set_ylabel(
                r'$\varepsilon$',dict(fontsize=18))
            self.ax_ep.set_xlabel(
                r'$w^{pl}$',dict(fontsize=18))
            #ax5 = gcf.add_subplot(133)  #Crss-distribution
            self.ax5 = self.fid.add_axes(
                (0.15,0.1,0.55,0.8),
                label='ax5') #Crss-distribution
            self.ax3.set_xlabel(r'$\varepsilon_{11}$',
                                dict(fontsize=18));
            self.ax4.set_xlabel(r'$\varepsilon_{11}$',
                                dict(fontsize=18))
            self.ax3.set_ylabel(r'R',dict(fontsize=18));
            self.ax4.set_ylabel(r'$\sigma_{11}$',
                                dict(fontsize=18))
            self.ax5.set_ylabel('Fraction', dict(fontsize=20))
            #ax5.set_xticks(
            #  np.arange(0.,bins,bins-1), ('low', 'high'))  
            plt.xticks(np.arange(0.,self.bins,self.bins-1),
                       ('Minimum','Maximum'))
            
            ## Accumulated E11 and von Mises strain
            # E11 = np.array(np.append([0,],rst['e11']))
            # evm = np.array(np.append([0,],rst['evm']))
            E11 = np.array(rst['e11'])
            evm = np.array(rst['evm'])
            svm = np.array(rst['svm'])
            pass
        else:
            # E11 = np.append([0],rst['e11'])+acc
            # evm = np.append([0],rst['evm'])+acc_vm
            E11 = np.array(rst['e11']) + self.acc
            evm = np.array(rst['evm']) + self.acc_vm
            svm = np.array(rst['svm']) #+ self.acc_svm
            pass
        
        self.acc = E11[-1]
        self.acc_vm = evm[-1]

        ## Accumulated E11 and eps_VM -------------------
        R = rst['R']
        ## hardwired flow stress along RD(11)
        flows = rst['s11'] - rst['s33']
            
        if max(R)>self.Rmax: self.Rmax = max(R)
        if max(flows)>self.Smax: self.Smax = max(flows)
        self.Rmax = max(R)
        self.Smax = max(flows)
        
        print 'len(E11)', len(E11)
        print 'len(R)', len(R)
            
        self.ax3.plot(E11, R)
        self.ax4.plot(E11, flows)
            
        self.ax_ep.plot(rst['e11'], label='e11')
        self.ax_ep.plot(rst['e22'], label='e22')
        self.ax_ep.plot(rst['e33'], label='e33')
        self.ax_ep.plot(rst['e12'], label='e12')
        self.ax_ep.plot(rst['e23'], label='e23')
        self.ax_ep.plot(rst['e13'], label='e13')
                       
        self.ax_sg.plot(rst['s11'], label='s11')
        self.ax_sg.plot(rst['s22'], label='s22')
        self.ax_sg.plot(rst['s33'], label='s33')
        self.ax_sg.plot(rst['s12'], label='s12')
        self.ax_sg.plot(rst['s23'], label='s23')
        self.ax_sg.plot(rst['s13'], label='s13')
        # crss distribution
        dist = np.array(dist)
        dist = dist/float(dist.sum())
        self.ax5.plot(dist); self.ax5.set_ylim(0.,1.0)
        #E11 = np.array([0.,])
        self.ax3.set_ylim(0.,self.Rmax*1.15)
        self.ax4.set_ylim(0.,self.Smax*1.15)
        ##------------------------------------------------------------
        
        ## irecover flag
        if i==0: self.irecover = 0 
        else: self.irecover = 1
        # rst.keys() #def tensile returns a dictionary
        # rst['e11']; rst['e22']; rst['e33']; rst['s11']
        # rst['d11']; rst['d22']; rst['d33']; rst['ag'] # aggregate
        ## current pole figure----------------------------------------

        #if fsx==None: poles=[[1,1,1]] #default
        #elif fsx!=None: poles=[[0,0,0,2]] #default for fsx!=None
        if icase==3:
            return rst['e11'], flows, R
            pass
        elif icase in [2,]:
            return rst['e11'], flows, R, rst['ag']
        elif icase in [4,]:
            return evm, flows, svm, rst['ag']
        else:
            return rst['ag']

    def __ysinvar__(self):
        """
        1) Caculates yield locus on the plane stress space,
        2) in-plane variation and,
        3) plots them.

        Add Sxy-Sxx  yield surface plotting (2011-Oct-21)
        """
        ##------------------------------------------------------------
        ## Yield stress from the saved state. 'irecover=1'
        ysx, ysy, ysz, d, ds = self.ys(
            irecover=1, stp=self.nprob,
            ftex=self.ftex, fsx=self.fsx,
            interaction=self.interaction, proj='3D',
            jobid=0, ihardlaw=self.ihardlaw)
        
        ysx = np.append(ysx, [ysx[0]]); ysy = np.append(ysy, [ysy[0]])
        ysz = np.append(ysz, [ysz[0]])

        ##------------------------------------------------------------
        ## in-plane variations of yield surface and R-value
        YSp, Rp, ang = self.var(
            irecover=1, stp=self.dang, ftex=self.ftex, ihardlaw=self.ihardlaw,
            fsx=self.fsx, interaction=self.interaction, jobid=0)
        
        ##-------------------------------------------------
        ## normalization of stress 
        ysxn = ysx/self.norm
        ysyn = ysy/self.norm
        yszn = ysz/self.norm
        yspp = YSp.copy()
        Rp = Rp/Rp[0]; YSp = YSp/YSp[0]
        ##-------------------------------------------------

        ## Yield locus + strain vector direction plotting 
        # ax1.plot(ysx,ysy,
        #       label=r'$\varepsilon_{11}^{pl}=$%s'%str(acc))
        self.ax1.plot(
            ysx,ysy,
            label=r'$\varepsilon_{VM}^{pl}=$%s'%str(self.acc_vm))
        #ax2.plot(ysxn,ysyn,
        #       label=r'$\varepsilon_{11}^{pl}=$%s'%str(acc))
        self.ax2.plot(
            ysxn, ysyn,
            label=r'$\varepsilon_{VM}^{pl}=$%s'%str(self.acc_vm))
        self.ax3ys.plot(
            ysx, ysy, ysz,
            label=r'$\varepsilon_{VM}^{pl}=$%s'%str(self.acc_vm))
        self.ax3ysn.plot(
            ysxn, ysyn, yszn,
            label=r'$\varepsilon_{VM}^{pl}=$%s'%str(self.acc_vm))

            
        ## strain rate around yield surface.
        for j in range(len(d)):
            color='gray'; alpha=0.5
            r = ysx[0]*0.1
            self.ax1.plot(
                [ysx[j], ysx[j] + r * math.cos(d[j])],
                [ysy[j], ysy[j] + r * math.sin(d[j])],
                color, alpha=alpha)
            self.ax2.plot(
                [ysxn[j], ysxn[j] + 0.1 * math.cos(d[j])],
                [ysyn[j], ysyn[j] + 0.1 * math.sin(d[j])],
                color, alpha=alpha)
            self.ax3ys.plot(
                [ysx[j], ysx[j] + r * math.cos(
                        ds[j]) * math.cos(d[j])],
                [ysy[j], ysy[j] + r * math.cos(
                        ds[j]) * math.sin(d[j])],
                [ysz[j], ysz[j] + r * math.sin(ds[j])],
                color, alpha=alpha)
            self.ax3ysn.plot(
                [ysxn[j],
                 ysxn[j] + 0.1 * math.cos(
                        ds[j]) * math.cos(d[j])],
                [ysyn[j],
                 ysyn[j] + 0.1 * math.cos(
                        ds[j]) * math.sin(d[j])],
                [yszn[j],
                 yszn[j] + 0.1 * math.sin(ds[j])],
                color, alpha=alpha)
            pass
        
        ## In-plane variations of R and YS
        self.ax6.plot(ang, Rp); self.ax7.plot(ang, YSp);self.ax8.plot(ang,yspp)
        pass
    
    def __detailingplots__(self):
        """
        Before saving the canvases,
        (a) Write down some information.
        (b) Make the plots more readable by desgnating axis labels.
        (c) and axis limits, and ticks.
        """
        import matplotlib.font_manager as fontm
        ### make planar yield locus in the square axes patch.
        self.ax1.plot([0],[0],'k+')#indicates a point origin on the stress space
        xl = self.ax1.get_xticks()[0]; xh = self.ax1.get_xticks()[-1]
        yl = self.ax1.get_yticks()[0]; yh = self.ax1.get_yticks()[-1]
        self.ax1.set_xlim(min(xl,yl), max(xh,yh))
        self.ax1.set_ylim(min(xl,yl), max(xh,yh))
        xl = self.ax2.get_xticks()[0]; xh = self.ax2.get_xticks()[-1]
        yl = self.ax2.get_yticks()[0]; yh = self.ax2.get_yticks()[-1]
        self.ax2.set_xlim(min(xl,yl), max(xh,yh))
        self.ax2.set_ylim(min(xl,yl), max(xh,yh))
        
        prop = fontm.FontProperties(size=8)
        self.ax1.legend(bbox_to_anchor=(1.05,1),loc=2, prop=prop)
        # self.ax2.set_ylim(-1.28, 1.28); self.ax2.set_xlim(-1.28, 1.28)
        # self.ax2.set_xticks(np.arange(-1.0, 1.01, 1))
        # self.ax2.set_yticks(np.arange(-1.0, 1.01, 1))
        self.ax2.grid()
        self.ax6.set_xlabel(r'$\theta$ from RD',dict(fontsize=20))# from RD')
        self.ax6.set_ylabel(r'$r / r_{0}$',dict(fontsize=20))
        self.ax6.set_xticks(np.arange(0.,180.001,45.0));self.ax6.grid()
        self.ax7.set_xlabel(r'$\theta$ from RD',dict(fontsize=20)) # from RD')
        self.ax7.set_ylabel(r'$\sigma / \sigma_{RD}$',dict(fontsize=20))
        self.ax7.set_xticks(np.arange(0.,180.001,45.0));self.ax7.grid()
        self.ax8.set_xticks(np.arange(0.,180.001,45.0));self.ax8.grid()

        self.ax3ys.set_xlabel(r'$\sigma_{x}$', dict(fontsize=20))
        self.ax3ys.set_ylabel(r'$\sigma_{y}$', dict(fontsize=20))
        self.ax3ys.set_zlabel(r'$\sigma_{z}$', dict(fontsize=20))
        self.ax3ysn.set_xlabel(r'$\sigma_{x}$', dict(fontsize=20))
        self.ax3ysn.set_ylabel(r'$\sigma_{y}$', dict(fontsize=20))
        self.ax3ysn.set_zlabel(r'$\sigma_{z}$', dict(fontsize=20))

        self.ax_ep.legend(loc='best')
        self.ax_sg.legend(loc='best')

        ## writes information!! ------------------------------------
        if   self.interaction==0: self.interaction='FC - Taylor'
        elif self.interaction==1: self.interaction='affine'
        elif self.interaction==2: self.interaction='secant'
        elif self.interaction==3: self.interaction=r'$n^{eff}=10$'
        elif self.interaction==4: self.interaction='tangent'
        elif self.interaction==5: self.interaction='$2^{nd}$ order'
        else: raise IOError

        self.fig.text(0.05, 0.15,r'interaction: %s '%(self.interaction))
        
        # self.fig.text(0.05, 0.05,r'$\tau_0$: %3.1f'%
        #               (self.tau0) + '  ' + r'$\tau_1$: %3.1f'%
        #               (self.tau1) + '  ' + r'$\theta_0$: %3.1f'%
        #               (self.thet0)+ '  ' + r'$\theta_1$: %3.1f'%
        #               (self.thet1))

        if self.ihardlaw==-2: #Rauch et al's Hardening
            self.fig.text(0.05,0.10,"Rauch el al.'s dislocation hardening law")
            self.fig.text(0.05,0.07,"ibau: %1i"%self.ibau)
            self.fig.text(0.05,0.04,"ilat: %1i"%self.ilat)
            pass

        ##
        x = 0.6; y= 0.22; dy = 0.03
        if self.iopsysx==0:temp='Unidirectional'
        elif self.iopsysx==1:temp='Bidirectional'
        self.fig.text(x,y,'iopsysx: %s'%temp, dict(fontsize=8))
        y = y - dy
        ##
        
        if self.dependency=='indv':
            ##
            self.fig.text(x,y,'Individual slip modes',
                          dict(fontsize=8))
            y = y - dy
            ##
            self.fig.text(x, y, r'$H_{S}$: %3.1f'%(self.hii),
                          dict(fontsize=8))
            y = y - dy
            ##
            self.fig.text(x, y, r'$H_{L}$: %3.1f'%(self.hij),
                          dict(fontsize=8))
            y = y - dy
            ##
            if self.iopsysx==0:
                """ only if the slip mode is unidirectional"""
                self.fig.text(x,y, r'$H_{B}$: %3.1f'%(self.hb),
                              dict(fontsize=8))
                y = y - dy
                pass
            if self.copl==True:
                self.fig.text(x,y, r'$H_{P}$: %3.1f'%(self.hp),
                              dict(fontsize=8))
                y = y - dy
                pass
            ##
            pass
        
        elif self.dependency=='iso':
            self.fig.text(
                x, y, 'Isotropic latent hardening',
                dict(fontsize=8)
                )
            y = y - dy
            ###---------------------------------------------------------
            pass
        pass
    
    def __sxmaker__(self):
        """
        define self.fsx
        """
        if self.fsx!=None:
            if not(os.path.isfile(self.fsx)): raise IOError, 'Wrong file name for fsx'
            return # exit
        
        # single crystal file generator
        # isotropic, Bauschinger, Latent hardening,
        # Bauschinger + Latent hardening        
        self.fsx ='temp.sx'
        tau0 = 1.07e2; tau1 =0.73e2
        thet0 = 2.80e02; thet1 = 1.45e02
            #hii = 1.0; hij = 2.0; hb = 1.0 #they are now arguments
            #dependency='indv'
        dependency = self.dependency

            #ssint = 'Bauschinger Latent';'Bauschinger'; 'Latent'
        ssint = self.sxinteraction #ios, indiv
        nrsx = 20  #iopsys is passed to def __init__ as argument
        if self.ihardlaw==0: shard ='voce'
        elif self.ihardlaw==-2: shard ='rauch'
        else: raise IOError, 'Unexpected ihardlaw'
        sx_maker.cubic(
            filename=self.fsx,
            
            #self-, latent-, bauschinger, co-planar coefficient
            hii=self.hii, hij=self.hij, hb=self.hb, hp=self.hp, copl=self.copl,
            
            ## voce hardening curve characterization
            tau0=tau0, tau1=tau1, thet0=thet0, thet1=thet1,
            header='**ys_parameteric study: ',
            
            ## dislocation hardening 
            ihardlaw=shard, #
            
            ## dicloation hardeing single crystal file ('rauch' case)
            idislmode=self.idislmode,  # idislmode for all is
            # act(1) or act-inactive mode(0)
            ilat=self.ilat,
            # ilat(1): latent hardening (1): no latent hardening(0)
            ibau=self.ibau,       # ibau(1):sudden drop,
            # ibau(0): no such drop
            pp=self.pp,           # used to drop the sudden density level.
            
            grsize=self.grsize,
            burgers=self.burgers,
            fdisln=self.fdisln,
            shearm=self.shearm,
            ftau0=self.ftau0,
            portion=self.portion,
            ftherm=self.ftherm,
            
            ## disloction hardening description.
            iopsysx=self.iopsysx, dependency=dependency,
            nrsx=nrsx, interaction=ssint,
            
            ## slip systems
            #Face-centered cubic
            b=[[1,1,0]], n=[[1,1,1]],    
            #Body-centered (pencile glide) cubic
            #b=[[1,1,1],[1,1,1],[1,1,1]],
            #n=[[1,1,0],[1,1,2],[1,2,3]], 
            )
            ## -----------------------------------------------
        pass

    def __initys__(self):
        """
        initial yield surface
        """
        print "\n*********************"
        print "Yield surface probing"
        print "*********************\n"
        #   ysx, ysy, d = self.ys(
        #   irecover=0, stp=nprob,
        #   axes=[[1.,1.,1.]], angle=[[0.,0.,0.]],
        #   ftex=self.ftex, fsx=self.fsx,
        #   interaction=interaction, proj='2D')        
        ysx, ysy, ysz, d, ds = self.ys(
            irecover=0, stp=self.nprob,
            ftex=self.ftex, fsx=self.fsx,
            interaction=self.interaction, proj='3D',
            ihardlaw=self.ihardlaw)
        
        ysx = np.append(ysx, [ysx[0]])
        ysy = np.append(ysy, [ysy[0]])
        ysz = np.append(ysz, [ysz[0]]) 
        
        ##------------------------------------------------
        ## This is just for getting to know the initial YS
        ## along hardwired RD direction
        norm = self.initialys(
            ftex=self.ftex , #ftex=ftex,
            fsx=self.fsx,# fsx will be further specified
            interaction=self.interaction, 
            mode='RD', #mode # VPSC class loading
            ihardlaw=self.ihardlaw #0:voce, -2:rauch
            )
        ##------------------------------------------------

        ## Yield locus normalization
        ysxn = ysx/norm; ysyn = ysy/norm; yszn = ysz/norm
                
        ## Yield locus + strain rate direction
        self.ax1.plot(ysx,ysy,label='Initial')
        self.ax2.plot(ysxn,ysyn, label='Initial')
        self.ax3ys.plot(ysx, ysy, ysz, label='Initial')
        self.ax3ysn.plot(ysxn, ysyn, yszn, label='Initial')
        for j in range(len(d)):
            r = ysx[0]/10.; color='black'; alpha=1.
            self.ax1.plot(
                [ysx[j], ysx[j] + r * math.cos(d[j])],
                [ysy[j], ysy[j] + r * math.sin(d[j])],
                color, alpha=alpha,)
            self.ax2.plot(
                [ysxn[j], ysxn[j] + 0.1 * math.cos(d[j])],
                [ysyn[j], ysyn[j] + 0.1 * math.sin(d[j])],
                color, alpha=alpha)
            self.ax3ys.plot(
                [ysx[j],
                 ysx[j] + r * math.cos(
                        ds[j]) * math.cos(d[j])],
                [ysy[j],
                 ysy[j] + r * math.cos(
                        ds[j]) * math.sin(d[j])],
                [ysz[j],
                 ysz[j] + r * math.sin(ds[j])],
                color, alpha=alpha)
            self.ax3ysn.plot(
                [ysxn[j],
                 ysxn[j] + 0.1 * math.cos(
                        ds[j]) * math.cos(d[j])],
                [ysyn[j],
                 ysyn[j] + 0.1 * math.cos(
                        ds[j]) * math.sin(d[j])],
                [yszn[j],
                 yszn[j] + 0.1 * math.sin(ds[j])],
                color, alpha=alpha) 
            pass
        ## end of ys and strain rate plotting. -----
                     
        #ax6: R    self.ax7: YS
        YSp, Rp, ang = self.var(
            irecover=0, stp=self.dang, ftex=self.ftex,
            ihardlaw=self.ihardlaw,
            fsx=self.fsx, interaction=self.interaction)
        self.ax8.plot(ang, YSp) #non normalized yield strength
        YSp = YSp/YSp[0]; Rp = Rp/Rp[0]
        self.ax6.plot(ang, Rp)
        self.ax7.plot(ang, YSp)
        
        ##----------------------------------------
        pass
        
        
    def __pf__(self, fsx=None):
        """
        returns mypf and initial pole figures
        """
        ## plotting the first pole figure--------------------
        ## arguments to upf.py
        if fsx==None:
            ## if fsx is not given, assumes sx as cubic
            self.csym = 'cubic'; self.cdim = [1., 1., 1.];
            self.cang = [90., 90., 90.] #Default single crystal
            self.poles = [[1, 1, 1]]
            pass
        else:
            ## SX is acceppted if either hexagonal or cubic
            FILE = open(fsx, 'r')
            info = FILE.readlines()[0:3]; FILE.close()
            if info[1].split()[0]=='HEXAGONAL':
                self.csym = 'hexag'
                self.poles = [[0, 0, 0, 2]]
                pass
            if info[1].split()[0]=='cubic':
                self.csym = 'cubic'
                self.poles = [[1, 1, 1]]
                pass
            else:
                print 'Unexpected crystal',
                print 'structure found in %s'%fsx
                raise IOError
            info[1] #crystal symmetry
            info[2] #cdim and cang
            self.cdim = map(float, info[2].split()[0:3])
            self.cang = map(float, info[2].split()[3:6])
            pass
            
        mypf = upf.polefigure(
            filename=self.ftex, csym=self.csym,
            cdim=self.cdim, cang=self.cang)
        ini = mypf.pf(
            pole=self.poles, mode=self.pfmode,# 'dot', 'contourf' 
            ifig=self.ifig+1, dm=self.dm, dn=self.dn,
            axes=[self.figpf.add_subplot(1,self.ncycle+1,1)])
        return mypf, ini
    
    def __figprep__(self): #pole figure preparation
        """
        Global figure preparation. All canvases!
        """
        from mpl_toolkits.mplot3d import axes3d, Axes3D

        ## prepare the canvas on which results are displayed
        # figures
        self.fig    = plt.figure(
            self.ifig, figsize=(14,5))  # Yield locus (2D)
        self.fig.clf()
        self.fig3ys = plt.figure(self.ifig+6) # Yield surface (3D)
        self.fig3ys.clf()
        self.fig3ysn = plt.figure(self.ifig+7)  # Yield surface noarmalized (3D)
        self.fig3ysn.clf()
        self.figpf  = plt.figure(
            self.ifig+1, figsize=(6*self.ncycle,3)) # pole figure
        self.figpf.clf()
        self.figinp = plt.figure(
            self.ifig+3, figsize=(20,5)) # in-plane variation(R, YS)
        self.figinp.clf()
        
        # axes 
        rectl = 0.07, 0.15, 0.38, 0.7  #l,b,w,h
        rectr = 0.58,  0.15, 0.38, 0.7
        # self.ax1   = self.fig.add_axes(rectl, label='ax1') #YS
        # self.ax2   = self.fig.add_axes(rectr, label='ax2') #Normalized YS-locus
        self.ax1 = self.fig.add_subplot(131, label='YS')
        self.ax2 = self.fig.add_subplot(133, label='normalized YS')
        self.ax3ys = Axes3D(self.fig3ys)                   # 3D yield surface
        self.ax3ysn= Axes3D(self.fig3ysn)

        # rectl = 0.07, 0.15, 0.38, 0.7  #l,b,w,h
        # rectr = 0.58, 0.15, 0.38, 0.7        
        # self.ax6 = self.figinp.add_axes(
        #     rectl, label='ax6') #Normalized in-plane R
        # self.ax7 = self.figinp.add_axes(
        #     rectr, label='ax7') #Normalized in-plane YS
        self.ax6 = self.figinp.add_subplot(131, label='ax6')
        self.ax7 = self.figinp.add_subplot(132, label='ax7')
        self.ax8 = self.figinp.add_subplot(133, label='ax8')
        self.figinp.subplots_adjust(
            bottom=0.3, top=0.7, wspace=0.5, hspace=0.3
            )
        ## -------------------------------------------------


        ## Uniaxial Tensile and YS cycle
        self.master = []
        self.ax1.set_aspect('equal'); self.ax2.set_aspect('equal')
        self.ax1.set_title('Yield surface evolution')

        
        self.ax2.set_title(
            'Normalized along the loading direction',
            dict(fontsize=8)
            )
        self.ax1.set_xlabel(
            r'$\sigma_{RD} $' + ' [MPa], ' + r'$\dot{\varepsilon}^{pl}_{RD}$',
            dict(fontsize=8))
        self.ax1.set_ylabel(
            r'$\sigma_{TD} $' + ' [MPa], ' + r'$\dot{\varepsilon}^{pl}_{TD}$',
            dict(fontsize=8))
        self.ax2.set_xlabel(
            r'$\sigma_{RD}$' + ', ' + r'$\dot{\varepsilon}^{pl}_{RD}$',
            dict(fontsize=8))
        self.ax2.set_ylabel(
            r'$\sigma_{TD}$' + ', ' + r'$\dot{\varepsilon}^{pl}_{TD}$',
            dict(fontsize=8))
        
        pass
        

    def initialys(self, ftex=None, fsx=None, interaction=None,
                  mode=None, proj='2D', jobid=0, ihardlaw=0):
        """
        Finds the initial uniaxial yield stress
        ---------
        Arguments
        ftex = None
        fsx  = None
        interaction = None
        mode = None: VPSC class' loading mode ('RD','BU','SS','inplanebiaxialfc')
        """
        rst = self.tensile(irecover=0, stp=1, isave=False,
                           axes=[[1.,1.,1.]], angle=[[0.,0.,0.]],
                           ftex=ftex, fsx=fsx,interaction=interaction,
                           mode=mode, jobid=jobid, ihardlaw=ihardlaw)
        return rst['s11'][0]-rst['s33'][0]

    def tensile(self, irecover=0, stp=None, eqc=None,
                axes= [[1.,1.,1.]], angle=[[0.,0.,0.]],
                iupdate=[1,1,1,0],
                ftex=None, fsx=None, interaction=3,
                ihardlaw=0,
                isave=False,bins=20, mode=None,
                u11=None,u22=None, jobid=0, ang=None):
        """
        Performs a uniaxial (or biaxial) tensile test. (axiality)
        Depending on the irecover and isave
        interacts with the postmorten files

        The simulated results are returned properly

        Must chose what to be returned among all the keys
        in the datamaster dictionary

        ## Here I chose,
        1. Uniaxial(or biaxial) tension stress tensor
        2. Uniaxial(or biaxial) tension strain tensor
           * All of the tensors will be that of 3x3 2nd order,
             and they are macroscopic entities.
        3. Polycrystalline aggregate
        4. CRSS distribution # May need be normalized...?
        5. Returns R-value!
        """
        import shutil
        if stp==None: raise IOError, 'stp is not given'
        if irecover==1:
            pass
        #     if os.path.isfile("POSTMORT.OUT"):
        #         os.rename("POSTMORT.OUT", "POSTMORT.IN")
        #     else: raise IOError
            # else:
            #     if os.path.isfile("POSTMORT.IN"): pass
            #     else: raise IOError
        else: pass # No recovery by reading the postmort.in file

        if isave==False: isave=0
        elif isave==True: isave=stp #Always the last step is to be saved.
        else: raise IOError, 'wrong isave argument'

        if stp==0: raise IOError

        #os.system('ls') ;raw_input('enter to proceed>>')

        ## decision on the mode
        # it can be "RD", "TD", "inplanebiaxialfc"
        job = vpsc(irecover=irecover, isave=stp, stp=stp, mode=mode, eqc=eqc,
                   texture=ftex, fsx=fsx, init_eld_rat=axes, init_eul_axe=angle,
                   interaction=interaction, ihardlaw=ihardlaw,
                   iupdate=iupdate, u11=u11, u22=u22, jobid=jobid, ang=ang)

        job.run(); job.pp()
        if isave!=0: #save the postmort output
            if os.name!='posix':
                print 'Unexpected os.system found'
                raw_input('Enter to proceed>>')
                raise IOError
            
            os.system("cp POSTMORT_%s.OUT POSTMORT_%s.IN"%(str(jobid).zfill(4),
                                                           str(jobid).zfill(4))
                      )
            os.system("rm -f POSTMORT_%s.OUT"%(str(jobid).zfill(4)))
            pass
        os.system('clear')
        #os.system('ls'); raw_input('enter to proceed >>>')
        if self.mode in ['ten_ang','TEN_ANG','TD','DD']:ipro = 1
        
        else: ipro =0
        sbar = job.datamaster['sbar'][ipro]  ## Macroscopic stress tensor
        dbar = job.datamaster['dbar'][ipro]  ## Macroscopic strain tensor
        eps  = job.datamaster['eps'][ipro]   ## Eps
        svm  = job.datamaster['svm'][0]  ## von Mises stress        
        evm  = job.datamaster['evm'][0]  ## von Mises strain
        agre = job.datamaster['lgr']      # the most recent aggregate's COD

        s11 = []; s22 = []; s33 = []
        s12 = []; s23 = []; s13 = []
        
        d11 = []; d22 = []; d33 = []
        d12 = []; d23 = []; d13 = []        
        
        e11 = []; e22 = []; e33 = []
        e12 = []; e23 = []; e13 = []

        for i in range(len(sbar)):
            s11.append(sbar[i][0,0])
            s22.append(sbar[i][1,1])
            s33.append(sbar[i][2,2])
            s12.append(sbar[i][0,1])
            s23.append(sbar[i][1,2])
            s13.append(sbar[i][0,2])
            
            d11.append(dbar[i][0,0])
            d22.append(dbar[i][1,1])
            d33.append(dbar[i][2,2])
            d12.append(dbar[i][0,1])
            d23.append(dbar[i][1,2])
            d13.append(dbar[i][0,2])

        for i in range(len(eps)):  #eps usually one iteration less than sbar or dbar
            e11.append(eps[i][0,0])
            e22.append(eps[i][1,1])
            e33.append(eps[i][2,2])
            e12.append(eps[i][0,1])
            e23.append(eps[i][1,2])
            e13.append(eps[i][0,2])

        master = {} # master data dictionary
        s11, s22, s33, d11, d22, d33, e11, e22, e33 = __makenp__(s11,s22,s33,
                                                                 d11,d22,d33,
                                                                 e11,e22,e33)

        e12, e23, e13, d12, d23, d13, e12, e23, e13 = __makenp__(e12,e23,e13,
                                                                 d12,d23,e13,
                                                                 e12,e23,e13)
        evm = np.array(evm)

        master['s11'] = s11; master['s22'] = s22; master['s33'] = s33 #stress
        master['s12'] = s12; master['s23'] = s23; master['s13'] = s13
        
        master['d11'] = d11; master['d22'] = d22; master['d33'] = d33 #strain rate
        master['d12'] = d12; master['d23'] = d23; master['d13'] = d13
        
        master['e11'] = e11; master['e22'] = e22; master['e33'] = e33 #strain history
        master['e12'] = e12; master['e23'] = e23; master['e13'] = e13
        
        master['R']   = d22/d33  #R-value (instantaneous)
        master['ag']  = agre
        master['evm'] = evm; master['svm'] = svm

        #print glob.glob('crss*.out');raw_input()
        
        master['cr'] = cr(bins=bins)  #cr method
        return master

    def var(self, irecover=0, stp=32, ihardlaw=0,
            ftex=None, fsx=None, interaction=None, jobid=0):
        """
        The in-plane variation of yield stress and R-value
        stp: incremental angle
        """
        if irecover==1:
            if os.path.isfile("POSTMORT_%s.IN"%str(jobid).zfill(4)): pass
            else: raise IOError
        job = vpsc(irecover=irecover, stp=stp, mode='lankf',
                   texture=ftex, fsx=fsx,
                   interaction=interaction,
                   ihardlaw=ihardlaw,
                   jobid=jobid                   
                   )

        job.run()
        job.pp()
        scau = job.datamaster['scau'][0]
        dbar = job.datamaster['dbar'][0]

        ang = np.arange(0.,90.00001, stp)
        ang = np.append(ang, ang+90.)

        Y = np.array(job.datamaster['YSprob'][0])
        Y = np.append(Y, job.datamaster['YSprob'][2])
        R = np.array(job.datamaster['Rprob'][0])
        R = np.append(R,job.datamaster['Rprob'][2])

        return Y, R, ang

    def ys(self, irecover=0, stp=32,
           axes=[[1.,1.,1.,]], angle=[[0.,0.,0.]],
           ftex=None, fsx=None, interaction=3,
           proj='2D', jobid=0,
           ihardlaw=0):
        """
        Yield surface probing

        1. Obtains the complete sbar, dbar.
        2. Post-processes sbar and dbar to get the yield surface plotting and
           strain direction at each point.
        3. Return stress contour together with its probing strain's free vector whose
           origin is coincident with the stress state.
        ** The values are in numpy array!!

        proj='2D' or '3D'
        """
        if irecover==1:
            if os.path.isfile("POSTMORT_%s.IN"%str(jobid).zfill(4)): pass
            else: raise IOError
            pass
        
        job = vpsc(irecover=irecover, stp=stp, mode='YS',
                   texture=ftex, fsx=fsx,
                   init_eld_rat=axes, init_eul_axe=angle,
                   interaction=interaction, #interaction (0:FC, .. 3:neff)
                   ihardlaw=ihardlaw, #0:voce, -2:rauch
                   jobid=jobid
                   )
        job.run()
        job.pp()

        sbar = job.datamaster['sbar'][0] ## Macroscopic stress tensor
        dbar = job.datamaster['dbar'][0] ## Macroscopic strain tensor

        # for i in range(len(sbar)):
        #     s11=sbar[i][0,0]; d11=dbar[i][0,0]
        #     s22=sbar[i][1,1]; d22=dbar[i][1,1]
        #     s33=sbar[i][2,2]; d33=dbar[i][2,2]
        #     s12=sbar[i][0,1]; d12=dbar[i][0,1]
        #     s23=sbar[i][1,2]; d23=dbar[i][1,2]
        #     s13=sbar[i][0,2]; d13=dbar[i][0,2]

        ## 1-2 components as x and y

        ## strains ----------------
        d   = np.zeros((len(dbar))) # returns the angle using atan2(y,x)
        phi = np.zeros((len(dbar))) #
        ssd = np.zeros((len(dbar))) # shear strain rate norm
        # shear strains
        e12 = np.zeros((len(sbar)))
        e23 = np.zeros((len(sbar)))
        e13 = np.zeros((len(sbar)))


        ## stresses --------------
        sx = np.zeros((len(sbar)))
        sy = np.zeros((len(sbar)))
        # shear stresses
        s12 = np.zeros((len(sbar)))
        s23 = np.zeros((len(sbar)))
        s13 = np.zeros((len(sbar)))
        sss = np.zeros((len(sbar)))  # shear stress norm
        
        for i in range(len(sbar)):
            #### stresses ---------------
            sx[i] = sbar[i][0,0]   #s11
            sy[i] = sbar[i][1,1]   #s22
            # shear stresses
            s12[i] = sbar[i][0,1]  #s12
            s23[i] = sbar[i][1,2]  #s23
            s13[i] = sbar[i][0,2]  #s13
            # shear stress norm
            sss[i] = np.sqrt(s12[i]**2 + s23[i]**2 + s13[i]**2)
            #### ------------------------

            #### strains ----------------
            ssd[i] = np.sqrt(dbar[i][0,1]**2 + dbar[i][1,2]**2 + dbar[i][0,2]**2)
            d[i]   = math.atan2(dbar[i][1,1], dbar[i][0,0]) #radian
            dn     = np.sqrt(dbar[i][1,1]**2 + dbar[i][0,0]**2) #
            #phi[i] = math.atan2(dn, ssd[i])
            phi[i] = math.atan2(ssd[i],dn)
            #### ------------------------
            pass

        if proj=='2d' or proj=='2D': return sx,sy,d
        #S11, S22, ShearStress, theta, phi
        elif proj=='3d' or proj=='3D': return sx, sy, sss, d, phi 

        
        else:
            raise IOError

    
## disposed method(s)
    # def __axe_ang_gr_from_ftex__(self, filename='TEX_PH1.OUT'):
    #     """
    #     Provided the file whose name is given is the VPSC
    #     convention of writing format, returns the
    #     axes, angles, and polycrystalline aggregate.
    #     """
    #     FILE = open(filename, 'r')
    #     lines = FILE.readlines()[1:3]
    #     axes   = map(float, lines[0].split()[0:3])
    #     angles = map(float, lines[1].split()[0:3])
    #     FILE.close()
    #     gr = np.genfromtxt(filename, skiprows=4) #
    #     return axes, angles, gr

## This class is to show the influence of choice on linearization scheme on
## VPSC prediction on initial YS surface with the rest of all parameters
## being the same.
## For each linearization scheme, tension along RD precedes, in order to
## perform normalization.


class Linearization():
    """
    comparison of the consequences upon choice of linearization schemes
    interpreted on the yield surface

    arguments:
    ifig
    stp=32
    texture=None
    fsx=None
    """


    def __init__(self, ifig=1, stp=32, texture=None, fsx=None):
        l = 0.15  #left blank
        b = 0.15  #bottom blank
        w = 0.75  #width
        h = 0.80  #height
        nax = 1.

        fig_width=7
        fig_height=5
        fig_rel_size= 0.9 

        fig_width, fig_height = np.array(
            [fig_width,fig_height]) * fig_rel_size

        fig = plt.figure(ifig  , [fig_width*nax, fig_height])
        fig1= plt.figure(ifig+1, [fig_width*nax, fig_height])

        el = l/nax
        ew = w/nax
        fullscale = 1.0
        incr = fullscale / nax

        ax1 = fig.add_axes((el,      b,ew,h), aspect='equal')
        ax2 = fig1.add_axes((el     ,b,ew,h), aspect='equal')

        markers = ['o','^','d','p','h','*','+',
                   'o','^','d','p','h','*','+']
        linestyle = ['-','--','-.',':','-','--',
                     '-.',':','-','--','-.',':',]    

        # ngrain=100
        # tex = re(ngrain=ngrain)
        # self.ftex = '%i.tex'%ngrain
        # self.fsx = None
        # tex.write(self.ftex)
        self.ftex=texture
        self.fsx=fsx

        #fig = plt.figure(ifig)
        # ax1 = fig.add_subplot(121,aspect='equal')
        # ax2 = fig.add_subplot(122,aspect='equal')

        interactions = [0,1,2,3,4]
        for i in interactions:
            uni = vpsc(stp=1, mode='rd', texture=self.ftex,
                       fsx=self.fsx, interaction=i)
            uni.run();uni.pp()
            s = uni.datamaster['sbar'][0][0]
            ys = s[0,0]-s[2,2]
            job = vpsc(stp=stp, mode='YS', texture=self.ftex,
                       fsx=self.fsx, interaction=i)
            job.run();job.pp()
            
            sbar=job.datamaster['sbar'][0]
            dbar=job.datamaster['dbar'][0]
            
            d = np.zeros((len(dbar)))
            sx = np.zeros((len(sbar)))
            sy = np.zeros((len(sbar)))
            for j in range(len(sbar)):
                sx[j] = sbar[j][0,0]
                sy[j] = sbar[j][1,1]
                d[j] = math.atan2(dbar[j][1,1], dbar[j][0,0])

            ## normalization by the uniaxial yield stress
            sx = np.append(sx, [sx[0]])
            sy = np.append(sy, [sy[0]])                
            sxn = sx.copy()/ys
            syn = sy.copy()/ys
            ax1.plot(sx,  sy,  label=i, color='black',
                     ls= linestyle[i],mfc='None',
                     marker=markers[i] )
            ax2.plot(sxn, syn, label=i, color='black',
                     ls= linestyle[i],mfc='None',
                     marker=markers[i] )
            pass
        
        ax1.set_xlabel(r'$\sigma_{x}$', dict(fontsize=20))
        ax1.set_ylabel(r'$\sigma_{y}$', dict(fontsize=20))

        ax2.set_xlim(-1.2, 1.2); ax2.set_ylim(-1.2,1.2)
        ax2.set_xticks(np.arange(-1, 1.001, 0.5))
        ax2.set_yticks(np.arange(-1, 1.001, 0.5))        
        ax2.set_xlabel(r'$\sigma_{x}$', dict(fontsize=20))
        ax2.set_ylabel(r'$\sigma_{y}$', dict(fontsize=20))
        
        ax2.grid()
        
        pass
                      
##
## Example 5
## Equivalent plastic work contour constructor, which is corresponding to
## that of in-plane biaxial tester.
## The original script is from vpsc_param.py
def angle_to_u(angle=None, nprob=None, init=0, fin=None):
    """
    Radian angle to vector (x,y)
    """
    if angle==None and nprob==None: raise IOError
    if angle!=None:
        angle = angle * math.pi /180.
        return math.cos(angle), math.sin(angle)
    elif nprob!=None:
        if fin ==None:fin = 2.*math.pi
        angle = np.linspace(init,fin,nprob)
        return np.cos(angle), np.sin(angle)

class work_contour:
    print 'See __init__.__doc__'
    def __init__(self, interaction=0, stp=10, ict=7, eqc= 0.005, nprob=32,
                 texture='texture/00500.cmb', fsx='sx/Neil.sx',
                 mode='inplanebiaxialfc', ifig=1, dstp=5, init=-45., fin=135.,
                 wwexp=None, jobid=33):
        """
        Work contour plotting.
        
        ---------
        Arguments:
        ---------
        interaction=0
        stp=10
        ict=7
        eqc=0.005
        nprob=32
        texture='texture/00500.cmb',
        fsx='sx/Neil.sx'
        mode='inplanebiaxialfc'
        ifig=1
        dstp=5     : # of incremental work levels
        init=-45
        fin =135
        wwexp=None : list of work levels around which the work is interpolated
        """
        self.dstp = dstp
        self.ifig = ifig
        self.jobs = []

        #----------------------------------------------
        vectors = angle_to_u(
            nprob=nprob,
            init=init*math.pi/180.,
            fin = fin*math.pi/180.)
        for nprb in range(len(vectors[0])):
            u11 = round(vectors[0][nprb],4)
            u22 = round(vectors[1][nprb],4)

            # Job submission
            currentjob = vpsc(mode=mode, stp=stp,
                              ict=ict, eqc=eqc,
                              u11=u11, u22=u22,
                              texture=texture, fsx=fsx,
                              interaction=interaction,
                              jobid=nprb+jobid)

            self.jobs.append(currentjob)
            self.jobs[nprb].run()
            #print self.jobs[nprb].histf[0];raw_input()
            #self.cnt = {} #self.cnt dictionary
            self.cnt = [] #self. cnt is now a list variable
            pass

        #post-process
        self.pp()
        #----------------------------------------------
        
        # Plotting the results
        self.plot(mode='evm',ifig=ifig, wwexp=wwexp)
        self.plot(mode='work',ifig=ifig+1, wwexp=wwexp)


        # self.w ; self.e : work-level and VM-strain-level
        ### perform uniaxial and estimate uniaxial stress at that level
        uni = vpsc(mode='RD', stp=stp,
                   ict=ict, eqc=eqc,
                   texture=texture, fsx=fsx,
                   interaction=interaction,
                   jobid=nprb+jobid+1)
        
        uni.run(); uni.pp()


        ## interpolate the uniaxial stress based on the work (self.w)
        scau = uni.datamaster['sbar'][0]  #sbar
        eps = uni.datamaster['eps'][0]    #eps
        self.s11 = []; self.e11 = []
        for i in range(len(scau)):
            self.s11.append(scau[i,0,0]-scau[i,2,2])
            self.e11.append(eps[i,0,0]-eps[i,2,2])
            pass
        
        self.uniwork = self.__work__(sbar=scau, dbar=eps)
        
        self.uni_s = []; self.uni_e = []
        for i in range(len(self.w)):
            x, y = __interpolate__(
                x=self.s11, y=self.e11,
                z=self.uniwork, value= self.w[i])
            self.uni_s.append(x)
        ## --------------------------------------------------------
        self.plot(mode='work', ifig=ifig+2, inorm=True, wwexp=wwexp)
        #Now, Let's plot the normalized work equivalent cotour!!

    def __work__(self, sbar, dbar):
        """
        Provided the sbar and dbar, calculates the total plastic work through
        cumulative trapezodial method
        """
        if len(sbar)!=len(dbar): raise IOError
        work = []
        for i in range(len(sbar)):
            temp = 0.
            for x in range(len(sbar[i])):
                for y in range(len(sbar[i,x])):
                    temp = temp + sbar[i,x,y]*dbar[i,x,y]
            work.append(temp)
        return np.array(work)

    def pp(self):
        """
        Post process.
        Saves the data to self.cnt : dictionary type variable
        """
        for i in range(len(self.jobs)):
            sig = self.jobs[i].datamaster['sbar'][0]
            eps = self.jobs[i].datamaster['eps'][0]
            evm = self.jobs[i].datamaster['evm'][0]
            svm = self.jobs[i].datamaster['svm'][0]

            ## strain --------------------------------
            # to = []
            # to.append(np.zeros((3,3)))
            # for j in range(len(eps)):
            #     to.append(eps[j])
            #     pass
            # eps = np.array(to)
            e11, e22 = [], []
            for j in range(len(eps)):
                e11.append(eps[j,0,0])
                e22.append(eps[j,1,1]); pass
            e11 = np.array(e11); e22 = np.array(e22)
            #-----------------------------------------

            ## stress --------------------------------
            s11, s22, s33 = [], [], []
            for ii in range(len(sig)):
                s11.append(sig[ii,0,0])
                s22.append(sig[ii,1,1])
                s33.append(sig[ii,2,2])
            s11 = np.array(s11); s22 = np.array(s22); s33 = np.array(s33)
            s11 = s11 - s33; s22 = s22 - s33
            #-----------------------------------------

            # calculates the work 
            wrk = self.__work__(sbar=sig, dbar=eps)
            
            temp = []
            for j in range(len(e11)):
                temp.append([s11[j],s22[j],e11[j],e22[j],wrk[j],evm[j],svm[j]])
            self.cnt.append(temp)
            #self.cnt['job%s'%str(i).zfill(3)] = temp

    def plot(self, mode='work',ifig=1, inorm=False, wwexp=None):
        """
        Plots and writes on to the files

        mode='work', 'evm', 'svm'
        """
        markers = ['+','*',',','.','1','2','3',
                   '4','<','>','D','H','^','_','h']
        linestyle = ['-','--','-.',':',' ','-','--',
                     '-.',':',' ','-','--','-.',':',' ']
        if mode=='work':
            if wwexp==None: w = self.__worklevels__(ini=0.02, step = self.dstp)
            else: w = wwexp

            fig = plt.figure(ifig)
            if len(fig.axes)!=0: ax = plt.gca()
            else: ax = fig.add_axes((0.12,0.15,0.65,0.80),
                                    label=r'$w^{pl}$'+' contour')
            for i in range(len(w)):
                #Interpolate the stress based on the iso work
                s11, s22 = self.__isolevels__(w=w[i])
                #Plot Stresses on the iso plastic work level
                if inorm==True:
                    ## noar
                    s11 = s11/ self.uni_s[i]; s22 = s22/self.uni_s[i]
                    ax.plot(s11, s22,
                            ls=linestyle[i],
                            marker=markers[i],
                            mfc='None', color='black',
                            label=r'w$^{pl}$= %s'%str(round(w[i], 2)))
                    pass
                else: ax.plot(s11, s22,
                              ls=linestyle[i],
                              marker=markers[i],
                              mfc='None', color='black',
                              label=r'w$^{pl}$= %s'%str(round(w[i], 2)))
                #writes on file
                fout('work_at%s'%str(w[i]).zfill(3),s11,s22)
                pass
            self.w = w
            pass

        elif mode=='evm':
            e = self.__evmlevels__(ini=0.001, step=10)
            fig = plt.figure(ifig)
            ax = fig.add_axes((0.12,0.15,0.60,0.80),
                              label=r'$\varepsilon^{VM}$'+' contour')
            for i in range(len(e)):
                #interpolates
                s11, s22 = self.__isolevels__(e=e[i])

                #plot stress on the iso evm level
                ax.plot(s11,s22,
                        ls=linestyle[i],
                        marker=markers[i],
                        color='black',
                        label=r'$\bar{\epsilon}^{VM}$= %s'%str(round(e[i],3)))

                #writes on file
                fout('evm_at%s'%str(e[i]).zfill(3),s11,s22)
                pass
            self.e = e
            pass

        ax.set_aspect('equal')
        
        ## axis right dimension
        ax.legend(bbox_to_anchor=(1.02,1), loc=2)
        ax.plot([0], [0], '+' ,color='black')
        ax.set_aspect('equal')        
        if inorm==True:
            ax.set_xlabel(r'$\sigma_{RD} / \bar{\sigma}$',
                          dict(fontsize=20))
            ax.set_ylabel(r'$\sigma_{TD} / \bar{\sigma}$',
                          dict(fontsize=20))
            ax.set_ylim(-0.1,);ax.set_xlim(-0.1,)
            pass
        else:
            ax.set_xlabel(r'$\sigma_{RD}$ [MPa]', dict(fontsize=20))
            ax.set_ylabel(r'$\sigma_{TD}$ [MPa]', dict(fontsize=20))
            ax.set_ylim(-10,);ax.set_xlim(-10)
            pass

            
        ## axis off ##
        #ax.set_axis_off()
        #ax.set_xticks(()); ax.set_yticks(())            
            
    def __isolevels__(self, w=None, e=None):
        """
        Makes the segment of the iso-level
        The argument can be either w or e
        """
        s11 = []
        s22 = []
        for i in range(len(self.jobs)):
            if w!=None:
                s1, s2 = self.__interpolate__(ijob=i, w=w)
            elif e!=None:
                s1, s2 = self.__interpolate__(ijob=i, e=e)
            if s1==0 and s2==0: pass
            else: s11.append(s1); s22.append(s2)
        return np.array(s11), np.array(s22)

    def __interpolate__(self, ijob=0, w=None, e=None):
        """
        Linearly interpolate thes s11 and s22
        """
        s11, s22, e11, e22, work, evm, svm = [],[],[],[],[],[],[]
        for i in range(len(self.cnt[ijob])):
            s11.append(self.cnt[ijob][i][0]) # s11
            s22.append(self.cnt[ijob][i][1]) # s22
            e11.append(self.cnt[ijob][i][2]) # e11
            e22.append(self.cnt[ijob][i][3]) # e22
            work.append(self.cnt[ijob][i][4]) # work
            evm.append(self.cnt[ijob][i][5])
            svm.append(self.cnt[ijob][i][6])
            
        s11, s22, e11, e22, work, evm, svm= __makenp__(s11, s22,
                                                       e11, e22,
                                                       work,
                                                       evm, svm)
        # debugging purpose -----------
        if w!=None and w>max(work):
            print 'ijob =', ijob,'\n'
            print 'w=', w, '\n'
            print 'work = ', work
            raw_input(); raise IOError
        if e!=None and e>max(evm):
            print 'ijob =', ijob,'\n'
            print 'e=', e, '\n'
            print 'evm = ', evm
            raw_input(); raise IOError
        # -----------------------------
        if w!=None:
            s1, s2 = __interpolate__(x=s11, y=s22,
                                     z=work, value=w)
        elif e!=None:
            try:
                s1, s2 = __interpolate__(x=s11, y=s22,
                                         z=evm, value=e)
            except:
                print 'evm:\n', evm
                print 'e: ', e, '\n'
                raw_input()
                raise IOError
        else: raise IOError
        
        return np.array(s1), np.array(s2)

    def __worklevels__(self, ini=0.01, step=10):
        """
        Finding the automatic work levels.

        Argument:
        ini = 0.01 : initial work level is fixed to be 0.01
        step = 10 : number of work level
        """
        mx = []
        for i in range(len(self.cnt)):
            cjob = np.array(self.cnt[i])
            cjob.copy().transpose()[0] #s11
            cjob.copy().transpose()[1] #s22
            cjob.copy().transpose()[2] #e11
            cjob.copy().transpose()[3] #e22
            plw = cjob.copy().transpose()[4] #work
            #print '\nplw=', plw; raw_input()
            mx.append(max(plw)) #the last work == the maximum work

        mx_work = min(mx) #The smallest work
        w = np.linspace(ini, mx_work*0.99 , step) #spacing the work
        return w

    def __evmlevels__(self, ini=0.01, step=10):
        """
        Finding the automatic von Mises strain levels

        ini = 0.01: initial VM strain is fixed to be 0.01
        step = 10: number of VM strain level
        """
        mx = []
        for i in range(len(self.cnt)):
            evm = np.array(self.cnt[i]).copy().transpose()[5]
            mx.append(max(evm))
        mx_evm = min(mx)
        e = np.linspace(ini, mx_evm*0.99, step) #spacing the VM strain
        return e
#
#
## Example 6 (ICOTOM16 -2)
## Simple-shear forward-backward loading (cyclic)
## Fitting procedure to the experimental counterpart
## 
class ss:
    def __init__(self, hl=1.0, hb=1.0, hp=1.0,
                 tau0=1., tau1=0.5, thet0=1., thet1=0.5,
                 interaction=3, ifig=1, ngrain=100,
                 ftex=None, fsx=None, iopsysx=0, eqc=0.01,
                 nneigh=0 ):
        """
        Simple shear test

        Things to do
        1. Selects the experimental data
           -- can be done in self.__exp__
        2. Predefines the strain
        3. Presets the axes labels and legend.


        Arguments:
          hl, hb, hp : Latent, Bauschinger and coplanar coefficient
          tau0, tau1, thet0, thet1 : Voce hardening parameters
          interaction
          ifig,
          ngrain
          ftex
          fsx
          iopsysx = 0 (non-directional)
          eqc=0.01
          nneigh : number of neighbouring (coupling) grains
        """
        if ftex!=None: self.ftex=ftex
        elif ftex==None:
            # random texture generation
            ngrain = ngrain
            tex = re(ngrain=ngrain)
            self.ftex = '%i.tex'%ngrain
            tex.write(self.ftex)
        if fsx!=None: self.fsx = fsx
        elif fsx==None:
            hl #latent hardening
            hb #Bauschinger hardening
            tau0=tau0; tau1=tau1; thet0=thet0; thet1=thet1 #voce- hardening
            self.fsx='temp_ss.sx'
            sx_maker.cubic(
                filename=self.fsx,
                hii=1.0, hij=hl, hb=hb, hp=hp, #copl=False,
                tau0=tau0, tau1=tau1,
                thet0=thet0, thet1=thet1, header='**ss_parameteric study: ',
                iopsysx=iopsysx, #dependency='Bauschinger latent',
                nrsx=20, b=[[1,1,0]], n=[[1,1,1]],    #Face-centered
                #b=[[1,1,1],[1,1,1],[1,1,1]],n=[[1,1,0],[1,1,2],[1,2,3]], #Body-centered
                )
            
        # cyclic simple shear
        eps_pl = [0.14, -0.22, 0.22, -0.22, 0.4]
        #eps_tot = abs(np.array(eps_pl)).sum() / 1.6
        eps_tot = 0.66

        ec = np.array(())#[0])
        sc = np.array(())#[0])

        for n in range(len(eps_pl)): #cyclic
            stp = abs(int(eps_pl[n]/eqc))
            if eps_pl[n]<0: mode='SSR'
            elif eps_pl[n]>0: mode='SS'
            else: raise IOError

            #print 'stp =', stp;raw_input()
            if n==0: irecover=0
            else: irecover=1

            ## cyclic shear
            ### parametric 

            if n>3: iupdate=[0,0,1,0]
            else: iupdate=[1,1,1,0]
            #iupdate=[1,1,1,0]
            rst = self.__ss__(irecover=irecover, stp=stp, eqc=eqc,
                              isave=True, interaction=interaction,
                              ftex=self.ftex, fsx=self.fsx, mode=mode,
                              nneigh=nneigh,iupdate=iupdate)

            if n==0:
                acc = 0.
                gcf = plt.figure(ifig, figsize=(9.6,4.8))
                #ax01 = gcf.add_subplot(111,label='ss')
                ax01 = gcf.add_axes((0.1,0.1,0.55,0.8), label='ss')
                pass
            if n==0:
                #E12 = np.array(np.append([0,], abs(rst['e12'])))
                E12 = np.array(abs(rst['e12']))
                evm = np.array(np.append([0,], rst['evm']))
            else:
                #E12 = np.append([0], abs(rst['e12'])) + acc
                E12 = abs(rst['e12']) + acc
                evm = np.append([0], rst['evm']) + acc_vm
                pass

            acc = abs(E12[-1])
            acc_vm = evm[-1]
            s12 = abs(np.array(rst['s12']))

            # print 's12 \n', s12
            # print 'E12 \n', E12
            # print 'len(s12) \n', len(s12)
            # print 'len(E12) \n', len(E12)
            # raw_input()            

            ec = np.append(ec, E12)
            sc = np.append(sc, s12)
            #ax01.plot(E12, s12, 'r.') #fixed color
            
            pass
            # run the vpsc (if initial ignores post-morten)
            # save the post-morten
            # save the flow stress and strain
        #print 'eq\n',ec,'\nsc',sc; raw_input()

        ax01.plot(ec, sc, 'ko', mfc='w', label='VPSC cyclic')
        ax01.text(1.5, 200.0,
                  r'$\tau_0$: %5.2f $\tau_1$: %5.2f'%(tau0, tau1),
                  dict(fontsize=15))
        ax01.text(1.5, 150.0,
                  r'$\theta_0$: %5.2f $\theta_1$: %5.2f'%(thet0, thet1),
                  dict(fontsize=15))
        ax01.text(1.5, 100,
                  r'$H_{L}$: %3.1f $H_{B}$: %3.1f'%(hl, hb),
                  dict(fontsize=15))
        # coplanar coefficient
        # ax01.text(1.5, 50,
        #           r'$H_{P}$: %3.1f'%(hp),
        #           dict(fontsize=15))

        stp = int(eps_tot/eqc)
        print 'stp:', stp, ' eqc=',eqc

        ## simple shear
        rst = self.__ss__(irecover=0, stp=stp, eqc=eqc, nneigh=nneigh,
                          isave=False, interaction=interaction,
                          ftex=self.ftex, fsx=self.fsx, mode='SS')
        
        #e12 = np.append([0,], rst['e12'])
        e12 = np.array(rst['e12'])
        s12 = np.array(rst['s12'])
        ax01.plot(e12, s12, marker='+',mec='gray',mfc='w',ls='None',
                  label='VPSC monotonic')
        
        ax01.set_ylim(0.,)
        ax01.set_ylabel(r'$\tau$'+' [MPa]', dict(fontsize=20))
        ax01.set_xlabel(r'$\gamma$', dict(fontsize=20))

        self.__exp__(ax=ax01, mode='ss')

        ## uniaxial tension
        rst = self.tensile(irecover=0, stp=20, eqc=0.01,
                           isave=False, interaction=interaction,
                           ftex=self.ftex, fsx=self.fsx, mode='RD',nneigh=nneigh)

        #e11 = np.append([0,], rst['e11'])
        e11 = np.array(rst['e11'])
        s11 = np.array(rst['s11']-rst['s33'])
        gf = plt.figure(ifig+1)
        ax02 = gf.add_axes((0.1,0.1,0.7,0.6), label='unix_')
        ax02.plot(e11,s11, label='VPSC uniaxial')
        self.__exp__(ax=ax02, mode='uni')
        
        ax01.legend(bbox_to_anchor=(1.05,1), loc=2)
        plt.show() #Show!
        pass

    def __exp__(self, ax=None, mode='ss'):
        """
        Experimental counterparts
        It is guided to manually revised upon situations.

        arguments
        ax = None
        mode = 'ss' or 'uni'
        """
        ## simple shear
        if mode=='ss':
            files = glob.glob('ref/bh/5_cy*.pp')  # BH STEEL
            # files = glob.glob('ref/fss/6_cy*.pp')   # FSS
            # #for i in range(len(files)):
            i=0
            data = np.loadtxt(files[i], skiprows=1)
            data = data.transpose()
            ax.plot(data[0], data[1], label='EXP cyclic',color='k')


            files = glob.glob('ref/bh/1_mono*.pp')  # BH steel
            
            #files = glob.glob('ref/fss/2_mono*.pp')  # FSS
            #for i in range(len(files)):
            i=0
            data = np.loadtxt(files[i], skiprows=1)
            data = data.transpose()
            ax.plot(data[0], data[1], label='EXP monotonic', color='gray')

            # prestrained
            # files = glob.glob('ref/fss/4_pre.pp')  # FSS
            # i=0
            # data = np.loadtxt(files[i], skiprows=1)
            # data = data.transpose()
            # ax.plot(data[0], data[1], label='EXP prestained', color='k') 

        # ## Uniaxial
        elif mode=='uni':
        #     files = glob.glob('ref/*RD1*')
        #     for i in files:
        #         data = np.loadtxt(i, skiprows=2,dtype='str')
        #         data = data.transpose()
        #         ax.plot(map(float,data[1]), map(float,data[2]),'--', label=i,color='r')
        #         pass
            pass
        else: raise IOError

    def __plot__(self):
        """
        Plots the simulated results
        """
        pass

    def __ss__(self, stp=None, irecover=0, eqc=None,
               ftex=None, fsx=None, interaction=3, isave=False,
               mode='SS', nneigh=0, iupdate=[1,1,1,0], jobid=0
               ):
        """
        Set off a simple shear simulation.

        arguments
        stp: number of the incremental step
        irecover: Read grain states (1) or not (0)
        eqc: incremental step size
        ftex: texture file name
        fsx: single crystal file name
        interacton: grain interaction(fc, affine, sec, neff, tangent, so)
        isave: True (save the grain states at step=stp) False (not doing so)
        mode:'SS' : simple shear mode indicator to the VPSC class
        nneigh: number of neighbouring grains
        iupdate: [1,1,1,0] update ori, grain shape, hardening, phase transf
        jobid: job id argument to the VPSC class
        """
        if isave==False: isave=0
        elif isave==True: isave=stp
        ### VPSC job submit
        try: 
            job = vpsc(irecover=irecover, isave=stp, stp=stp, mode=mode, eqc=eqc,
                       texture=ftex, fsx=fsx, interaction=interaction, nneigh=nneigh,
                       iupdate=iupdate, jobid=jobid)
        except: raise IOError
        ### ---------------
        
        job.run(); job.pp()
        if isave!=0:
            shutil.move(
                'POSTMORT_$s.OUT POSTMORT_%s.IN'%
                (str(jobid).zfill(4),
                 str(jobid).zfill(4))
                )
            pass
        os.system('clear')
        sbar = job.datamaster['sbar'][0]
        dbar = job.datamaster['dbar'][0]
        eps  = job.datamaster['eps'][0]
        svm  = job.datamaster['svm'][0]
        evm  = job.datamaster['evm'][0]
        agre = job.datamaster['lgr']

        s11 = []; s22 = []; s33 = []
        s12 = []; s23 = []; s13 = []
        
        d11 = []; d22 = []; d33 = []
        d12 = []; d23 = []; d13 = []        
        
        e11 = []; e22 = []; e33 = []
        e12 = []; e23 = []; e13 = []

        for i in range(len(sbar)):
            s11.append(sbar[i][0,0])
            s22.append(sbar[i][1,1])
            s33.append(sbar[i][2,2])
            s12.append(sbar[i][0,1])
            s23.append(sbar[i][1,2])
            s13.append(sbar[i][0,2])
            
            d11.append(dbar[i][0,0])
            d22.append(dbar[i][1,1])
            d33.append(dbar[i][2,2])
            d12.append(dbar[i][0,1])
            d23.append(dbar[i][1,2])
            d13.append(dbar[i][0,2])

        for i in range(len(eps)):  #eps usually one iteration less than sbar or dbar
            e11.append(eps[i][0,0])
            e22.append(eps[i][1,1])
            e33.append(eps[i][2,2])
            e12.append(eps[i][0,1])
            e23.append(eps[i][1,2])
            e13.append(eps[i][0,2])

        master = {} # master data dictionary
        s11, s22, s33, d11, d22, d33, e11, e22, e33 = __makenp__(s11,s22,s33,
                                                                 d11,d22,d33,
                                                                 e11,e22,e33)

        s12, s23, s13, d12, d23, d13, e12, e23, e13 = __makenp__(s12,s23,s13,
                                                                 d12,d23,e13,
                                                                 e12,e23,e13)
        evm = np.array(evm)

        master['s11'] = s11; master['s22'] = s22; master['s33'] = s33 #stress
        master['s12'] = s12; master['s23'] = s23; master['s13'] = s13
        
        master['d11'] = d11; master['d22'] = d22; master['d33'] = d33 #strain rate
        master['d12'] = d12; master['d23'] = d23; master['d13'] = d13
        
        master['e11'] = e11; master['e22'] = e22; master['e33'] = e33 #strain history
        master['e12'] = e12; master['e23'] = e23; master['e13'] = e13

        master['evm'] = evm; master['svm'] = svm
        #master['cr'] = cr(bins=bins)  #cr method
        return master
    
    def tensile(self, irecover=0, stp=None, eqc=None,
                axes= [[1.,1.,1.]], angle=[[0.,0.,0.]],
                iupdate=[1,1,1,0],
                ftex=None, fsx=None, interaction=3,
                isave=False,bins=20, mode=None,
                u11=None,u22=None, nneigh=0):
        """
        Performs a uniaxial (or biaxial) tensile test. (axiality)
        Depending on the irecover and isave
        interacts with the postmorten files

        The simulated results are returned properly

        Must chose what to be returned among all the keys
        in the datamaster dictionary

        ## Here I chose,
        1. Uniaxial(or biaxial) tension stress tensor
        2. Uniaxial(or biaxial) tension strain tensor
           * All of the tensors will be that of 3x3 2nd order,
             and they are macroscopic entities.
        3. Polycrystalline aggregate
        4. CRSS distribution # May need be normalized...?
        5. Returns R-value!
        """
        if stp==None: raise IOError
        if irecover==1:
            pass
        #     if os.path.isfile("POSTMORT.OUT"):
        #         os.rename("POSTMORT.OUT", "POSTMORT.IN")
        #     else: raise IOError
            # else:
            #     if os.path.isfile("POSTMORT.IN"): pass
            #     else: raise IOError
        else: pass # No recovery by reading the postmort.in file

        if isave==False: isave=0
        elif isave==True: isave=stp #Always the last step is to be saved.
        else: raise IOError

        if stp==0: raise IOError

        #os.system('ls') ;raw_input('enter to proceed>>')

        ## decision on the mode
        # it can be "RD", "TD", "inplanebiaxialfc"

        job = vpsc(irecover=irecover, isave=stp, stp=stp, mode=mode, eqc=eqc,
                   texture=ftex, fsx=fsx, init_eld_rat=axes, init_eul_axe=angle,
                   interaction=interaction, iupdate=iupdate, u11=u11, u22=u22,
                   nneigh=nneigh)

        job.run(); job.pp()

        if isave!=0:
            if os.name!='posix':
                print 'Unexpected os.system found'
                raw_input('Enter to proceed>>')
                raise IOError
            os.system("cp POSTMORT.OUT POSTMORT.IN")
            os.system("rm -f POSTMORT.OUT")
        
        os.system('clear')
        #os.system('ls'); raw_input('enter to proceed >>>')        
        sbar = job.datamaster['sbar'][0]  ## Macroscopic stress tensor
        dbar = job.datamaster['dbar'][0]  ## Macroscopic strain tensor
        eps  = job.datamaster['eps'][0]   ## Eps
        svm  = job.datamaster['svm'][0]  ## von Mises stress        
        evm  = job.datamaster['evm'][0]  ## von Mises strain
        agre = job.datamaster['lgr']      # the most recent aggregate's COD

        s11 = []; s22 = []; s33 = []
        s12 = []; s23 = []; s13 = []
        
        d11 = []; d22 = []; d33 = []
        d12 = []; d23 = []; d13 = []        
        
        e11 = []; e22 = []; e33 = []
        e12 = []; e23 = []; e13 = []

        for i in range(len(sbar)):
            s11.append(sbar[i][0,0])
            s22.append(sbar[i][1,1])
            s33.append(sbar[i][2,2])
            s12.append(sbar[i][0,1])
            s23.append(sbar[i][1,2])
            s13.append(sbar[i][0,2])
            
            d11.append(dbar[i][0,0])
            d22.append(dbar[i][1,1])
            d33.append(dbar[i][2,2])
            d12.append(dbar[i][0,1])
            d23.append(dbar[i][1,2])
            d13.append(dbar[i][0,2])

        for i in range(len(eps)):  #eps usually one iteration less than sbar or dbar
            e11.append(eps[i][0,0])
            e22.append(eps[i][1,1])
            e33.append(eps[i][2,2])
            e12.append(eps[i][0,1])
            e23.append(eps[i][1,2])
            e13.append(eps[i][0,2])

        master = {} # master data dictionary
        s11, s22, s33, d11, d22, d33, e11, e22, e33 = __makenp__(s11,s22,s33,
                                                                 d11,d22,d33,
                                                                 e11,e22,e33)

        e12, e23, e13, d12, d23, d13, e12, e23, e13 = __makenp__(e12,e23,e13,
                                                                 d12,d23,e13,
                                                                 e12,e23,e13)
        evm = np.array(evm)

        master['s11'] = s11; master['s22'] = s22; master['s33'] = s33 #stress
        master['s12'] = s12; master['s23'] = s23; master['s13'] = s13
        
        master['d11'] = d11; master['d22'] = d22; master['d33'] = d33 #strain rate
        master['d12'] = d12; master['d23'] = d23; master['d13'] = d13
        
        master['e11'] = e11; master['e22'] = e22; master['e33'] = e33 #strain history
        master['e12'] = e12; master['e23'] = e23; master['e13'] = e13
        
        master['R']   = d22/d33  #R-value (instantaneous)
        master['ag']  = agre
        master['evm'] = evm; master['svm'] = svm

        #print glob.glob('crss*.out');raw_input()
        
        master['cr'] = cr(bins=bins)  #cr method
        return master

##
## Example 7 YS surface using parallel python package
## 
class YieldSurface:
    def __init__(self):
        """
        """
        pass
    
    def __vpscset__(self):
        """
        VPSC setting and assigns a job id
        """
        pass

    def __run__(self):
        """
        runs!
        """
        pass
    pass


def vp_run(mode='RD', jobid=0, stp=100, eqc=0.001,
           interaction=3, ftex='texture/08000.cmb'):
    import vp_f2py
    reload(vp_f2py) #debugging purpose on vp_f2py which has been simultaneous going
    job = vp_f2py.vpsc(mode=mode, jobid=jobid,
                       stp=100, eqc=0.001,
                       interaction=3,
                       texture='texture/08000.cmb')
    data = job.run()
    return data['svm'][0][-1]


try: import pp; reload(pp)
except: print 'parellel python package is not found'


def ppex():
    """
    VPSC - PP examplary run
    """
    ppservers=()
    #job_server = pp.Server(ppservers=ppservers)
    ncpus=3
    job_server = pp.Server(ncpus, ppservers=ppservers)

    jobs = []

    results = []
    for i in range(4):
        f = job_server.submit(
            vp_run, #function to be excuted
            ( #args to excuted function
                ['RD','BU','DD','SS'][i],
                i,
                3,
                0.01,
                3,
                'texture/01000.cmb',
                )
            )
        print 'f:', f
        results.append(f())
    print 'results=', results
    job_server.wait()
    job_server.print_stats()
    
    #return results
    pass


#    for x,f in results:
#        val = f()

    #    job_server = pp.Server(ppservers=ppservers)        
        

    
########################################################
# Example 8 Optimization of the hardening parameters  #
#                using scipy optimization package      #
#                                                      #
# The list of definitions as below                     #
#                                                      #
# def _interpolate_(x,y,x0)                            #
# def func(parameters=[], fsx, ftex)                   #
# def voce_optimization()                              #
########################################################

"""
Steps to follow to optimize the hardening parameters.

First of all, the experimental data file must be prepared to which
the hardening curve is optimized.

Then, just follow below.


>>> import vpsc_examples
>>> vo = vpsc_examples.voce_optimization

         ### tau0, tau1, thet0, and thet1 are
         ### the initial guess on the Voce hardening parameters.
         
>>> xopt = vo(tau0=?, tau1=?, thet0=?, thet1=?, eqc=0.01,
              expfile = 'path and file name',

              ftex = "Give texture file otherwise iso random file",
              ## or one can provide the curve directly as below.
              ssexp = default is 'None'
              
              loadingmode='Loading mode: RD, TD ..',
              interaction=3,
              ifig = 1)

              
*** One may add or subtract the parameter
*** by tweaking the original code.

"""
def _interpolate_(x, y, x0):
    """
    Estimates new y by linearly interpolate the x-y curve
    at the given point, p
    x0: 1-d array
    """
    y0 = 0
    if x0> max(x) or x0< min(x):
        print 'The given x0 is out of the allowed range'
        print 'The given x0: %f'%x0
        print 'max(x): %f'%max(x), 'min(x): %f'%min(x)
        raise IOError
    
    for i in range(len(x)):
        if x[i]>x0: break
        else: pass
        pass
    
    x1 = x[i-1]; x2 = x[i]
    y1 = y[i-1]; y2 = y[i]
    ## slope
    slope = (y2-y1) / (x2-x1)
    
    ## linear interpolation of interest y0 point.
    y0 = slope*(x0-x1) + y1
    return y0

def _interpolate_2(x,y,x0):
    """
    Estimates new y be linearly interpolate the x-y curve
    at the given point x0, even if x0 is out of the range.
    Whenever x0 is out of the x range, extrapolate the data.
    """
    if x0 > max(x):
        x1 = x[-2]; x2 = x[-1]
        y1 = y[-2]; y2 = y[-1]
        ## slope
        slope = (y2 - y1) / (x2 - x1)

        ## linear extrapolate
        y0 = slope * (x0 - x1) + y1
        return y0
    elif x0 < min(x):
        x1 = x[0]; x2 = x[1]
        y1 = y[0]; y2 = y[1]
        ## slope
        slope = (y2 - y1) / (x2 - x1)
        y0 = slope * (x0 - x1) + y1
        return y0
    else:
        # call _interpolate_
        return _interpolate_(x,y,x0)

def voce_func(
    parameters=[1., 0.5, 1.0, 0.5], eqc=0.01,
    ftex='100.tex', mode='RD', ifig=1,
    expfile = 'expfile', ngrain = 100,
    rstflag='fit', ssexp = None, jobid=0,
    
    #0:FC, 1:affine, 2:sec, 3:neff=10, 4:tangent, 5:SO
    interaction = None,
    ):
    """
    Preliminary func for fitting the curve.
    Acts as an non-linear function to be minimized to find
    the designated hardening parameters.

    if rstflag=='fit': Works as a fitting function with
    returning the standard deviation
    """
    if interaction==None: raise IOError,'interaction should be given'
    ## Experimental results ------------------
    if ssexp==None:
        ## number of rows to be skipped is fixed to be 1!!
        SE = np.loadtxt(expfile, skiprows=0)
        SE = SE.transpose()
        Eexp = SE[0]#np.array([0.01,   0.1,   0.2,  0.3])
        Sexp = SE[1]#np.array([3, 3.2, 3.3, 3.35])
    else:
        Eexp = ssexp[0]
        Sexp = ssexp[1]
    ## ---------------------------------------

    ## hardening parameters and single crystal file making
    tau0 = parameters[0]; tau1 = parameters[1]
    thet0 = parameters[2]; thet1 = parameters[3]
    
    fsx = 'temp%s.sx'%str(jobid).zfill(3)
    sx_maker.cubic(filename=fsx, tau0=tau0, tau1=tau1,    
                   thet0=thet0, thet1=thet1, iopsysx=1,
                   dependency='iso')

    ## ----------------------------------------------------
    ## Texture making if ftex = None
    if ftex!=None: pass
    elif ftex==None:
        # random texture generation
        ngrain = ngrain
        tex = re(ngrain=ngrain)
        ftex = '%i_%s.tex'%(ngrain,str(jobid).zfill(3))
        tex.write(ftex)
        pass

    ## Simulate in accordance wit the vpsc mode ---------------
    stp = int(max(Eexp)*1.20/eqc) # 20% overshoot!
    job = vpsc(
        mode=mode, fsx=fsx, texture= ftex,
        stp=stp, eqc=eqc, jobid=jobid,
        interaction=interaction
        )
    job.run(); job.pp()
    ## --------------------------------------------------------

    ## Stress strain curve ---------
    eps = job.datamaster['eps'][0]    
    if mode=='RD' or mode=='rd':
        if interaction==0:
            scau = job.datamaster['sdev'][0]
            pass
        else:
            scau = job.datamaster['scau'][0]
    elif mode=='SS': scau = job.datamaster['sdev'][0]
    else:
        print "Beside 'RD' and 'SS', loadingode is not prepared"
        raise IOError,'Unexpected mode'
    
    s11 = []; e11 = []
    if mode=='RD' or mod=='rd': xi = 0; yi = 0
    elif mode=='SS' : xi = 1; yi = 0
    else: print "Inappropriate loading mode"; raise IOError
    for i in range(len(scau)):
        #print '\n', scau[i], eps[i], '\n\n\n'
        if mode=='RD' or mode=='rd':
            s11.append(scau[i,xi,yi] - scau[i,2,2])
        elif mode=='SS':
            s11.append(scau[i,xi,yi])
        else:
            print 'given mode=', mode
            raise IOError, 'Unexpected mode'
        e11.append(eps[i,xi,yi])
        #print 's11=', s11, '\ne11=', e11, '\n'
        #raw_input()
        pass
    ## -----------------------------
    
    ## Interpolate the experimental results --
    # for example, as follow

    ## fits the simulated results along the experimental strain axis
    ys = []
    for i in range(len(Eexp)):
        ys.append(_interpolate_(x=e11, y=s11, x0=Eexp[i]))
        pass
    ## -------------------------------------------------------------
    
    ## Estimate the error to be returned ---
    ys = np.array(ys)
    diff = abs(ys - Sexp)
    stdev = (diff**2).sum()
    stdev = math.sqrt(stdev)
    ## -------------------------------------
    if rstflag=='fit': return stdev
    elif rstflag=='regular': return e11, s11
    else: raise IOError

def voce_print(parameters = []):
    """
    A call back function
    parameters = []
    """
    if os.name=='posix':
        os.system('clear')
    print '\n\n\nparameters =', parameters, '\n\n\n'

    ### ifig
    # fig = plt.gcf()
    # ax = plt.gca()
    pass

def voce_optimization(
    xtol=0.000001, maxiter=400, ftex=None,
    tau0=1., tau1=0.5, thet0=1.2, thet1=0.5, eqc=0.01,
    algorithm='simplex',
    expfile=None, loadingmode=None,
    ifig=1, ssexp = None, jobid=0, interaction=None,
    ax=None,
    ):
    """
    Optimization over the modified Voce hardening parameters
    (tau0, tau1, thet0, thet1) over a particular loading.

    ## arguments
    xtol = 0.000001 :tolerance
    maxiter = 4000  :Maximum # of iteration
    ftex = None     :Texture file
    algorithm='simplex', 'bfgs'
    expfile=None
    loadingmode='None' -> 'RD'
    ifig=1,
    ssexp=None (complementary to expfile, array of stress-stain curve)
    jobid
    interaction
    """
    from scipy.optimize import fmin
    from scipy.optimize import fmin_bfgs
    #tau0 = 1; tau1 = 0.5; thet0 = 1.2; thet1 = 0.5
    parameters = np.array([tau0, tau1, thet0, thet1])
    
    ## Optimize the parameters
    if ftex==None: ftex = '100.tex'
    if loadingmode==None:
        raise IOError, 'loadingmode is missing'

    if expfile!=None and ssexp!=None:
        raise IOError, 'Both expfile and ssexp are given'
    
    if interaction==None: raise IOError, 'input interaction!'

    ### Optimization ------------------------------------------------
    #expfile = 'expfile'
    ngrain = None
    #loadingmode='RD'
    ### Nelder-Mead simplex algorithm (fmin)
    if algorithm=='simplex':
        xopt = fmin(voce_func,  parameters,
                    args=(eqc, ftex, loadingmode, ifig,
                          expfile, ngrain, 'fit', ssexp, jobid, interaction),
                    xtol=xtol,
                    maxiter=maxiter,
                    maxfun=maxiter,
                    callback=voce_print,
                    retall=True,# full_output=True #flags
                    )
        pass
    elif algorithm=='bfgs':
        
        print '##############################################################'
        print '    Broyden-Fletcher-Goldfarb-Shanno algorithm (fmin_bgfs)'
        print '   It is still unstable. End-user is guided to use simplex'
        print '##############################################################'
        return -1
        xopt = fmin_bfgs(
                    voce_func, parameters,
                    args=(eqc, ftex, loadingmode, ifig,
                          expfile, ngrain, 'fit', ssexp, jobid, interaction),
                    fprime= None,
                    maxiter=maxiter,
                    callback=voce_print,
                    gtol = 10e-8,
                    epsilon= 10e-5,
                    retall=True,# full_output=True #flags
                    )
        ## currently no monitoring nor post-process on bfgs
        ## write in xopt's dimension for reference, with which
        ## I can't remember!!! and get confused all the time!
        return xopt
    else: raise IOError


    # xopt = fmin_bfgs(
    #     ipbiaxial,[x0],
    #     args=(eqincr, ftex, fsx, jobid, theta0),fprime=None,
    #     maxiter=maxiter,
    #     gtol=10e-8,
    #     epsilon=10e-5,
    #     retall=True
    #     )

    
    ### -------------------------------------------------------------
    
    #return xopt
    # print 'Initial parameters are as below'
    # print parameters
    # print 'Final parameters are as below'
    # print xopt[1][-1]
    # print 'total parameters are as below'
    # print xopt

    ## experimental results

    ## -----------------------------------------------------------------------
    ## figure template
    l = 0.15  #left blank
    b = 0.15  #bottom blank
    w = 0.75  #width
    h = 0.80  #height
    fullscale = 1.0
    fig_width=7
    fig_height=5
    fig_rel_size= 1.0
    nax = 1
    nax = float(nax)
    
    el = l/nax; ew = w/nax
    incr = fullscale / nax
    fig_width, fig_height = np.array([fig_width,fig_height]) * fig_rel_size
    
    fig0 = plt.figure(ifig, figsize=(fig_width*nax,fig_height) )
    fig1 = plt.figure(ifig+1, figsize=(fig_width*nax,fig_height) )
    ax0 = fig0.add_axes((el, b,ew,h))
    ax1 = fig1.add_axes((el, b,ew,h))
    ## -----------------------------------------------------------------------    
    
    if ssexp==None:
        curve = np.loadtxt(expfile, skiprows=1)
        curve = curve.transpose()
    else: curve=ssexp
    
    ax0.plot(curve[0], curve[1], 'x', label='exp', color='black')

    ## initial
    e11, s11 = voce_func(xopt[1][0], eqc, ftex, loadingmode,
                         ifig, expfile, ngrain,
                         'regular', ssexp, interaction=interaction)
    ax0.plot(e11, s11, 'o', label='initial', color='black', mfc='None')

    ## optimized
    print "\n\nxopt[1][-1]\n" , xopt[1][-1], '\n'
    ee11, ss11 = voce_func(xopt[1][-1], eqc, ftex, loadingmode,
                           ifig, expfile, ngrain,
                           'regular', ssexp, interaction=interaction)
    ax0.plot(ee11, ss11, 'd', label='optimized',color='gray', mfc='None')

    # if there's passed-in axis plot on it too
    if ax!=None: ax.plot(ee11, ss11, 'd',
                         color='gray', mfc='None')
        

    ## hardening parameters
    pa = np.array(xopt[1]).transpose()

    labels = [r'$\tau_0$', r'$\tau_1$', r'$\theta_0$', r'$\theta_1$']


    for i in range(4): ax1.plot(
        pa[i], ['x','o','d','--'][i],
        label=labels[i], color='black'
        )
    
    ax0.set_ylim(0.,)
    ax0.legend(loc='best');ax1.legend(loc='best')
    ax0.set_ylabel(r'$\sigma$ [MPa]', dict(fontsize=20))
    ax0.set_xlabel(r'$\varepsilon$', dict(fontsize=20))
    ax1.set_ylabel('[MPa]', dict(fontsize=20))
    ax1.set_xlabel('iteration', dict(fontsize=17))
    
    
    print'pa=', pa
    return xopt #, e11, s11, ee11, ss11

##############################################################################
# Example 9 Checking the dislocation density based hardening implementation #
#                                                                            #
# After Rauch et al's dislocation density based hardening model              #
#       - 2011-07-21                                                         #
##############################################################################
def ddh2(fsx='disl.sx', stp=10, eqc=0.005):
    """
    Test the dilsocation hardening incorporated in VPSC
    with some hardwired situations.
    def ddh is called for each situations
    """
    from cmb import random as rgr
    import matplotlib.pyplot as plt
    plt.ioff()
    gr0 = rgr(ngrain=100)
    gr1 = np.loadtxt('100_0.5.cmb', skiprows=4).transpose()
    ddh(mode='ssfr', fsx=fsx,ifig=1,gr=gr0)
    ddh(mode='ssfr', fsx=fsx,ifig=1,gr=gr1)
    ddh(mode='tc', fsx=fsx,ifig=1,gr=gr0)
    ddh(mode='tc', fsx=fsx,ifig=1,gr=gr1)
    fig = plt.gcf()
    ax = plt.gca()        
    pass

def ddh(mode='ssfr', fsx='disl.sx',stp=10, eqc=0.005, ifig=1, gr=None):
    """
    Test the dislocation density hardening model (DDH)
    on the loading modes of 'ssfr' or 'tc' of the python-wrapped vp_f2py.vpsc
    """
    if gr==None: print 'Input gr'; return 0

    import matplotlib.pyplot as plt
    import vp_f2py
    import os
    plt.ioff()  #plt interactivity is turned off due
    myjob = vp_f2py.vpsc(
        mode=mode, gr=gr, fsx=fsx, ihardlaw=-2, stp=stp, eqc=eqc, jobid=1)
    data = myjob.run()
    print 'VPSC run finished'
    sbar = data['sbar']; ebar = data['eps']
    evm = data['evm']
    x = []; y = []
    if mode=='ssfr':
        for i in range(2):
            for j in range(stp+1):
                #x.append(abs(ebar[i][j][0,1])) #e12
                x.append(abs(evm[i][j]))
                y.append(abs(sbar[i][j][0,1])) #S12
                pass
            pass
        title ='Forward and reverse shear'
        ylabel = '$\sigma_{12}$'; xlabel = '$\epsilon_{VM}$'
        pass
    elif mode=='tc':
        for i in range(2):
            for j in range(stp+1):
                #x.append(abs(ebar[i][j][0,0])) #e11
                x.append(abs(evm[i][j]))
                y.append(abs(sbar[i][j][0,0])) #S11
                pass
            pass
        title = 'Tension and compression'
        ylabel = '$\sigma_{11}$'; xlabel = '$\epsilon_{VM}$'
        pass

    ########################################################
    "plot details"
    fig = plt.figure(ifig)
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.plot(x,y,ms='o')
    ax.set_ylim(0.,)
    ax.set_ylabel(r'%s'%ylabel); ax.set_xlabel(r'%s'%xlabel)
    ########################################################
    
    fname = 'temp000.pdf'
    i = 0
    while True:
        if os.path.isfile(fname):
            i = i + 1
            fname = '%s%s%s'%('temp',str(i).zfill(3),'.pdf')
            pass
        else: break
        pass
    fig.savefig(fname)
    fig.clf()
    print 'the figure has been saved as %s'%fname
    os.system('rm -f *.out')
    os.system('rm -f *.OUT')
    pass


def rauch_test0(ngr=100, eqc=0.005, fsx='disl.sx'):
    """
    For isotropic random texture,
    with changing the prestrain
    observe how the flow curve after reloading
    will be like under simple shear?

    Arguments:
      ngr, eqc, fsx
    """
    from cmb import random as rgr
    import matplotlib.pyplot as plt
    import vp_f2py
    plt.ioff()
    gr = rgr(ngrain=ngr) #100 gr should be sufficient for tests

    ## run the simulations with different prestrain level
    prestr = np.arange(10, 50.01,20) # prestrain steps
    jobs = []
    fig = plt.figure();
    fig.clf()
    ax = fig.add_subplot(111)
    for i in range(len(prestr)):
        stp0 = int(prestr[i])
        stp1 = 60
        print 'stp0: %i stp1: %i ' %(stp0, stp1)
        jobs.append(
            vp_f2py.vpsc(
                mode='ssfr1', stp0=stp0,
                stp1=stp1, stp=0, gr=gr, 
                ihardlaw=-2, fsx=fsx)
            )
        data = jobs[i].run()
        sbar = data['sbar']
        eps = data['eps']
        s = [] ; e = []
        for j in range(len(sbar[0])):
            s.append(sbar[0][j][0][1])
            e.append(j*eqc)
            pass
        for j in range(len(sbar[1])):
            s.append(sbar[1][j][0][1])
            e.append(j*eqc)
            pass
        s = np.abs(np.array(s))
        e = np.array(e)
        ax.plot(e, s, '.',
                label=r'$\varepsilon=%s$'%(
                str(round(prestr[i]*eqc,4))))
        pass
    ax.legend(loc='best')
    title='ss_prestr.pdf'
    fig.savefig(title)
    print 'The plot has been saved to: %s'%title
    pass

class Rauch:
    """
    Class Rauch after the Rauch et al.'s dislocation
    density based hardening

    arguments to __init__
    ----------------------------------------
    mode=None, gr=None,
    stp=None, eqc=0.005, ang=0,
    fsx='fcc.sx', iupdate=[1,1,1,0]
    figname = 'temp.pdf'
    ----------------------------------------
    """
    import matplotlib.pyplot as plt
    import gam
    def __init__(self, mode=None, gr=None,
                 stp=None, eqc=0.005, ang=0,
                 fsx='fcc.sx',
                 iupdate=[1,1,1,0],#ori, shape, hardening, itransph
                 figname='temp.pdf', ngr=10, nsm=10,
                 ): 
        self.mode = mode
        self.gr = gr
        self.eqc = eqc
        self.ang = ang
        self.job = vp_f2py.vpsc(
            mode=mode, gr=self.gr, stp=stp,
            eqc=eqc, ang=ang, iupdate=iupdate)
        self.rst = self.job.run()
        
        self.eps, self.sig, self.ep, self.ec = self.__pp__(
            stp=stp)
        if None in [self.ep, self.ec]:
            self.theta = None
            print 'No orthogonality of the loading',
            print 'path change(theta) is available'
        else:
            self.theta = self.__orthogonal__(
                eps=self.ec, epsp = self.ep)
            print 'Theta: ', self.theta
            pass
        if self.theta==None: figname = '%s_%s.pdf'%(
            figname.split('.pdf')[0],'None')
        else: figname = '%s_%s.pdf'%(
            figname.split('.pdf')[0],
            str(round(self.theta,4)))
        self.__gam__(ngr=ngr, nsm=nsm,
                     relative=False,
                     figname=figname) ## plots randomly selected igr ism
        # self.gam is assigned as gam.gr_plot's return,
        # which is GamDot class
        
        # note below:
        # GamDot.grain(igr,sm=None) returns grain igr's gamdot
        # GamDot.masterdata[igr,ism,istp] is the master data
        self.gam; self.sgr; self.ssm; #self.gam.masterdata[igr,ism,istp]

        nrv = 0; totnrv = 0.
        if self.theta==None: pass
        else:
            # Loading path change for the selected slip modes&grains
            for g in self.sgr:
                for s in self.ssm:
                    psign = np.sign(self.gam.masterdata[g,s,int(
                                self.stp/2-1)])
                    csign = np.sign(self.gam.masterdata[g,s,int(
                                self.stp/2)])
                    if psign!=csign: nrv = nrv + 1
                    pass
                pass
            # Loading path reversed for the all slip modes&grains
            for g in range(len(self.gam.masterdata)):
                for s in range(len(self.gam.masterdata[g])):
                    psign = np.sign(self.gam.masterdata[g,s,int(
                                self.stp/2-1)])
                    csign = np.sign(self.gam.masterdata[g,s,int(
                                self.stp/2)])
                    if psign!=csign: totnrv = totnrv + 1.
                    pass
                pass
            pass
        self.nrv = nrv; self.totnrv = totnrv
        # index depicting the probability of
        # reversing the shear direction: self.avgnrv
        self.avgnrv = self.totnrv / len(self.gr) / len(self.gam.masterdata[0])
        print 'avgnrv: %f'%self.avgnrv
        print 'nrv: %i'%self.nrv
        print 'totnrv: %i'%self.totnrv
        plt.ioff()
        plt.gcf().clf()
        pass
    
    def __pp__(self, stp):
        if self.mode=='RD':
            self.stp = stp
            sbar = self.rst['sbar'][0]
            ebar = self.rst['eps'][0]
            s = [] ; e = []
            for i in range(len(sbar)):
                s.append(sbar[i][0,0])
                e.append(ebar[i][0,0])
                pass
            s = np.array(s); e = np.array(e)
            ep = None
            ec = None 
            pass
        elif self.mode=='ssfr':
            sbar = self.rst['sbar'][0]
            s = []
            e0 = np.arange(0,(stp+1)*self.eqc,self.eqc)
            e1 = np.arange(0,(stp+1)*self.eqc,self.eqc)
            e1 = e1 + e0[-1]
            e = np.array([e0,e1]).flatten()
            for i in range(len(sbar)):
                s.append(sbar[i][0,1])
                pass
            sbar = self.rst['sbar'][1]
            for i in range(len(sbar)):
                s.append(abs(sbar[i][0,1]))
                pass
            s = np.array(s)
            ep = self.rst['dbar'][0][-1]
            ec = self.rst['dbar'][1][0]
            self.stp = (stp+1)*2
            pass
        elif self.mode=='tc':
            sbar = self.rst['sbar'][0]
            s = []
            e0 = np.arange(
                0,(stp+1)*self.eqc,self.eqc)
            e1 = np.arange(
                0,(stp+1)*self.eqc,self.eqc)
            e1 = e1 + e0[-1]
            e = np.array([e0,e1]).flatten()
            for i in range(len(sbar)):
                s.append(sbar[i][0,0])
                pass
            sbar = self.rst['sbar'][1]
            for i in range(len(sbar)):
                s.append(abs(sbar[i][0,0]))
                pass
            s = np.array(s); e = np.array(e)
            ep = self.rst['dbar'][0][-1]
            ec = self.rst['dbar'][1][0]
            self.stp = (stp+1)*2
            pass
        elif self.mode=='ten_ang2':
            sbar = self.rst['sbar'][0]
            s = []
            e0 = np.arange(
                0,(stp+1) * self.eqc,self.eqc)
            e1 = np.arange(
                0,(stp+1) * self.eqc,self.eqc)
            e1 = e1 + e0[-1]
            e = np.array([e0,e1]).flatten()
            for i in range(len(sbar)):
                s.append(sbar[i][0,0])
                pass
            sbar = self.rst['sbar'][2]
            for i in range(len(sbar)):
                s.append(sbar[i][0,0])
                pass
            s = np.array(s); e = np.array(e)
            ep = self.rst['dbar'][0][-1]
            ec = self.rst['dbar'][2][0]

            ## Rotate the strain back ##
            th = self.ang * np.pi/180.
            sth = np.sin(th) ; cth=np.cos(th)
            rmat = np.array(
                [[cth,sth,0],[-sth,cth,0],[0,0,1]])
            rmat = rmat.T #inverse of the rmat
            ec = np.dot(np.dot(rmat,ec),rmat.T)
            ## end of epsilon rotation ##
            self.stp = (stp+1) * 2
            pass
        elif self.mode=='ss_ang':
            sbar = self.rst['sbar'][0]
            s = []
            e0 = np.arange(
                0,(stp+1) * self.eqc,self.eqc)
            e1 = np.arange(
                0,(stp+1) * self.eqc,self.eqc)
            e1 = e1 + e0[-1]
            e = np.array([e0,e1]).flatten()
            for i in range(len(sbar)):
                s.append(sbar[i][0,1])
                pass
            sbar = self.rst['sbar'][2]
            for i in range(len(sbar)):
                s.append(sbar[i][0,1])
                pass
            s = np.array(s); e = np.array(e)
            ep = self.rst['dbar'][0][-1]
            ec = self.rst['dbar'][2][0]

            ## Rotate the strain back ##
            th = self.ang * np.pi/180.
            sth = np.sin(th) ; cth = np.cos(th)
            rmat = np.array(
                [[cth,sth,0],[-sth,cth,0],[0,0,1]])
            rmat = rmat.T
            ec = np.dot(np.dot(rmat, ec), rmat.T)
            ## end of epsilon rotation ##
            self.stp = (stp+1) * 2
            pass
        else:
            print "Given mode <%s> is not available."%self.mode
            print "please add it in the source code"
            raise IOError
        return e, s, ep, ec
    
    def __gam__(self, igr=None, ngr=1, relative=False,
                ifig=1, figname='temp.pdf',
                ism=None, nsm=None):
        import gam
        self.gam, self.sgr, self.ssm = gam.gr_plot(
            igr=igr, ngr=ngr, relative=relative,
            ifig=ifig, figname=figname, ism=ism, nsm=nsm)
        # self.gam.masterdata[igr,ism,istp]
        # gr and sm is the list of selected gr and their slip modes
        pass

    def __orthogonal__(self, eps, epsp):
        """
        Calculates the orthogonality (theta)
        thet = ep:e/sqrt(ep:ep)/sqrt(e:e)
        """
        thet = np.tensordot(epsp, eps)
        thet = thet / np.sqrt(
            np.tensordot(epsp,epsp)) / np.sqrt(
            np.tensordot(eps,eps))
        return thet
    pass


def Rauch_ang(dang=15.,stp=20., eqc=0.005, mode='ss_ang', ngr=100):
    """
    With given texture sample (population of grains),
    if mode=='ss_ang'
      Do shear-rot-shear test with different rotation angles
    if mode=='ten_ang2'
      Do tension-rot-tension with different rotation angles
    """
    import cmb
    import matplotlib.pyplot as plt
    angles = np.arange(0.,180.+ dang/1000., dang)
    gr = cmb.random(ngrain=ngr)
    nrv = []; avgnrv = []; theta = []
    for i in range(len(angles)):
        job = Rauch(mode=mode, gr=gr, fsx='fcc.sx',
                    figname='%s.pdf'%mode, ang=angles[i],
                    stp=stp, eqc=eqc)
        nrv.append(job.nrv)
        avgnrv.append(job.avgnrv)
        theta.append(job.theta)
        pass
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(angles, theta, label=r'$\theta$ factor',color='red')
    axt = ax.twinx()
    #axt.plot(angles, nrv, label=r'$N_{rv}$')
    axt.plot(angles, avgnrv, label=r'$\bar{P}_{rv}$',color='blue')
    ax.set_xlabel('Rotation Angle')
    ax.set_ylabel(r'$\theta$')
    axt.set_ylabel('Probability of shear direction reversal')
    axt.legend(loc='best') ; ax.legend(loc='best')
    axt.set_ylim(0.,1.)
    fig.savefig('%s.pdf'%mode)
    return angles, theta, nrv, avgnrv
    pass

class RauchDen:
    """
    The Ultimate class for debugging and checking the
    implemented Rauch et al's dislocation density based
    hardening.

    2011-Aug-03

    1. Specific loading mode
       a. Monotonic uniaxial ('rd','dd', and 'td')
       b. Monotonic simple shear ('ss')
       c. Forward-Reverse simple shear  ('ssfr')
       c. Simple shear (RD) and simple shear after rotation ('ten_ang2')
       
    2. Commond post-processing
       a. Pbar (Probability to change its shear direction)
       b. Average dislocation density evolution (denf, denr, dent)
       c. Flow stress curve plotting depending on loading modes
       d. Pole figure (111 pole)
       e. If loading path is changed please calculate theta (orthogonality)

    3. Possible parametric tweaks
       a. Influence of different hardening parameters
       b. Different texture (ftex)
       c. Different prestrain (stp0)
       d. Different rotation angle (via ang incase of 'ten_ang2' and 'ss_ang' modes)

    ** Guide to parametric tweaks
    Write an independent defs for parametric tweaks.
    That'll be less confusing and helpful for controlling
    parameters systematically.

    Refer to the below defs.
    """
    def __init__(self, stp0=None, stp1=None, stp=None, ang=0.,
                 eqc=0.005, fsx=None,
                 poles=[[1,1,1]], ngr=100, ftex=None, gr=None,
                 ifig=1,

                 mode='ten_ang2', #'ss_ang','rd','dd','td','tc', 'ssfr'...
                 
                 # single crystal parameters
                 grsize=3.0e-05, burgers=2.46e-10,
                 fdisln=180., shearm=8.5e+04,
                 ftau0=38., portion=0.8, ftherm=2.8, pp = 0.8,
                 ibau=0, idislmode=1,
                 ):
        ##
        import matplotlib.pyplot as plt
        import numpy as np
        from cmb import random as rgr
        import sx_maker
        plt.ioff()
        ##

        #-------------------------------------------------------#
        # ***  figures to plot                                  #
        # pfs (before, after)                                   #
        # flow curve                                            #
        # shear gamma change selected                           #
        # density evolution average, density evolution selected #
        #-------------------------------------------------------#
        
        ## global mode to run
        self.mode=mode
        self.stp = stp
        self.stp0 = stp0
        self.stp1 = stp1
        self.eqc = eqc
        
        ## texture
        if ftex==None and gr==None and ngr!=None: self.gr = rgr(ngrain=ngr)
        elif ftex!=None and gr==None:
            self.gr = np.loadtxt(ftex, skiprows=4)
            pass
        elif ftex==None and gr!=None:
            self.gr = gr
            pass
        else: raise IOError,'Unexpected case on gr and ftex'

        ## single crystal
        if fsx==None:
            fsx ='temp_disl.sx'
            sx_maker.cubic(
                filename=fsx,
                ihardlaw='rauch', dependency='iso', grsize=grsize,
                burgers=burgers, fdisln=fdisln, shearm=shearm,
                ftau0=ftau0, portion=portion, ftherm=ftherm, iopsysx=1,
                idislmode=idislmode, ibau=ibau, pp=pp
                )
            pass
        
        ## run and save data
        self.myjobs = vpsc(mode=self.mode, fsx=fsx, gr=self.gr, ang=ang,
                           stp0=stp0, stp1=stp1, ihardlaw=-2, eqc=self.eqc)
        self.mydata = self.myjobs.run()
        ##----------------------------------------------------

        ## orthogonality change
        if mode=='ss_ang' or mode=='ten_ang2':
            ep = self.mydata['dbar'][0][-1]
            ec = self.mydata['dbar'][2][0]
            ## rotation of the ec ##
            th = ang * np.pi/180.
            sth = np.sin(th) ; cth=np.cos(th)
            rmat = np.array(
                [[cth,sth,0],[-sth,cth,0],[0,0,1]])
            rmat = rmat.T #inverse of the rmat
            ec = np.dot(np.dot(rmat,ec),rmat.T)
            ## end of rotation --- 
            pass
        elif mode=='tc' or mode=='ssfr':
            ep = self.mydata['dbar'][0][-1]
            ec = self.mydata['dbar'][1][0]
            pass
        elif mode=='rd' or mode=='td' or mode=='dd' or mode=='ss':
            ep = None; ec = None
            pass
        else: raise IOError,'Unexpected mode'

        if None in [ep, ec]: self.thet=None
        else:
            thet = np.tensordot(ep, ec)
            thet = thet / np.sqrt(
                np.tensordot(ep,ep))
            thet = thet / np.sqrt(
                np.tensordot(ec,ec))
            self.thet = thet
            pass
        ## End of orthogonality cal. -----------
        self.__pp__(
            figname='%s'%mode,
            ifig = ifig,
            poles = poles
            )

        ## erase of *.out ##
        # fremove = glob.glob('*.out')
        # for f in fremove: os.remove(f)
        ##
        
        pass
    
    def __pp__(self, figname, ifig, poles):
        """
        Post-processing and plotting
        Deal with denf_, denr_, dent_, gam_, and crss_ files

        arguments:
           figname, ifig, poles
        """
        import upf
        import den
        import gam
        # Average density change
        plt.figure(ifig).clf()
        self.df, self.dr, self.dt = den.gr_plotavg(
            figname='denavg_%s.pdf'%figname, ifig=ifig
            )
        # gammadot plot
        plt.figure(ifig).clf()
        self.gam, self.sgr, self.ssm = gam.gr_plot(
            ngr=10, ifig=ifig,
            figname='selected_gamdot_%s.pdf'%figname, nsm=2)
        
        if self.thet==None:
            print 'theta is not available'
            raw_input()
            self.nrv=None; self.totnrv=None
            self.avgnrv=None
            pass
        else:
            nrv = 0; totnrv = 0.
            # Loading path change for the selected slip modes&grains
            for g in self.sgr:
                for s in self.ssm:
                    psign = np.sign(self.gam.masterdata[g,s,int(
                                self.stp0)])
                    csign = np.sign(self.gam.masterdata[g,s,int(
                                self.stp0+1)])
                    if psign!=csign: nrv = nrv + 1
                    pass
                pass
            # Loading path reversed for the all slip modes&grains
            for g in range(len(self.gam.masterdata)):
                for s in range(len(self.gam.masterdata[g])):
                    psign = np.sign(self.gam.masterdata[g,s,int(
                                self.stp0)])
                    csign = np.sign(self.gam.masterdata[g,s,int(
                                self.stp0+1)])
                    if psign!=csign: totnrv = totnrv + 1.
                    pass
                pass
            self.nrv = nrv ; self.totnrv = totnrv
            self.avgnrv = self.totnrv / len(self.gr)/ len(self.gam.masterdata[0])
            pass
        
        # flow curves
        if self.mode=='rd' or self.mode=='td' or self.mode=='dd':
            sbar = self.mydata['sbar'][0]
            eps = self.mydata['eps'][0]
            s = []
            e = []
            for i in range(len(sbar)):
                s.append(sbar[i][0,0])
                e.append(eps[i][0,0])
                pass
            self.s = np.array(s); self.e = np.array(e)
            pass
        elif self.mode=='ten_ang2' or  self.mode=='ss_ang':
            sbar1 = self.mydata['sbar'][0]
            sbar2 = self.mydata['sbar'][2]
            s = []
            e0 = np.arange(0, (self.stp0+1) * self.eqc, self.eqc)
            e1 = np.arange(0, (self.stp1+1) * self.eqc, self.eqc)
            ee = []
            for i in range(len(e0)): ee.append(e0[i])
            for i in range(len(e1)): ee.append(e1[i])
            e1 = e1.copy() + e0[-1]
            e = []
            for i in range(len(e0)): e.append(e0[i])
            for i in range(len(e1)): e.append(e1[i])
            e = np.array(e)
            if self.mode=='ten_ang2': ix, iy = 0, 0
            elif self.mode=='ss_ang': ix, iy = 0, 1
            else: raise IOError, 'Unexpected self.mode'
            for i in range(len(sbar1)):
                s.append(sbar1[i][ix,iy])
                pass
            for i in range(len(sbar2)):
                s.append(sbar2[i][ix,iy])
                pass
            self.e = np.array(e); self.s = np.array(s)
            self.ee = np.array(ee)
            pass
        elif self.mode=='tc' or self.mode=='ssfr':
            sbar1 = self.mydata['sbar'][0]
            sbar2 = self.mydata['sbar'][1]
            s = []            
            e0 = np.arange(0, (self.stp0+1) * self.eqc, self.eqc)
            e1 = np.arange(0, (self.stp1+1) * self.eqc, self.eqc)
            ee = []
            for i in range(len(e0)): ee.append(e0[i])
            for i in range(len(e1)): ee.append(e1[i])            
            e1 = e1.copy() + e0[-1]
            e = []
            for i in range(len(e0)): e.append(e0[i])
            for i in range(len(e1)): e.append(e1[i])
            e = np.array(e)
            if self.mode=='tc': ix, iy = 0, 0
            elif self.mode=='ssfr': ix, iy = 0, 1
            else: raise IOError, 'Unexpected self.mode'
            for i in range(len(sbar1)):
                s.append(sbar1[i][ix,iy])
                pass
            for i in range(len(sbar2)):
                s.append(abs(sbar2[i][ix,iy]))
                pass
            self.e = np.array(e); self.s = np.array(s)
            self.ee = np.array(ee)
            pass
        else: raise IOError, 'Unexpected self.mode'
        # Texture
        self.lgr = self.mydata['lgr']
        #self.gr
        pass
    pass # end of class

def rauchloadings(stp0, stp1, stp, ang,
                  eqc=0.005, fsx='disl.sx'):
    """
    Conducts different loading path
    """
    modes = ['rd','dd','td','ss','ssfr','tc','ten_angs','ss_ang']
    pass


def __rotate(ang=0., fn='temp.rot'):
    th = ang * np.pi / 180.
    sth = np.sin(th) ; cth = np.cos(th)
    rmat = np.array([[cth, sth, 0.],[-sth, cth, 0.],[0., 0., 1.]])
            
    FILE = open(fn, 'w')
    FILE.writelines('Rotation matrix for polycrystalline aggregate\n')
    for i in range(3):
        for j in range(3):
            FILE.writelines('%13.9f  '%rmat[i][j])
            pass
        FILE.writelines('\n')
        pass
    FILE.close()


def rauchbench(fsx=None, ftex=None, gr=None,

               # hardening control.
               idislmode=1,
               ibau=0,
               pp=0.0,

               # single crystal parameters
               grsize=3.0e-05, burgers=2.46e-10,
               fdisln=180., shearm=8.5e+04,
               ftau0=38., portion=0.8, ftherm=2.8, eqincr=0.005):
    """
    Benchmark problem done by Dr. Tome and Koshiro
    process: tension 10% RD, tension 5% 45 degree, and tension 5% RD
    """
    if ftex==None: ftex='IRON_5A.TEX'
    from sx_maker import cubic as cubsx
    histloader = vpsc(mode='RD', gr=[[0,0,0,0]],
                      stp=1, eqc=1,
                      jobid=100).__loading__
    rot = __rotate
    ##
    nstep =int( 0.10/ eqincr)
    histloader(histfile='rd0', mode='unitension',
               eqincr=eqincr, nstep=nstep)
    rot(ang=45, fn='rbr45')
    nstep = int(0.05 / eqincr)
    histloader(histfile='rd1', mode='unitension',
               eqincr=eqincr, nstep=nstep)
    rot(ang=-45, fn='rbr-45')
    nstep =int( 0.05 / eqincr)
    histloader(histfile='rd2', mode='unitension',
               eqincr=eqincr, nstep=nstep)
    del histloader
    del rot
    
    ##
    if fsx==None:
        fsx = 'rauchbench.sx'
        cubsx(
            filename=fsx,
            ihardlaw='rauch', dependency='iso', grsize=grsize,
            burgers=burgers, fdisln=fdisln, shearm=shearm,
            ftau0=ftau0, portion=portion, ftherm=ftherm, iopsysx=1,
            idislmode=idislmode, ibau=ibau, pp=pp
            )
        pass
    ##
    prcs = ['0', 'rd0', '4', 'rbr45', '0',
            'rd1', '4', 'rbr-45', '0', 'rd2']
    job = vpsc(prcs=prcs, texture=ftex, gr=gr,
               fsx=fsx, ihardlaw=-2, jobid=2)
    job.run()

    ## deletion of gr files. --------
    files = glob.glob('gr_*_*')     #
    for f in files: os.remove(f)    #
    ## ------------------------------

    ## pyplot setting ##
    plt.ioff()
    fig0 = plt.figure(1); ax0 = fig0.add_subplot(111) #str-str
    fig1 = plt.figure(2); ax1 = fig1.add_subplot(111) #
    ##

    evm = np.append(job.datamaster['evm'][0],job.datamaster['evm'][1])
    evm = np.append(evm, job.datamaster['evm'][2])
    svm = np.append(job.datamaster['svm'][0],job.datamaster['svm'][1])
    svm = np.append(svm, job.datamaster['svm'][2])

    ax0.plot(evm,svm,'o',mfc='None'); ax0.set_ylim(0.,)
    ax0.set_ylabel(r'$\sigma^{VM}$ [MPa]', fontsize=20)
    ax0.set_xlabel(r'$\varepsilon^{VM}$', fontsize=20)
    
    fig0.savefig('str-str_rauch_bench.pdf')
    fig0.clf();
    
    return job

def rauchprestrain(nstp=3, inistp=10, finstp=30, stp=60,
                   mode='ssfr', eqc=0.005,
                   fsx=None, ftex=None, ngr=None, gr=None,
                   ang=0., ifig=1,
                   
                   # single crystal parameters
                   idislmode=1, #0:active-nonactive, 1:all active
                   ibau=0, #0 for no bauschinger 1: pp protion deletion
                   pp=0.8,

                   # dislocation and microstructure parameters
                   grsize=3.0e-05, burgers=2.46e-10,
                   fdisln=180., shearm=8.5e+04,
                   ftau0=38., portion=0.8, ftherm=2.8
                   ):
    """
    Different amount of prestrains

    Arguments:
     nstp=3 : number of prestrain levels
     inistp=10 : initial prestrain steps
     finstp=30 : final prestrain steps
     stp=60: number of strain steps after reloading
     mode='ssfr': vpsc mode
     eqc=0.005
     fsx='disl.sx'
     ftex=None
     ngr=100
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from cmb import random as rgr
    plt.ioff()
    
    if ftex==None and ngr!=None and gr==None:
        gr = rgr(ngrain=ngr)
        ngr = None
        pass
    elif ftex!=None and ngr==None and gr==None:
        gr = np.loadtxt(ftex, skiprows=4)
        pass
    elif ftex==None and ngr==None and gr!=None:
        pass
    else: raise IOError, 'Unexpected case'

    jobs = []
    theta = []
    s = []; e = []; ee = [];
    df, dr, dt = [], [], []
    lgr = []; avgnrv = []
    
    steps = np.array(np.linspace(inistp, finstp, nstp),
                     dtype='int')
    
    for i in range(len(steps)):
        jobs.append(
            RauchDen(
                stp0=steps[i], stp1=stp, mode=mode, ang=ang,
                eqc=eqc, fsx=fsx, ftex=ftex, gr=gr, ibau=ibau,
                idislmode=idislmode, pp=pp
                )
            )
        theta.append(jobs[i].thet)
        s.append(jobs[i].s)
        e.append(jobs[i].e)
        ee.append(jobs[i].ee)
        df.append(jobs[i].df)
        dr.append(jobs[i].dr)
        dt.append(jobs[i].dt)
        lgr.append(jobs[i].lgr)
        avgnrv.append(jobs[i].avgnrv)
        pass

    #----------------------------------------------------------#
    #   ** plotting                                            #
    # Visualization in this case is very important             #
    # It is as much important as the idea of the model itself. #
    # Therefore, plot the data clear enough and always try to  #
    # find the best way of plotting and impove it furthermore. #
    #----------------------------------------------------------#

    #--------------------------------------------------------------#
    # * Two ways of flow stress curves                             #
    #   - Regular type in which the strain is accummulated one     #
    #   - dragging the last segment to the origin to facilitate    #
    #     the comparison between the flow curve before path change #
    #     and after path change.                                   #    
    # * Two ways of dislocation density tracker                    #
    #   - Log scale                                                #
    #   - Normal scale                                             #    
    # * Pole figure                                                #
    #   - initial-final pole figures                               #
    # * All in one type                                            #
    #   - Flow curve, density, theta are in the same graph         #
    #--------------------------------------------------------------#

    ## - flow curves
    fig0 = plt.figure(ifig); fig0.clf()     ## flow curves
    fig1 = plt.figure(ifig+1,figsize=(16,12));
    fig1.clf() ## fig1: average dislocation density tracker
    fig2 = plt.figure(ifig+2); fig2.clf()   ## theta, avgnrv
    ax00 = fig0.add_subplot(211) #Regular flow curve
    ax01 = fig0.add_subplot(212) #Rauch's flow curve style
    ax02 = fig1.add_subplot(221) #Log sacle average density tracker
    ax03 = fig1.add_subplot(223) #Normal scale average density tracker
    ax04 = fig2.add_subplot(111) #Theta and average 

    markers=['o','.','+','x','d','*','o','.','+','x','d','*']
    colors=['k','r','b','g','gray','m','k','r','b','g','gray','m']
    for i in range(len(steps)):
        ax00.plot(
            e[i], s[i],
            label=r'$\varepsilon_p = %s$'%str(
                round(steps[i] * eqc, 3))
            )
        ax01.plot(
            ee[i], s[i], markers[i],
            label=r'$\varepsilon_p = %s$'%str(
                round(steps[i] * eqc, 3))
            )
        ax02.semilogy(
            e[i], df[i], 'o', color= colors[i],
            label=r'$\rho_f$ at $\varepsilon_p = %s$'%str(
                round(steps[i] * eqc, 3))
            )
        ax02.semilogy(
            e[i], dr[i], 'x', color=colors[i],
            label=r'$\rho_r$ at $\varepsilon_p = %s$'%str(
                round(steps[i] * eqc, 3))
            )
        # ax02.semilogy(
        #     e[i], dt[i], color=colors[i],
        #     label=r'$\tilde{\rho}$ at $\varepsilon_p = %s$'%str(
        #         round(steps[i] * eqc, 3))
        #     )        
        ax03.plot(
            e[i], df[i], 'o', color= colors[i],
            label=r'$\rho_f$ at $\varepsilon_p = %s$'%str(
                round(steps[i] * eqc, 3))
            )
        ax03.plot(
            e[i], dr[i], 'x', color= colors[i],
            label=r'$\rho_r$ at $\varepsilon_p = %s$'%str(
                round(steps[i] * eqc, 3))
            )
        # ax03.plot(
        #     e[i], dt[i], color=colors[i],
        #     label=r'$\tilde{\rho}$ at $\varepsilon_p = %s$'%str(
        #         round(steps[i] * eqc, 3))
        #     )        
        pass
    
    ax00.legend(loc='best'); ax01.legend(loc='best')
    ax02.legend(loc=2, bbox_to_anchor=(1.05,0.8))
    ax00.grid('--'); ax01.grid('--'); ax02.grid('--'); ax03.grid('--')    
    ax00.set_ylim(0.,); ax01.set_ylim(0.,)
    ax00.set_ylabel(r'$\sigma$'); ax00.set_xlabel(r'$\epsilon$')
    ax01.set_ylabel(r'$\sigma$'); ax01.set_xlabel(r'$\epsilon$')
    #ax03.legend(loc=2, bbox_to_anchor=(1.05,1))
    
    ## save figures
    fig0.savefig('fc.pdf');fig1.savefig('den.pdf')
    fig0.clf();fig1.clf()
    ##
    pass


def rauchportion():
    """
    Different portion level
    """
    pass
    
def rauchangles():
    """
    Conducts series of different angle rotation
    """
    pass

def rauchtex():
    """
    Conducts a series of different textures
    """
    pass





"""
 example 10 ! MSMSE2011 paper to be submitted revision
    Systematic assess of grain number nad linearization 
   sensitivity on the prediction capability.

 ------------------------------------------------------------
  The response to the MSMSE2011 paper reviews
  Which linearization scheme upon a certain number of grains works
  the best?

  outline
 
  loop - number of grains
    loop - making RVE (several RVEs)
       loop - over linearization (Taylor, Secant, Tangent, Affine, Neff=10)
          1. find the parameters
          2. do the tensile tests
          3. analyze how much they are well-fitted!(standard deviation)
              --> local statistical analysis
       loop
    end
        --> macro statistical analysis
  end

  suggest the best case! (The best combincation of number
  of grains and  linearization scheme)
 -------------------------------------------------------------
"""

class BestFitting:
    """
    Find the best combination of number of grains and
    linearization scheme

    Assess the capability of the given self-consistent scheme
    based on experimental stress-strain curve.

    plan:
       Expand the experimental reference on which the capability
      is estimated to R-value curve.

    Modified strategy: (it's rather a recipe than a strategy though)
      1. Find suitable parameters for 10 gr with 20 iterations
         then start with the optimized parameter set
      2. Then, try with 20 gr with 20 iterations
      3. then, 100 gr with 20 iterations
      4. then, 500 gr with 2 iterations

    Do this with the globally designated linearization methods

    --> 1. Pre testing
        2. Save the parameters (to a file as well for latter use)
           for each of linearization
        3. do the simulations as globablly given
           by making full use of obtained (optimized parameters)
    """
    def __init__(
        self, cod='304_corrected(surface).TXT',
        
        #304 STS ref str-eps curve along RD
        sscurve='eps-sig-rd-meancurve',
        #304 STS ref str-eps curve along TD
        ss_cf = 'eps-sig-td-meancurve',
        #304 STS cf r-eps curve along RD,
        rcurve='eps-R-rd-meancurve',
        #304 STS cf r-eps curve along TD,
        r_cf = 'eps-R-td-meancurve',

        ## number of grains
        ngr=[10, 50], #, 1000, 2000, 4000, 8000],
        
        # Linearization : [0,1,2,3,4] ------------------------
        # 0:FC, 1:Affine, 2:secant, 3:neff=10, 4:tangent, 5:SO
        # if 0:FC, comparision using R-value is not possible
        linearization = [0,1], #,2,3,4],
        # ----------------------------------------------------

        # number of sampling of RVE
        nsample = 1, #5
        
        # initial guess for Voce hardening
        tau0=98, tau1=14.82, thet0=448, thet1=243.,
        #tau0=88., tau1=20., thet0=200., thet1=100.,

        ## xtolerance
        #xtol=0.001, maxiter=10, (replaced by recipe)

        ## pretesting recipe
        recipe = 5, # 0: is the most time-consuming.
        ifig=30,
        ):
        
        import glob
        print "Best Fitting Class will erase all *.pdf, *.out, *each_stat* files"
        temp = '*.pdf'; r = glob.glob(temp)
        for fn in r: os.remove(fn)
        temp = '*.eps'; r = glob.glob(temp)
        for fn in r: os.remove(fn)        
        temp = '*.out'; r = glob.glob(temp)
        for fn in r: os.remove(fn)
        temp = '*.OUT'; r = glob.glob(temp)
        for fn in r: os.remove(fn)
        temp = 'VPSC*_*.in*'; r = glob.glob(temp)
        for fn in r: os.remove(fn)        
        temp = '*each_stat*'; r = glob.glob(temp)
        for fn in r: os.remove(fn)
        
        import matplotlib.pyplot as plt
        import sys
        plt.ioff()
        import cmb, time    
        
        if cod==None: raise IOError,'No cod file'
        if not(os.path.isfile(cod)):
            raise IOError,'no such cod file found'
        
        ## presettings --------------------------------------------
        self.cod = cod

        ## experimental data settings -----------------------------
        # exp stress strain curve
        self.eps_ref,  self.sig_ref,dum = np.loadtxt(sscurve).T
        self.ssexp_ref = np.array([self.eps_ref, self.sig_ref])   #flow curve
        self.eps_cf,   self.sig_cf, dum = np.loadtxt(ss_cf).T 
        self.ssexp_cf  = np.array([self.eps_cf, self.sig_cf])     #flow curve       
        # exp r-value strain curve
        self.reps_ref, self.rval_ref,  dum = np.loadtxt(rcurve).T
        self.rexp_ref  = np.array([self.reps_ref, self.rval_ref]) #r-values
        self.reps_cf,  self.rval_cf,   dum = np.loadtxt(r_cf).T
        self.rexp_cf   = np.array([self.reps_cf, self.rval_cf])   #r-values
        ## experimental data settings -----------------------------

        self.linearization = linearization
        self.ngr = ngr
        self.nsample = nsample
        ## --------------------------------------------------------
        
        ## log file -----------------------------------------------
        i = 0
        while True:
            logfilename = 'bestfitting_log%s'%str(i).zfill(3)
            if not(os.path.isfile): i = i + 1
            else:
                log = open(logfilename, 'w') # log file naming
                break
            if i>100:
                print 'too many iterations in setting log file nam'
                raise IOError
            pass
        ## --------------------------------------------------------

        ### master data file -----------------------------------------
        masterfilename = 'bestfitting_master%s'%str(i).zfill(3)
        masterd = open(masterfilename, 'w')
        masterd.writelines(' *** BESTFITTING MASTER DATA FILE *** \n')
        masterd.writelines(' %s\n'%time.ctime())
        masterd.close(); masterd = open(masterfilename,'a')
        ### ----------------------------------------------------------
        """
        Basic strategy
         1. Standard deviation of each of linearization method upon each grain
            1-1. Average the standard deviation of each RVE sampling.
         2. Based on this, suggest the best linearization method for
            each grain number
        """
        ## pretesting
        time_0 = time.time()
        st = 'Parameter presetting starts'
        print st
        log.write(st+'\n');log.close(); log = open(logfilename, 'a')

        ## file descriptor is redirected within self.__pretesting__
        print '__pretesting__ starts'
        self.voce, totgr = self.__pretesting__( 
            #linearization=self.linearization, default: self.linearization
            tau0=tau0, tau1=tau1, thet0=thet0, thet1=thet1,
            # use the default recipe
            ssexp_ref = self.ssexp_ref,
            recipe=recipe
            )
        
        # totgr[len(linearization:prescribed in __init__)+ len(testngr: given from the __recipe__)]
        
        # Note that self.voce is in the shape of recipe.##
        
        tt, u = self.__timecnvt__(time = (time.time()-time_0))
        print '__pretesting__ ends  elapsed time:%4.2f [%s]'%(
            tt, u)
        
        testngr = self.recipe[0]
        testmaxiter = self.recipe[1]
        print 'The recipe that has used: '
        print '# grains:', testngr,'\nmaxiter:', testmaxiter
        ## --------------------------------------------------------

        import matplotlib.pyplot as plt
        import matplotlib.font_manager
        import numpy.random as ra

        rand = ra.rand
        plt.ioff()
        ## figures setting ---------------------------------------- ##
        #fac = 30.
        # width, height = len(self.linearization), self.nsample
        # norm = np.sqrt(width**2.+height**2.)
        # width, height = width/norm * fac, height/norm * fac

        h, w = plt.figaspect(
            rand(len(self.linearization),
                 self.nsample))
        
        figs = []; gax = [] #total grid axes for stress
        figsr = []; gsxr =[] #total grid axes for r-value
        print 'figure declaration'

        for i in range(len(testngr)):
            ## ss curve
            figs.append(
                plt.figure(
                    ifig+i*2,
                    figsize=(w, h)
                    )
                )
            figs[i].clf()
            figs[i].subplots_adjust(
                bottom=0.3, top=0.7, wspace=0.5, hspace=0.3
                )
            ## r curve
            figsr.append(
                plt.figure(
                    ifig+i*2+1,
                    figsize=(w, h)
                    )
                )
            figsr[i].clf()
            figsr[i].subplots_adjust(
                bottom=0.3, top=0.7, wspace=0.5, hspace=0.3
                )
            pass

        for j in range(len(testngr)):
            for i in range(len(self.linearization)*nsample):
                ## ss curve
                figs[j].add_subplot(
                    nsample, len(self.linearization), i+1 #6
                    )
                ## r curve
                figsr[j].add_subplot(
                    nsample, len(self.linearization), i+1 #6
                    )
                pass
            pass
        print '\n%i figures are to be plotted'%(len(self.ngr))
        
        ## figures setting ends ----------------------------------- ##

        ## stress deviation
        sdv_ref_tot = np.zeros(
            (len(self.linearization), len(testngr)) #,len(nsample)
             )
        sdv_cf_tot = sdv_ref_tot.copy()

        ## r-value deviation
        rdv_ref_tot = np.zeros(
            (len(self.linearization), len(testngr)) #,len(nsample)
             )
        rdv_cf_tot = rdv_ref_tot.copy()
        


        print '\n *** Check the fitted curves ***\n'

        

        ithgtot = 0
        for i in range(len(self.linearization)):
            print '%i - testngr : %i grains'%(j, testngr[j])
            for j in range(len(testngr)): #deletion of the loop -> revival
                print '%i - linearization'%i
                temp = 0
                for k in range(nsample):
                    ###------------------------------------------###
                    ## test ngrain is set to be 4000
                    #tngr = 4000

                    print '%i - sample'%k
                    #self.__nullstdout__(opt=0) #start nulling stdout
                    
                    if nsample==1:
                        """
                        if nsample equals 1, the previously used rve is
                        again used...
                        """
                        
                        self.rve2ftex(filename='check.cmb', rve=totgr[ithgtot])
                        pass
                    elif nsample>1:
                        """
                        if nsample is larger than one, RVEs are re-made
                        """
                        a = cmb.RVE(
                            ngrain=testngr[j], odf=self.cod,
                            cmbfile='check.cmb'
                            )
                        pass

                    tau0, tau1, thet0, thet1 = self.voce[i][j] # use max ngr result
                    gax  = figs[ j].axes[i+ (k)*len(self.linearization)]
                    rgax = figsr[j].axes[i+ (k)*len(self.linearization)]

                    ## caution! if self.linearization[i]==0: Taylor
                    sdv_ref, rdv_ref = self.each_stat(
                        param = self.voce[i][j], #self.voce[i,j] (i:lin, y:ngr)
                        exp = self.ssexp_ref,
                        rexp = self.rexp_ref,
                        ngr = testngr[j], #this is only for label purpose 
                        ftex='check.cmb',
                        nsample=1, #
                        loadingmode='RD',
                        linearization=self.linearization[i],
                        ifig=1,
                        figname ='check_prefit_lin%s'%str(self.linearization[i]),
                        sgax=gax,
                        rgax=rgax,
                        )
                    #print 'sdv_ref, rdv_ref'
                    #print sdv_ref, rdv_ref
                    #raw_input()
                    
                    sdv_cf, rdv_cf = self.each_stat(
                        param = self.voce[i][j], #self.voce[i,j] (i:lin, y:ngr)
                        exp = self.ssexp_cf,
                        rexp = self.rexp_cf,
                        ngr = testngr[j], #this is only for label purpose 
                        ftex='check.cmb',
                        nsample=1, #
                        loadingmode='TD',
                        linearization=self.linearization[i],
                        ifig=1,
                        figname ='check_prefit_lin%s'%str(self.linearization[i]),
                        sgax=gax,
                        rgax=rgax,
                        )
                    
                    gax.set_xlabel(r'$\varepsilon$')
                    gax.set_ylabel(r'$\sigma$[MPa]')
                    rgax.set_xlabel(r'$\varepsilon$')
                    rgax.set_ylabel('R-value')

                    
                    #self.__nullstdout__(opt=1) #stop nulling stdout
                    temp = temp + sdv_ref
                    ###------------------------------------------###
                    pass
                ithgtot = ithgtot + 1
                #average if necessary (effective only if nsample!=1)
                sdv_ref_tot[i][j] = sdv_ref / nsample
                sdv_cf_tot[i][j] = sdv_cf / nsample
                rdv_ref_tot[i][j] = rdv_ref / nsample
                rdv_cf_tot[i][j] = rdv_cf / nsample
                pass
            pass

        for i in range(len(testngr)):
            figs[i].savefig('grid_ss_%s.pdf'%str(testngr[i]).zfill(5))
            figs[i].savefig('grid_ss_%s.eps'%str(testngr[i]).zfill(5))
            pass
        for i in range(len(testngr)):
            figsr[i].savefig('grid_re_%s.pdf'%str(testngr[i]).zfill(5))
            figsr[i].savefig('grid_re_%s.eps'%str(testngr[i]).zfill(5))
            pass        
        
        #np.savetxt('voce', self.voce)
        log.write(
            "voce parameters from pretesting has been saved to 'voce' file\n")
        ## ------------------------------------------------------
        elapsed = time.time() - time_0
        elapsed, u = self.__timecnvt__(time=elapsed)
        st0 = 'parameter presetting ends'
        st1 = 'elasped time: %4.2f [%s]'%(elapsed, u)
        st2 = 'check out the results'
        st = st0 + '\n' + st1 + '\n' + st2
        
        print st
        
        log.write(st)
        log.write('preset parameters\n')
        log.write(
            '%3s %8s %13s %13s %13s %13s %13s %13s\n'%(
                'lin','ngr','tau0','tau1','thet0','thet1', 'sdv_ref', 'sdv_cf')
            )
        
        for i in range(len(self.linearization)):
            for j in range(len(testngr)): 
                log.write(
                    '%3i %8i %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n'%(
                        self.linearization[i],
                        testngr[j], #tngr, #ngr[j],
                        self.voce[i][j][0], # use result from the maximum (the last) 
                        self.voce[i][j][1], # pretesting voce parammeters
                        self.voce[i][j][2],
                        self.voce[i][j][3],
                        sdv_ref_tot[i][j],
                        sdv_cf_tot[i][j]
                        ))
                pass
            pass
        log.close() #log = open(logfilename, 'a')
        masterd.close() #masterd = open(masterfilename, 'a')
        self.sdv_ref_tot = sdv_ref_tot
        self.sdv_cf_tot = sdv_cf_tot

        
        
        #--- sdv figure template #1  ---#
        fig = plt.figure()
        fig.clf()
        sdvr = sdv_ref_tot.T
        sdvc = sdv_cf_tot.T
        for i in range(len(testngr)):
            plt.plot(sdvr[i], 'x', label='RD %i'%testngr[i])
            plt.plot(sdvc[i], 'o', label='TD %i'%testngr[i])
            pass
        sint = []
        for j in range(len(self.linearization)):
            sint.append(__lins__(lin=self.linearization[j]))
            pass

        plt.gca().set_xlim(-1, len(linearization))
        plt.xticks(np.arange(len(self.linearization)), sint)
        prop = matplotlib.font_manager.FontProperties(
            size=8)        
        plt.legend(loc='best')
        plt.gcf().savefig('sdvt.pdf')
        plt.gcf().savefig('sdvt.eps')
        
        #--- sdv figure template #2  ---#
        fig = plt.figure()
        fig.clf()
        colors=['r','g','b','k','y']
        for i in range(len(self.linearization)):
            plt.plot(testngr, sdv_ref_tot[i], '--',color=colors[i],
                     label='RD %s'%__lins__(lin=self.linearization[i])
                     )
            plt.plot(testngr,sdv_cf_tot[i], 'o',mec=colors[i], mfc='none',
                     label='TD %s'%__lins__(lin=self.linearization[i])
                     )
            pass
        xl, xh = plt.gca().get_xlim()
        df = xh - xl 
        plt.gca().set_xlim(xl - df/len(testngr), xh + df/len(testngr))
        plt.xticks(testngr, map(str, testngr))
        plt.gca().set_xlabel('Number of grains in RVE')
        plt.gca().set_ylabel(
            #r'$\frac{\sqrt{\sum^n_i{(S^{simul}_i-S^{exp}_i)^2}}}{n-1}$',
            'deviation',
            dict(fontsize=15),
            )
        

        plt.legend(loc='best', prop=prop )
        plt.gcf().savefig('sdv.pdf')
        plt.gcf().savefig('sdv.eps')




        #--- rdv figure based on the above template #2  ---#
        fig = plt.figure()
        fig.clf()
        colors=['r','g','b','k','y']
        for i in range(len(self.linearization)):
            if self.linearization[i]==0: pass #skip Taylor
            else:
                plt.plot(testngr, rdv_ref_tot[i], '--',color=colors[i],
                         label='RD %s'%__lins__(lin=self.linearization[i])
                         )
                plt.plot(testngr, rdv_cf_tot[i], 'o',mec=colors[i], mfc='none',
                         label='TD %s'%__lins__(lin=self.linearization[i])
                         )
                pass
            pass
        
        xl, xh = plt.gca().get_xlim()
        df = xh - xl 
        plt.gca().set_xlim(xl - df/len(testngr), xh + df/len(testngr))
        plt.xticks(testngr, map(str, testngr))
        plt.gca().set_xlabel('Number of grains in RVE')
        plt.gca().set_ylabel(
            'devation',
            #r'$\frac{\sqrt{\sum^n_i{(R^{simul}_i-R^{exp}_i)^2}}}{n-1}$',
            dict(fontsize=15)
            )
        

        plt.legend(loc='best', prop=prop )
        plt.gcf().savefig('rdv.pdf')
        plt.gcf().savefig('rdv.eps')            
        pass

        
    
    def test1(self, nsample=1,
              maxiter=10,
              xtol = 10.,
              fn = 'test1',
              ifig=10,
              ifigdum = 112,
              ):
        """
        use self.voce -> for each of linearization

        resulting standard deviation is saved in the file
        which named as fn

        --
        nsample : Number of repeatition upon the same condition
        upon different RVE sampling.
        maxiter: maxiteration for fine tuning upon each sampling again.
        xtol=10. tolerance
        fn='test1', testfile save name
        """
        ### master data file -----------------------------------------
        filename = fn
        masterd = open(filename, 'w')
        masterd.writelines(' *** BESTFITTING MASTER DATA FILE *** \n')
        ### ----------------------------------------------------------
        
        import time, cmb
        import matplotlib.pyplot as plt
        #import mpl_toolkits.axes_grid as ag #axes_grid
        plt.ioff()

        ## figures setting
        fac = 7
        width, height = fac * len(self.linearization), fac
        figs = []; gax_tot = [] #total grid axes
        for i in range(len(self.ngr)):
            figs.append(plt.figure(ifig+i,figsize=(width, height))) #default=10
            figs[i].clf()
            figs[i].subplots_adjust(
                bottom=0.2, top=0.8, wspace=0.3, hspace=0.2
                )
            pass
        
        for i in range(len(self.ngr)):        
            for j in range(len(self.linearization)*nsample):
                figs[i].add_subplot(
                    nsample, len(self.linearization), j+1 #6
                    )
                pass
            pass
        print '%s figs are to be plotted'%str(len(self.ngr)).zfill(2)
        ## figure setting over

        ### 
        t_0 = time.time()
        dum_ref = np.zeros(
            (len(self.ngr),
             len(self.linearization),
             nsample))
        dum_cf = dum_ref.copy()
        dum_voce = np.zeros(
            (len(self.ngr),
             len(self.linearization),
             nsample,
             4))
        ###
        
        ## Main loop -----------------------------------------------
        for i in range(len(self.ngr)):
            ithax=0
            for j in range(len(self.linearization)):
                for k in range(nsample):
                    ###------------------------------------------###
                    self.__nullstdout__(opt=0) #start nulling stdout
                    ###------------------------------------------###
                    ftex = 'temp.cmb'
                    a = cmb.RVE(
                        ngrain = self.ngr[i],
                        odf= self.cod,
                        cmbfile=ftex)
                    ###------------------------------------------###
                    self.__nullstdout__(opt=1) #start nulling stdout
                    ###------------------------------------------###

                    ifig0 = plt.figure(ifigdum)
                    ifig1 = plt.figure(ifigdum + 1)
                    ifig0.clf()
                    ifig1.clf()

                    #only use the last pretested voce parameter set one
                    # in belief that the more number of grains, the
                    # better the result will be.
                    tau0, tau1, thet0, thet1 = self.voce[j,-1]
                    
                    ###------------------------------------------###
                    self.__nullstdout__(opt=0) #start nulling stdout
                    ###------------------------------------------###
                    # fine tuning 
                    voce_param = self.param(
                        eqc=0.005,
                        
                        #voce parameters
                        tau0=tau0, tau1=tau1,
                        thet0=thet0, thet1=thet1,
                        ##
                        
                        ssexp=self.ssexp_ref,
                        loadingmode='RD',#more general argument is needed
                        ftex=ftex,
                        linearization=self.linearization[j],
                        xtol=xtol,
                        maxiter=maxiter,
                        ifig=ifigdum,
                        #ax = gax_tot[i].axes_row[j][k]#[row][col]
                        )
                    ###------------------------------------------###
                    self.__nullstdout__(opt=1) #stop nulling 
                    ###------------------------------------------###

                    #print 'voce_param\n',voce_param
                    voce_param = voce_param[-1]
                    dum_voce[i,j,k] = np.array(
                        (voce_param[0], voce_param[1],
                         voce_param[2], voce_param[3])
                        )
                    plot0 = 'strstrc_ngr%s_lin%s_sample_%s.pdf'%(
                        str(self.ngr[i]).zfill(5),
                        str(self.linearization[j]).zfill(1),
                        str(k).zfill(3)
                        )
                    plot1= 'param_trend_ngr%s_lin%s_sample_%s.pdf'%(
                        str(self.ngr[i]).zfill(5),
                        str(self.linearization[j]).zfill(1),
                        str(k).zfill(3)
                        )
                    ifig0.savefig(plot0)
                    ifig1.savefig(plot1)
                    ifig0.clf(); ifig1.clf()
                    
                    # Assess the predictive capability based on
                    # uniaxial loading along the 'TD' 
                    ###------------------------------------------###
                    self.__nullstdout__(opt=0) #start nulling 
                    ###------------------------------------------###
                    #print gax_tot[i].axes_row[j][k]
                    #gax_tot[i].axes_row[j][k].plot([0,0.12,0.4],[-100,200,0],'--')
                    #figs[i].savefig('temp__.pdf')
                    #raw_input(">>> temp__.pdf has been saved check it out.")
                    grid_cax = plt.figure(ifig+i).axes[ithax]
                    
                    sdv_ref = self.each_stat(
                        param = voce_param, exp = self.ssexp_ref,
                        ngr = self.ngr[i],
                        linearization = self.linearization[j],
                        ftex = ftex,
                        loadingmode = 'RD', nsample = k,
                        ifig = ifigdum + 1,
                        sgax = grid_cax#gax_tot[i].axes_row[j][k]#[row][col]
                        )
                    
                    sdv_cf = self.each_stat(
                        param = voce_param, exp = self.ssexp_cf,
                        ngr = self.ngr[i],
                        linearization = self.linearization[j],
                        ftex = ftex,
                        loadingmode = 'TD', nsample = k,
                        ifig = ifigdum,
                        sgax = grid_cax#gax_tot[i].axes_row[j][k]#[row][col]
                        )
                    
                    grid_cax.set_xlabel(r'$\varepsilon$')
                    grid_cax.set_ylabel(r'$\sigma$[MPa]')
                    #grid_cax.set_xlim(0.,0.3)#
                    ###------------------------------------------###
                    self.__nullstdout__(opt=1) #stop nulling
                    ###------------------------------------------###
                    
                    dum_ref[i,j,k] = sdv_ref
                    dum_cf[i,j,k] = sdv_cf

                    ithax = ithax + 1 #axis index
                    pass                    
                pass
            pass
        self.sdv_ref = dum_ref
        self.sdv_cf = dum_cf
        elapsed = time.time() - t_0
        elpased, u = self.__timecnvt__(time=elapsed)
        print 'Elasped time: %4.2f [%s]'%(elapsed, u)
        # Save data onto masterd
        masterd.write(
            '%8s %8s %8s %8s %8s %8s %8s %8s %8s\n'%(
                'ngr', 'lin', 'nsample','sdv_ref','sdv_cf',
                'tau0', 'tau1', 'thet0','thet1'
                )
            )
        print '%8s %8s %8s %8s %8s %8s %8s %8s %8s\n'%(
            'ngr', 'lin', 'nsample','sdv_ref','sdv_cf',
            'tau0', 'tau1', 'thet0','thet1'
            )
        for i in range(len(self.ngr)):
            for j in range(len(self.linearization)):
                for k in range(nsample):
                    print '%8i %8i %8i %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n'%(
                        self.ngr[i],
                        self.linearization[j],
                        k,
                        self.sdv_ref[i][j][k],
                        self.sdv_cf[i][j][k],
                        dum_voce[i,j,k][0],
                        dum_voce[i,j,k][1],
                        dum_voce[i,j,k][2],
                        dum_voce[i,j,k][3],
                        )
                    masterd.write(
                        '%8i %8i %8i %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n'%(
                            self.ngr[i],
                            self.linearization[j],
                            k,
                            self.sdv_ref[i][j][k],
                            self.sdv_cf[i][j][k],
                            dum_voce[i,j,k][0],
                            dum_voce[i,j,k][1],
                            dum_voce[i,j,k][2],
                            dum_voce[i,j,k][3],                            
                            )
                        )
                    pass
                pass
            pass
        
        masterd.close()

        for i in range(len(self.ngr)):
            figs[i].savefig(
                'grid_ngr_%s.pdf'%(str(self.ngr[i]).zfill(5))
                )
            pass
        if os.name=='posix':
            os.system('cat %s'%filename)
            pass
        print '\n%s is as below: '%filename
        print open(filename, 'r').read() 



        print '********************************'
        print '    Summary of the last run\n'
        print ' grains: '
        print self.ngr
        print ' Linearization schems'
        print self.linearization
        print ' Number of sampling'
        print nsample
        print '********************************\n'
        
        print '\n********************************'
        print ' sdv_ref and sdv_cf is saved to %s'%fn
        print ' You have figures from id: %i ~ %i'%(ifig, ifig+ithax)
        print ' Tweak further as you wish'
        print '********************************\n'
        pass
        
    def __nullstdout__(self, opt=None):
        """
        opt==0 starts nulling stdout
        opt==1 stops nulling stdout
        """
        # open 2 fds
        import os
        if opt==0:
            try:
                self.null_fds = [os.open(
                        os.devnull,
                        os.O_RDWR
                        ) for x in xrange(2)]
                # save the current file descriptors to a tuple
                self.save = os.dup(1), os.dup(2)
                # put /dev/null fds on 1 and 2
                os.dup2(self.null_fds[0], 1)
                os.dup2(self.null_fds[1], 2)
                pass
            except:
                print'Could not null the fds'
                pass
            pass
                
        elif opt==1:
            # restore file descriptors
            # so I can print the results
            try:
                os.dup2(self.save[0], 1)
                os.dup2(self.save[1], 2)
                
                # close the temporary fds
                os.close(self.null_fds[0])
                os.close(self.null_fds[1])
                pass
            except:
                raise IOError
            pass
        elif opt==None:
            pass
        else:
            return -1
        pass
        
    def __timecnvt__(self, time):
        """
        Given time, return the proper unit of time
        """
        if time <= 60:
            time = time /1.
            unit='s'
            pass
        elif time>60 and time<= 3600:
            time = time / 60.
            unit='min'
            pass
        elif time>3600:
            time = time / 3600.
            unit ='h'
            pass
        return time, unit

    def __recipe__(self, opt=None):
        """
        recipe storage.
        """
        if opt==-1:
            # Recipe -1--------------------------------------------
            self.testngr =     [100,500,1000, 2000, 4000]
            self.testmaxiter = [200, 20, 20,  20, 20]
            pass        
        elif opt==0:
            # Recipe 0---------------------------------------------
            self.testngr =     [10,20,40,100,200,500,1000,2000,4000]
            self.testmaxiter = [20,20,20,200,100,100, 20,   20,  20]
            pass
        elif opt==1:
            # Recipe 1---------------------------------------------
            self.testngr =     [10,  20, 40, 100, 200, 500,1000, 2000, 4000]
            self.testmaxiter = [20, 20, 20,  200, 100,  100,  15,   10,  10]
            pass
        elif opt==2:
            #Recipe 2----------------------------------------------
            self.testngr =     [10,  20, 40, 100, 200, 500]
            self.testmaxiter = [20, 20, 20,  200, 100,  50]
            pass
        elif opt==3:
            # Recipe 3---------------------------------------------
            self.testngr=    [ 10, 100, 200,500]
            self.testmaxiter=[ 20,  20, 100, 20]
            pass
        elif opt==4:
            #Recipe 4----------------------------------------------
            self.testngr=    [ 10, 20]
            self.testmaxiter=[ 20, 20]
            pass
        elif opt==5:
            # Recipe 5---------------------------------------------
            self.testngr=    [ 10, 20]
            self.testmaxiter=[ 10, 10]
            pass
        elif opt==6:
            # Recipe 6---------------------------------------------
            self.testngr=    [ 10, 20]
            self.testmaxiter=[ 10,  5]
            pass
        elif opt==7:
            # Recipe 7---------------------------------------------
            self.testngr=    [ 10, 20, 30, 40, 60]
            self.testmaxiter=[ 10,  5,  5,  5,  5]
            pass
        elif opt==8:
            # Recipe 8---------------------------------------------
            self.testngr=    [ 100, 500, 1000, 2000]
            self.testmaxiter=[ 10,     5,    5,  5]
            pass                
        else:
            raise IOError,'input suitable opt for recipe selection'
            pass
        return self.testngr, self.testmaxiter
        
    def __pretesting__(
        self,
        ssexp_ref,
        ## Voce hardening parameters
        tau0, tau1, thet0, thet1,
        #linearization, #given as array: use self.linearization
        recipe=5, 
        loadingmode='RD',
        ):
        
        ## recipe
        testngr, testmaxiter = self.__recipe__(opt=recipe)
        self.recipe = np.array([testngr, testmaxiter])
        """
        Given the recipe, preliminarily optimize the
        tau0, tau1, thet0, thet1, for each of linearizations
        that the end-user wants.
        """
        import time, sys

        
        
        logfilename = '__pretesting__elog'
        log = open(logfilename, 'a')
        log.write('  new pretesting starts at %s\n'%time.asctime())
        testparam=np.zeros(
            (len(self.linearization), len(testngr), 4))

        totgr = []
        
        for i in range(len(self.linearization)): ## linearization loop
            log.write(' ** linearization %i\n'%self.linearization[i])
            print ' ** linearization %i'%self.linearization[i]
            ithax = 0
            for j in range(len(testngr)): ## ngr loop
                log.write('%i grains\n'%(testngr[j]))
                print '%i grains'%(testngr[j])
                print 'maxiter: %i '%testmaxiter[j]

                ###------------------------------------------###                
                self.__nullstdout__(opt=0)
                itt0 = time.time()
                param, cgr = self.main(
                    ngr=testngr[j],
                    linearization=self.linearization[i],
                    codf=self.cod,
                    tau0=tau0, tau1=tau1,
                    thet0=thet0, thet1=thet1,
                    ssexp_ref = ssexp_ref,
                    loadingmode=loadingmode,
                    maxiter=testmaxiter[j],
                    xtol=0.001
                    )
                self.__nullstdout__(opt=1)
                ###------------------------------------------###
                totgr.append(cgr) #save the RVE!
                

                elapsed = time.time() - itt0
                elapsed, u = self.__timecnvt__(time=elapsed)
                print 'Elasped time: %4.2f [%s]'%(elapsed, u)
                log.write('Elasped time: %4.2f [%s]\n'%(elapsed, u))                
                ## ----------------------------------------------------

                log.close()
                log = open(logfilename, 'a')
                tau0, tau1, thet0, thet1 = param
                log.write('parameters \n')
                log.write(
                    '%5.2f %5.2f %5.2f %5.2f \n\n'%(
                        tau0, tau1, thet0, thet1))
                log.close()
                log = open(logfilename, 'a')                
                        
                #raw_input('ddd')
                testparam[i,j][0] = tau0
                testparam[i,j][1] = tau1
                testparam[i,j][2] = thet0
                testparam[i,j][3] = thet1
                pass
            #present parameters for each of linearization
            pass
        log.close()
        return testparam, totgr

    def main(
        self, ngr, linearization, codf,
        tau0, tau1, thet0, thet1,
        ssexp_ref,
        maxiter, xtol, loadingmode
        ):
        """
        Initially, the __init__ defunc was way too large,
        unnecessarily.
        Therefore, further modulization is necessary. (2011-Sep-26)

        ngr
        linearization
           - 0:FC, 1:Affine, 2:Secant, 3:Neff
             4:Tangent, 5:Secondorder
        codf: Discrete orientation distribution file (Labotex format)
        """
        import time
        import cmb
        t_0 = time.time()

        ## RVE into temp.cmb file
        cgr = cmb.RVE(ngrain=ngr, odf=codf, cmbfile='temp.cmb')
        cgr = cgr.rve
        
        ifig = 1
        ifig0 = plt.figure(ifig)
        ifig1 = plt.figure(ifig + 1)
        #raw_input()
        
        voce_param = self.param(
            eqc=0.005,
            tau0=tau0, tau1=tau1, thet0=thet0, thet1=thet1,
            ssexp=ssexp_ref, loadingmode=loadingmode,
            ftex='temp.cmb',
            linearization=linearization,
            xtol=xtol, ifig=ifig,
            maxiter=maxiter
            )
                          
        ## Elapsed time calculation 
        t_1 = time.time()        
        elapsed = t_1 - t_0
        if elapsed>60 and elapsed<3600:
            unit='mins'
            elapsed = elapsed / 60.
            pass
        elif elapsed>=3600:
            elapsed = elapsed / 3600.
            unit='hours'
            pass
        else:
            unit='secs'
            pass
        print 'Elapsed time %5.3f (%s)'%(
            elapsed, unit)
        ## ---------------------------------------------
        return voce_param[-1], cgr
        
    def rve2ftex(self, filename=None, rve=None):
        """
        Write representative volume element to ftex
        """
        f = open(filename, 'w')
        for i in range(3):
            f.writelines('dummy\n')
            pass
        f.writelines('B  %i\n'%len(rve))
        for n in range(len(rve)):
            f.writelines(
                '%8.3f %8.3f %8.3f %10.7e\n'%(
                    rve[n][0], rve[n][1],
                    rve[n][2], rve[n][3]))
            pass
        pass

    def macro_stat(self):
        """
        macro statistical analysis
        """
        return
    
    def param(
        self, eqc,
        tau0, tau1, thet0, thet1,
        ssexp, loadingmode, ftex,
        linearization, #linearization scheme
        xtol, maxiter,
        ifig, #ifig (ifig: curve comparision, ifig+1: parameters trend)
        ax=None, #passing axes
        ):
        """
        Find the voce hardening parameters.
        Makes use of example 8 (optimization using Simplex)
        """
        #vo = voce_optimization
        #plt.ioff()
        xopt = voce_optimization(
            # hardeing parameter guess
            tau0=tau0, tau1=tau1, thet0=thet0, thet1=thet1,
            #
            eqc=0.005,
            ssexp = ssexp, # experimental str-str curve
            loadingmode='RD', # default is along RD
            # figure
            ifig=ifig,
            ax=ax,
            
            ftex=ftex,

            #convergency
            xtol = xtol, maxiter = maxiter,


            # jobid = 0
            interaction = linearization
            )
        
        #plt.gcf().clf()
        return xopt[-1]

    def local_stat(self):
        """
        local statistical analysis
        """
        return

    def each_stat(
        self, param,  ngr, 
        ftex, nsample, linearization,
        loadingmode='TD', eqc = 0.005,

        ifig=None, sgax=None, rgax=None, figname='dum',
        exp=None,  #str-str experiment
        rexp=None, #r-value experiment
        ):
        """
        statistical analysis upon the single comparison
        between a simulated stress-strain curve and ex-
        perimental stress-stress curve. Conducts simul-
        ation with given enviornment and compare to the
        designated experimental counterpart, which is
        followed by calculation of standard deviation
        of delta (difference in stress level) stress'
        standard deviation. 
        """
        import matplotlib.pyplot as plt
        import matplotlib.font_manager

        tau0, tau1, thet0, thet1 = param

        ## single crystal file maker ----------------------- ##
        fsx = 'sx_for_each_stat_ngr-%s_linear-%s_id-%s.sx'%(
            str(ngr).zfill(5), str(linearization),
            str(nsample).zfill(8))
                

        sx_maker.cubic(
            filename=fsx,
            tau0=tau0, tau1=tau1, thet0=thet0, thet1=thet1,
            iopsysx=1, dependency='iso')
        ## ------------------------------------------------- ##


        ## experimental data ------------------------------- ##
        Eexp = exp[0]
        Sexp = exp[1]
        if rexp!=None:    # r-value exp data
            #Ee = rexp[0]
            Rexp = rexp[1]
            pass
        ## ------------------------------------------------- ##
            

        ## ------------------------------------------------- ##
        stp = int( max(Eexp) * 1.05 / eqc ) # 5% overshoot
        job = vpsc(
            mode=loadingmode, fsx=fsx, texture=ftex,
            stp=stp, eqc=eqc, jobid=315, #jobid is fixed to be 315
            interaction=linearization
            )
        rst = job.run() #masterdata
        ## ------------------------------------------------- ##
        

        if loadingmode=='td' or loadingmode=='TD':
            eps = rst['eps'][1]
            sig = rst['sdev'][1]
            dbar = rst['dbar'][1]
            pass
        elif loadingmode=='rd' or loadingmode=='RD':
            eps = rst['eps'][0]
            sig = rst['sdev'][0]
            dbar = rst['dbar'][0]
            pass
        else: raise IOError,'unexpected loading mode'
        

        s11 = []; e11 = []

        
        if loadingmode=='RD' or loadingmode=='rd':
            xi = 0; yi = 0
        elif loadingmode=='TD' or loadingmode=='td':
            # mode='TD' first rotates the polycrystalline aggregate
            # then tensile. The polycrystalline aggregate is rotat-
            # ed back again. Thus, [0,0] indices should be used
            # for having the right component in the case of
            # uniaxial tension
            xi = 0; yi = 0
        else: raise IOError, 'unexpected loading mode'


        l = loadingmode
        if l in ['td','TD','RD','rd']:
            r = []
        for i in range(len(sig)):
            s11.append(sig[i][xi][yi]-sig[i][2][2]) # ensure to subtract s33
            e11.append(eps[i][xi][yi])
            if l in ['td','TD','RD','rd']:
                r.append(dbar[i][1,1]/dbar[i][2,2])
                pass
            else: raise IOError,'unexpected case'
            pass


        ## ---------------------------------------------------- ##
        ## interpolates the experimental data in order to ----- ##
        ## compare the stress level directly
        ys = [] # flow stress
        yr = [] # r-value
        diff = [] #
        diffr = [] #r-value difference
        
        for i in range(len(e11)):
            if e11[i]< Eexp[-1]:
                """
                Interpolate to find the stress level at the each of
                simulated deformation increment..
                """
                ## stress interpolation ----------------------- ##
                temps = _interpolate_2(x=e11, y=s11, x0=Eexp[i])
                ys.append(temps)  #interpolated simulated value
                diff.append(
                    abs(ys[i] - Sexp[i])
                    )

                ## r-value interpolation 
                if rexp!=None:
                    tempr = _interpolate_2(x=e11, y=r, x0=Eexp[i])
                    yr.append(tempr)
                    diffr.append(
                        abs(
                            yr[i] - Rexp[i]
                            )
                        )
                ## -------------------------------------------- ##
                pass
            
            else: pass
            
            pass
        ## interpolation ends --------------------------------- ##
        ## ---------------------------------------------------- ##
        
        
        if ifig!=None:

            ## flow curve
            fig = plt.figure(ifig)
            try:fig.clf()
            except: pass
            ax = fig.add_subplot(111)
            ax.plot(e11, s11, 'o', label='sim')
            ax.plot(Eexp, Sexp, 'x', label='exp')
            fig.savefig('%s_fcurv.pdf'%figname)
            fig.clf()

            if rexp!=None:
                #fig = plt.figure(ifig)
                ax = fig.add_subplot(111)
                ax.plot(e11, r, 'o', label='sim')
                ax.plot(Eexp, Rexp, 'x', label='exp')
                fig.savefig('%s_rval.pdf'%figname)
                fig.clf()
                pass
            
        if sgax!=None:
            if loadingmode=='RD' or loadingmode=='rd':
                color = 'red'
            elif loadingmode=='TD' or loadingmode=='td':
                color = 'blue'
            else: raise IOError, 'unexpected loadingmode'
            lins = __lins__(lin=linearization)
            sgax.plot(
                Eexp, Sexp, 'x', color=color,
                label='EXP %s'%loadingmode
                )
            label = '%s, %s, %i grains'%(lins, loadingmode, ngr)
            sgax.plot(e11, s11, 'o', label=label, mec = color, mfc = 'None') #color=color)
            prop = matplotlib.font_manager.FontProperties(
                size=5.#/ len(self.linearization) #/ self.nsample
                )

            ## x ticks setting : 3 ticks only ##
            xl = sgax.get_xlim()[0]
            xh = sgax.get_xlim()[1]
            xh = (int(xh / 0.1) + 1.0)* 0.1            
            sgax.legend(loc='lower right', prop=prop)
            sgax.set_xticks(
                np.linspace(
                    xl,
                    xh,
                    3
                    )
                )

            ## y ticks setting : 4 ticks only ##
            yl = sgax.get_ylim()[0]
            yh = sgax.get_ylim()[1]
            yh = (int(yh / 100) + 1) * 100
            sgax.legend(loc='lower right', prop=prop)
            sgax.set_yticks(
                np.linspace(
                    0, #yl,
                    yh,
                    3
                    )
                )            
            pass

        if rgax!=None and rexp!=None:
            if loadingmode=='RD' or loadingmode=='rd':
                color = 'red'
            elif loadingmode=='TD' or loadingmode=='td':
                color = 'blue'
            else: raise IOError, 'unexpected loadingmode'

            ## plotting
            lins = __lins__(lin=linearization)
            rgax.plot(
                Eexp, Rexp, 'x', color=color,
                label='EXP %s'%loadingmode
                )
            
            label = '%s, %s, %i grains'%(lins, loadingmode, ngr)
            rgax.plot(
                e11, r, 'o', label=label,
                mec = color, mfc = 'None'
                )
            ##

            
            prop = matplotlib.font_manager.FontProperties(
                size=5.#/ len(self.linearization) #/ self.nsample
                )

            if True:
            ## x ticks setting : 3 ticks only ##
                xl = rgax.get_xlim()[0]
                xh = rgax.get_xlim()[1]
                xh = (int(xh / 0.1) + 1.0) * 0.1
                rgax.legend(loc='lower right', prop=prop)
                rgax.set_xticks(
                    np.linspace(
                        xl,
                        xh,
                        3
                        )
                    )
                
            ## y ticks setting : 4 ticks only ##
                yl = rgax.get_ylim()[0]
                yh = rgax.get_ylim()[1]
                yh = (int(yh / 0.5) + 1) * 0.5
                rgax.legend(loc='lower right', prop=prop)
                rgax.set_yticks(
                    np.linspace(
                        0., #yl,
                        yh,
                        3
                        )
                    )                        
                pass
            pass



        diff = np.array(diff)
        stdev = (diff**2).sum()
        stdev = math.sqrt(stdev) / (len(diff)-1)

        if rexp!=None:
            diffr = np.array(diffr)
            stdev_r = (diffr**2).sum()
            stdev_r = math.sqrt(stdev_r) / (len(diffr)-1)
            pass

        if rexp==None and exp!=None:
            return stdev #standard deviation
        elif rexp!=None and exp!=None:
            return stdev, stdev_r
        else:
            raise IOError, 'Unexpected case occured'

    
def __lins__(lin=None):
    if lin==0: return 'Taylor'
    elif lin==1: return 'Affine'
    elif lin==2: return 'Secant'
    elif lin==3: return r'$n^{eff}=10$'
    elif lin==4: return 'Tangent'
    elif lin==5: return 'Second-Order'
    else: raise IOError
    return -1

    
