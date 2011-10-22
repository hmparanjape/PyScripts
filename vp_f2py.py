"""
FORTRAN OBJECT wrapped by the f2py module
As of now, it is known to be only compatible with Linux.

Set-off date: 2011-02-06

comments: 
  (2011-02-07)
   It's much faster than vpsc_param.vpsc class. Is 
  this due to compiler? Here in this project, for
  unknown reasons, only linux fortran compiler is
  working than g95. No further specification for
  the compiler type is necessary or should be banned.

   'VPSC7.for' is modified to make consequent processes 
  smooth. In 'vpsc7.dim', two parameters are added,
  NMXSTP and NMXPRCS. They denote the maximum number
  of steps for each process and the maximum number of 
  process.

   Then in the down below, a common block is added,
  named 'f2py'. In it, DBAR_ AND SBAR_ is assigned.
  (2011-02-23)
  
   vpsc7, vpsc7_gfortran, vpsc7_g95 are the possible
  wrapped *.so files. Each of which has a appending
  compiler names.
"""
#------------------------------------------------------------
## proper vpsc7?
ivpsc7, ivpsc7_gfortran, ivpsc7_g95 = False, False, False

#1 VPSC7
try: reload(vpsc7)
except: 
    try: import vpsc7 ; reload(vpsc7)
    except: print " No vpsc7 module loaded"
    else: print "vpsc7 was loaded"
else:
    print "vpsc7 was loaded"
    ivpsc7 = True
    
#2 vpsc7_gfortran
try: reload(vpsc7_gfortran)
except:
    try: import vpsc7_gfortran; reload(vpsc7_gfortran)
    except: print " No vpsc7_gfortran"
    else: print "vpsc7_gfortran was loaded"
else:
    print "vpsc7_gfortran was loaded "
    ivpsc7_gfortran=True

#3 vpsc7_g95 
try: reload(vpsc7_g95)
except:
    try: import vpsc7_g95; reload(vpsc7_g95)
    except: print " No vpsc7_g95 loaded"
    else: print "vpsc7_g95 was loaded"
else:
    print "vpsc7_g95 was loaded"
    ivpsc7_g95 = True
#------------------------------------------------------------

import vpsc_in, sx_maker, chg_basis
import os, glob
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math, crss, euler, shutil, pf 
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

def __rmall__():
    """ Deletes all the existing files """
    osname = os.name
    if osname=='nt':
        print "This module is not compatible to nt system."
        print "You mus use posix system"
    elif osname=='posix':
        df = glob.glob('*.out')+ glob.glob('*.OUT')
        for i in df: os.remove(i)

def fin(filename,nhead,ix,*args ):
    """
    with number of header given
    Things are returned as strings
    """
    FILE = open(filename,'r')
    lines = FILE.readlines()
    FILE.close() ## FILE CLOSED
    LINE = [ ]
    lines = lines[nhead:]
    for i in range(len(lines)):
        cl = lines[i].split()
        LINE.append(cl)
    LINE = np.array(LINE)
    LINE = LINE.transpose()
    ycolumns =[]
    for i in range(len(args)):
        ycolumns.append(LINE[args[i]])
    return LINE[ix], ycolumns

def recompile(compiler='gfortran'):
    """
    Handy recompile module for fortran object vpsc
    from signature file to compilation.

    Try to compile the code with the gfortran compiler,
    whenever possible.
    """
    try:
        os.remove('vpsc7.so')
        os.remove('vpsc7.pyf')
        os.remove('vpsc7_g95.so')
        os.remove('vpsc7_gfortran.so')
    except: pass
    # """
    # #signature file making
    # iflag = os.system('f2py VPSC7.for -m vpsc7 -h vpsc7.pyf --overwrite-signature')
    # if iflag!=0:
    #     print "Something's come up during signature file making"
    #     return -1
    # #raw_input('Enter to proceed to wrap the code with signature file >>')
    # """

    # """
    # #Wrapping it through signature file
    # iflag = os.system("f2py -c --fcompiler=%s --f77flags='%s' vpsc7.pyf VPSC7.for"%
    #                   #('g95','-freal=zero -r8') #coupled with g95 compiler
    #                   ('gnu95','-fdefault-real-8 -fdefault-double-8  -fno-align-commons') #coupled with gfortranf77 flags
    #                   #-fdefault-real-8
    #                   #-fdefault-double-8
    #                   #-fdefault-integer-8
    #                   #-fno-align-commons
    #         )
    # """

    ## Wrapping without signature file making
    ###
    if compiler=='g95':
        cmppath = os.popen('which g95').readline().split('\n')[0]
        iflag = os.system(
            "f2py -c %s %s %s -m %s  VPSC7.for -DF2PY_REPORT_ON_ARRAY_COPY=1>f2pytemp"%
            (
                #'--fcompiler=gnu95', #compiler options
                "--f90exec='%s' --f77exec='%s'"%(cmppath, cmppath), #compiler path
                "--f77flags='-d8 -r8 -fzero'",#"--f77flags='-fdefault-double-8 -fdefault-real-8'",
                "--f90flags='-d8 -r8 -fzero'", #"--f90flags='-fdefault-double-8 -fdefault-real-8'",
                #'-fdefault-real-8 -fdefault-double-8 -fno-align-commons', #f77 flags
                "vpsc7_g95" #module name
                )
            )
        pass

    elif compiler=='gfortran':
        cmppath = os.popen('which gfortran').readline().split('\n')[0]
        iflag = os.popen(
            "f2py -c %s %s %s -m %s VPSC7.for -DF2PY_REPORT_ON_ARRAY_COPY=1 > f2pytemp"%
            (
                ## location of gfortran bin
                "--f90exec='%s' --f77exec='%s'"%(cmppath, cmppath),
                ## flags                
                "--f77flags='%s %s %s %s %s %s'"%(
                    '-fdefault-real-8', '-fdefault-double-8',
                    '-fno-align-commons', '-finit-local-zero',
                    '-finit-real=zero', '-finit-integer=0'),
                "--f90flags='%s %s %s %s %s %s'"%(
                    '-fdefault-real-8', '-fdefault-double-8',
                    '-fno-align-commons', '-finit-local-zero',
                    '-finit-real=zero', '-finit-integer=0'),
                ## end of flags
                "vpsc7_gfortran"
                )
            )
        pass
    elif compiler=='gfortran_mac':
        cmppath = os.popen('which gfortran').readline().split('\n')[0]
        iflag = os.system(
            "f2py -c %s %s %s -m %s VPSC7.for -DF2PY_REPORT_ON_ARRAY_COPY=1 > f2pytemp"%
            (
                "--f90exec='%s' --f77exec='%s'"%(cmppath, cmppath),
                "--f77flags='-fdefault-real-8 -fdefault-double-8'",
                "--f90flags='-fdefault-real-8 -fdefault-double-8'",
                "vpsc7_gfortran"
                )
            )        
        pass
    else: return -1
    return 0

class vpsc:
    """
    The OBJETIVE VISCO-PLASTIC SELF-CONSISTENT script

    It will be designed for a single run.
    However, the data and conditions are saved so that
    post-execution access to data is feasible.

    Data are saved internally not to files.
    Having files as an access point of store is not 
    efficient in point of management and post-processing

    MPI??? (I'm not sure... yet!!!)

    comments
      2011-02-06: f2py module is used now for wrapping
      the fortran code under the project name of vpsc_f2py
      
      2011-02-06: This is the f2py wrapper version of
      the original vpsc class embedded in vpsc_param.py
      script. Further modification will be added to this.

      2011-03-15: The inverse pole figure representation
      of stabilizing rotation of tensile axis is shown.
      In the meantime, a model has been written for objective
      plotting of pole figure of which I'm very proud!

      2011-04-27: The VPSC7.for code is further modified
      for file I/O not to be duplicated. Doing so requires
      to pass the jobid argument.

      2011-07: class vpsc has been used in many applicational
      scripts under many different circumstances. Please find
      the vpsc_examples.py useful and explore it.
      
    Author: Youngung Jeong
            Materials Mechanical Laboratory
            Graduate Institute of Ferrous Technology
            Pohang University of Science and Technology

            Some minor modification of the script is performed
            during 2011-July 2013-Jan
            National Institute of Standards and Technology
            Metallurgy Devision
    """
    class data:
        pass
    class VPSCin:
        nph = 0
        pass
    class ph:
        pass
    def manual(self):
        """ manual in order to use class vpsc.so """
        print " Still empty!! "
    def __init__(self, iwait=True, 
                 ## grain morphology contorl
                 ## angles and principals lengths of 
                 ## axes of the ellipsoid

                 ## loading conditions
                 u11=None, u22=None ,
                 stp=None, eqc=None, ict=None,
                 stp0=None, stp1=None,
                 mode=None, 
                 ang=None,
                 shear_ang=None,  #active when mode='ysangle'
                 ysxc=1,ysyc=2, # yscomponent mode

                 #ang is useful when mode='ten_ang'
                 #for tension along after rotating as much as ang

                 ##texture file
                 texture=None,
                 gr = None, #ex: [[45.0, 30.0, 15.0, 0.4],  [0., 10., 75.0, 0.6]]

                 ## single crystal file
                 fsx=None,

                 ## Enviroments
                 del_crss=True,

                 ## misc - options  (passed to __vpsc_in__ module)
                 irecover=0, isave=0, icubcomp=0, 
                 nwrite=0, ihardlaw=0,
                 iflu=0, iratesens=0, interaction=3, 
                 iupdate=[1,1,1,0],  # ori, shape, hardening, itran
                 nneigh=0,

                 ### precision settings for convergence prodcedure
                 err=0.001, errd=0.001, 
                 errm=0.001, errso=0.01,

                 nmxiter=100, exiter=100, initer=100,
                 irsvar=0, xrsini=2, xrsfin=10, xrstep=2,

                 ### manual prcs 
                 prcs=None,

                 ### inclusion property (angle and axes)
                 init_eld_rat=[[ 1., 1., 1.], [ 1., 1., 1.], [ 1., 1., 1.]],
                 init_eul_axe=[[90.,90.,90.], [90.,90.,90.], [90.,90.,90.]],

                 ## job id
                 jobid=0,
                 ):   #stp can be given (in general it's step)
        """
        Comments 
        (2011-01-30)
          Is it okay for me to include all the possible
          arguments into def __init__? It is not quite sure
          but due to python's ability to have defaulted 
          arguments in modules, more arguments does not really
          mean more confusions. However, having more arguments
          is always making me wonder it is really necessary.

        (2011-02-07)
          This is now in the vp_f2py.py. The original vpsc 
          class is in vpsc_param.py script. 

        Arguments:
          mode = 'ten_ang=ang_ten', 'RD=rd','DD=dd','TD=td',
                 'BULGE=bulge=bu=BU','SIMPLESHEAR=SS=SSRD',
                 'SSTD=sstd','SSDD=ssdd',
                 'lankf=LANKFORD=lankford=lnkf',
                 'inplanebiaxial=ipb',
                 'inplanebiaxialfc=ipbfc',
                 'ys=YS',
                 'ev_pcys'
                 'tc=tensioncompression'
                 'ssfr=SSFR' Forward and reverse simple shear
                 
          u11, u22
          stp, eqc
          ict, ang(effective only when mode =='ten_ang')
          stp0, stp1 : for 'ssfr1' mode
          
          texture, fsx, gr 
          del_crss=True (other wise crss*.out file are stacked)
          
          irecover=0, isave=0, icubcomp=0, nwrite=0, ihardlaw =0,
          iflu=0, iratesens=0, interaction=3(neff=10), iupdate=[1,1,1,0],
          nneigh=0

          ### precision settings for convergence prodcedure
          err=0.000001, errd=0.000001, 
          errm=0.000001, errso=0.00001,

          nmxiter=1000, exiter=1000, initer=1000,
          irsvar=0, xrsini=2, xrsfin=10, xrstep=2,

          ### manual prcs
          prcs= None,
          
        """
        ### The old style of allocating jobs into several pieces of 
        ### prescribed directories. It turned out that it takes too
        ### much time to mkdir and copy files. Thus here in vp_f2py
        ### no more such jobs will be done.

        ## job id is fixed
        if jobid > 9999:
            print 'jobid must not exceed 9999'
            raise IOError
        self.jobid = str(jobid).zfill(4)
        
        self.cwd = os.getcwd()   #full path name of current working directory 
        self.osname = os.name
        
        self.histf = [] ## initialized here but written in the __historyfileprocess
        self.vpscsx = []
        
        #--------------------------------------------------------
        #   VPSC input (texture, morphology, phase)
        #     number of phase, phase volume fraction, grain shape
        #     ftx, fsx, fax, 
        #--------------------------------------------------------

        self.VPSCin.nph = 1
        self.VPSCin.ph = []
        self.VPSCin.wph = [1.0, 0.9][0:self.VPSCin.nph]
        self.VPSCin.ishape = [0,0][0:self.VPSCin.nph]

        #TEXTURE FILE
        if texture!=None:  
            if type(texture).__name__=='str':
                self.ftx = [texture]  #when the texture is given as a string
                #print 'The given texture file = %s'%texture
            elif type(texture).__name__=='list':
                self.ftx = texture
        else:
            if gr!=None:
                temp_filename='temp_vpsc_gr.tex'
                FILE = open(temp_filename, 'w')
                FILE.writelines(
                    'dummy\ndummy\ndummy\n%s %i\n'%('B',len(gr)))
                for i in range(len(gr)):
                    FILE.writelines('%8.3f %8.3f %8.3f %11.3e\n'%(
                            gr[i][0], gr[i][1], gr[i][2], gr[i][3]
                            ))
                FILE.close()
                self.ftx=[temp_filename]
                pass
            else:
                print " You have not input any texture file to class vpsc"
                print " You must assign a file name"
                raise IOError
            pass

        #SINGLE CRYSTAL FILE
        self.fsx = [ ]
        if fsx!=None: 
            if type(fsx).__name__=='list': pass
            elif type(fsx).__name__=='str': fsx=[fsx]
            else: raise IOError
            self.fsx=fsx        
        elif fsx==None:
            for i in range(self.VPSCin.nph):
                sxfilename = 'temp_%s.sx'%str(i).zfill(3)
                self.fsx.append(sxfilename)
                sx_maker.cubic(
                    filename=sxfilename, #self.VPSCin.ph[i].fsx,
                    #Hardening matrix
                    hii=[1.0,1.0][i], 
                    hij=[1.0,1.0][i], 
                    hb=[1.0,1.0][i],
                    hp=[1.0,1.0][i],
                    #Voce parameters
                    tau0=[1.,1.][i], tau1=[0.,0.][i], 
                    thet0=[0.,0.][i], thet1=[0.,0.][i],
                    
                    ####
                    dependency=['iso','iso'][i],  #'iso', 'indv'
                    #dependency=['indv','indv'][i],  #'iso', 'indv'
                    
                    #slip-twiningg system
                    n=[ [[1,1,1]], [[1,1,0]] ][i], 
                    b=[ [[1,1,0]], [[1,1,1]] ][i],
                    iopsysx=[1,1][i],
                    #iopsysx=[0,0][i],
                    )            
            ## In this case, an isotropic aggregate is made
            # self.fsx = ['iso.sx',]
            # sx_maker.cubic(filename=self.fsx[0], hii=1., hij=1., hb=1., hp=1.,
            #                tau0=1.0, tau1=0, thet0=0, thet1=0,
            #                iopsysx=1 , dependency='iso' )
            # os.system('cat iso.sx')
            
            # self.fsx = ['sx%sph01.sx'%os.sep,
            #             'sx%sph02'%os.sep][0:self.VPSCin.nph]
        else: raise IOError

        #AXES FILE
        self.fax = ['shape1.100'][0:self.VPSCin.nph]
        for i in range(self.VPSCin.nph): 
            self.VPSCin.ph.append(self.ph)

        ANG = []  #Euler angles of inclusion principal axes 
                  #(for multiple inclusions)

        AXE = []  #Principal lengths of the inclusion
                  #(for multiple inclusions)

        for i in range(self.VPSCin.nph):
            self.VPSCin.ph[i].ftx = self.ftx[i]
            self.VPSCin.ph[i].ngr = self.__seekngr__(filename=self.ftx[i])
            self.VPSCin.ph[i].fsx = self.fsx[i]
            self.VPSCin.ph[i].fax = self.fax[i]
            self.VPSCin.ph[i].grainshapectrl = 0
            self.VPSCin.ph[i].fragmentn = 0
            self.VPSCin.ph[i].criticalaspect = 0
            self.VPSCin.ph[i].init_elpsd_rat=init_eld_rat[i]
            self.VPSCin.ph[i].init_Eul_axe  =init_eul_axe[i]
            ANG.append(self.VPSCin.ph[i].init_Eul_axe)   #angle
            AXE.append(self.VPSCin.ph[i].init_elpsd_rat) #axis

        """
        -- initialization
             prcs, step, eqincr, nhist
        """
        utemp = [u11,u22]
        if any(utemp[i]==None for i in range(2)): pass
        else:
            u33= -(u11+u22)  #INCOMPRESSIBILITY
            u12=0; u13=0; u21=0; u23=0; u32=0; u31=0

        if prcs==None: self.prcs = []
        else: self.prcs=prcs

        self.step = 0; self.eqincr = []
        self.nhist = 0; self.nrot = 0
        
        """
        -- PROCESSES (LOADING, PCYS, LANK .. )
        """
        #### general parameters over historyfileprocesses
        if stp==None: stp = 5 
        if eqc==None: eqc = 0.005 
        if ict==None: ict = 7
        if ang==None: ang = 0.

        #### Global mode variable
        self.mode = mode

        if mode==None: 
            # the below block according to their needs.
            if prcs!=None: pass
            else:
                print "Prcs and mode is complementary to each other"
                print "You must assign either of the two."
                raise IOError
            """
            self.__historyfileprocess(step=80, eqincr=0.005, ictrl=7, 
            mode='inplanebiaxial',u11=u11,u22=u22)
            self.__rotate(ang=-45)
            self.__lankfordprocess(step=5)
            self.__historyfileprocess(step=20, eqincr=0.005,
                                      ictrl=7, mode='unitension')
            self.__pcysplprocess(step=9, shear_ang=90.)
            self.__rotate(ang=45)
            """
            pass
        elif mode=='ten_ang' or mode=='TEN_ANG' or mode=='ang_ten' or mode=='ANG_TEN':
            ## tension along a particular direction
            self.__rotate(ang=ang)
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, mode='unitension',
                                      jobid=self.jobid)
            self.__rotate(ang=-ang)
            pass
        elif mode=='ten_ang2' or mode=='TEN_ANG2':
            """
            A sequence of uniaxial loadings after a planar rigid body rotation
            """
            if None in [stp0, stp1]: raise IOError,'stp0 or stp1 are missing'
            self.__historyfileprocess(step=stp0, eqincr=eqc,
                                      ictrl=ict, mode='unitension',
                                      jobid=self.jobid)
            self.__rotate(ang=ang)
            self.__historyfileprocess(step=stp1, eqincr=eqc,
                                      ictrl=ict, mode='unitension',
                                      jobid=self.jobid)
            pass
        elif mode=='RD' or mode=='rd':
            #RD tension
            if interaction==0:lmode='unitensionfc'
            else: lmode='unitension'
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode=lmode,
                                      jobid=self.jobid)
            pass
        elif mode=='-RD' or mode=='-rd':
            #RD compression
            if interaction==0:lmode='unicompressionfc'
            else: lmode='unicompression'
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode=lmode,
                                      jobid=self.jobid)
            pass            
        elif mode=='DD' or mode=='dd':
            #DD tension
            if interaction==0: lmode='unitensionfc'
            else: lmode='unitension'
            self.__rotate(ang=45)
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode=lmode,
                                      jobid=self.jobid)
            self.__rotate(ang=-45)
            pass
        elif mode in['-td','-TD']:
            ## TD uniaxial compression
            self.__rotate(ang=90.)
            if interaction==0:lmode='unicompressionfc'
            else: lmode='unicompression'
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode=lmode,
                                      jobid=self.jobid)
            self.__rotate(ang=-90.)
            pass                        
        elif mode=='TD' or mode=='td':
            #TD tension
            # The rotation is now excluded
            self.__rotate(ang=90)
            if interaction==0: lmode='unitensionfc'
            else: lmode='unitension'
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode=lmode,
                                      jobid=self.jobid)
            self.__rotate(ang=-90)     
            pass
        elif mode=='td_without_rotation':
            ## tensile along e22 without causing any regid body rotatin
            if interaction==0: raise IOError, 'add unitensiony_fc'
            else: lmode='unitensiony'
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode=lmode,
                                      jobid=self.jobid)
            pass
        elif mode=='BULGE' or mode=='bulge' or mode=='bu' or mode=='BU':
            #BULGE
            self.__historyfileprocess(
                step=stp, eqincr=eqc, 
                ictrl=ict, mode='bulge', jobid=self.jobid)
            pass
        elif mode=='SIMPLESHEAR' or mode=='SS' or mode=='SSRD':
            #SIMPLESHEAR
            self.__historyfileprocess(
                step=stp, eqincr=eqc, 
                ictrl=ict, mode='simpleshearfc', jobid=self.jobid)
            pass
        elif mode in ['-ss','-ssrd','-simpleshear',
                      '-SS','-SSRD','-SIMPLESHEAR']:
            mode=='SSR'
            pass
        elif mode=='SIMPLESHEARREVERSE' or mode=='SSR' or mode=='SSRDR':
            #REVERSED SIMPLESHEAR
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, mode='SSR', jobid=self.jobid)
            pass
        elif mode=='SIMPLESHEARREVERSEFC' or mode=='SSRFC' or mode=='SSRDRFC':
            #REVERSED SIMPLESHEAR FULLY CONSTRAINT
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, mode='SSRFC', jobid=self.jobid)
            pass        
        elif mode=='SSTD' or mode=='sstd':
            #SIMPLESHEAR ALONG TRANSVERSE DIRECTION
            self.__rotate(ang=90)
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, mode='simpleshearfc',
                                      jobid=self.jobid)
            self.__rotate(ang=-90)
            pass
        elif mode=='SSDD' or mode=='ssdd':
            #SIMPLE SHEAR ALONG DIAGONAL DIRECTION
            self.__rotate(ang=45)
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, mode='simpleshearfc',
                                      jobid= self.jobid)
            self.__rotate(ang=-45)
            pass
        elif mode=='ssfr' or mode=='SSFR':
            #Forward and Reverse simple shear
            self.__historyfileprocess(
                step=stp0, eqincr=eqc, ictrl=ict,
                mode='simpleshear', jobid=self.jobid)
            self.__historyfileprocess(
                step=stp1, eqincr=eqc, ictrl=ict,
                mode='simpleshearreverse', jobid=self.jobid)
            pass
        elif mode=='ssfr1' or mode=='SSFR1':
            if None in [stp0,stp1]:
                raise IOError,'stp0 and stp1 should be given'
            self.__historyfileprocess(
                step=stp0, eqincr=eqc, ictrl=ict,
                mode='simpleshear', jobid=self.jobid)
            self.__historyfileprocess(
                step=stp1, eqincr=eqc, ictrl=ict,
                mode='simpleshearreverse', jobid=self.jobid)
            pass
        elif mode=='ss_ang':
            #simple shear, a rotation, and simple shear
            self.__historyfileprocess(
                step=stp, eqincr=eqc, ictrl=ict,
                mode='simpleshear', jobid=self.jobid)
            self.__rotate(ang=ang)
            self.__historyfileprocess(
                step=stp, eqincr=eqc, ictrl=ict,
                mode='simpleshear', jobid=self.jobid)
            pass
        elif mode=='lankf'or mode=='LANKFORD'or mode=='lankford'or mode=='lnkf':
            self.__lankfordprocess(step=stp)
            self.__rotate(ang=90.)
            self.__lankfordprocess(step=stp)
            self.__rotate(ang=-90.)
            pass
        elif mode=='inplanebiaxial' or mode=='ipb':
            self.__historyfileprocess(
                step=stp, eqincr=eqc,u11=u11,u22=u22,
                ictrl=ict, mode='inplanebiaxial',jobid= self.jobid)
            pass
        elif mode=='inplanebiaxialfc' or mode=='ipbfc':
            self.__historyfileprocess(step=stp, eqincr=eqc,u11=u11,u22=u22,
                                      ictrl=ict, mode='inplanebiaxialfc',
                                      jobid= self.jobid)
            pass

        elif mode=='ys' or mode=='YS' or mode=='YIELDSURFACE' or mode=='pcyspl':
            self.__pcysplprocess(step=stp, shear_ang=90.) #step == stp
            pass
        elif mode in ['yscomponent','ys_components']:
            self.__pcysplprocess(step=stp, shear_ang=90, xcomp=ysxc, ycomp=ysyc)
            pass
        elif mode=='ysangle':
            self.__pcysplprocess(step=stp, shear_ang=shear_ang) #
            pass
        elif mode=='ev_pcys':
            # ntot = 3 ## INTERMIDEATE STOPPING NUMBER TO PERFORM PCYS :hardwired
            # totalstrain = stp*eqc
            # eachstep_inc = totalstrain/ntot
            # eachstep_inc = eachstep_inc/eqc
            # eachstep_inc = int(eachstep_inc)
            # print "\n loading is %i times in Total"%ntot
            # print " Total Strain imposed is %6.3f"%totalstrain
            # print " Each loading has %i"%eachstep_inc,
            # print " number of incremental step"
            # flag = raw_input('type q to raise IError >>')
            # if flag=='q' or flag =='Q':
            #     raise IOError
            #     pass

            nprobing = 32 #number of probings
            self.__pcysplprocess(
                step=nprobing, shear_ang=90.,
                jobid=self.jobid) ##initial PCYS-pl
            
            for i in range(ntot):
                # self.__historyfileprocess(step=eachstep_inc, eqincr=eqc,
                #                           u11=u11,u22=u22,
                #                           ictrl=ict, mode='inplanebiaxialfc')
                self.__historyfileprocess(step=eachstep_inc, eqincr=eqc,
                                           mode='unitension', jobid= self.jobid)
                self.__pcysplprocess(step=nprobing, shear_ang=90.)
            
        elif mode=='eqeps': #Equivalent strain component case
            """
            Equivalent in-plane strain (balanced strain case)
            """
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, u11=0.5, u22=0.5,
                                      mode='inplanebiaxial',
                                      jobid= self.jobid)
            pass
        elif mode=='tc' or mode=='tensioncompression': # tension-compression
                         # (checking for dislocation density hardening model
                         # (Implementation of Rauch et al.'s scheme)
            # RD tension-compression
            if interaction==0:
                lmodet='unitensionfc'; lmodec='unicompressionfc'
            else: lmodet='unitension'; lmodec='unicompression'
            self.__historyfileprocess(step=stp0, eqincr=eqc,
                                      mode=lmodet, jobid=self.jobid)
            self.__historyfileprocess(step=stp1, eqincr=eqc,
                                      mode=lmodec, jobid=self.jobid)
        else:
            print "The passed mode %s to class vpsc is not among "%mode,
            print "available mode options.\n Please check the code!!"
            raise IOError



        #VPSC_INPUTFILE
        self.__vpsc_in__(
            nelem=1, nph=self.VPSCin.nph,
            wph=self.VPSCin.wph, ishape=self.VPSCin.ishape,
            #wph=[1.0], ishape=[0],

            tfile=self.ftx,     #tex
            sxfile=self.fsx,    #sx 
            fileaxes=self.fax,  #shape

            #####  precision settings for convergence prodcedure #####
            err=err, errd=errd,
            errm=errm, errso=errso,

            nmxiter=nmxiter, exiter=exiter, initer=initer,
            irsvar=irsvar, xrsini=xrsini, xrsfin=xrsfin, xrstep=xrstep,
            ####-------------------------------------------------------

            ### Ellipsoid angle and principal lengths
            ellipangl = ANG,
            ellipaxes = AXE,

            ## Controls
            irecover=irecover, isave=isave, icubcomp=icubcomp, nwrite=nwrite,
            ihardlaw=ihardlaw, iflu=iflu, iratesens=iratesens, 
            interaction=interaction, iupdate=iupdate, nneigh=0, 

            prcs=self.prcs, npro=len(self.prcs)/2
            )

        FILE = open('VPSC7_%s.in'%(self.jobid),'r')
        lines = FILE.read()
        self.vpscin = lines        
        FILE.close()
        for i in range(len(self.fsx)):
            FILE = open(self.fsx[i],'r')
            lines = FILE.read()
            self.vpscsx.append(lines)
        FILE.close()

        self.datamaster = {} #INITIALIZATION UPON self.datamaster

    def VPSC_process(self, u11=None, u22=None):
        """
        Aims at further automation of processes in VPSC7.in file
        Still need to be worked out, though...
        Arguments: u11, u22
        """
        if any([u11,u22][i]==None for i in range(2)): pass
        else:
            u33 = -(u11+u22)
            u12 = 0; u13 = 0; u21 = 0; u23 = 0; u32 = 0; u31 = 0
        self.prcs = []; self.step = 0; self.eqincr = []
        self.nhist = 0
        # Yield surface projection
        self.__pcysplprocess(step = 9, shear_ang=90.)        
        # Loading through history file #1 -unitension
        self.__historyfileprocess(step=20, eqincr=0.005, ictrl=7, mode='unitension')
        # Lankford coefficient probing
        self.__lankfordprocess(step = 2)
        # Loading through history file #2 -inplane biaxial
        self.__historyfileprocess(step=1, eqincr=0.005, 
                                  ictrl=7, mode='inplanebiaxial', 
                                  u11=0.5, u22=0.3)

    def __historyfileprocess(self, step=1, eqincr=0.005, mode=None,
                             ictrl=7, u11=None, u22=None, jobid=None):
        ####
        ####  HISTORY FILE PROCESS 
        ####
        """
        LIST OF MODES : 'unitension', 'unitensionfc',
            'inplanebiaxial', 'bulge', 'simpleshear', 
            'simpleshearfc', 'unitensiony', unitensionz',
            'SSR',
            'unicompression','unicompressionfc',
            

        This module is the key control over history file input 
        system of VPSC. The monotonic boundary conditions are
        passed through this module. All the experimentally 
        obtainable boundary conditions are imposable. Refer to 
        the above the list of modes. 

        These modes are again passed to __loading__ module. There, 
        history file for the given mode is made. VPSC class does 
        not  accepts the exisiting file. 

        ---------
        Arguments
         step = 1, # of strain increment
         eqincr = 0.005 size of the strain increment
         mode 
         ictrl
         u11
         u22
         jobid

        log:
        mode 'unicompression' is added
        """
        if mode==None:
            print " You must input a proper mode to __historyfileprocess"
            print "\navailable modes are below"
            print " 'unitension', 'unitensionfc', 'inplanebiaxial'"
            print " 'bulge', 'simpleshear', 'simpleshearfc'"
            raise IOError

        histfile = 'hist%stemp%i_%s'%(os.sep, self.nhist, jobid)
        
        self.prcs=self.prcs+self.__loading__(histfile=histfile,
                                             mode=mode,
                                             nstep=step, 
                                             eqincr=eqincr,
                                             ictrl=7,
                                             u11=u11, 
                                             u22=u22)

        FILE = open(histfile, 'r')
        lines = FILE.read()
        self.histf.append(lines)
        FILE.close()

        self.eqincr.append(eqincr)
        self.step = self.step+step
        self.nhist = self.nhist + 1

    def __lankfordprocess(self, step = 30, ):
        ####
        ####  LANKFORD coefficient
        ####
        remainder = 90.%step 
        if remainder < 0.0001 : pass
        else:
            print " Your angle increment for lankford probing",
            print " is not properly given"
        self.prcs = self.prcs + self.__lankf__(step = step)
        self.step = self.step + 90./step + 1
        pass
    
    def __pcysplprocess(self, step=72, shear_ang=90., xcomp=1, ycomp=2):
        ####
        ####  POLYCRYSTALLINE YIELD SURFACE on plane stress space
        ####
        step = step 
        shear_ang = shear_ang #angle between E12 and E11-E22 plane
        self.prcs = self.prcs + self.__pcyspl__(
            nstep=step, shear=shear_ang,
            xcomp=xcomp, ycomp=ycomp
            )
        self.step = self.step + step
        pass
    
    def __rotate(self, ang=1, fn=None):
        """ make rotate file """
        th = ang*np.pi/180.
        sth = np.sin(th) ; cth = np.cos(th)
        rmat = np.array([[cth,sth,0],[-sth,cth,0],[0,0,1]])
        if fn==None: filename = '%s%s%s%i'%(
            'rot', os.sep, 'r', self.nrot)
        else: filename = fn
        
        FILE = open(filename,'w')
        FILE.writelines('Rotation matrix for polycrystalline aggregate\n')
        for i in range(3):
            for j in range(3):
                FILE.writelines('%9.4f  '%rmat[i][j])
            FILE.writelines('\n')
        FILE.close()
        self.prcs = self.prcs + ['4', filename]
        self.nrot = self.nrot + 1
    
    def maxstep(self):
        """
        suggest best guess for process number and maximum step number
        """
        mx = 0
        for i in range(len(self.prcs)/2):
            cpr = self.prcs[i*2] 
            if cpr=='0':
                FILE = open('prcs[i*2+1]','r')
                cl = FILE.readline(); FILE.close()
                nstep = cl.split()[0] + 1 
                if nstep>mx: mx = nstep
                pass #hist
            elif cpr=='1':
                print 'process #1 is not ready '
                raise IOError
                pass #file load
            elif cpr=='2': 
                print "process #2' max step number",
                print " is fixed\n to be 108"
                if 108>mx: mx = 108
                pass #yield surface
            elif cpr=='3':
                arg = self.prcs[i*2+1]
                dang = int(arg)
                nstep = len(np.arange(0.,180.001,dang))
                if nstep>mx: mx = nstep
                pass #lankford
            elif cpr=='4': 
                if mx>1: pass
                else: mx = 2
                pass #rbr
            elif cpr=='5': 
                arg = self.prcs[i*2+1]
                n = int(arg.split()[2])
                if n>mx: mx=n
                pass #PCYS-pl
            else:
                print "self.maxstep is not ready for this"
                raise IOError
                
    def run(self,  #Np=None,  --> replaced to be len(self.prcs)/2+1
            Ni=200):
        """
        Runs the vpsc7.vpsc_deluxe_ss()
           * Now I started to have two different f2py-wrapped VPSC module.
           * One is vpsc7_g95, the other is vpsc7_gfortran
           * Which are vpsc7 compiled by g95 and gfrotran compilers, respectively
        ----
         Arguments:
          Np: # of maximum processes
          Ni: # of maximum steps
        """
        try:
            import vpsc7_gfortran as vpsc7; reload(vpsc7)
        except:
            import vpsc7_g95 as vpsc7; reload(vpsc7)
            
            
        ## number of processes
        Np = len(self.prcs)/2 + 1
        
        # vpsc7.vpsc_deluxe_ss()
        dfile = glob.glob('*.OUT') + glob.glob('*.out')
        for i in dfile: os.remove(i)

        
        Ng = 0
        for itx in self.ftx:
            Ng = Ng + self.__seekngr__(filename = itx)
        Nphase = len(self.ftx)

        ## global products initialization
        self.dbar   = np.zeros((Np,Ni,5))
        self.sbar   = np.zeros((Np,Ni,5))
        self.epstot = np.zeros((Np,Ni,3,3))
        self.scau   = np.zeros((Np,Ni,3,3))
        self.sdev   = np.zeros((Np,Ni,3,3))
        self.dsim   = np.zeros((Np,Ni,3,3))
        self.svm    = np.zeros((Np,Ni))
        self.dvm    = np.zeros((Np,Ni))
        self.epsvm  = np.zeros((Np,Ni))
        self.gr     = np.zeros((Nphase,Np,Ni,Ng,4)) #nph/proc/step/ngr/4

        ## add the initial texture to self.gr
        for iph in range(len(self.ftx)):
            cgr = np.genfromtxt(self.ftx[iph], skiprows=4)
            bgr = cgr.transpose()
            tot_wgt = bgr[3].sum()
            bgr[3] = bgr[3] / tot_wgt
            cgr = bgr.transpose()
            self.gr[iph,0,0] = cgr
        

        ## delete existing gr_0_0000_0000 files
        grfiles = glob.glob('gr_%s*'%self.jobid)
        for i in grfiles: os.remove(i)


        ################################################
        # Visco Plastic Self-Consistent tweaked solver #
        # Resulting data                               #
        #   dbar                                       #
        #   sbar                                       #
        #   epstot                                     #
        #   scau                                       #
        #   sdev                                       #
        #   dsim                                       #
        #   svm                                        #
        #   dvm                                        #
        #   epsvm                                      #
        ################################################
        ## ----------------------------------------------------------------------------------
        self.dbar,self.sbar,self.epstot,self.scau,self.sdev,self.dsim,self.svm,self.dvm,self.epsvm  = vpsc7.vpsc_deluxe_ss(
            self.dbar,self.sbar,self.epstot,self.scau,self.sdev,self.dsim,self.svm,self.dvm,self.epsvm, self.jobid
            )
        ## ----------------------------------------------------------------------------------

        ## Reads gr_0_0000_0000 files
        grfiles = glob.glob('gr_%s*'%self.jobid)
        for i in range(len(grfiles)):
            nph = int(grfiles[i].split('_')[2])
            ipr = int(grfiles[i].split('_')[3])
            ist = int(grfiles[i].split('_')[4])
            try: cgr = np.genfromtxt(grfiles[i])#current poly-xtal aggregate
            except:
                print grfiles[i]
                print i
                raise IOError, "could not find the file"
            try: self.gr[
                nph-1][ipr][ist]=cgr  #nph and ist start from 1...
            except: pass
            else: pass
            
            pass
        
        ## Deletes crss_*.out
        for i in grfiles: os.remove(i) # del?
        self.datamaster['lgr'] = cgr #save the last texture

        ## Reads modulus_0 files
        modulus_file_name = 'modulus_%s'%self.jobid
        FILE = open(modulus_file_name); FILE.close()
        try: moduli = np.loadtxt(
            modulus_file_name, skiprows=2).transpose()
        except:
            print '%s is not found.'%modulus_file_name
            os.listdir(os.getcwd())
            #raw_input()
        
        else:
            "moduli[0]--> deformation step "
            self.datamaster['ssc11'] = moduli[1]
            self.datamaster['ssc22'] = moduli[2]
            self.datamaster['ssc33'] = moduli[3]
            os.remove(modulus_file_name) #deletes
            pass

        """  #old style of returning variables (2011-02-10) """
        self.pp()
        return self.datamaster#self.sbar,self.dbar


    ### POLYCRYSTALLINE AGGREGATE'S CONDITIONS ###
    "- PHASE INFORMATION"
    "- SINGLE CRYSTAL INFO"
    "- FILEAXES(GRAIN MORPHOLOGY) "

    ### MODELLING CONDITIONS FOR THE RUN ###


    ### RUNNING PROCESSES
    
    #--------------   VPSC7.in    ----------------
    #  using v_in
    #  import vpsc_in.py   v_in = vpsc_in.vpsc_in
    #---------------------------------------------

    def __vpsc_in__(self,
                    nelem=1, nph=None, wph=None, 
                    ishape=None, fragmentn=None, crit=None,
                    tfile=None, sxfile=None, fileaxes=None,
                    irecover=0, isave=0, icubcomp=0, nwrite=0, ihardlaw=0,
                    iflu=0, iratesens=0, interaction=None, iupdate=None,
                    nneigh = 0, npro=None, prcs=[],
                    ### precision settings for convergence prodcedure
                    err=0.001, errd=0.001, errm=0.001, errso=0.01,
                    nmxiter=100, exiter=100, initer=100,
                    irsvar=0, xrsini=2, xrsfin=10, xrstep=2,
                    ####----------------------------------------------
                    ellipangl = None, #[ [0.,0.,0.],[0.,0.,0.]],
                    ellipaxes = None, #[ [1.,1.,1.],[1.,1.,1.]]
                    ):

        """ vpsc input file manipulator """
        vpscinfilename = 'VPSC7_%s.in'%self.jobid
        v_in(
            nelem=nelem, nph=nph, wph =wph, ishape=ishape,
            vin_name = vpscinfilename,
            tfile= tfile,
            sxfile=sxfile,
            fileaxes=fileaxes,

            #Characterizing the ellipsoid
            ellipangl = ellipangl[0:nph],
            ellipaxes = ellipaxes[0:nph],

            ### precision settings for convergence prodcedure
            err=err, errd=errd, errm=errm, errso=errso,
            nmxiter=nmxiter, exiter=exiter, initer=initer,
            irsvar=irsvar, xrsini=xrsini, xrsfin=xrsfin, xrstep=xrstep,
            ####----------------------------------------------

            #Controls 
            irecover = irecover, isave=isave,
            icubcomp = icubcomp, nwrite = nwrite, ihardlaw =ihardlaw, 
            
            iflu = iflu,

            iratesens = iratesens, #rate insensitive
            interaction = interaction, #(0:FC,1:aff,2:sec,3:neff=10,4:tngt,5:SO)
            iupdate = iupdate, #(orient, grain shape, hardening, itran)
            nneigh = nneigh, #(0 for no neighbors, 1 for pairs, etc.)
            npro = npro,
            prcs = prcs    #processes
            )

    ### SAVE DATA FROM FILES ##
    def _db_master_(self):
        ## HISTFILE
        self.datamaster['dbar']=[]
        self.datamaster['sbar']=[]
        self.datamaster['eps']=[]
        self.datamaster['scau']=[]
        self.datamaster['sdev']=[]
        self.datamaster['dsim']=[]
        self.datamaster['svm'] = []
        self.datamaster['dvm'] = []
        self.datamaster['evm'] = []

        ## LANKF
        self.datamaster['YSprob']=[]
        self.datamaster['Rprob']=[]
        self.datamaster['ang']=[]

        ## orientation
        self.datamaster['amat']=[] #rotation matrix
        
    def pp(self):
        """
        Based on self.prcs, do what is necessary to do
        """
        self.amat = eul(ph=0,th=0,tm=0,echo=False)
        ##
        # Strategy
        # 1. Loop over prcs.
        # 2. With the each loop's process id number (flag)
        #    and its argument passed, run a proper sub-pp def.
        # 3. Establishes self.dbar, self.sbar which is accessible later on.
        # 4. Use chg_basis to map the 5-D variable into 6 components-variable
        self._db_master_()  #initialize the dictionary in which all keys are included
        nprc = 0

        for i in range(len(self.prcs)/2):
            #self.prcs(i*2)   --> process id
            #self.prcs(i*2+1) --> corresponding arguments

            if self.prcs[i*2]=='0':
                ## history file loading
                self.__ppHIST__(
                    arg=self.prcs[i*2+1],
                    dbar=self.dbar[nprc],
                    sbar=self.sbar[nprc],
                    epstot=self.epstot[nprc],
                    scau=self.scau[nprc],
                    sdev=self.sdev[nprc],
                    svm=self.svm[nprc],
                    dvm=self.dvm[nprc],
                    evm=self.epsvm[nprc],
                    n=nprc)

            elif self.prcs[i*2]=='1':
                ## subroutine loading
                self.__ppLOADING__(arg=self.prcs[i*2+1],
                                   v1=self.dbar[nprc],
                                   v2=self.sbar[nprc],
                                   n=nprc)
                print " PP for subroutine loading is not ready"

            elif self.prcs[i*2]=='2':
                ## pcys 
                self.__ppPCYS__(arg=self.prcs[i*2+1],
                                v1=self.dbar[nprc],
                                v2=self.sbar[nprc],
                                n=nprc)

            elif self.prcs[i*2]=='3':
                ## lankford probing
                """
                It does not have neither scau nor sdev written
                """
                self.__ppLANK__(arg=self.prcs[i*2+1],
                                dbar=self.dbar[nprc],
                                sbar=self.sbar[nprc],
                                scau=self.scau[nprc],
                                dsim=self.dsim[nprc],
                                n=nprc)

            elif self.prcs[i*2]=='4':
                self.__ppRBR__(arg=self.prcs[i*2+1])
                ## rotation of the polycrystal
                ## should be tracked in self.amat
                ## as a rotation matrix
                """
                arg=self.prcs[i*2+1],
                v1=self.dbar[nprc],
                v2=self.sbar[nprc],
                n=nprc)
                """
                
            elif self.prcs[i*2]=='5':
                self.__ppPCYSPL__(arg=self.prcs[i*2+1],
                                  dbar=self.dbar[nprc],
                                  sbar=self.sbar[nprc],
                                  n=nprc)
            nprc = nprc + 1
        self.datamaster['eul'] = []
        for i in self.datamaster['amat']:
            self.datamaster['eul'].append(eul(a=i,echo=False))
        pass

    def __ppHIST__(self,arg,dbar,sbar,epstot,scau,sdev,svm,dvm,evm,n):
        """
        post-process history file
        Open the argument
        dbar, sbar, epstot, scau, sdev, n
        """
        #c4 = chg_basis(ce2,iopt=1)
        nstep, eqc = self.__ppHIST2__(arg=arg) #returning the # step and eqc
        #print 'nstep=',nstep
        dbar = dbar[0:nstep+1]
        sbar = sbar[0:nstep+1]
        epstot = epstot[0:nstep]
        scau = scau[0:nstep+1]
        sdev = sdev[0:nstep+1]
        svm  = svm[0:nstep+1]
        dvm  = dvm[0:nstep+1]
        evm  = evm[0:nstep+1]
        
        V1 = []; V2 = []; V3 = []; V4 = []; V5 = []; V6 = []; V7 = []; V8 = []
        
        for i in range(len(dbar)): V1.append(chg_basis(ce2=dbar[i], iopt=1))
        for i in range(len(sbar)): V2.append(chg_basis(ce2=sbar[i], iopt=1))
        for i in range(len(epstot)): V3.append(epstot[i]) #eps is in the 3x3 components
        for i in range(len(scau)): V4.append(scau[i])
        for i in range(len(sdev)): V5.append(sdev[i])
        for i in range(len(svm)): V6.append(svm[i])
        for i in range(len(dvm)): V7.append(dvm[i])
        for i in range(len(evm)): V8.append(evm[i])

        ## plastic work calculation
        """
        It is believed to be no good to provide plastic work in vp_f2py.
        Instead, user is advised to do so in each front-end modules.
        """

        # make lists as np array types.
        V1, V2, V3, V4, V5, V6, V7, V8 = __makenp__(V1,V2,V3,V4,V5,V6,V7,V8)

        # print 'V3', V3
        # raw_input()
        V3 = np.append([np.zeros((3,3))], V3, axis=0)

        self.datamaster['dbar'].append(V1)
        self.datamaster['sbar'].append(V2)
        self.datamaster['eps'].append(V3)
        self.datamaster['scau'].append(V4)
        self.datamaster['sdev'].append(V5)
        self.datamaster['svm'].append(V6)
        self.datamaster['dvm'].append(V7)
        self.datamaster['evm'].append(V8)
        
        self.datamaster['dsim'].append(None)
        self.datamaster['YSprob'].append(None)
        self.datamaster['Rprob'].append(None)
        self.datamaster['ang'].append(None)
        ## ORIENTATION 
        self.datamaster['amat'].append(self.amat)
        pass

    def __ppHIST2__(self,arg):
        """ Returns the nstep of the history file"""
        FILE = open(arg,'r')
        line = FILE.readline(); FILE.close()
        return int(line.split()[0]), int(line.split()[1])

    def __ppLOADING__(self,arg,v1,v2,n):
        """
        post-process history subroutine
        """
        ## HISTFILE
        self.datamaster['dbar'].append(None)
        self.datamaster['sbar'].append(None)
        self.datamaster['eps'].append(None)
        self.datamaster['scau'].append(None)
        self.datamaster['sdev'].append(None)
        self.datamaster['dsim'].append(None)
        ## LANKF
        self.datamaster['YSprob'].append(None)
        self.datamaster['Rprob'].append(None)
        self.datamaster['ang'].append(None)
        ## ORIENTATION 
        self.datamaster['amat'].append(self.amat)
        return None

    def __ppPCYS__(self,arg,v1,v2,n):
        """
        post-process pcys 
        """
        ## HISTFILE
        self.datamaster['dbar'].append(None)
        self.datamaster['sbar'].append(None)
        self.datamaster['eps'].append(None)
        self.datamaster['scau'].append(None)
        self.datamaster['sdev'].append(None)
        self.datamaster['dsim'].append(None)
        ## LANKF        
        self.datamaster['YSprob'].append(None)
        self.datamaster['Rprob'].append(None)
        self.datamaster['ang'].append(None)
        ## ORIENTATION 
        self.datamaster['amat'].append(self.amat)
        return None

    def __ppLANK__(self,arg,dbar,sbar,scau,dsim,n):
        """
        post-process lankford
        """
        dang = arg
        if int(float(dang)) - float(dang)==0.0:
            dang = float(dang)
            nstp = int(90./dang+1)
            ang = np.arange(0.,90.0001,dang)
        else: 
            print "Check the angular increment ",
            print "for LANKF for #%i-th prcs"%n
            print "prcs id starting from 0"
        e = []; s = []; scauchy = []; d = []
        r = []
        ys = []
        for i in range(nstp):
            e.append(chg_basis(ce2=dbar[i], iopt=1))
            s.append(chg_basis(ce2=sbar[i], iopt=1))
            scauchy.append(scau[i])
            d.append(dsim[i])
            r.append(dsim[i][1][1] / dsim[i][2][2])
            ys.append(scau[i][0,0])

        #HISTFILE
        self.datamaster['dbar'].append(e)
        self.datamaster['sbar'].append(s)
        self.datamaster['scau'].append(scauchy)
        self.datamaster['dsim'].append(d)
        self.datamaster['eps'].append(None)
        self.datamaster['sdev'].append(None)
        ##LANKF
        self.datamaster['YSprob'].append(ys)
        self.datamaster['Rprob'].append(r)
        self.datamaster['ang'].append(ang)
        ## ORIENTATION 
        self.datamaster['amat'].append(self.amat)        
        return None

    def __ppRBR__(self,arg):
        """
        post-process RBR
        Global strain and stress fields are fetched to be zero.
        However, at the orientation with respect to laboratory
        must be saved.
        """
        FILE = open(arg,'r')
        lines = FILE.readlines()
        lines = lines[1:4]
        a = np.zeros((3,3))
        a[0] = map(float, lines[0].split())
        a[1] = map(float, lines[1].split())
        a[2] = map(float, lines[2].split())
        # Rotation matrix updates
        eul(a=a,echo=False)
        self.amat = np.dot(self.amat,a)

        ## HISTFILE
        self.datamaster['dbar'].append(None)
        self.datamaster['sbar'].append(None)
        self.datamaster['eps'].append(None)
        self.datamaster['scau'].append(None)
        self.datamaster['sdev'].append(None)
        self.datamaster['dsim'].append(None)
        ## LANKF
        self.datamaster['YSprob'].append(None)
        self.datamaster['Rprob'].append(None)
        self.datamaster['ang'].append(None)
        ## ORIENTATION 
        self.datamaster['amat'].append(self.amat)
        return None
    
    def __ppPCYSPL__(self,arg,dbar,sbar,n):
        """ Post-process polycrystal yield surface

        Refer to the deviatoric space that VPSC uses.
        1 = (s22 - s11)/sqrt(2)
        2 = (2*s33-s22-s11)/sqrt(6)
        3 = s23*sqrt(2)
        4 = s13*sqrt(2)
        5 = s12*sqrt(2)
        """
        ix, iy, nstep, dummy = map(int, arg.split())
        ang_shear = float(arg.split()[-1])

        if ix==1 and iy==2 :
            pass
        elif ix==1 and iy==4:
            #s13 - (s22-s11)/sqrt(2)
            pass
        elif ix==1 and iy==5:
            #S12 - (
            pass
        else:
            print "Unexpected argument for axis of pcyspl projection"
            print "Default x and y axis are 1 and 2 "
            raise IOError
        
        s = []; d = []
        for i in range(nstep):
            sdum = chg_basis(ce2=sbar[i], iopt=1)
            ddum = chg_basis(ce2=dbar[i], iopt=1)
            for j in range(3):
                sdum[j][j] = sdum[j][j] - sdum[2][2]
                ddum[j][j] = ddum[j][j] #- ddum[2][2] # rather controversial...
            s.append(sdum);d.append(ddum)
            pass

        ## HISTFILE
        self.datamaster['dbar'].append(d)
        self.datamaster['sbar'].append(s)
        self.datamaster['eps'].append(None)
        self.datamaster['scau'].append(None)
        self.datamaster['sdev'].append(None)
        self.datamaster['dsim'].append(None)
        ## LANKF
        self.datamaster['YSprob'].append(None)
        self.datamaster['Rprob'].append(None)
        self.datamaster['ang'].append(None)
        ## ORIENTATION 
        self.datamaster['amat'].append(self.amat)
        return None

    #----------------------------
    ###
    ### prcs MAKERS 
    ###
    #----------------------------

    def __pcyspl__(self,nstep=None,shear=None,xcomp=1,ycomp=2):
        """
        subroutine pcys_pl
        
        nstep=None
        shear=None : angle from the shear axis.
        """
        if os.name=='nt': os.system('cls')
        elif os.name=='posix': os.system('clear')
        print " *** POLYCRYSTALLINE YIELD SURFACE",
        print " ON PLANE STRESS SPACE *** "
        pcrs = ['5','%i %i %i %i'%(xcomp, ycomp, nstep,shear)]
        return pcrs

    def __pcys_tau__(self, nstep, ix, iy):
        """
        subroutine pcys_pl to draw a shear axis included
        yield locus
        2011-Aug-18
        """
       
        return

    def __lankf__(self,step):
        """
        lankford prcs returning definition
        """
        if os.name=='nt': os.system('cls')
        elif os.name=='posix': os.system('clear')
        print "*** LANKFORD PROBING ***"
        print "angle increment : %3i"%(step)
        prcs = ['3', str(step)]
        return prcs
        
    def __loading__(self, histfile=None, mode= 'unitension', 
                    nstep=1, eqincr=0.005,
                    ictrl=7,
                    u11=None,u22=None):
        """ 
        uniaxial or multiaxial loadings  (uniaxial compression: 2011-07)
        returns prcs : process list
        available mode : unitension, unitensionfc, bulge, 
                         inplanebixial, simpleshear, SSR,
                         SSRFC,
                         unicompression, unicompressionfc
        unitension : prcs = ['0',histfile]
        bitension : prcs = ['0', histfile]
        bb : prcs = ['0', histfile]

        Since the process flag, e.g. '0' for loading, decided to be
        integer not the string. During post-execution process,
        which processes been through is detected by this flag. To 
        tell this integer flag from other number in 'self.prc', it
        is better to have a different type. Since this module is not
        the best place to perform, it is done after execution.
        Change in flags are done in the 'self.pp()' module.
        """
        if histfile==None:
            print "histfile is hardwired to be 'hist%stemp'%os.sep"
        elif histfile!=None: pass
        prcs = ['0', histfile]
        if mode=='unitension':
            #---- uniaxial tension along 1-axis (SC) ----#
            iudot=    [[ 1,   0,   0 ],
                       [ 1,   0,   0 ],
                       [ 1,   1,   0 ]]
            udot =    [[ 1 , 0  , 0   ],
                       [ 0 ,-0.5, 0   ],
                       [ 0 , 0  ,-0.5 ]]
            iscau=    [0,1,1,1,1,1]
            scauchy = [0,0,0,0,0,0]
        elif mode =='unicompression':
            #---- uniaxial compression along 1-axis (SC) ----#
            iudot=    [[ 1,   0,   0 ],
                       [ 1,   0,   0 ],
                       [ 1,   1,   0 ]]
            udot =    [[-1 , 0  , 0   ],
                       [ 0 , 0.5, 0   ],
                       [ 0 , 0  , 0.5 ]]
            iscau=    [0,1,1,1,1,1]
            scauchy = [0,0,0,0,0,0]            
        elif mode=='unitensiony': 
            iudot =   [[ 0,   0,   0],
                       [ 1,   1,   0],
                       [ 1,   1,   0]]
            udot  =   [[-0.5,   0,   0],
                       [0   ,   1,   0],
                       [0   ,   0,-0.5]]
            iscau =   [1,0,1,1,1,1]
            scauchy = [0,0,0,0,0,0]
        elif mode=='unitensionz': 
            iudot =   [[ 0,   0,   0],
                       [ 1,   0,   0],
                       [ 1,   1,   1]]
            udot  =   [[-0.5,   0,   0],
                       [0   ,-0.5,   0],
                       [0   ,   0,   1]]
            iscau =   [1,1,0,1,1,1]
            scauchy = [0,0,0,0,0,0]

        elif mode=='unitensionfc':
            #---- uniaxial tension along 1-axis (FC) ----#
            iudot=    [[ 1,   1,   1 ],
                       [ 1,   1,   1 ],
                       [ 1,   1,   1 ]]
            udot =    [[ 1 , 0  , 0   ],
                       [ 0 ,-0.5, 0   ],
                       [ 0 , 0  ,-0.5 ]]
            iscau=    [0,0,0,0,0,0]
            scauchy = [0,0,0,0,0,0]
        elif mode=='unicompressionfc':
            #---- uniaxial compression along 1-axis (FC) ----#
            iudot=    [[ 1,   1,   1 ],
                       [ 1,   1,   1 ],
                       [ 1,   1,   1 ]]
            udot =    [[-1 , 0  , 0   ],
                       [ 0 , 0.5, 0   ],
                       [ 0 , 0  , 0.5 ]]
            iscau=    [0,0,0,0,0,0]
            scauchy = [0,0,0,0,0,0]            
        elif mode=='inplanebiaxial':
            #---- in-plane biaxial tension ----#
            u11 = round(u11,6); u22 = round(u22,6)
            u33 = -u11-u22; u33 = round(u33,6)
            if u33 + u11 + u22 !=0: 
                """ 
                This may be an indication of roundoff error,
                mostly not after modification of the code slightly
                """
                #print " alert, sum of diagoanl terms is not zero"
                #print " u11,u22,u33 = %f %f %f"%(u11,u22,u33)
                #raw_input("press enter to proceed anyway >> ")
            u = [u11,u22]; 
            if any(u[i]==None for i in range(2)):
                print "for inplanebiaxial mode, one has to give",
                print "u11 and u22"
                raise IOError
            iudot =   [[ 1,  0,  0 ],
                       [ 1,  1,  0 ],
                       [ 1,  1,  0 ]]
            udot =    [[u11, 0,  0 ],
                       [ 0, u22, 0 ],
                       [0,   0, -u11-u22]]
            iscau=     [0,0,1,1,1,1]
            scauchy =  [0,0,0,0,0,0]

        elif mode=='inplanebiaxialfc':
            #---- in-plane biaxial tension (full constraint) ----#
            u11 = round(u11,6); u22= round(u22,6)
            u33 = -u11-u22; u33=round(u33,6)
            if u33 + u11 + u22 !=0: 
                """ 
                This may be an indication of roundoff error,
                mostly not after modification of the code slightly
                """
                #print " alert, sum of diagoanl terms is not zero"
                #print " u11,u22,u33 = %f %f %f"%(u11,u22,u33)
                #raw_input("press enter to proceed anyway >> ")
            u = [u11,u22]; 
            if any(u[i]==None for i in range(2)):
                print "for inplanebiaxial mode, one has to give",
                print "u11 and u22"
                raise IOError
            iudot =   [[ 1,  1,  1 ],
                       [ 1,  1,  1 ],
                       [ 1,  1,  1 ]]
            udot =    [[u11, 0,  0 ],
                       [ 0, u22, 0 ],
                       [0,   0, -u11-u22]]
            iscau =    [0,0,0,0,0,0]
            scauchy =  [0,0,0,0,0,0]            

        elif mode=='bulge': 
            #---- disc-compression ----#
            iudot =   [[ 0,   0,   0 ],
                       [ 1,   0,   0 ],
                       [ 1,   1,   1 ]]
            udot =    [[ 0.5, 0  ,  0 ],
                       [ 0  , 0.5,  0 ],
                       [ 0  , 0  , -1.0]]
            iscau=     [1,1,0,1,1,1]
            scauchy =  [0,0,0,0,0,0]

        elif mode=='simpleshear':  
            #---- simpleshear mode ----#
            iudot =   [[ 0,   1,   1 ],
                       [ 1,   0,   1 ],
                       [ 1,   1,   0 ]]
            udot =    [[ 0  , 2  ,  0 ],
                       [ 0  , 0  ,  0 ],
                       [ 0  , 0  ,  0 ]]
            iscau=     [1,1,1,0,0,0]
            scauchy =  [0,0,0,0,0,0]
            
        elif mode=='SSR' or mode=='simpleshearreverse':
            #---- Reverse simpleshear mode ----#
            iudot =   [[ 0,   1,   1 ],
                       [ 1,   0,   1 ],
                       [ 1,   1,   0 ]]            
            udot =    [[ 0  ,-2  ,  0 ],
                       [ 0  , 0  ,  0 ],
                       [ 0  , 0  ,  0 ]]
            iscau=     [1,1,1,0,0,0]
            scauchy =  [0,0,0,0,0,0]
            pass
        elif mode=='SSRFC':
            #---- Reverse simpleshear FC mode ----#
            iudot =   [[ 1,   1,   1 ],
                       [ 1,   1,   1 ],
                       [ 1,   1,   1 ]]            
            udot =    [[ 0  ,-2  ,  0 ],
                       [ 0  , 0  ,  0 ],
                       [ 0  , 0  ,  0 ]]
            iscau=     [0,0,0,0,0,0]
            scauchy =  [0,0,0,0,0,0]
            pass
                      
        elif mode=='simpleshearfc':
            #---- simpleshear fc mode: all zeroes for diagonal terms ----#
            iudot =   [[ 1,   1,   1 ],
                       [ 1,   1,   1 ],
                       [ 1,   1,   1 ]]
            udot =    [[ 0  , 2  ,  0 ],
                       [ 0  , 0  ,  0 ],
                       [ 0  , 0  ,  0 ]]
            iscau=     [0,0,0,0,0,0]
            scauchy =  [0,0,0,0,0,0]

        else:
            print 'Unexpected mode given to __hist__ method.'
            print 'Given mode is', mode
            raise IOError

        ## hist file making
        vpsc_in.histmaker(
            filename = histfile,
            nstep = nstep, eqincr = eqincr, ictrl = ictrl,
            iudot = iudot, udot = udot, iscau = iscau,
            scauchy = scauchy)
        return prcs

    def __showfile__(self,filename):
        sys = os.name
        if sys=='nt':
            cmd = 'type'
            try: os.system('%s %s%s%s'%(cmd, self.cwd,os.sep,filename))
            except: 
                print "Could not open the file"
                print "Tried command is as follow"
                print '%s %s%s%s'%(cmd, self.cwd,os.sep,filename)
                return -1
        if sys=='posix':
            cmd = 'cat'
            try:
                os.system('%s .%s%s%s'%(cmd, self.cwd.split(os.getcwd())[1],
                                        os.sep,filename))
            except:
                print "Could not open the file"
                print "Tried command is as follow"
                print '%s .%s%s%s'%(cmd, self.cwd.split(os.getcwd())[1],
                                        os.sep,filename)
                return -1
        pass

    def show_lankford(self):
        """  """
        self.__showfile__('LANKFORD.OUT')
        pass

    def show_pcys(self):
        """  """
        self.__showfile__('PCYS.OUT')
        pass

    def show_pcyspl(self):
        """  """
        self.__showfile__('pcys_pl.out')
        pass
    
    def show_vpscin(self):
        """
        Shows 'VPSC7.in'
        """
        if os.name=='nt': 
            cmd = 'type'
            os.system('%s %s%s%s'%(cmd, self.cwd, os.sep,
                                   'VPSC7_%s.in'%self.jobid
                                   ))
        elif os.name=='posix':
            cmd ='cat'
            os.system('%s .%s%s%s'%(cmd, 
                                    self.cwd.split(os.getcwd())[1],
                                    os.sep,
                                    'VPSC7_%s.in'%self.jobid
                                    ))

    def show_sx(self):
        """
        Shows single crystal files
        """
        if os.name=='nt': cmd = 'type'
        elif os.name=='posix': cmd = 'cat'
        for i in range(self.VPSCin.nph):
            print "\n"
            print " **********************"
            print " single crystal file #%1i"%i
            print " **********************\n"
            if os.name=='nt':
                os.system("%s %s%s%s"%(cmd,self.cwd,os.sep,self.fsx[i]))
            elif os.name=='posix':
                os.system('%s .%s%s%s'%(cmd, self.cwd.split(os.getcwd())[1],
                                        os.sep,self.fsx[i]))                
            raw_input('\n Press enter ')

    def show_hist(self):
        """
        Shows history file(s)
        """
        if os.name=='nt': cmd='type'
        elif os.name=='posix': cmd='cat'
        nhist = 0
        for i in range(len(self.prcs)/2): 
            if self.prcs[i*2]=='0':
                if i>0 :  raw_input('\npress enter >>');print'\n'
                print " \n ***************************"
                print " LOADING FILES (histfile) #%i"%nhist
                print " ***************************\n"
                if os.name=='nt':
                    os.system("%s %s%s%s%s%s%i"%(cmd, self.cwd,os.sep,
                                                 'hist',os.sep, 
                                                 'temp',nhist))
                elif os.name=='posix':
                    os.system('%s .%s%s%s%s%s%i'%(cmd, self.cwd.split(os.getcwd())[1],
                                            os.sep,'hist',os.sep,'temp',nhist))
                    """
                    os.system('%s .%s%s%s'%(cmd, self.cwd.split(os.getcwd())[1],
                                            os.sep,self.fsx[i]))                
                    os.system("%s .%s%s%s%s%s%i"%(cmd, self.cwd, os.sep,
                                                 'hist',os.sep, 
                                                 'temp',nhist))
                    """
                nhist = nhist + 1
            pass

    ### SOME HANDY MODULES EXPECTED ###
    def __seekngr__(self, filename):
        """ finds ngrain in the textfile """
        FILE = open(filename, 'r')
        lines = FILE.readlines()[3]
        return int(lines.split()[1]) #number of grain
        #self.ngrain = self.__seekngr__(self.ftex)

    
    ### preliminaries ###
    " default file making"
    " ODF management "  #done (related libraries: cod_section, cmb, )
    " random grains "   #done (cmb.random)
    " RBR file making " #done
    
    ### to be continued upon issues like ... ###
    " POSTMORT management " #somehow, it's been done
    " fileaxes maker (shape file)" #done



