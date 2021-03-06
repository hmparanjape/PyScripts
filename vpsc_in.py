
import os
import numpy as np
try: plt
except: 
    try: 
        import matplotlib.pyplot as plt
    except:
        print " no plotting is possible under the current environment"

def vpsc_in(
    ## # elem, ph and vol fraction of phases.
    ## and ishape for flag for morphology updates
    nelem=1, nph=1, wph =[1.0, 0.0], ishape=[0,],

    #framgmentation options and critical value limit 
    # in distortion of the inclusion 
    # as well as ellipsoid description: orientation and 
    # principal lengths of'em
    fragmentn = [0,0], crit = [25,25],
    ellipangl = [[0., 0., 0.],[0.,0.,0.]],
    ellipaxes = [[1., 1., 1.],[1.,1.,1.]],
    
    # vpsc7.in file name (fixed)
    # texture, singlecrystal and shape file names in list format
    vin_name = 'VPSC7.in',
    tfile= ['texture/500.tex',],
    sxfile=['sx/hijhiihb.sx',],
    fileaxes=['shape1.100', 'shape1.100'],

    ###  precision settings for convergence prodcedure
    ##
    err=0.001, errd=0.001, errm=0.001, errso=0.01,
    nmxiter=100, exiter=100, initer=100,
    irsvar=0, xrsini=2, xrsfin=10, xrstep=2,

    ### post-mortem flags, texture component calc. flag.
    ## 
    irecover = 0, isave=0, icubcomp = 0, nwrite = 0, 

    ## Single crystal hardening rule (-2: dislocation, 0:Voce, 1: MTS, 2: composite)
    ihardlaw = 0, 

    ## fluctuation flag
    iflu = 0,

    ##  rate insensitive, interaction :(FC,aff,2:sec,3:neff,4:tan,5:SO)
    ##  update flags for orient, grain shape, hardening, itran(youngung)
    ##  nneigh (0 for no neighbors, 1 for pairs, etc.)
    iratesens = 0, interaction = 3, iupdate = [1,1,1,0], nneigh = 0, 

    ## proce number and processes
    npro = None,
    prcs = [],  #processes
    
    ishow=False
    ):
    """
    VPSC INPUT FILE MANIPULATOR
    Creats a vpsc code input file (a.k.a 'VPSC7.in')

    ---------
    ARGUMENTS:

    ## # elem, ph and vol fraction of phases.
    ## and ishape for flag for morphology updates
    nelem=1, nph=1, wph =[1.0, 0.0], ishape=[0,],

    #framgmentation options and critical value limit 
    # in distortion of the inclusion 
    # as well as ellipsoid description: orientation and 
    # principal lengths of'em
    fragmentn = [0,0], crit = [25,25],
    ellipangl = [[0., 0., 0.],[0.,0.,0.]],
    ellipaxes = [[1., 1., 1.],[1.,1.,1.]],
    
    # vpsc7.in file name (fixed)
    # texture, singlecrystal and shape file names in list format
    vin_name = 'VPSC7.in',
    tfile= ['texture/500.tex',],
    sxfile=['sx/hijhiihb.sx',],
    fileaxes=['shape1.100', 'shape1.100'],

    ###  precision settings for convergence prodcedure
    ##
    err=0.001, errd=0.001, errm=0.001, errso=0.01,
    nmxiter=100, exiter=100, initer=100,
    irsvar=0, xrsini=2, xrsfin=10, xrstep=2,

    ### post-mortem flags, texture component calc. flag.
    ## 
    irecover = 0, isave=0, icubcomp = 0, nwrite = 0, 

    ## Single crystal hardening rule (0:Voce, 1: MTS, 2: composite)
    ihardlaw =0, 

    ## fluctuation flag
    iflu = 0,

    ##  rate insensitive, interaction :(FC,aff,2:sec,3:neff,4:tan,5:SO)
    ##  update flags for orient, grain shape, hardening, itran(youngung)
    ##  nneigh (0 for no neighbors, 1 for pairs, etc.)
    iratesens = 0, interaction = 3, iupdate = [1,1,1,0], nneigh = 0, 

    ## proce number and processes
    npro = None,
    prcs = []  #processes


    """
    try: os
    except: import os
    else: pass
    cwd = os.getcwd()

    # READ EXISITING 'vpsc7.in'
    if os.path.isfile(vin_name): #if there already exists the file!
        f_vin = file(vin_name, 'r')
        lines = f_vin.read()
        lines = lines.split('\n')
        f_temp = file('%s.tmp'%vin_name, 'w')
        for i in lines:
            f_temp.writelines(i)
            f_temp.writelines('\n')
            pass
        f_temp.close()
        f_vin.close()
    """
    # DEL EXISITING 'vpsc7.in'
    f_vin.close()
    os.remove(vin_name)
    """
    # MAKE NEW 'vpsc7.in'
    f_vin = file(vin_name, 'w')
    f_vin.writelines('%i   nelem\n'%(nelem))
    f_vin.writelines('%i   nph\n'%(nph))
    for i in range(nph):
        f_vin.writelines('%f  '%(wph[i]))
    f_vin.writelines('  wph[i]\n')
    for i in range(nph): 
        f_vin.writelines('***  INFORMATION ABOUT PHASE #%i\n'%(i+1))
        f_vin.writelines('%i %i %i'%(ishape[i], fragmentn[i], crit[i]))
        f_vin.writelines('        ishape, fragmentn, crit aspect\n')
        f_vin.writelines('%f %f %f'%(ellipaxes[i][0], 
                                     ellipaxes[i][1], 
                                     ellipaxes[i][2]))
        f_vin.writelines('        initial ellipsoid ratios (dummy if ishape=4)\n')
        f_vin.writelines('%f %f %f'%(ellipangl[i][0], 
                                     ellipangl[i][1], 
                                     ellipangl[i][2]))
        f_vin.writelines('        init Eul ang ellips axes(dummy if ishape=3,4)\n')
        # WRITES TEXTURE FILE
        f_vin.writelines('-----------   filetext\n')
        f_vin.writelines(tfile[i]+'\n')
        f_vin.writelines('-----------   sxfile\n')
        # WRITES SINGLE CRYSTAL FILE
        f_vin.writelines(sxfile[i]+'\n')
        # WRITES SHAPE FILE
        f_vin.writelines('-----------   fileaxes (dummy if ishape=0) \n')
        f_vin.writelines('%s \n'%(fileaxes[i]))


    #### precision control block  ###
    f_vin.writelines('***  PRECISION SETTINGS FOR')
    f_vin.writelines(' CONVERGENCE PROCEDURES (default values)\n')
    f_vin.writelines('%f %f %f %f'%(err, #err
                                    errd, 
                                    errm,
                                    errso))
    f_vin.writelines('       errs,errd,errm,errso\n')
    f_vin.writelines('%i %i %i %25s'%( nmxiter, exiter, initer,' '))
    f_vin.writelines('itmax: max # of iter, external, internal and SO loops\n')
    f_vin.writelines('%i %i %i %i %30s'%(irsvar, xrsini, xrsfin, xrstep,' '))
    f_vin.writelines('irsvar & xrsini, xrsfin, xrstep (dummy if irsvar = 0) \n')
    f_vin.writelines("%i       ibcinv (0: don't use <Bc>**-1, 1: use \n"%(1))
    f_vin.writelines('***  INPUT/OUTPUT SETTINGS FOR THE RUN(default is zero)\n')
    f_vin.writelines('%i       irecover:read grain stats from POSTMORT.IN(1) OR NOT(0)\n'%(irecover))
    f_vin.writelines("%i       isave: write grain states in postmor.out at step 'isave'?\n"%(isave))
    f_vin.writelines("%i       icubcompL calcuate fcc rolling components?\n"%(icubcomp))
    f_vin.writelines("%i       nwrite (frequency of texture downloads)\n"%(nwrite))
    f_vin.writelines("***  MODELING CONDITIONS FOR THE RUN \n")
    f_vin.writelines("%i        ihardlaw(0:VOCE, 1: MTS, 2: composite grain\n"%(ihardlaw))
    f_vin.writelines('%1i' %(iratesens))
    f_vin.writelines('         iratesens (0:rate insensitive, 1:rate sensitive) \n')
    f_vin.writelines('%1i' %(interaction))
    f_vin.writelines('         interaction (0:FC,1:affine,2:secant,3:neff=10,4:tangent,5:SO) \n')
    f_vin.writelines('%1i%3i%3i%3i'%(iupdate[0],iupdate[1],iupdate[2],iupdate[3]))
    f_vin.writelines('        iupdate: update ori, grain shape, hardening, itran \n')
    f_vin.writelines('%1i'%(nneigh))
    f_vin.writelines('        nneigh (0 for no neighbors, 1 for pairs, etc.)\n')
    f_vin.writelines('%i'%(iflu))
    f_vin.writelines("        iflu(0:don't calc, 1: calc fluctuations\n")

    f_vin.writelines("*NUMBER OF PROCESSES\n")
    f_vin.writelines("%i  \n"%(npro))
    f_vin.writelines("**  (0,1: loadings, 2:pcys, 3:lankf, 4: RBR,")
    f_vin.writelines("5:pcys_pl, 6: pysc_pl 6d(ix,iy, nprob, ang")
    f_vin.writelines("( or if negative behave like pcys))) \n")
    for i in range(len(prcs)):
        f_vin.writelines(prcs[i]+'\n')
    
    # Closure of the file
    f_vin.close() 
    # -------------------

    # Printing-out the just-made-vpsc7.in file
    if ishow==True:
        if os.name=='nt' or os.name=='posix':
            print "****** VPSC input file ******"
            print "*  VPSC7 input parameters   "
            print "* has been written down to  "
            print "* '%8s' as below"%vin_name
            print "*****************************\n"
        if os.name=='nt':
            os.system('%s %s'%('type',vin_name))
        elif os.name=='posix':
            os.system('%s %s'%('cat',vin_name))
    else: pass

def histmaker(
              #filename
              filename='hist/monotonic',
              # increment number, control flag, increment, temperature
              nstep=None, ictrl=None, eqincr=None, temp=278.,
              #velocity gradient control flag (1:known, 0:unknown)
              iudot=[[1,1,1],[1,0,1],[1,1,0]],
              #velocity gradient tensor
              udot=[[1.,0.,0.],[0.,-0.5,0.],[0.,0.,-0.5]],
              #Stress control flag (1:known, 0:unknown)
              iscau=[0,1,1,0,0,0],
              #Cauchy stress tensor
              scauchy=[0.,0.,0.,0.,0.,0.,], mode=None
              ):
    """
    History file maker as an input to VPSC7 core code.

    --------
    Arguments

        #filename
        filename='hist/monotonic',

        # increment number, control flag, increment, temperature
        nstep=None, ictrl=None, eqincr=None, temp=278.,

        #velocity gradient control flag (1:known, 0:unknown)
        iudot=[[1,1,1],[1,0,1],[1,1,0]],

        #velocity gradient tensor
        udot=[[1.,0.,0.],[0.,-0.5,0.],[0.,0.,-0.5]],

        #Stress control flag (1:known, 0:unknown)
        iscau=[0,1,1,0,0,0],

        #Cauchy stress tensor
        scauchy=[0.,0.,0.,0.,0.,0.,]

        #genuine
    """
    
    FILE = open(filename,'w')
    FILE.writelines('%i %i %f %f '%(nstep, ictrl, eqincr, temp))
    FILE.writelines('  nsteps ictrl  eqincr  temp\n')
    FILE.writelines('* boundary conditions *\n')
    for i in range(3):
        for j in range(3):
            FILE.writelines('%8i  '%(iudot[i][j]))
        FILE.writelines('\n')
    FILE.writelines('\n')

    for i in range(3):
        for j in range(3):
            FILE.writelines('%15.9f  '%(udot[i][j]))
        FILE.writelines('\n')
    FILE.writelines('\n')

    
    FILE.writelines('%15i  %15i  %15i\n'%(iscau[0],iscau[5],iscau[4]))
    FILE.writelines('%15s  %15i  %15i\n'%(' ',     iscau[1],iscau[3]))
    FILE.writelines('%15s  %15s  %15i\n\n'%(' ',   ' ',     iscau[2]))
    
    FILE.writelines('%15i  %15i  %15i\n'%(scauchy[0], scauchy[5], scauchy[4]))
    FILE.writelines('%15s  %15i  %15i\n'%(' ',        scauchy[1], scauchy[3]))
    FILE.writelines('%15s  %15s  %15i\n'%(' ',        ' ',        scauchy[2]))


def incr_modifier(fname='hist/monotonic', 
                  nstep=1, 
                  inc=0.005, ictrl=7):
    """
    documents to be filled in

    This module is not actively being used anymore.(2011-02-06)
    """
    try: f = file(fname, 'r')
    except:
        print 
        print '#############################################'
        print ' File : ', fname, ' was not found'
        print ' Please check the filename or directory again'
        print '#############################################'
        print
        raw_input('IOError will be raised. Please press any key')
        f.close()
        raise IOError
    
    lines = f.read()
    lines = lines.split('\n')
    f.close()
    
    f = file(fname,'w') #open the file again for over-writing
    
    row1 = lines[0].split()
    row1[0] = str(nstep)
    row1[1] = str(ictrl)
    row1[2] = str(inc)

    for i in range(len(lines)):
        if i == 0 :
            for j in range(len(row1)):
                f.writelines(row1[j] + '   ')
            f.writelines('\n')
        else:
            f.writelines(lines[i])
            f.writelines('\n')
    
    f.close()


def itran(alpha=2.1, beta=0.368, n=1.8,iplot=False):
    """
    transformation file writer
    """
    FILE = open('transformation.in','w')
    FILE.writelines('* PHASE TRANSFORMATION SETTINGS FOR SOME AUSTENITIC STEELS\n')
    FILE.writelines('%2i %2i   nph(AUST) to nph (MART) currently no effect \n'%(1,2))
    FILE.writelines('* VARIANT SELECTIOIN CRITERIA')
    FILE.writelines('(strain:1 or stress:2 or strain*stress:3) \n')
    FILE.writelines('%2i \n'%(1))
    FILE.writelines('* VARIANTS INFORMATION WHERE MARTENSITE NUCLEATE AND GROW\n')
    FILE.writelines('sx/invariants/gamma_ass.SX   !currently no effect\n')
    FILE.writelines('sx/invariants/alpha_ass.SX   !currently no effect\n')
    FILE.writelines('* PHASE EVOLUTION PARAMETERS(alpha,beta,n)\n')
    FILE.writelines('%7.3f  %7.3f  %7.3f '%(alpha, beta, n))

    if iplot==False: pass
    elif iplot==True:
        fig = plt.figure()
        print " It's still under construction!"
        print " You ought to see if this will be helpful"
        print " If so, please try to fill this in"

    else: 
        print 'iplot to itran must be either True or False'
        raise IOError

def rotfile(filename='r', ang=1):
    """ make in-plane rotate file under cwd's 'rot' directory """
    ##checking if the directory exists
    if os.path.exists('rot') and os.path.isdir('rot'): pass
    else: os.mkdir('rot')

    ##In-plane rotation matrix based on the angle
    th = ang*np.pi/180.
    sth = np.sin(th) ; cth = np.cos(th)
    rmat = np.array([[cth,sth,0],[-sth,cth,0],[0,0,1]])

    filename = '%s%s%s'%('rot', os.sep, filename)
    FILE = open(filename,'w')
    FILE.writelines('rotation matrix for polycrystalline aggregate\n')
    for i in range(3):
        for j in range(3):
            FILE.writelines('%9.4f  '%rmat[i][j])
        FILE.writelines('\n')
    FILE.close()
