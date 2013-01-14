"""
Provides VPSC parametric study frame with its applications.

Author : Youngung Jeong
         Materials Mechanics Laboratory
         Graduate Institute of Ferrous Technology
         Pohang University of Science and Technology

important modules and class
  --
  def vp
  class wrk_cnt  --> ex#1 using class vpsc for work contour plotting
  class multipath --> ex#2 for multiple number of typical monotonic loading paths
                      A particular loadings are normalized so that compared to 
                      experimental counterparts which are available in MML.

necessary libraries to operate with this program
 * Third party libraries (available in internet)
  - numpy             python numerical library 
  - scipy             python scientific library
  - matplotlib        python plotting library

 * libraries of developer's own
  - sx_maker.py       author's single crystal file maker
  - vpsc_in.py        author's VPSC7.in file maker
  - pf.py             author's pole figure plotting script incorporated with 
                                pole8.exe program provided by Tome and Lebensohn

 * VPSC
  - VPSC fortran code executable(s).
  - VPSC7.exe (for NT system) or VPSC (for POSIX system)
 * pole figure contour plotter
  - pole8 program that should be compile for the given system
     (either for NT:Windows or POSIX: Linux)
     
Environment:
  - Ipython (recommended) : ipython -pylab
     under either windows or linux
  - Author's own method to maximize the performance of this script
   is that folders of the local computer are shared to the linux/GNU
   workstation in MML to make use of the 8 CPUed computing.
   Through 'putty' a small SSH terminal, one can access to that linux
   provided that he/she has an account there. 
   
**** USEFUL EXAMPLES
  -- wrk_cnt:
        * Equivalent plastic work contour class.
        * Mimics in-plane biaxial tester (Kokusai)
        * Controls over strain rate ratio (u11/u22) to span a specific areas
          of the resuling two major stress components
  -- multipath
        * Simulate typical mechanical tests
        * Writes and plots the results easily and altogether.
  -- r_ngr
        * Lankford coefficient probing for different texture file.
          If one has different texture file with each's own characteristic,
          One can compare its effects upon lankford coefficient.
  -- r_text
        * Lankford coefficient (R-value) probing from RD to -RD (0~pi)
        * To see if the material's anisotropic axes are well-alligned.
  -- ev_pcys
        * Impose Polycrystalline aggregate of Monotonic proportional loading
        * At the same time, yield surface on each status is observed.
  -- R_prob 
        * LANKFORD coefficient probing spanning the given directions between
          +RD to -RD (0~pi). Unlike r_text, polycrystalline aggregate is all 
          the way tensiled upto a certain amount of deformation. Then its
          transverse and longitudinal strains are linearly fitted to form a
          straight line. The slope of the fitted line is returned. Noticeably
          this is equivalent to experimental way to obtain R-value, except here
          everything is plastic. Polycrystal code including elastic is being
          developed by other people.
        * I checked the R-value measurement using two different method here 
          available, one by subroutine and one by manunal tension. The result
          is included in the MST presentation. 
        * the initial values of eps is fixed to be zero.
  -- practical scripting for specific problem.
        * What F. Barlat asked me to do.#1 
  
  -- FLD (planned)
  
"""

#------------------------------
#  IMPORTS important LIBRARIES 
import matplotlib.pyplot as plt
import glob
import os
import math
import crss
import euler
eul = euler.euler  #Euler module imported as eul.

if os.name=='nt': import winsound  #This is actually for fun

try: import scipy.integrate as integrate
except: print 'Could not import sciy.integrate'; raise IOError
try: import subprocess ; Popen = subprocess.Popen
except: print 'Could not import subprocess.Popen'; raise IOError
try: import pp
except: 
    print 'Could not import pp ' #; raise IOError

try: reload(np)
except: import numpy as np
try: reload(vpsc_in)
except: import vpsc_in
try: reload(sx_maker)
except: import sx_maker
try: reload(pf)
except: import pf
try: reload(shutil)
except: import shutil

pf = pf.pf
v_in = vpsc_in.vpsc_in

def __diff__(x,y):
    """ Differentiates x and y """
    xdiff = np.diff(x)
    ydiff = np.diff(y)
    delta = []
    for i in range(len(xdiff)):
        delta.append(ydiff[i]/xdiff[i])
    return delta

def fin(filename,nhead,*args):
    """ 
    given the file, returns the colums
    of whose indices are given

    ---------
    Arguments
    filename, nhead, *args-->column id
    """
    
    FILE = open(filename,'r')
    lines = FILE.readlines()
    iline = 0
    ## get rid of head
    lines = lines[nhead:len(lines)]
    LINE = [ ]
    for i in range(len(lines)):
        try:
            LINE.append(map(float,lines[i].split()))
        except: pass
    LINE = np.array(LINE)
    LINE = LINE.transpose()
    value = []
    for i in args: value.append(LINE[i])
    return value

def __slope__(x,y):
    """ Calculates x-y's slope\n arg: x, y """
    z = np.polyfit(x,y,1)
    return z[0]

def __diff2__(x,y,delt):
    """
    x = [0,1,2,3,4,5,6... 99]
    len(x)? -> 100
    
    """
    z=[]
    for i in range(len(x)):
        if i-delt>1 and i+delt< len(x)-2: 
            z.append(__slope__(x[i-delt:i+delt], y[i-delt:i+delt]))
        else: z.append(0.,)
    z=np.array(z)
    return z

def __rmtree__():
    """ Remove vpsc_cwd__* directories  """
    files = glob.glob('vpsc_cwd__*')
    print files
    raw_input('Is that Okay? >> ')
    for i in files: shutil.rmtree(i)

def __rmany__(ask=True,*args):
    """
    Remove either directory or file
    *args 
    ask=True  
          - If True ask if to proceed, 
           otherwise passes by it without asking.
    """
    failedFILES = []
    files = []
    for i in range(len(args)):
        files.append(args[i])
        print files
    if ask == True: pass
    else:
        raw_input('Is that Okay? >> ')
    try: ##try to delete a directory if it fails think it as a file
        for i in range(len(args)): shutil.rmtree(args[i])
    except:
        try: 
            for j in range(len(args)): os.remove(args[j])
        except:
            print "Fail to delete %s"%args[j]; 
            failedFILES.append(args[j])
            #raw_input("Enter to proceed >>> ")
    if len(failedFILES)>0:
        return failedFILES

def __rmall__():
    """
    Delete all resulting files generated by vpsc_param.py
    
    ** If you want to add a bunch of files to delete,
    just addd its wildcard to the 'wildcards' list.
    """
    wildcards = ['vpsc_cwd__*', 'R_*', '_R_*', 'plwrk_*',
                 'wrk_cntr_*', 'norm_work*', 'Rprob*', 'YSprob*',
                 'flowwrk_*', 'interwrk_*', 'flowstress_*',
                 'workings*','loadings*','straining*',
                 'sig_*.png','eps_*.png','se_*.png',
                 'YS_*']
    FAILED = []
    for card in wildcards:
        cfiles = glob.glob(card)
        for cfile in cfiles: FAILED.append(__rmany__(True, cfile))

    if any(FAILED[i] for i in range(len(FAILED)))!=None and len(FAILED)!=0:
        print "\nFailed files are as below"
        for i in FAILED:
            if i!=None: print i
            else: pass

    """
    files = glob.glob('vpsc_cwd__*')
    for i in range(len(files)): __rmany__(True, files[i])
    files = glob.glob('R_*')
    for i in range(len(files)): __rmany__(True, files[i])
    files = glob.glob('_R_*')
    for i in range(len(files)): __rmany__(True, files[i])    
    files = glob.glob('plwrk_*')
    for i in range(len(files)): __rmany__(True, files[i])
    files = glob.glob('wrk_cntr_*')
    for i in range(len(files)): __rmany__(True, files[i])
    files = glob.glob('norm_work*')
    for i in range(len(files)): __rmany__(True, files[i])
    files = glob.glob('Rprob*')
    for i in range(len(files)): __rmany__(True, files[i])   
    files = glob.glob('YSprob*')
    for i in range(len(files)): __rmany__(True, files[i])     
    files = glob.glob('flowwrk_*')
    for i in range(len(files)): __rmany__(True, files[i])    
    files = glob.glob('interwrk_*')
    for i in range(len(files)): __rmany__(True, files[i])    
    files = glob.glob('flowstress_*')
    for i in range(len(files)): __rmany__(True, files[i])        
    """

def __makenp__(*args):
    temp = []
    for i in range(len(args)):
        temp.append(np.array(args[i]))
    return temp

def u11u22(step):
    """
    """
    theta = np.linspace(0.,np.pi/2.,step)
    for i in range(len(theta)):
        x = np.cos(theta)
        y = np.sin(theta)
    return x,y

def _rot_(filename=None,ang=0.):
    """rotate the given texture file"""
    rmat = _rotf_(ang=ang)
    rmat = np.array(rmat)
    FILE = open(filename,'r')
    lines = FILE.readlines()
    FILE.close()
    header = lines[0:4]
    grains = []
    for i in range(len(lines)):
        if len(lines[i])< 3: pass
        else:
            try: map(float,lines[i].split())
            except: pass
            else: grains.append(map(float,lines[i].split()))

    for i in range(len(grains)):
        a = eul(ph=grains[i][0], th=grains[i][1], tm=grains[i][2],echo=False)
        a = np.array(a)
        temp = np.dot(a,rmat)
        ph1,ph,ph2=eul(a=temp,echo=False)
        #print 'grains[i]=',grains[i]
        grains[i]=[ph1,ph,ph2,grains[i][-1]]
        #print '\ngrains[i]=',grains[i]
        
    FILE = open(filename,'w')
    for i in range(4): FILE.writelines('%s'%header[i])
    for i in range(len(grains)):
        FILE.writelines('%8.2f%8.2f%8.2f%13.5e\n'%(grains[i][0],
                                                      grains[i][1],
                                                      grains[i][2],
                                                      grains[i][3]))
    pass
    
def _rotf_(filename=None, ang=0):
    """
    Make a rotation matrix file for  polycrystalline aggregate
    """
    th = ang*np.pi/180.
    sth = np.sin(th) ; cth = np.cos(th)
    rmat = np.sin(th) ; cth = np.cos(th)
    rmat = np.array([[cth,sth,0],[-sth,cth,0],[0,0,1]])
    if filename==None: return rmat
    FILE = open(filename,'w')
    FILE.writelines('*** rotation maxtrix for polycrysalline aggregate \n')
    for i in range(3):
        for j in range(3):
            FILE.writelines('%9.4f  '%rmat[i][j])
        FILE.writelines('\n')
    FILE.close()


def __interpolate__(x,y,z,value):
    """
    list x,y,z.
    Estimates x and y for given z-value
    """
    zmn = min(z); zmx = max(z)
    print ' zmn, zmx = ', zmn, zmx
    #raw_input("Enter >> ")
    
    if value < zmn:   # zmn is always 0
        print "Given value is below the minimum of z"
        #raw_input()
        #return 0,0
    elif value > zmx: # safe
        print 'Given value is higher than the maximum of z'
        #raw_input()
        return 0,0

    item = 0
    while True:
        if z[item]> value: break
        else:item += 1
    
    z0 = z[item-1]
    z1 = z[item]
    
    x0 = x[item-1]
    x1 = x[item]
    
    y0 = y[item-1]
    y1 = y[item]
    
    xslope = (x1-x0)/(z1-z0)
    yslope = (y1-y0)/(z1-z0)
    
    newx = xslope * (value - z0) + x0
    newy = yslope * (value - z0) + y0
    
    return newx, newy
    
def __texf__(filename=None):
    """ 
    Given the filename of textfile analyze upon
    1. Ellipsoid axes
    2. Euler angles (deg)
    3. Texture nomenclature (Bunge.. hopefully)
    4. number of grain

    Note that in a texture file, there can be several blocks in it.
    It is recommended that a filename argument includes its full
    path-filled name.
    
    *** UNDER CONSTRUCTION TO BE CONTINUED AFTER MY GATE REVIEW
    """
    class grain:
        def __init__(self,ngr,odfs):
            self.ODF = np.resize((),(ngr,4))
            self.ODF = odfs
            """
            for i in range(len(ngr)):
                ODF[i] = [
            """
    class GRAIN:
        def __init__(self, ngr):
            self.ODF = np.resize((),(ngr,4))
                    
    class info:
        def __init__(self,ib):
            self.nomen ='BUNGE'
            self.ellipax = []
            self.ellipag = []
            self.eps = []
    # DATA CLASS TREE
    #text - blocks - grain
    #              - info
    class text:
        def __init__(self,ngr,nb,odfs):
            class blocks: 
                """ 
                Under a blocks , one will have gr, the list.
                a gr list is actually a block of texture.
                Under gr[i] you will have ODF.
                That ODF is a list consisted of three euler angles 
                and ODF intensity.
                """
                def __init__(self):
                    self.gr = []

            self.blocks = blocks()
            self.blocks.info= info(ib=nb)
            print 'nb= ',nb
            for i in range(nb):
                self.blocks.gr.append(grain(ngr[i],odfs[i]))

    FILE = open(filename, 'r')
    contents = FILE.read(); FILE.close()
    lines = contents.split('\n')
    ib = 0;
    iline = 0
    ngr = []; nomen = []; ellipax = []; ellipag = []; eps = []
    texblocks = []

    #print '\n'
    while True: #LOOP OVER BLOCKS
        try: lines[iline]
        except IndexError: break
        if len(lines[iline]) < 2: break

        print '#%i block'%(ib)
        nerr = 0
        while True:               #### LOOP OVER ONE BLOCKS' HEADER
            try: map(float, lines[iline].split())
            except ValueError: nerr += 1; iline += 1
            else: break
        ngr.append(int(lines[iline-1].split()[1]))
        nomen.append(lines[iline-1].split()[0])
        ellipax.append(map(float,lines[iline-3].split()[0:3]))
        ellipag.append(map(float,lines[iline-2].split()[0:3]))
        eps.append(float(lines[iline-4].split('STRAIN =')[1]))
        """        
        print 'ngr =' , ngr[ib]
        print 'nomen = ', nomen[ib]
        print 'ellipax = ' , ellipax[ib]
        print 'ellipag = ' , ellipag[ib]
        """

        texblocks.append(GRAIN(ngr[ib]))

        print 'ngr[ib] = ', ngr[ib]
        for i in range(ngr[ib]):  #### LOOP OVER ONE BLOCKS' GRAINS
            "[ph1, phi, phi2, ODF]"
            "print map(float, lines[iline].split())"
            "raw_input()"
            texblocks[ib].ODF[i] = map(float,lines[iline].split())
            iline += 1
        ib += 1
    
    odfs = []
    for i in range(ib):
        odfs.append(texblocks[i].ODF)
        
    mytext = text(ngr=ngr, nb = ib, odfs=odfs)
    print "\n*************************************************"
    print "... %15s"%(filename[len(filename)-15:len(filename)])
    print "file has %i blocks of data in it"%(ib)
    print "*************************************************"

    
    return mytext


def __slc_col__(ind, FILE=None, filename=None, 
                fout = None, header =' ** HEADER'):
    """
    Provided the file object, returns colums of the given index
        Arguments:
            FILE : a file object
            ind  : a list object consisting of integer elements
                   for index of column to be returned
            fout= None : The name of output file
            header = ' ** HEADER'
    """
    if type(FILE).__name__!='file': 
        if filename==None:
            print "Error: Given FILE argument is not an actual 'file' object"
            print "You must at least give either FILE or filename"
            raise IOError
        else:FILE=open(filename,'r')

    if type(ind).__name__!='list':
        print "Error: Given ind argument is not an actual 'list' object"
        raise IOError

    rst = []; iline = 0

    lines = FILE.readlines()
    for i in range(len(lines)):
        if len(lines[i])< 3: break  # This is to see if e.o.f is reached
        try: 
            map(float,lines[i].split())
        except ValueError: pass
        else:
            cwd = map(float, lines[i].split())
            rst.append([])
            for j in ind:
                try: rst[iline].append(cwd[j])
                except: pass
            iline = iline + 1
    
    if fout == None:  return rst
    else:
        if type(fout).__name__=='str':
            outf = file(fout, 'w')
            outf.writelines(header); outf.writelines('\n')
            for i in range(len(rst)):
                for j in range(len(rst[i])):
                    try: outf.writelines('%18.9e'.rjust(11)%(rst[i][j]))
                    except: pass
                outf.writelines('\n')
        else:
            print 'Error: Inappropriate fout argument as a file name'
            raise IOError

#                   Pole figure plotting 
#-----------------------------------------------------------------------
def __pf_plot__(switch, osname, pf_name, ftex, ishow, idot, sx_file):
    """
    Provided the arguments as follow, plots the relevant pole figure
    """
    if switch == False: return -1
    if os.name == 'nt' : delfiles = glob.glob('*.dat')
    elif os.name == 'posix' : delfiles = glob.glob('*.dat') + glob.glob('*.DAT')
    else: print 'Error: Unexpected os.name'; raise IOError
    for i in delfiles: os.remove(i)
    try:
        pf(pf_name=pf_name, ftex=ftex, 
           ishow=ishow, idot=idot, sx_file = sx_file)
    except IOError: pass


#                   MODIFIED VOCE HARDENING CURVE 
#-----------------------------------------------------------------------
def voce(x, tau0=1., tau1=0.5, thet0=0.8, thet1=0.2, ifig=None):
    """
    Returns modified voce hardening flow curve for single slip system
    Arguments
        x : Supposed to be x linear spaced values
        tau0= 1, tau1 = 0.5, thet0 = 0.8, thet1 = 0.2
        ifig = None (can be a integer for figure numbering)
    """
    temp = []

    for i in range(len(x)):
        temp.append(tau0+(tau1+thet1*x[i])*(1-math.exp(-x[i]*abs(thet0/tau1))))
    if ifig==None: pass
    else:
        f = plt.figure(ifig)
        ax=f.add_subplot(111)
        ax.plot(x,temp)
        plt.ylim(ymin=0, ymax=max(temp)*1.4)
    return temp
#-----------------------------------------------------------------------



#             FILE PLOT FUNCTION fplot declaration block
#-----------------------------------------------------------------------
def fplot(ix, iy, filename, 
          ifig = None, nhead=1, ABS=False, 
          idot=False, label=None, marker='o' ):
    """
    Provided filename as well as indices of x and y column,
    plot and returns x vs y only if ifig is given as None

    Arguments:
        ix, iy : indicies of columns
        filename : file name
        ifig : figure id
        nhead=1  : number of head lines in the given file object
        ABS = False : If true, apply abs(values) to the columns of data
    """
    try:
        FILE = open(filename,'r')
    except IOError:
        print 'Could not open ', filename
        print 'fplot method returns -1'
        return -1,-1
    try: plt
    except:
        import matplotlib.pyplot as plt
    #c   X and Y initilization
    x, y = [], []
    
    #c  examining the type of File
    if FILE.mode =='r':
        pass
    else:
        print 'Wrong file mode'
        print 'function column_plot returns -1'
        return -1

    if type(FILE).__name__ == 'file':
        pass
    else:
        print '******************************************'
        print 'Error: wrong type variable given to FILE'
        print "The argument file should be type of 'file'"
        print '******************************************'
        print 
        raise ValueError
    source = FILE.read()
    FILE.close()

    lines = source.split('\n')  #c  data are splitted to list of lines
    iline = 0
    iline = iline + nhead  #This can be risky, however, is being used for getting exact lines of header
    idum = 0
    y.append([]); x.append([])
    #c  LOOP OVER LINES TO PUT THEM INTO X and Y
    while True:
        #c  breaks if E.O.F
        try: lines[iline]
        except IndexError: break

        if len(lines[iline]) < 5:
            break
        
        try:
            float(lines[iline].split()[0])
        except ValueError:
            x.append([]); y.append([])
            idum = idum + 1   #This dummy index is for telling the block of data.
            #print 'lines[iline].split()[0] = ', lines[iline].split()[0]
            #raw_input('...')
            pass
        else:
            if ABS == False:
                try:
                    #print lines; raw_input()
                    float(lines[iline].split()[ix])
                    float(lines[iline].split()[iy])
                except ValueError: pass
                else:
                    x[idum].append(float(lines[iline].split()[ix]))
                    y[idum].append(float(lines[iline].split()[iy]))
                
            elif ABS == True:
                try:
                    float(lines[iline].split()[ix])
                    float(lines[iline].split()[ix])
                except ValueError: pass
                else:
                    x[idum].append(abs(float(lines[iline].split()[ix])))
                    y[idum].append(abs(float(lines[iline].split()[iy])))
            else: print 'Error in ABS argument'; raise IOError

        iline = iline + 1

    #c  plotting
    if ifig==None: return x,y
    else: 
        fig = plt.figure(ifig)
        ax = fig.add_subplot(111)
        if idot==True:
            for i in range(len(x)): plt.plot(x[i],y[i],marker)
        else:
            for i in range(len(x)): plt.plot(x[i],y[i])
        return x, y
    

#-----------------------------------------------------------------------
#             CYCLIC SIMPLE SHEAR STRESS-STRAIN OUTPUT
#-----------------------------------------------------------------------
def ss(filename='ss_curv.out',
       __fnstr__='str_str.out',
       header='Header not given',
       indx=7, indy=13,
       iplot=False,
       ifig = None, 
       ABS = False
       ):
       #incr = 0.0025):
       #c  --> increment size is found from str_str.out
    """
    Provided the result file listing columns of variables,
    x and y columns are plotted and saved to a file whose name is given.
   
    This method is optimized to cyclic simple shear test.
    General solution can be found from fplot method.

    Arguments:
         filename ='ss_curv.out'
         __fnstr__ ='str_str.out'
         header ='header not given'
         indx = 7, indy = 13
         iplot =False
         ifig = None
    """
    try: f = file(__fnstr__,'r')
    except IOError: 
        print 'No such a file is found. IOError is raised in vpsc_param.ss function'
        raise IOError
    lines =f.read()
    lines = lines.split('\n')
    f.close()

    x, y = [], []
    #c checking if header of string
    iline = 0
    istep = 0
    while True:
        if len(lines[iline]) < 3 : break
        #print 'lines[iline].split() \n', lines[iline].split()
        try: float(lines[iline].split()[0])
        except: iline = iline + 1
        else:
            col = map(float, lines[iline].split()[0:max(indx,indy)+1])
            try:
                xtmp, ytmp = col[indx], col[indy]
            except:
                print indx, indy
                raw_input()
            
            if istep == 0:
                x.append(xtmp)
                y.append(abs(ytmp))
                pass
            else:
                #c  determining xstrain increment
                if istep == 1: 
                    inc = xtmp - x[0]
                #c  to see if loading path has changed
                if abs(x[istep - 1] - xtmp)< 1.0e-6 :
                    x.append(xtmp)
                    y.append(abs(ytmp))
                else:
                    x.append(x[istep - 1] + inc)
                    y.append(abs(ytmp))
                    pass #c  case that loading path is not changed
                
            istep = istep + 1 #c  not decided if this is the right place
            iline = iline + 1

    try: reload(plt)
    except: import matplotlib.pyplot as plt

    if ifig==None: fig = plt.figure()
    else: fig=plt.figure(ifig)
    ax = fig.add_subplot(111)
    ax.plot(x,y)
    f = file(filename,'w')
    f.writelines(header+'\n')
    for i in range(len(x)):
        f.writelines('%8.5f %8.5f \n'%(x[i],y[i]))
    f.close()
#-----------------------------------------------------------------------

def vpscrun():
    os.system('a')
    ss()

def __files_move__(iloop = 0):
    """
    make current working directory : cwd_dir
    into which everything ending with the extension of 
    '.out' is put
    """
    prep = 0
    while True:
        cwd_dir = 'it'+str(prep).zfill(3)+'_'+str(iloop).zfill(5)
        if len(glob.glob(cwd_dir))==0: break
        else: 
            prep = prep + 1
            pass
        
    try:
        os.mkdir(cwd_dir)
    except:
        print ;print 'Unexpected Error in making a directory';print
        raise IOError
    
    if os.name=='nt': files = glob.glob('*.out')
    elif os.name=='posix': files = glob.glob('*.out') + glob.glob('*.OUT')

    for i in files: shutil.move(i,cwd_dir)

def __plot__(plotlist=[], iloop=0):
    """
    Provided the list of plot indicies,
    decide which to proceed to plot certain graphs 

    Arguments:
        plotlist = ['True','False' ... ] 
        iloop : # of looping in case of multi runs
    """
    #polefigure
    if plotlist[0] == True:
        #-----------------  POLE FIGURES  ---------------
        ftex = glob.glob('%s%s%s'%(self.cwd+os.sep+'TEX_PH*.OUT'))
        for itexf in ftex:
            __pf_plot__(switch = True, osname=os.name, pf_name='aaaaa',
                        ftex=itexf, ishow=True, idot=False,
                        sx_file='sx/hijhiihb.sx')


    #uniaxial stress-strain curve
    if plotlist[1] == True:
        #---------------  STRESS-STRAIN CURVES ------------
        """
        " Cyclic Stress strain curve plotting"
        ss(filename='ss_curv_' + str(iloop) + '.ss',
           __fnstr__='STR_STR.OUT',
           header='Header not given',
           indx=2, indy=8,      #E11-S11
           indx=7, indy=13,     #E12-S12
           iplot=True, ifig=20,
           )
        """
        dumx, dumy = fplot(filename='STR_STR.OUT', 
                           ix=2, 
                           iy=8, 
                           ABS=True, 
                           idot=True)#ifig = 20
        fig = plt.figure(20)
        plt.plot(dumx, dumy, label=[None,'RD -VPSC','TD -VPSC'][iloop])

        et, ew = fplot(filename='STR_STR.OUT',
                       ix=4,
                       iy=3,
                       ABS=False)

        """ Backward instantaneous R-value measure """
        R=[]
        x=[]
        for i in range(len(et)):
            for j in range(len(et[i])):
                x.append(dumx[i][j])
                if j==0: R.append(None); 
                else:    
                    R.append((ew[i][j]-ew[i][j-1])/(et[i][j]-et[i][j-1]))
                    
        R[0] = R[1]
        
        fig = plt.figure(21)
        ax = fig.add_subplot(111)
        ax.plot(x, R, marker='o',label = [None, 'RD -VPSC','TD -VPSC'][iloop])
        plt.ylim(0,2)

    #pcys
    if plotlist[2] == True:
        #------------  YIELD SURFACE PLOTS -------------

        #cwf = open('pcys_pl.out', 'r')
        cwf = 'pcys_pl.out'
        cwx_, cwy_ =[], []
        #YIELD SURFACE made by D1-D5 probings
        cwx, cwy = fplot(ix=0, iy=2, 
                         filename=cwf, ifig=40)   
        
        #cwx, cwy = fplot(ix=0, iy=1,filename=cwf, ifig=40)   #YIELD SURFACE made by D1-D2 probings
        cwf.close()
        cwf = open('YS_'+str(iloop)+'.ys', 'w')
        cwf.writelines('*** Yield surface probings \n x  y \n')
        for i in range(len(cwx)):
            for j in range(len(cwx[i])):
                try: cwf.writelines('%8.4e  %8.4e \n'%(cwx[i][j], cwy[i][j]))
                except: raw_input("Error: what's wrong here?")
                cwf.writelines('%8s %8s\n'%('--','--'))
        cwf.close()
        
    #pcys in 3d
    if plotlist[3]==True:
        #------------- YIELD SURFACE PLOTS --------------
        #---------- AS IN A FORM OF POINT CLOUD ---------
        ysfn = 'ex_'+str(iloop)+'.ys'
        __slc_col__(FILE=file('pcys_pl.out','r'), 
                    ind=[0,1,2],
                    fout=ysfn)
        irst = os.system('qhull <'+ysfn+' -G > '+ysfn.split('.')[0]+'.qhl')
        if irst!= 0 :
            print 'Excuting qhull possibly resulted in an error'
            print 'Please check this out'
            raise IOError

    if plotlist[4]==True:
        #------------- LANKFORD COEFFICIENT PROBING -----
        #
        cwx, cwy = fplot(ix=0, iy=2, filename='LANKFORD.OUT',ifig=90)
        #cwx, cwy = fplot(ix=0, iy=2, filename=cwf, ifig=40)   
        #YIELD SURFACE made by D1-D5 probings

    #pcyses processed with qhull and moving to geomview



    # BULGE PLOT _EFFECTIVE
    if plotlist[5]==True:
        cwx,cwy = fplot(ix=4, iy=10, filename='STR_STR.OUT', ABS=True)
        plt.plot(cwx,cwy, label = ['bulge -VPSC',None,None][iloop])
        

def cod(FILE, phi2=50.):
    """
    Getting a particular slice from DIOR.exe

    2010-12-23 : Displosed now, since different way of getting to the slice
                 of the Euler space is more feasible. I'm thinking that maybe
                 it is more plausible to get the tex_ph*.out file into grided.

    2011-06-14 : Please refer to the cod_section.py.
    """
    lines = FILE.read()
    lines = lines.split('\n')
    i = 0
    ir = 0
    rst = []
    switch=False
    while True:
        "if all(len(lines[j]) < 4 j=1,2 ): break"

        try: lines[i]
        except:
            print 'lines[i]', lines[i]
            print ' i = ', i
            return lines[i]

        if len(lines[i])<3:
            if (switch==True): return cwb
            else: iblock_started=False
            if len(lines[i+1])<3:
                print 'EOF'
                return -1
        else: pass

        try: map(float,lines[i].split())
        except: i = i + 1;  iblock_started=True; cwb=[]

        else:
            if i == 2:
                header = lines[i-1]
                try:
                    phi2_ = float(header.split('phi2=')[1])
                    #phi1 = float(header.split('phi1=')[1])
                    #phi = float(header.split('phi=')[1])
                    phi1i=float(header[5:10])
                    phi1f=float(header[10:15])
                    phii =float(header[15:20])
                    phif =float(header[20:25])
                except: print 'Error: Unexpected header type'; raise IOError
            print 'phi2=', phi2
            print 'phi2_=', phi2_

            if abs(phi2_-phi2)<0.1 : switch=True
            
            cwl = lines[i]
            cwl = cwl[1:len(cwl)]
            cwc = []
            step = len(cwl)/4
            for j in range(step):
                cwc.append(int(cwl[0+j*4:0+(j+1)*4]))
            cwb.append(cwc)
            #cwb=np.array((cwb[0:-1],cwc]))
            #print 'cwc = ', cwc
            #raw_input()
            #if i == 10: return cwc
            i = i + 1
            ir = ir + 1


# note that this should be changed into class so that access to data
# can be easily handled aftermath.


#            VPSC WRAPPING SCRIPTS FOR PARAMETRIC STUDY
#-----------------------------------------------------------------------
def vp(del_crss=False, norun=False, nloop = None):
    """
    VPSC PARAMETERIC STUDY MANIPULATOR.
    ENHANCES AUTOMATIZATION OF VPSC RUN AND PLOTTING THE RESULTS.

    IT DOES NOT PROVIDE THE ARGUMENTS THAT AN END-USER CAN PLAY WITH.
    INSTEAD, THE SCRIPTS NEED DIRECT MODIFICATION.


    (2011-01-30) It is now almost disposed. I'm currently into class vpsc.
    class vpsc is thought to be better for making a repeated stuffs. 
    Using class vpsc, some derivatives are readily made. Some examples are
    appended. And some more are planned to be done.
    """
    #
    # COMMENTS:
    # REFER TO THE DOCUMENT FOR GENERAL PURPOSE OF THIS FUNCTION!
    #
    
    #### PREPARE FOR THE RUN ####
    
    """
    flags in VPSC.in file

      update flags (0: off, 1:on)
      1. orientation, 2. grain shape, 3. hardening

      interaction(linearization scheme):
      FC, affine, secant, neff, tanget, second-order

      rate-sensitive (0: insentitive 1:sensitive)

      nneigh (0: no 1 for pairs, and so on)

      texture file
      single crystal file
    """
    

    if nloop == None:
        print 'Type the number of loop (nloop)'
        nloop = int(raw_input(' >> '))



    'TEXTURE FILE INPUT : POTENTIAL TEXTURE FILES'

    ngr = [500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
    ngr = map(str,ngr)
    fcortex, fucotex=[],[]
    for i in range(len(ngr)): fcortex.append('304_cor_'  +ngr[i]+'g.tex')
    for i in range(len(ngr)): fucotex.append('304_uncor_'+ngr[i]+'g.tex')


    ##### ITERATION STARTS #####

    for iloop in range(nloop):
        """
        LOOP OVER SOME PARAMETRIC VARIATIONS
        PROCESSES ARE CATEGORIZED IN ACCORDANCE WITH THEIR ROLES.

        AN END USER IS SUPPOSED TO APPLY HIS OR HER MODIFICATIONs TO
        EACH BLOCK IN THE LOOP DIRECTLY.
    
        DUE TO EXISTENCE OF TOO MANY PARAMETERS, MODIFICATIONS ARE APPLIED TO CODE
        DIRECTLY WITHOUT HAVING TO GO THROUGH ARUGMENTS MANIPULATIONS
        """

        """
        common variables:
           can be sx file
           can be hist file
           can be txt file
           and so on...
        """

        #------ VPSC SINGLE CRYSTAL FILE MODIFICATION -----
        #crystal #1
        iopsysx = 0
        sx_maker.cubic(
                       hii=1.0, hij=1.0, hb=-0.4,
                               # 304 1phase / 500 grain:
                               #"tau0=1.02e2, tau1=1.5e2, thet0=3.8e2,thet1=0.99e2,"

                       #304 1phase / 1000grains (DF corrected)
                       tau0 = 1.051e2, tau1=1.42e2, thet0=3.66e2, thet1=1.05e2,  
                             
                       #tau0=1.0, tau1=0.5, thet0=1.2, thet1=0.2,
                       filename='sx/austenite.sx',
                       dependency='indv',
                       n=[[1,1,1]], b=[[1,1,0]],
                       iopsysx=iopsysx 
                       )
        #crystal #2
        sx_maker.cubic(
                       hii= 1.0, hij=1.0, hb=-0.4,
                       tau0=3.0, tau1=0.0, thet0=0, thet1=0,
                       filename ='sx/martensite.sx',
                       dependency = 'iso',
                       n=[[1,1,0]], b=[[1,1,1]],
                       iopsysx=iopsysx 
                       )


        #-----  VPSC HISTORY FILE MANIPUALTIONS  ------
        """ VPSC HISTORY FILE MANIPUALTIONS """
        """
        vpsc_in.incr_modifier(fname='hist/simple_shear_1',
                              nstep=10, inc=0.01, ictrl=7)
        vpsc_in.incr_modifier(fname='hist/simple_shear_2', 
                              nstep=40, inc=0.01, ictrl=7)
        vpsc_in.incr_modifier(fname='hist/simple_shear_3', 
                              nstep=60, inc=0.01, ictrl=7)
        """
        vpsc_in.incr_modifier(fname='hist/axial_x',        
                              nstep=42, inc=0.005, ictrl=7)
        vpsc_in.histmaker(filename='hist/bulge',
                          nstep = 42, eqincr= 0.005, ictrl = 7,
                          iudot = [[  0, 0, 0],[1,   0, 0],[ 1, 1,  1]], 
                          #Cannot release both off-diagonal term at the same time
                          udot = [[0.5, 0, 0],[0, 0.5, 0],[ 0, 0, -1]],
                          iscau = [1,1,0,1,1,1],
                          scauchy = [0,0,0,0,0,0])
        """
        vpsc_in.incr_modifier(fname='hist/axial_x_comp',  
                              nstep=40, inc=0.01, ictrl=7)
        vpsc_in.incr_modifier(fname='hist/axial_x_1', 
                              nstep=10, inc=0.002, ictrl=7)
        vpsc_in.incr_modifier(fname='hist/axial_x_2',
                              nstep=10, inc=0.01, ictrl=7)
        vpsc_in.incr_modifier(fname='hist/axial_x_3',
                              nstep=10, inc=0.01, ictrl=7)
        vpsc_in.incr_modifier(fname='hist/axial_x_comp_1',
                              nstep=10, inc=0.01, ictrl=7)
        vpsc_in.incr_modifier(fname='hist/axial_x_comp_2',
                              nstep=10, inc=0.01, ictrl=7)
        """
        #-----------------------------------------------
 


        #-------- VPSC INPUT FILE MANIPULATION --------- 
        """
        VPSC INPUT FILE MANIPULATIONS
        #----------------------------------------------------------------------#
        # v_in is a loaded function, which originally aims at                  #
        # manipulationg vpsc7.in, the input file, using python script.         #
        # Therefore, for further modification of this function is strongly     #
        # recommend to be done through adjusting the arguments of the vpsc_in  #
        # Or, more for the complete exploitation of this scripts, one is       #
        # expected to manipulate the source script directly                    #
        #----------------------------------------------------------------------#
        """
        
        "0,1 : loading, 2: pcys, 3:LANKF, 4:RBR, 5:pcys_pl, 6:pcys_pl 6d"

        prcs  = ['3','15']
        tenRD = ['0','hist/axial_x']
        tenDD = ['4','rot/rotation.45','0','hist/axial_x','4','rot/rotation.-45']
        tenTD = ['4','rot/rotation.90','0','hist/axial_x','4','rot/rotation.-90']
        bb    = ['0','hist/bulge']
        prcs =[bb,tenRD,tenTD]
        pr = prcs[iloop]
        # multiphase case accounted
        tfile = [['texture/304_cor_5000g.tex',   'texture/rand_500.tex'],
                 ['texture/304_uncor_2000g.tex', 'texture/rand_2000.tex']]
        
        iorient     = 1
        igrainshape = 1
        ihardening  = 1
        itrans      = 0
        update_ctrl = [iorient, igrainshape, ihardening, itrans]

        v_in(
            nph = 1, wph=[1., 0.001],

            ishape=[0,0], fragmentn =[0,0], crit = [25,25],
            ellipangl = [[0.,0.,0],[0.,0.,0.]],
            ellipaxes = [[1.,1.,1.],[1.,1.,1]],

            tfile = tfile[1],
            sxfile = ['sx/austenite.sx', 'sx/martensite.sx'],   
            iratesens = 0, #rate insensitive
            interaction = 3, #(0:FC,1:affine,2:secant,3:neff=10,4:tangent,5:SO)
            iupdate = update_ctrl, #(orient, grain shape, hardening)
            nneigh = 0, #(0 for no neighbors, 1 for pairs, etc.)
            npro = len(pr)/2,
            prcs = pr
            )

        # PHASE TRANSFORMATION 
        #vpsc_in.itran(alpha=2.1, beta=0.368, n=3.8, iplot=False)
        #parent phase
        
        """
        sx_maker.invariant(
                            filename='sx/invariants/gamma.sx',
                            header = '** VARIANTS - gamma',
                            crysym = 'cubic',
                            b=[[1,1,0]], n=[[1,1,1]],
                            iopsysx = iopsysx
                          )
        #child phase
        sx_maker.invariant(
                            filename='sx/invariants/alpha.sx',
                            header = '** VARIANTS - alpha',
                            crysym = 'cubic',
                            b=[[1,1,1]], n=[[1,1,0]],
                            iopsysx = -1
                           )
        """
        #-----------------------------------------------


        #-----------------  VPSC RUN  ------------------
        # I hope I can change this lousy commands
        is_good = []
        try:
            if norun==True: pass
            else:
                #For Linux!
                if os.name=='posix': is_good.append(os.system('./VPSC7'))
                #For Windows!
                elif os.name=='nt' :   is_good.append(os.system('VPSC7.exe'))  
                # The VPSC7.for is supposed to be compiled by g95 compiler
                else: 
                    print "\n\n Unexpected operating system"
                    raise IOError

        except: pass
        if all(is_good[i] != 0 for i in range(len(is_good))):
            print 'Something fishy happened during VPSC7.exe'
            print 'Thus IOError is raised'
            raise IOError
        #-----------------------------------------------


        #-----------------------------------------------
        #            SIMULTANEOUS PLOTTINGS
        #-----------------------------------------------
        """
        Plot simultaneously if it is thought to be handy
          List of possible plots:
              PCYS in 2D
              POLE FIGURE
              STR_STR
              LNKF AND SO ON
        """
        __plot__(
                 plotlist=[
                           False,                     #pole figure
                           [False,False,False][iloop],  #stress-strain(E11-S11)+ R-VALUE
                           False,                     #Yield-surface type1
                           False,                     #Yield-surface type2
                           False,                     #LANKF
                           [False,False, False][iloop] #bulge (E33-S33)
                           ],
                 iloop = iloop
                 )

        if itrans == 1:
            " plot volume fraction evolution "
            cwx, cwy = fplot(ix=0, iy=1, filename='VOL_FRC.OUT', nhead=0, ifig = 10)
            plt.ylim(0,1)

        #-----------------------------------------------
        #   move *.out files into a directory
        #-----------------------------------------------
        __files_move__(iloop = iloop)


        #-----------------------------------------------
        #         notification of end of loop
        #-----------------------------------------------
        if os.name=='nt':
            if len(prcs)/2 > 20:
                winsound.Beep(400, 500)
                winsound.Beep(500, 500)
                winsound.Beep(600, 500)
                winsound.Beep(1000,500)
        else: pass
    ###### ITERATION ENDS ######










class vpsc:
    """
    The OBJETIVE VISCO-PLASTIC SELF-CONSISTENT model

    It will be designed for a single run.
    However, the data and conditions are saved so that
    post-execution access to data is feasible.

    data are saved internally not to files.
    Having files as an access point of store is not 
    efficient in point of management and post-processing

    MPI??? (I'm not sure...)


    It starts from the idea that 
    "
    To save data for a run, better to have a class
    to consecutively access to the data. (2011-Jan-01)
    ".

    VPSC executable compiled by g95 fortran compiler
    is located into a directory together with all necessary files
    copied. This is done to avoid possible interruption between two 

    The executable file is called through "subprocess.Popen" method.
    subprocess.Popen is gloablly imported as Popen

    

    comments
      2011-02-06: f2py module is used now for wrapping
      the fortran code under the project name of vpsc_f2py
      



    #### Bug(s)
    2011-01-31
      It was found that VPSC7 executables compiled g95 for 
      NT and POSIX produced different behavior. One for 
      Linux(POSIX) tends to often raise NEWTON RAPHSON 
      failure, while, with all being the same, that for
      windows does not. Is this a compileer problem or any
      possible bug in the legacy code? I have no idea.
      Ever since I started using Linux for computation,
      compiled executable for Linux produced a reasonable
      results. Despite of the suspicious failure raised 
      in Newton-Raphson loop, all the lankford, stress-
      strain curves looked reasonable. I probably need 
      to get to this point later. Otherwise, I have to
      double check with the results from Windows system.

    2011-02-07
      It turned out that in Linux/GNU system, double preci-
      sion is not default when compiling fortran code. 
      It was informed by Dr. Lebensohn. He kindly suggested
      to invoke double precision when compiling. This seems
      to solve most of the problem(not perfectly since I
      still observe the the singularities.. )


    Author: Youngung Jeong
            Materials Mechanical Laboratory
            Graduate Institute of Ferrous Technology
            Pohang University of Science and Technology
    """ 
    class data:
        pass
    class VPSCin:
        nph = 0
        pass
    class TEXTURE:
        pass
    class ph:
        pass
    class loading:
        class uniaxial:pass
        class biaxial:pass
    def manual(self):
        """ manual in order to use class vpsc """
        print " initiation "
        print " import vpsc_param"
        print " myvp = vpsc_param.vpsc()"
        print " myvp.pp()  --> post-processing "
        print " myvp.datamaster()",
        print "---> a data class including all the resulting data"
        print " myvp.getback() ---> "

    def __init__(self, iwait=True, 

                 ## grain morphology contorl
                 ## angles and principals lengths of 
                 ## axes of the ellipsoid

                  
                 ## loading conditions
                 u11=None, u22=None ,
                 stp=None, eqc=None, ict=None,
                 mode=None, 
                 ang = None, 
                            #ang is useful when mode='ten_ang'
                            #for tension along after rotating as much as ang

                 ##texture file
                 texture=None,

                 ## single crystal file
                 fsx=None,

                 ## Enviroments
                 del_crss=True,

                 ## misc - options  (passed to __vpsc_in__ module)
                 irecover=0, isave=0, icubcomp=0, nwrite=0, ihardlaw=0,
                 iflu=0, iratesens=0, interaction=3, iupdate=[1,1,1,0],
                 nneigh=0,


                 ### precision settings for convergence prodcedure
                 err=0.001, errd=0.001, 
                 errm=0.001, errso=0.01,

                 nmxiter=100, exiter=100, initer=100,
                 irsvar=0, xrsini=2, xrsfin=10, xrstep=2,

                 ### manual prcs
                 prcs=None,

                 ):   #stp can be given (in general it's step)
        """
        Comments 
        (2011-01-30)
          Is it okay for me to include all the possible arguments 
          into def __init__? It is not quite sure but due to 
          python's ability to have defaulted arguments in modules,
          more arguments does not really mean more confusions.
          However, having more arguments is always making me wonder
          it is really necessary.

        Arguments:
          mode = 'ten_ang=ang_ten', 'RD=rd','DD=dd','TD=td',
                 'BULGE=bulge=bu=BU','SIMPLESHEAR=SS=SSRD',
                 'SSTD=sstd','SSDD=ssdd','lankf=LANKFORD=lankford=lnkf',
                 'inplanebiaxial=ipb','inplanebiaxialfc=ipbfc','ev_pcys'
          u11, u22
          stp, eqc
          ict, ang(effective only when mode =='ten_ang')
          

          texture, fsx,
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
          prcs=None,          
        """


        #self.manual()  --> showing manual 
        self.del_crss = del_crss
        self.procwait = iwait
        self.cwd_old = os.getcwd()
        self.cwd = self.__getnewdirectory__() #suggests a non-exisiting 
                                              #directory
    
        #----------------------------------------------------------------------
        self.__mkcopy__(dst = self.cwd) #makes a directory and copies stuffs in
        #----------------------------------------------------------------------
        os.chdir(self.cwd)       #move on to the working directory
        self.cwd = os.getcwd()   #full path name of current working directory 
        self.osname = os.name

        #--------------------------------------------------------
        #-- VPSC input (texture, morphology, phase)
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
                """
                print "The give texture file is below"
                for tex in texture:
                    print "%s"%tex
                """
        else: 
            print " You have not input any texture file to class vpsc"
            print " You must assign a file name"
            raise IOError
            self.ftx = ['texture%s304_cor_1000g.tex'%os.sep,
                        ][0:self.VPSCin.nph]

        #SINGLE CRYSTAL FILE
        if fsx==None:
            self.fsx = ['sx%sph01.sx'%os.sep,
                        'sx%sph02'%os.sep][0:self.VPSCin.nph]
        else:
            if type(fsx).__name__=='list': pass
            elif type(fsx).__name__=='str': fsx=[fsx]
            else: raise IOError
            self.fsx=fsx

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
            self.VPSCin.ph[i].init_elpsd_rat=[[1.,1.,1.],
                                              [1.,1.,1.]][i]
            self.VPSCin.ph[i].init_Eul_axe=  [[0.,0.,0.],
                                              [0.,0.,0.]][i]
            ANG.append(self.VPSCin.ph[i].init_Eul_axe)
            AXE.append(self.VPSCin.ph[i].init_elpsd_rat)

        """
        -- initialization
             prcs, step, eqincr, nhist
        """
        utemp = [u11,u22]
        if any(utemp[i]==None for i in range(2)):
            pass
        else:
            u33= -(u11+u22)  #INCOMPRESSIBILITY
            u12=0; u13=0; u21=0; u23=0; u32=0; u31=0

        if prcs!=None: 
            self.prcs = prcs
            if mode!=None:
                print "Only one of prcs or mode should be given"
                raise IOError
        else: self.prcs = []

        self.step = 0; self.eqincr = []
        self.nhist = 0; self.nrot = 0
        
        """
        -- PROCESSES (LOADING, PCYS, LANK .. )
        """
        # general parameters over historyfileprocesses
        if stp==None: stp =60 
        if eqc==None: eqc = 0.005 
        if ict==None: ict = 7
        if ang==None: ang = 0.

        #### Global mode variable
        self.mode = mode
        ###

        if mode==None: 
            if prcs==None:
                print "Either prcs or mode must be given"
                raise IOError
            # If mode is not given, end-user must modify 
            #the below block according to their needs.
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
        elif mode=='ten_ang' or mode=='TEN_ANG' or mode=='ang_ten' or mode=='ANG_TEN':
            ## tension along a particular direction
            self.__rotate(ang=ang)
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                       ictrl=ict, mode='unitension')
            self.__rotate(ang=-ang)
        elif mode=='RD' or mode=='rd':
            #RD tension
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode='unitension')
        elif mode=='DD' or mode=='dd':
            #DD tension
            self.__rotate(ang=45)
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode='unitension')
            self.__rotate(ang=-45)
        elif mode=='TD' or mode=='td':
            #TD tension
            # The rotation is now excluded
            #self.__rotate(ang=90)
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode='unitensiony')
            #self.__rotate(ang=-90)     
        elif mode=='BULGE' or mode=='bulge' or mode=='bu' or mode=='BU':
            #BULGE
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode='bulge')
        elif mode=='SIMPLESHEAR' or mode=='SS' or mode=='SSRD':
            #SIMPLESHEAR
            self.__historyfileprocess(step=stp, eqincr=eqc, 
                                      ictrl=ict, mode='simpleshearfc')
        elif mode=='SSTD' or mode=='sstd':
            #SIMPLESHEAR ALONG TRANSVERSE DIRECTION
            self.__rotate(ang=90)
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, mode='simpleshearfc')
            self.__rotate(ang=-90)
        elif mode=='SSDD' or mode=='ssdd':
            #SIMPLE SHEAR ALONG DIAGONAL DIRECTION
            self.__rotate(ang=45)
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, mode='simpleshearfc')
            self.__rotate(ang=-45)
        elif mode=='lankf' or mode=='LANKFORD' or mode=='lankford' or mode=='lnkf':
            self.__lankfordprocess(step=stp)
            self.__rotate(ang=90.)
            self.__lankfordprocess(step=stp)
            self.__rotate(ang=-90.)
        elif mode=='inplanebiaxial' or mode=='ipb':
            self.__historyfileprocess(step=stp, eqincr=eqc,u11=u11,u22=u22,
                                      ictrl=ict, mode='inplanebiaxial')
        elif mode=='inplanebiaxialfc' or mode=='ipbfc':
            self.__historyfileprocess(step=stp, eqincr=eqc,u11=u11,u22=u22,
                                      ictrl=ict, mode='inplanebiaxialfc')
            
        elif mode=='ys' or mode=='YS' or mode=='YIELDSURFACE':
            self.__pcysplprocess(step=72, shear_ang=90.)
        elif mode=='ev_pcys':
            ntot = 8 ## INTERMIDEATE STOPPING NUMBER TO PERFORM PCYS 
            totalstrain = stp * eqc
            eachstep_inc = totalstrain/ntot
            eachstep_inc = eachstep_inc/eqc
            eachstep_inc = int(eachstep_inc)
            print " loading is %i times in Total"%ntot
            print " Total Strain imposed is %6.3f"%totalstrain
            print " Each loading has %i"%eachstep_inc,
            print " number of incremental step"
            if raw_input(' Type q to raise IOError >>')=='q': raise IOError
            ##initial PCYS-pl
            nprobing = 32 #number of probings
            self.__pcysplprocess(step=nprobing, shear_ang=90.)
            for i in range(ntot):
                ## LOADING AND PCYS-pl
                self.__historyfileprocess(step=eachstep_inc, eqincr=eqc,
                                          u11=u11,u22=u22,
                                          ictrl=ict, mode='inplanebiaxialfc')
                self.__pcysplprocess(step=nprobing, shear_ang=90.)
            
        elif mode=='eqeps': #Equivalent strain component case
            """
            Equivalent in-plane strain (balanced strain case)
            """
            self.__historyfileprocess(step=stp, eqincr=eqc,
                                      ictrl=ict, u11=0.5, u22=0.5,
                                      mode='inplanebiaxial')
            pass 
        else:
            print "The passed mode to class vpsc is not among ",
            print "options.\n Please ckeck the code"
            raise IOError

        # single crystal file maker or not
        if fsx!=None: pass
        elif fsx==None:
            for i in range(self.VPSCin.nph):
                sx_maker.cubic(
                    filename=self.VPSCin.ph[i].fsx,
                    #Hardening matrix
                    hii=[1.0,1.0][i], 
                    hij=[1.0,1.0][i], 
                    hb=[1.0,1.0][i],
                    #Voce parameters
                    tau0=[1.02e2,1.02e2][i], tau1=[1.5e2,1.5e2][i], 
                    thet0=[3.8e2,3.8e2][i], thet1=[0.99e2,0.99e2][i],
                    
                    #
                    dependency=['indv','iso'][i],  #'iso', 'indv'

                    #slip-twiningg system
                    n=[ [[1,1,1]], [[1,1,0]] ][i], 
                    b=[ [[1,1,0]], [[1,1,1]] ][i],
                    iopsysx=[0,0][i],
                    )

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

        #Execute the VPSC7.exe file
        #self.run() #self.runflag  :Popen that is called in self.run() method
        os.chdir(self.cwd_old)

        self.datamaster = {} #INITIALIZATION UPON self.datamaster

        #post-process
        #self.pp() 
            #post-execution process:
            #manage data files into data object
        os.chdir(self.cwd_old)
        #plot post-process data 


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
                             ictrl=7, u11=None, u22=None):
        ####
        ####  HISTORY FILE PROCESS 
        ####
        """
        LIST OF MODES : 'unitension', 'unitensionfc',
            'inplanebiaxial', 'bulge', 'simpleshear', 
            'simpleshearfc', 'unitensiony', unitensionz'

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
        
        """
        if mode==None:
            print " You must input a proper mode to __historyfileprocess"
            print "\navailable modes are below"
            print " 'unitension', 'unitensionfc', 'inplanebiaxial'"
            print " 'bulge', 'simpleshear', 'simpleshearfc'"
        self.prcs=self.prcs+self.__loading__(histfile='hist%stemp%i'%(os.sep,
                                                                      self.nhist),
                                             mode=mode,
                                             nstep=step, 
                                             eqincr=eqincr,
                                             ictrl=7,
                                             u11=u11, 
                                             u22=u22)
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
        
    def __pcysplprocess(self, step = 72, shear_ang=90.):
        ####
        ####  POLYCRYSTALLINE YIELD SURFACE on plane stress space
        ####
        step = step 
        shear_ang = shear_ang #angle between E12 and E11-E22 plane
        self.prcs = self.prcs + self.__pcyspl__(nstep = step, shear = shear_ang)
        self.step = self.step + step

    def __rotate(self, ang=1):
        """ make rotate file """
        th = ang*np.pi/180.
        sth = np.sin(th) ; cth = np.cos(th)
        rmat = np.array([[cth,sth,0],[-sth,cth,0],[0,0,1]])
        filename = '%s%s%s%i'%('rot', os.sep, 'r', self.nrot)
        FILE = open(filename,'w')
        FILE.writelines('rotation matrix for polycrystalline aggregate\n')
        for i in range(3):
            for j in range(3):
                FILE.writelines('%9.4f  '%rmat[i][j])
            FILE.writelines('\n')
        FILE.close()
        self.prcs = self.prcs + ['4', filename]
        self.nrot = self.nrot + 1
        
    def getback(self):
        print "\n\n"
        print "***********************************************"
        print "You have now got back to the original directory"
        print "***********************************************\n\n"
        os.chdir(self.cwd_old)
        
    def rmd(self):
        self.getback()
        print " **************************"
        print " Delete the saved directory"
        print "  ~ %9s"%(self.cwd[len(self.cwd)-15:15])
        print " **************************\n"
        shutil.rmtree(self.cwd)

    def run(self, ntrun = 'VPSC7.exe', posixrun='./VPSC7'):
        """
        Runs the VPSC executable.
        Wait until it ends only when globaly flag self.procwait is True
        """
        os.chdir(self.cwd)
        #print "self.cwd = ", self.cwd 
        if self.osname=='nt': p = Popen("%s"%ntrun)    #p = Popen(self.gwd+ntrun)
        elif self.osname=='posix': p = Popen("%s"%posixrun) #p = Popen(posixrun)
        #if self.procwait==True:
        else: pass
        p.wait()
        #self.pp()

        ### deletes crss file
        if self.del_crss==True:
            files = glob.glob('crss*')
            for i in files:os.remove(i)
        os.chdir(self.cwd_old)
        return p


    def running_conditions(self):
        """ governs the running condition """
    
    def write(self):
        """ writiting activities """

    def enviroment(self):
        """ check which is loaded and not """
        """ e.g. plotting is possible? """


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
        v_in(
            nelem=nelem, nph=nph, wph =wph, ishape=ishape,
            vin_name = 'VPSC7.in',
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
    def pp(self):
        """
        Calling post-processings.
        This is done after a modification of the flags of 
        each processes. Flags are given as string. Here,
        they are converted into integer to avoid possible
        overlap with other 'number' which originally 
        intends to be descriptions for each loading.
        
        example)
        Following block is for monotonic loading and its 
        subsequent rigid body rotation.
        -----------
         0
         hist/temp0
         4
         rot/rot0
        -----------
        For this, self.prcs'd be ['0','hist/temp0','4','rot/rot0'].
        Now, it is changed, here in this module, to
        [0, 'hist/temp0', 4 ,'rot/rot0']
        
        This is particularly benficial in the case that you
        have mixed processes such as ['3','5','0','hist/temp0']
        This is origianlly for lankford coefficient with an 
        incremental angle of 5.0 degree. However, this can be
        interpreted as a rigid body rotation.
        
        """

        print "\n\n\n*********************"
        print " POST-PROCESS MODULE "
        print "*********************\n\n"
        print " processes performed listed below "
        for i in range(len(self.prcs)/2):
            print "%1s  %s"%(self.prcs[i*2],self.prcs[i*2+1])


        if any(self.prcs[i*2]=='0' for i in range(len(self.prcs)/2)):
            self.__hist__()
            self.__texture__()
            self.__str_wrk_R__()
        if any(self.prcs[i*2]=='1' for i in range(len(self.prcs)/2)):
            print " You need to work for post-processing for",
            print " hardwired loading using subroutine "
            raw_input( "  type enter >> ")
        if any(self.prcs[i*2]=='2' for i in range(len(self.prcs)/2)):
            self.__pcys__()
        if any(self.prcs[i*2]=='3' for i in range(len(self.prcs)/2)):
            self.__lank__()
        if any(self.prcs[i*2]=='4' for i in range(len(self.prcs)/2)):
            self.__rbr__()
        if any(self.prcs[i*2]=='5' for i in range(len(self.prcs)/2)):
            self.__pcys_plane_stress__()



    def __hist__(self):
        """ hist - post __loading__ boundary condition process """
        print " __hist__ post-process"
    def __activity__(self):
        """ acitivity """
        print " __activity__ post-process"
    def __rbr__(self):
        """ rigid body rotation """
        print " __rbr__ post-process"
    def __texture__(self):
        """
        to be continued after my gate review preparation ends (2011-01-05)
        Now I'm back into this (2011-01-11)
        
        texout is a list containing mytext out from def __texf__
        
        **  texout[i].blocks.gr[j].ODF
        Stands for in the i-th phase's texture file's j-th block's ODF
        An ODF object is a list containing each consistituent grain's 
        three euler angles and OD intensity
        """
        print " in the __texture__ post-process"
        texout = []
        for i in range(self.VPSCin.nph):
            FILE = self.cwd+os.sep+'TEX_PH%i.OUT'%(i+1)
            texout.append(__texf__(filename = FILE))
        self.datamaster['texout'] = texout

        # Guide below 
        # """
        # print "        myvpsc = vpsc_param.vpsc()"
        # print "        myvpsc.pp()"
        # print "        myvpsc.datamaster['texout'][i].blocks.gr[j].ODF"
        # print "            : i phase"
        # print "            : j blocks in i-phase's texture file"
        # """
        
    def __lank__(self):
        """ subroutine pcys """
        print " __lank__ post-process"
        ang, lank = fplot(filename = self.cwd+os.sep+'LANKFORD.OUT',ix=0, iy=2)
        ang, lank = __makenp__(ang, lank)
        self.datamaster['lank'] = np.array([ang,lank])

    def __pcys__(self):
        """ polycrystal yield surface file post-processing """
        print " __pcys__ post-process"
        Sx,Sy = fplot(filename=self.cwd + os.sep+'PCYS.OUT', ix=0, iy=1 )
        self.datamaster['pcys'] = np.array([Sx,Sy])
        pass
    
    def __pcys_plane_stress__(self):
        """ polycrystal yield surface on plane stress plane post-processing """
        print " __pcys_plane_stress__ post-process"
        Sx, Sy = fplot(filename=self.cwd+os.sep+'pcys_pl.out', 
                       ix=0, iy=1, nhead=1)
        """
        Sx=np.array(Sx);Sy=np.array(Sy)
        """
        self.datamaster['pcyspl'] = ([Sx,Sy])
        pass

    def __str_wrk_R__(self):
        """ 
        1. uniaxial stress-strain curve 
        2. and R-value (instantaneouse R-value by differentiation 
        3. plastic- work vs stress

        note that fplot returns values
        assuming that separating block can exist, like the way the
        texture file is saved quite frequently.
        
        AN IMPORTANT REMARK: 
        R valus is one element less than strain or stress list.
        Work was also one element less when it is coming out from
        cummulative trapzoidal method for integration. Later(2010-01-28)
        it is modified to have zero as a prepending element to the rest
        """
        print " __str_wrk_R__ post-process"
        EVM, SVM = fplot(filename=self.cwd+os.sep+'STR_STR.OUT', ix=0, iy=1)
        E11, S11 = fplot(filename=self.cwd+os.sep+'STR_STR.OUT', ix=2, iy=8)
        E22, S22 = fplot(filename=self.cwd+os.sep+'STR_STR.OUT', ix=3, iy=9)
        E33, S33 = fplot(filename=self.cwd+os.sep+'STR_STR.OUT', ix=4, iy=10)
        E23, S23 = fplot(filename=self.cwd+os.sep+'STR_STR.OUT', ix=5, iy=11)
        E13, S13 = fplot(filename=self.cwd+os.sep+'STR_STR.OUT', ix=6, iy=12)
        E12, S12 = fplot(filename=self.cwd+os.sep+'STR_STR.OUT', ix=7, iy=13)
        
        R = __diff__(x=E33, y=E22)

        #wrk = integrate.cumtrapz(y=S11, x=E11)
        #wrk = np.resize((),(6))
        w11 = integrate.cumtrapz(y=S11, x=E11).flatten(1)
        w22 = integrate.cumtrapz(y=S22, x=E22).flatten(1)
        w33 = integrate.cumtrapz(y=S33, x=E33).flatten(1)
        w23 = integrate.cumtrapz(y=S23, x=E23).flatten(1) * 2.
        w13 = integrate.cumtrapz(y=S13, x=S13).flatten(1) * 2.
        w12 = integrate.cumtrapz(y=S12, x=E12).flatten(1) * 2.
        w11 = np.insert(w11,[0],[0])
        w22 = np.insert(w22,[0],[0])
        w33 = np.insert(w33,[0],[0])
        w23 = np.insert(w23,[0],[0])
        w13 = np.insert(w13,[0],[0])
        w12 = np.insert(w12,[0],[0])

        wrk = []

        for i in range(len(w11)):
            wrk.append(w11[i]+w22[i]+w33[i]+w23[i]+w13[i]+w12[i])
        wrk = np.array(wrk)
        """
        tm = [0]
        for i in range(len(wrk)): tm.append(wrk[i])
        wrk = tm; del tm
        """

        EVM,SVM,R = __makenp__(EVM,SVM,R)
        E11,E22,E33,E23,E13,E12 = __makenp__(E11,E22,E33,E23,E13,E12)
        S11,S22,S33,S23,S13,S12 = __makenp__(S11,S22,S33,S23,S13,S12)

        self.datamaster['R'] = R; 
        self.datamaster['wrk'] = wrk; self.datamaster['wrk11'] = w11
        self.datamaster['EVM'] = EVM; self.datamaster['SVM'] = SVM
        self.datamaster['E11'] = E11; self.datamaster['E22'] = E22
        self.datamaster['E33'] = E33; self.datamaster['E23'] = E23
        self.datamaster['E13'] = E13; self.datamaster['E12'] = E12
        self.datamaster['S11'] = S11; self.datamaster['S22'] = S22
        self.datamaster['S33'] = S33; self.datamaster['S23'] = S23
        self.datamaster['S13'] = S13; self.datamaster['S12'] = S12

        print "\n\n"
        print "**************************************************"
        print " Sig-Eps, R and wrk post-process info\n"
        print "    ----    "
        print "   Data has %i number of block(s)"%E11.shape[0]
        print "   You have %i number of histfile processes(s)"%self.nhist
        print "    ----    "
        print "**************************************************\n"


    #----------------------------
    ###
    ### DATAMASTER WRITES
    ###
    #----------------------------
    def __ppwrite__(self): pass

    #----------------------------
    ###
    ### prcs MAKERS 
    ###
    #----------------------------
    def __pcyspl__(self,nstep,shear):
        """ subroutine pcys_pl """
        if os.name=='nt': os.system('cls')
        elif os.name=='posix': os.system('clear')
        print " *** POLYCRYSTALLINE YIELD SURFACE",
        print " ON PLANE STRESS SPACE *** "
        pcrs = ['5','1 2 %i %i'%(nstep,shear)]
        return pcrs

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
        uniaxial or multiaxial loadings 
        returns prcs : process list
        available mode : unitension, unitensionfc, bulge, 
                         inplanebixial, simpleshear
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

        elif mode=='unitensiony': 
            iudot =   [[ 0,   0,   0],
                       [ 1,   1,   0],
                       [ 1,   1,   0]]
            udot  =   [[-0.5,   0,   0],
                       [0   ,   1,   0],
                       [0   ,   0,-0.5]]
            iscau = [1,0,1,1,1,1]
            scauchy = [0,0,0,0,0,0]
        elif mode=='unitensionz': 
            iudot =   [[ 0,   0,   0],
                       [ 1,   0,   0],
                       [ 1,   1,   1]]
            udot  =   [[-0.5,   0,   0],
                       [0   ,-0.5,   0],
                       [0   ,   0,   1]]
            iscau = [1,1,0,1,1,1]
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
            raise IOError

        ## hist file making
        vpsc_in.histmaker(
            filename = histfile,
            nstep = nstep, eqincr = eqincr, ictrl = ictrl,
            iudot = iudot, udot = udot, iscau = iscau,
            scauchy = scauchy)
        return prcs

    def __histvar__(self):
        """ loading through a subroutine """

    ### RUN RESULTS ###
    def crss(self):
        """ infor on crss evolution """

    def phases(self):
        """ info on phases """

    def grainmorph(self):
        """ grain morphology """

    def texture(self):
        """ Governing texture """

    def sx(self):
        """ Governing single crystal """    

    def show_str_str(self):
        """
        Shows 'STR_STR.OUT' file
        """
        if os.name=='nt': 
            cmd = 'type'
            os.system('%s %s%s%s'%(cmd, self.cwd,os.sep,'STR_STR.OUT'))
        elif os.name=='posix':
            cmd = 'cat'
            os.system('%s .%s%s%s'%(cmd, self.cwd.split(os.getcwd())[1],
                                   os.sep,'STR_STR.OUT'))
        pass

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

    def show_activity(self):
        """  """
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
            os.system('%s %s%s%s'%(cmd, self.cwd, os.sep, 'VPSC7.in'))
        elif os.name=='posix':
            cmd ='cat'
            os.system('%s .%s%s%s'%(cmd, 
                                    self.cwd.split(os.getcwd())[1],
                                    os.sep, 'VPSC7.in'))

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

    ### POST-EXECUTABLE PROCESS ###
    def print_info(self):
        """ Print information of current VPSC ojective"""
    
    def plot(self):
        """ 
        Plot a data in a suitable form 
    
        Lankford coefficients
        """

    ### SOME HANDY MODULES EXPECTED ###
    def __seekngr__(self, filename):
        """ finds ngrain in the textfile """
        FILE = open(filename, 'r')
        lines = FILE.readlines()[3]
        return int(lines.split()[1]) #number of grain
        #self.ngrain = self.__seekngr__(self.ftex)

    def record(self):
        """ write details into a global log file """
    
    def add(self):
        """  add an element into a data object  """
    
    def subtract(self):
        """ Remove an element from an data data """
    
    def directory_control(self):
        """ current working directory control """
    
    def __getnewdirectory__(self):
        """
        suggests and returns a non-existing directory name
        """
        i = 0 #trial number
    
        while True:
            temp = 'vpsc_cwd__'+('%i'%(i)).zfill(5)
            if len(glob.glob(temp))==0:
                return temp
            else: i = i + 1
    
    def __mkcopy__(self, dst=None):
        """ 
        Copies everything into the self.cwd,
        global current-working-directory.

        Arguments
          dst = None
          directories and filenames that are copied and moved are hardwired 
        """
        os.mkdir(dst)
        ## I have to use make directory rather than copy them
        ## and consider to move only what is necessary!! commented : (2011-01-27)
        directories = ['sx', 'texture', 'rot', 'hist']
        F = ['VPSC7','VPSC7.exe','p.exe','VPSC7.in']
        #dst = os.getcwd()+'/'+dst
        for dic in directories: 
            if dic=='hist':
                files=glob.glob(dic+'%s*'%os.sep)
            else: files = glob.glob(dic+'%s*.*'%os.sep)
            #making direcotries
            os.mkdir(dst+os.sep+dic)
            for i in range(len(files)):
                shutil.copyfile(files[i], dst+os.sep+files[i])
        for name in F:
            shutil.copyfile(name, dst+os.sep+name)



    def append(self,filename, contents):
        """ append contents to the file with the filename """
        FILE = open(filename,'r')
        lines = FILE.readlines()
        FILE.close()
        FILE = open(filename,'w')
        for i in lines:
            if len(lines)<3: pass
            else: FILE.writelines(lines[i]+'\n')
        #append the contents !
        FILE.write(contents)
        FILE.close()


    ### preliminaries ###
    " default file making"
    " ODF management "
    " random grains "
    " RBR file making "


    ### to be continued upon issues like ... ###
    " POSTMORT  management "
    " fileaxes maker (shape file)"


def angle_to_u(angle):
    angle = angle * 3.141592 /180.
    return math.cos(angle), math.sin(angle)

def __bb__(nstep = 0):
    nstep = nstep
    ini_ang = -45.
    fil_ang =  135.
    span_angle = fil_ang - ini_ang
    ang_incr = span_angle / (nstep-1)
    angles = []
    for i in range(nstep):
        curangle = ini_ang + ang_incr * i
        u = angle_to_u(curangle)
        currne_vpsc_run = vpsc(u11=u[0],u22=u[1])
    return curren_vpsc_run



##
## Example #1
## pseudo yield surface (based on plastic work) -- monotonic proportional loading
##
class wrk_cnt:
    """
    Plastic work equivalent contour class
    -->  A run is consisted of many inplane biaxials 
        having u11 u22 non zero among other
        velocity gradient components. In this way, 
        it mimics the in-plane biaxial experiments.
         One running also include an uniaxial for 
         such case that normalization contour is
         necessarily to be plotted.

    >>> import vpsc_param (importing vpsc_param.py python script)
    >>> mywrk = vspc_param.wrk_cnt(ini_ang=-30, fil_ang=120, 
                                   nstep=100, ict=7, eqc=0.005, stp=50,
                                   #nstep == number of loadings
                                   texture='texture/00500.cmb',
                                   interaction = 3,
                                   #mode can be probing or RDTD(rdtd)
                                   mode='probing', 
                                   )

            # mode='rdtd' is moved into a new whole class
    >>> mywrk = vpsc_param.wrk_cnt(ict=7, eqc=0.005, stp=50,
                                   texture='texture/00500.cmb',
                                   fsx='sx/Neil.sx',
                                   mode='rdtd',
                                   )


        ---> This may remove the existing files so that 
             overlapping can be prevented

    >>> mywrk.run()
    >>> mywrk.pp()  
           -- Now mywrk.work_contour() is in 'pp' module
                #You'll have wrk_cntr_%i type files
    >>> mywrk.work_contour()
    >>> mywrk.plwrk(self, *args)  --> mywrk.plwrk(1.,2.,4., ... ) will work!!
              >>> mywrk.__plwrk__(self, wrk = 0.1)
              >>> mywrk.__plwrk__(self, wrk = 0.2) .. and so on

    >>> vpsc_param.plt_plwrk(ifig)  --> plot the written plwrk_* files

    >>> vpsc_param.norm_plwrk(inorm=0, ifig=None):
         --> Normalize the resulting friles from plwrk
            by the loading along jobs[inorm]
    """
    def __init__(self, ini_ang=-180., fil_ang=180., nstep=2,  #angular terms
                 stp=100, ict=7, eqc=0.0025,  #Strain increment..
                 
                 ##Texture and single crystal file
                 texture=None, fsx=None, 
                 mode='inplanebiaxialfc', 
                 interaction=3,#neff
                 norm=True,
                 ):
        """
        jobs are usually biaxial strain rate governed boundary conditions.
        Special case is only consisted of uniaxial and transverse
        """
        self.jobs = []
        ##----------
        ## RDTD mode 
        ##----------
        if mode=='rdtd' or mode=='RDTD':
            """
            The last job will be producing the reference values during the 
            post-execution process so as to normalize stress values.
            """
            #for TD loading
            print "***********"
            print " TD tension"
            cjob = vpsc(mode='TD', stp=stp, ict=ict, eqc=eqc,
                        texture=texture, fsx=fsx,
                        interaction=interaction)
            self.jobs.append(cjob)
            #for RD loading
            print " RD tension"
            cjob = vpsc(mode='RD', stp=stp, ict=ict, eqc=eqc, 
                        texture=texture, fsx=fsx,
                        interaction=interaction)
            self.jobs.append(cjob)
            print "***********"
            #Some room left behind for following addition.

        ##--------------------
        ## inplanebiaxial mode
        ##--------------------
        else:
            span_angle = fil_ang - ini_ang
            ang_incr = span_angle / (nstep-1)
            angles = []
            for i in range(nstep):
                if i==0:  #for the first run
                    print "**********"
                    print "Job number"
                print "  %s/%i"%(str(i+1).zfill(3),nstep)
                # for the last run
                if i==nstep-1: print "**********"
                curangle = ini_ang + ang_incr * i
                u = angle_to_u(curangle)
                ## It is advised to put texture argument as an 
                ## end-user wants. Most probably, the texture file 
                ## gives a hardwired one if not given.
                cjob = vpsc(u11=u[0], u22=u[1], 
                            mode=mode,
                            texture=texture,  #texture=['500.cmb'] ... 
                            ict=ict, eqc=eqc,stp=stp,fsx=fsx,
                            interaction=interaction)
                self.jobs.append(cjob)

            self.norm = norm
            if norm==True:
                ## The last job will be uniaxial along RD direction
                ## from which the contours will be normalized
                cjob = vpsc(mode='RD', stp=stp, ict=ict, eqc=eqc, 
                            texture=texture,fsx=fsx,
                            interaction=interaction)
                self.jobs.append(cjob)

    def run(self):
        """ Execute the VPSC jobs .. """
        for i in range(len(self.jobs)):
            self.jobs[i].run()

    def pp(self):
        """
        post-execution process. It incluse self.work_contour
        """
        for i in range(len(self.jobs)):
            self.jobs[i].pp()
        self.work_contour()

    def work_contour(self):
        """
        Writes down the work_contour files, 'plwrk_???' files
        """
        for ijob in range(len(self.jobs)): ## Note the last is reference uni-RD
            if ijob ==len(self.jobs)-1 and self.norm==True: 
                FILE = open('wrk_ref','w')
                FILE.writelines('** plastic work contour job #%i ** \n'%(ijob))
            else: 
                FILE = open('wrk_cntr_%s'%str(ijob).zfill(3),'w')  # 
                FILE.writelines('** plastic wrk_ cntr reference loading **\n')
            # S11, S22, work
            FILE.writelines(' %8s  %8s  %8s  %8s\n'%('S11','S22',
                                                     'wrk11','wrk'))
            e11 = self.jobs[ijob].datamaster['E11'].flatten(1)
            e22 = self.jobs[ijob].datamaster['E22'].flatten(1)
            s11 = self.jobs[ijob].datamaster['S11'].flatten(1)
            s22 = self.jobs[ijob].datamaster['S22'].flatten(1)            
            wrkDAT    = self.jobs[ijob].datamaster['wrk'].flatten(1)
            wrk11DAT  = self.jobs[ijob].datamaster['wrk11'].flatten(1)
            FILE.writelines('--  0 th black \n')
            for i in range(len(e11)):
                FILE.writelines(' %8.4e  '%s11[i])
                FILE.writelines(' %8.4e  '%s22[i])
                try:
                    FILE.writelines(' %8.4e  '%wrk11DAT[i])
                except:
                    print "Error at job #%i"%ijob
                    print "wrk11DAT"
                    print "len(wrk11DAT):%i"%len(wrk11DAT)
                    print "***********"
                    print wrk11DAT
                    print "***********\n"
                    print "e11"
                    print "len(e11)=%i"%len(e11)
                    print e11
                    raw_input()

                    
                
                FILE.writelines(' %8.4e\n'%wrkDAT[i])
            FILE.close()
                
            
            # """  ## this is before work as flattened as it is created.
            #      ## Then it was modified to flatten.
            # eps_block = self.jobs[ijob].datamaster['E11']
            # E11DAT    = self.jobs[ijob].datamaster['E11']
            # S11DAT    = self.jobs[ijob].datamaster['S11']
            # S22DAT    = self.jobs[ijob].datamaster['S22']
            # wrkDAT    = self.jobs[ijob].datamaster['wrk']
            # wrk11DAT  = self.jobs[ijob].datamaster['wrk11']
            

            # for ib in range(len(eps_block)):  #loop over block
            #     FILE.writelines('---  %s th block\n'%(ib))
            #     for IS in range(len(eps_block[ib])): #if IS ==0, no wrk
            #         if IS==0:
            #             FILE.writelines(' %8.4e  '%(S11DAT[ib][IS]))
            #             FILE.writelines(' %8.4e  '%(S22DAT[ib][IS]))
            #             FILE.writelines(' %8.4e  '%(0.))
            #             FILE.writelines(' %8.4e\n'%(0.))
            #         else:
            #             File.writelines(' %8.4e  '%(S11DAT[ib][IS]))
            #             FILE.writelines(' %8.4e  '%(S22DAT[ib][IS]))
            #             FILE.writelines(' %8.4e  '%(wrk11DAT[ib][IS-1]))
            #             FILE.writelines(' %8.4e\n'%(wrkDAT[ib][IS-1]))
            # """
            FILE.close()


        # """ ## Old way of appending uniaxial but now it is a part of the above loop
        # ## The uniaxial loading(the reference) case will have difference name
        # FILE = open('wrk_ref','w')
        # FILE.writelines('** plastic work contour reference loading ** \n')
        # FILE.writelines(' %8s  %8s  %8s  %8s\n'%('S11','S22','wrk11','wrk'))
        # eps_block = self.jobs[-1].datamaster['E11']
        # E11DAT    = self.jobs[-1].datamaster['E11']
        # S11DAT    = self.jobs[-1].datamaster['S11']
        # S22DAT    = self.jobs[-1].datamaster['S22']
        # wrkDAT    = self.jobs[-1].datamaster['wrk']
        # wrk11DAT  = self.jobs[-1].datamaster['wrk11']
        # for ib in range(len(eps_block)):  #loop over block
        #     FILE.writelines('---  %s th block\n'%(ib))
        #     for IS in range(len(eps_block[ib])): #if IS ==0, no wrk
        #         if IS==0:
        #             FILE.writelines(' %8.4e  '%(S11DAT[ib][IS]))
        #             FILE.writelines(' %8.4e  '%(S22DAT[ib][IS]))
        #             FILE.writelines(' %8.4e  '%(0.))
        #             FILE.writelines(' %8.4e\n'%(0.))
        #         else:
        #             FILE.writelines(' %8.4e  '%(S11DAT[ib][IS]))
        #             FILE.writelines(' %8.4e  '%(S22DAT[ib][IS]))
        #             FILE.writelines(' %8.4e  '%(wrk11DAT[ib][IS-1]))
        #             FILE.writelines(' %8.4e\n'%(wrkDAT[ib][IS-1]))
        # FILE.close()
        # """

    def __rm__(self):
        """
        Remove exisiting self class related files "plwrk_*" .. and so on
        """
        FILES = glob.glob('wrk_cntr_*')
        print 'FILES ARE BELOW'
        print FILES
        raw_input('Press enter to proceed >>> ')
        for i in FILES: os.remove(i)
        FILES = glob.glob('plwrk_*')
        print 'FILES are below'
        print FILES
        raw_input('Press enter to proceed >>> ')
        for i in FILES: os.remove(i)
        FILES = glob.glob('norm_works*')
        print 'FILES are below'
        print FILES
        raw_input('Press enter to proceed >>> ')
        for i in FILES: os.remove(i)
    
    def plwrk(self, *args):
        """ 
        Calls self.__plwrk__ as many times as the number of *args 
        """
        for i in range(len(args)):
            self.__plwrk__(wrk=float(args[i]))

    def __plwrk__(self, wrk):
        """
        Work contour levels, given as an argument wrk, output into a file
        """
        if wrk==0 :
            print "Work level should be larger than 0."
            return -1
        FILE = open('plwrk_%s'%str(wrk).zfill(6), 'w')
        FILE.writelines(' ** plastic work of %8.3f\n'%(wrk))
        FILE.writelines(' %8s  %9s  %9s\n'%('JOB#','S11','S22'))
        FILES = glob.glob('wrk_cntr_*')
        for i in range(len(self.jobs)): #for wrk_ref addition
            if i==len(self.jobs)-1 and self.norm==True: 
                filename = 'wrk_ref'     # The last job is reference uniaxial loading along RD
                print "Now you are doing for wrk_ref"
            else: filename = 'wrk_cntr_%s'%str(i).zfill(3) 

            S11, S22 = fplot(filename=filename,
                             ix=0, iy=1, nhead=3)
            W11, W   = fplot(filename=filename,
                             ix=2, iy=3, nhead=3)
            S11, S22, W11, W = __makenp__(S11,S22,W11,W)
            try: len(S11); len(S22); len(W11)
            except TypeError: pass
            else:
                for ib in range(len(S11)):
                    try:
                        S11_, S22_ = __interpolate__(x=S11[ib],y=S22[ib],
                                                     z=W11[ib], value=wrk)
                        S11_new, S22_new = __interpolate__(x=S11[ib],y=S22[ib],
                                                           z=W[ib], value=wrk)
                    except ValueError: 
                        print " Value Error occured"
                        raw_input("press Enter >>")
                        pass
                    else:
                        if [S11_new, S22_new ] == [0,0]:
                            print 'not proper results'
                            raw_input('press enter to proceed >>')
                        FILE.writelines(' %8i  %9.3f  %9.3f\n'%(i,S11_new, S22_new))

def plt_plwrk(ifig):
    """
    Plots all glob.glob('plwrk_*') files in the working directory
    """
    files = glob.glob('plwrk_*') 
    fig = plt.figure(ifig)
    ax = fig.add_subplot(111)
    for i in range(len(files)):
        x, y = fplot(filename=files[i], ix=1, iy=2,
                     ifig=ifig, nhead=2, idot=True, marker='--')

def norm_plwrk(inorm=0,ifig=None):
    """
    Taking "plwrk_*" files and calculate the normalized plastic work contour
    normalized by a particular strain path given whose index is given by inorm

    >>> vpsc_param.norm_plwrk(inorm=0,  ifig=1)
    >>> vpsc_param.norm_plwrk(inorm=-1, ifig=1) 
                      --> Normalized by the last job submitted (supposed to be uniaxial loading)
    """
    files = glob.glob('plwrk_*')
    #if os.path.exists('wrk_ref'): files.append('wrk_ref')
    print 'files are as below \n', files
    xmast = []; ymast = []
    for i in range(len(files)):
        x, y = fplot(filename = files[i], ix=1, iy=2,
                     ifig=None, nhead=2)
        x = x[0]; y = y[0]
        fact = x[inorm]
        if inorm==-1:  # inorm==-1 is 
                       # the element in which uniaxial tension along RD is done
            x=x[0:len(x)-1]; y=y[0:len(y)-1]
        x = np.array(x); y = np.array(y)
        x = x/fact; y = y/fact

        fout('norm_work%s'%str(i).zfill(3), x, y)
        xmast.append(x); ymast.append(y)
        
    if ifig==None: return xmast, ymast
    else:
        fig = plt.figure(ifig)
        ax=fig.add_subplot(111)
        for i in range(len(xmast)):
            ax.plot(xmast[i],ymast[i])
    ax.set_aspect('equal')
    ax.set_ylim(0.,);ax.set_xlim(0.,)
    ax.set_xticks((0.,0.5,1.0));ax.set_yticks((0.,0.5,1.0))
    ax.grid(True)





###
###  EXAMPLE #2  multipath at once.
###
###  This compares the flow stresses of different monotonic loading paths.
###

class multipath:
    """
    Uniaxial loading along RD, TD, DD as well as bulge(BU) 
    and simple shear(SS) that are all, so far my knowledge is 
    correct, experimentally measurable.

    These are simulated simulataneously so that how the modeled 
    polycrystalline shows its anisotropy in terms of its 
    respective loadings all together. Because these are all
    experimentally obtainable basic mechanical tests, 
    predictions can be easily corelated to experiments.
    This may provide further insights upon characteristitics 
    of given polycrystalline materials. All the kinds of tests
    are now available in my lab, Materials Mechanics Lab, GIFT, 
    POSTECH, except simple shear (2011-Jan).

    Hyun Jin is now working on the simple shear device, which
    is expected to be soon operative.

    Again, it is also done for effective stress vs work curves
    (2011-Jan-28). Then all loadings are normalized by the
    result of RD to see how equivalent work contour is evolving.
    The equivalent work contour is more effective than yield
    surface since the evolution of yield surface is not 
    measurable data as far as polycrystal's yield surface is
    path dependent. (This may be observable from Mg alloys
    which goes through drastic anisotropic changes due to
    mechanical twinning during deformation -  Refer to the work
    with Dr. Steglich during his stay in MML,
    2009 spring to 2010 Autumn)

    ** Multipaths available
    -Uniaxial tension along RD, DD, and TD (which corresponds 
      respectively to 1,2 and 45 angle from 1 on the 1-2 plane)
    -Simple shear along RD
    -Bulge state (actually its uniaxial compression along ND,
                  3 axis)
    
    ** Typical way to use
    >>> import vpsc_param
    >>> mypath = vpsc_param.multipath(texture='texture/00500.cmb')
    >>> mypath.run()
           . . .
    >>> mypath.pp()
    >>> mypath.flowstress()
    >>> mypath.eq_work(*args) --> Equivalent-work contour
         ********************************************
         one can do the below or not. Once it is done 
         the output is the result after normalization 
         by one particular loading's stress
         ********************************************
         >>> mypath.eq_work_norm()  --> Normalization of the work

    ** Output filenames are 
      1. flowstress_*    from flowstress()
      2. flowwork_*      from eq_work()
         The wild card symbol '*' indicates the kinds of the loading 
                                  : (RD, DD, TD, SS, BU)

    """
    def __init__(self,texture='texture/00500.cmb',
                 stp=3, eqc=0.005, ict=7 ):
        """
        Assigns multiple jobs into myjobs object.
        End-user can add more jobs to the down below block.
        Arugments
        texture ='texture/00500.cmb'
        * args  : modes for class vpsc  'RD','DD','TD','SS' ...
        """
        #flag indicating whether the pp was performed
        self.waspp = False 
        
        ## hardwired modes are used.
        self.modes = ['RD','DD','TD','SS','BU'] 
        self.myjobs = []

        ## Shoots Information on to the screen
        print "\n***************************************************"
        print "You have total %i number",
        print " of different loading paths"%len(self.modes)
        print "%3s %15s"%('id','mode name')
        for i in range(len(self.modes)):
            print "%3i %15s"%(i,self.modes[i])
        print "***************************************************"
        
        for mode in self.modes:
            self.myjobs.append(vpsc(mode=mode, 
                                    texture=texture, 
                                    stp=stp, eqc=eqc, ict=ict
                                    )
                               )
        """
        print "\n Guide .. "
        print self.__doc__
        print "\n"
        """
        pass
        
    def run(self):
        """ Deform the polycrystalline aggregate!!! """
        for i in range(len(self.myjobs)):
            self.myjobs[i].run()
        pass
            
    def pp(self):
        """ post-execution process """
        for i in range(len(self.myjobs)):
            self.myjobs[i].pp()
        self.waspp = True
        pass

    def flowstress(self):
        """ 
        1. Extract flow stress level 
        2. Writes them down to flowstress_* file where '*'
            wild card indicates the kind of loading, such 
            as RD or TD.
        
        If your polycrystal were uniaxialy deformed 
        along other directions RD, you have to consider 
        you have rotated the aggregate before so doing. 

        """
        if self.waspp ==False:
            print "Post-Execution process was not performed"
            return -1

        if self.modes !=['RD','DD','TD','SS','BU']: 
            print "self.flowstress is only valid when the modes"
            print "are as default, which is",
            print " ['RD','DD','TD','SS','BU']"
            return -1
        
        ## this is under assumption that your loading were 
        ## RD DD TD SS BU, which is the given default

        RD = self.myjobs[0]; DD = self.myjobs[1]
        TD = self.myjobs[2]; SS = self.myjobs[3]
        BU = self.myjobs[4]

        ### ASSIGNS GLOBAL VARIABLES for further processes later on.
        ### BULGE and simple shear DATA are absolutes
        self.E = []; self.S = []
        self.E.append(RD.datamaster['E11'][0])
        self.E.append(DD.datamaster['E11'][0])
        self.E.append(TD.datamaster['E11'][0])
        self.E.append(abs(SS.datamaster['E12'][0]))
        self.E.append(abs(BU.datamaster['E33'][0]))

        self.S.append(RD.datamaster['S11'][0])
        self.S.append(DD.datamaster['S11'][0])
        self.S.append(TD.datamaster['S11'][0])
        self.S.append(abs(SS.datamaster['S12'][0]))
        self.S.append(abs(BU.datamaster['S33'][0]))

        ## Writes the effective flow stress-strain curves 
        self.fout(filename='flowstress_RD',
                  x=RD.datamaster['E11'][0],
                  y=RD.datamaster['S11'][0])
        self.fout(filename='flowstress_DD',
                  x=DD.datamaster['E11'][0],
                  y=DD.datamaster['S11'][0])
        self.fout(filename='flowstress_TD',
                  x=TD.datamaster['E11'][0],
                  y=TD.datamaster['S11'][0])
        self.fout(filename='flowstress_SS',
                  x=abs(SS.datamaster['E12'][0]),
                  y=abs(SS.datamaster['S12'][0]))
        self.fout(filename='flowstress_BU',
                  x=abs(BU.datamaster['E33'][0]),
                  y=abs(BU.datamaster['S33'][0]))

    def eq_work(self,*args): #norma
        """ 
        Calculates equivalent work contours. One can keep 
        proceeding to normalizaed data by invoking 
        self.eq_work_norm() 
        ('self' here refers to your own class name)

        *args refers to the plastic work levels that one'd
        like to get using interpolating the incremental 
        work results for each loading

        This module is designed to be compatible with 
        general self.modes case meaning that any arbitrary
        combination of modes is possible.
        """
        """
        print "\nWork level given is as follow"
        for i in range(len(args)):
            print "%5.2f "%args[i],
        """
            
        ## Plastic work level global variable
        self.wrk_levels = map(float,args)
        ##

        if self.waspp==False:
            print "Post-Execution process was not performed"
            return -1

        RD = self.myjobs[0]; DD = self.myjobs[1]
        TD = self.myjobs[2]; SS = self.myjobs[3]
        BU = self.myjobs[4]

        #List in which all loadings' work is saved.
        self.wrk = []
        self.wrk.append(RD.datamaster['wrk'].flatten())
        self.wrk.append(DD.datamaster['wrk'].flatten())
        self.wrk.append(TD.datamaster['wrk'].flatten())
        self.wrk.append(SS.datamaster['wrk'].flatten())
        self.wrk.append(BU.datamaster['wrk'].flatten())

        """
        print "self.S is as below"
        print self.S
        raw_input()
        """
        
        #effective stress vs. work will be the curve's axes
        #self.S #---> This is the effective stress
               # and this is assumed already obtained
               # in flowstress module
        self.sv = np.zeros((len(args), len(self.modes)))

        ### interpolation for each 
        for iwrk in range(len(args)): # Loop over working contour
            for ind in range(len(self.modes)): 
                # Loop over loading(RD, DD, TD, SS, BU)
                s = self.__interpolate__(s=self.S[ind], 
                                         w=self.wrk[ind], 
                                         wrk=args[iwrk])
                #print 'interpolated s=', s
                self.sv[iwrk][ind] = s #This is stress value
                # RD, DD, TD, SS, ND(bulge or BU) compress

        self.eq_work_norm(ind=0) 
        #hardwired to have RD as reference

        ## writing 
        #self.__writes_wrk_cnt__()

    def eq_work_norm(self, ind=0): #ind=0,1,2,3,4 (RD,DD,TD,SS,ND)
        """
        Normalizes equivalent work contours based on one particular 
        loading among others. RD is preferred as a convention.
        
        *Argument: ind=0
        """
        for iwrk in range(len(self.sv)):
            ref = self.sv[iwrk][ind]
            for iloading in range(len(self.sv[iwrk])):
                self.sv[iwrk][iloading] = self.sv[iwrk][iloading]/ref
        pass

    def __writes_wrk_str_curves__(self,prep='flowwrk'):
        """
        Writes each modes' wrk-effective stress curves
        """
        filenames = []

        for i in range(len(self.S)):
            ctem = "%s_%s"%(prep,self.modes[i]) 
            filenames.append(ctem)
        for i in range(len(self.S)):
            cfile = open(filenames[i], 'w')
            cfile.writelines('** work - effective stress **\n')
            cfile.writelines('%5s  %8s\n'%('wrk','stress'))
            for j in range(len(self.S[i])):
                cfile.writelines('%5f  %8f\n'%(self.wrk[i][j],
                                               self.S[i][j]
                                               )
                                 )
        pass

    def __writes_interpolated_wrk_cnt__(self,prep='interwrk'):
        """
        Writes self.sv
        Argument : prep : prependix to the output filenames
        """
        ## To check the possibility to become buggy again
        if self.modes!=['RD','DD','TD','SS','BU']:
            print "The self.modes are as default"
            print "__writes_wrk_cnt__ will return -1"
            raw_input("Press enter to proceed >>> ")

        filenames = []
        for i in range(len(self.wrk_levels)):
            # File for each work_levels
            ctem = "%s_%s"%(prep, 
                            str(round(self.wrk_levels[i], 3)
                                ).zfill(5))
            filenames.append(ctem)
        FILE = open(prep,'w')
        FILE.writelines('*** work_contour output ***\n')
        for i in range(len(self.modes)):
            FILE.writelines('%10s  '%self.modes[i])
        FILE.writelines('\n')
        for i in range(len(self.wrk_levels)):
            # File for each work levels
            cfile = open(filenames[i], 'w')
            cfile.writelines('***   wrk_cnt output   ***\n')
            cfile.writelines('%5s  %8s\n'%('modes','Stress'))
            for j in range(len(self.modes)):
                cfile.writelines('%5s  %8f\n'%(self.modes[j],
                                               self.sv[i][j])
                                 )
                FILE.writelines('%10s  '%self.sv[i][j])
            FILE.writelines('\n')


        """
        for i in range(len(self.sv)):  # Work contour
            for j in range(len(self.sv[i])):  # Kinds of loading
                fout(filename=filenames[j],)
                pass
        """
        pass

    def __interpolate__(self,s,w,wrk):
        """ 
        Interpolates the stress based on the work level given 
        as wrk arugment
        Arguments: s(list), w(list), wrk(constant)
        """
        
        #print 'Given s and w has the shape of below'
        zmn = min(w); zmx = max(w)
        #print ' zmn, zmx = ', zmn, zmx
        if wrk < zmn: 
            print "Given value is below the minimum of work"
        elif wrk > zmx: 
            print "Given value is higher than the maxmum of work"
            return 0.,0.
        item = 0
        while True:
            if w[item] > wrk : break
            else: item += 1
        w0 = w[item - 1]
        w1 = w[item]
        s0 = s[item - 1]
        s1 = s[item]
        slope = (s1 - s0)/(w1 - w0)
        """
        print 'type of slope = ', type(slope)
        print 'slope =', slope
        print 'type of w = ', type(w)
        print 'w = ', w
        print 'w0, w1, s0,s1 = ', w0,w1,s0,s1
        """
        #-------------------------------------
        #Returns the interpolated stress level
        return slope * (wrk - w0) + s0 
        #-------------------------------------    
        pass

    def fout(self,filename, x, y):
        """
        Given the filename, x and y is written\n
        Arguments: filename, x, y
        """
        FILE = open(filename,'w')
        FILE.writelines('** FLOW STRESS OF TENSION ALONG ROLLING DIRECTION **\n')
        for i in range(len(x)):
            FILE.writelines('%11.3e  %11.3e \n'%(x[i],y[i]))
        FILE.close()
        
    def plot(self,ifig=1):
        """ 
        Plots the stress-strain curves in a figure
        Argument: ifig = 1 
        """
        fig = plt.figure(ifig)
        ax = fig.add_subplot(111)
        for i in range(len(self.E)):
            plt.plot(self.E[i], self.S[i],'o')

def fout(filename, *args):
    """ 
    filename and *args
    *args
      Each argument will form a column in the file 
      whose name is the given filename argument
    vpsc_param.fout(filename='temp',[3,2,3],[2,1,3])

    *** FOUT is now being used many other modules. This turns out to be
       the most useful modules among others

    This module is now being used many places in vpsc_param.py script,
    which was not the intention when it was written. 
    
    """
    FILE = open(filename,'w')
    FILE.writelines("** writes values in order, vpsc_param.fout module **\n")
    """
    for i in range(len(args)):
        for j in range(len(args[i])):
            FILE.writelines('%11.3e  '%args[i][j])
        FILE.writelines('\n')
    """
    # find the maximum len return from *arguments
    mx = []
    for i in range(len(args)):
        mx.append(len(args[i]))
    mx = max(mx)

    for j in range(mx):
        for i in range(len(args)):
            try:
                FILE.writelines('%11.3e  '%args[i][j])
            except:
                FILE.writelines('%10s   '%'')
        FILE.writelines('\n')
    FILE.close()
                
def plt_flow(ifig=1):
    """
    Plot the flowstress_* files using fplot from vpsc_param.py
    """
    files = glob.glob('flowstress_*')
    print '\nfiles:',
    print files
    for i in files:
        fplot(filename=i, ix=0, iy=1, ifig=ifig, idot=True)

###
###  EXAMPLE #3  LANKFORD PROBING PROFILE with respect to # OF GRAINS
###

class r_ngr:
    """
    LANKFORD PROBING WITH RESPECT TO THE # of grains (DIFFERENT TEXTURE FILES)

    >>> myr_ngr = vpsc_param.r_ngr (ftex='texture/00500.cmb', dang=1.5)
    >>> myr_ngr.run()
    >>> myr_ngr.pp()
    >>> myr_ngr.write()

    -- where you can plot the output files results type below --
    vpsc_param.lankfplot(ifig=5)
    """
    def __init__(self, ftex, dang=1):
        self.dang = dang
        self.myjobs = []
        self.R = []; self.x = []
        for i in range(len(ftex)):
            print 'current textfile: %s'%ftex[i]
            self.myjobs.append(vpsc(texture=ftex[i],
                                    mode='lankf', stp=dang))
    def run(self):
        for i in range(len(self.myjobs)):
            self.myjobs[i].run()
    def pp(self):
        for i in range(len(self.myjobs)):
            self.myjobs[i].pp()
        self.prof()
        self.write()
    def prof(self):
        x = np.linspace(0,90,90./self.dang+1)
        x0 = np.linspace(90.,0.,90/self.dang+1)
        x = np.array(( x, x0))
        x = np.array(x)
        total_step = x.shape[0] * x.shape[1]
        x = x.reshape((total_step,))
        for i in range(len(self.myjobs)):
            R = self.myjobs[i].datamaster['lank'][1]
            R = np.array(R)
            self.x.append(x)
            self.R.append(R.reshape((total_step,)))
    def write(self):
        for i in range(len(self.myjobs)):
            FILENAME = 'LANK%s'%str(i).zfill(2)
            print "LANKFORD PROBING: '%s' IS BEING WRITTEN DOWN"%(FILENAME)
            fout(FILENAME, self.x[i],self.R[i])

def lankfplot(ifig=None):
    """
    lankford plot
    arguments : ifig = None
    """
    FILES = glob.glob('LANK*')
    ANGLES = []; R = []
    for i in FILES:
        ang, l = fplot(filename=i, ix=0,iy=1,nhead=1)
        ang = np.array(ang[0]); l = np.array(l[0])
        ANGLES.append(ang); R.append(l)
    if ifig==None:
        return ANGLES, R
    else:
        fig = plt.figure(ifig)
        ax = fig.add_subplot(1,1,1)
        for i in range(len(ANGLES)):
            plt.plot(ANGLES[i], R[i])

###
###  EXAMPLE #4  LANKFORD PROBING PROFILE with respect to # OF GRAINS
###
###  This compares the flow stresses of different monotonic loading paths.
###
class r_text:
    """
    Rotate texture file to fit the texture in lankford probing symmetry
    myr = r_text(ftex='texture/00500.cmb', dang=45, rang=5.)
            myr.run()--> It's in __init__ now
    myr.pp()
    myr.plot()

    myr.rotate(range=??)  --> rotate the texture and run the simulation
    myr.pp()
    myr.plot()
    myr.write(filename='temp')
      .
      .
      .

    """
    def __init__(self,ftex,fsx=None,dang=1,rang=0.,ifig=3,interaction=3):
        self.ifig=ifig
        self.dang=dang
        self.myjob = vpsc(texture=ftex,mode='lankf',
                          stp=dang,fsx=fsx,interaction =interation )
        self.history=rang
        self.rotate(rang=rang)
        #plt.figure(ifig)
        #plt.plot([45,45],[-0.1,4],'--',color='gray')
        self.run()

    def rotate(self,rang):
        """ Rotate the texture using _rot_ module """
        ftx = "%s%s%s"%(self.myjob.cwd,os.sep,self.myjob.ftx[0])
        #print "texture file that you are rotating is '%s'"%ftx
        _rot_(filename=ftx,ang=rang)
        self.history = self.history + rang
    def run(self):
        self.myjob.run()
    def pp(self):
        """ 
        Perform post-execution rocess for myjob 
        as well as for save lank data
        """
        self.myjob.pp()
        self.prof()
        print "\n**************************************"
        print "currently you have rotated %f degree"%self.history
        print "**************************************"
    def prof(self):
        """ This is a post-process for each job's 'lank' datamaster """
        self.x=[]
        self.R=[]
        x = np.linspace(0,90,90./self.dang+1)
        x0 = np.linspace(90.,0.,90/self.dang+1)
        x = np.array(( x, x0))
        x = np.array(x)
        total_step = x.shape[0] * x.shape[1]
        x = x.reshape((total_step,))
        R = self.myjob.datamaster['lank'][1]
        R = np.array(R)
        self.x=x
        self.R=R.reshape((total_step,))
    def plot(self,ifig=None,marker='.',color=None):
        if ifigure==None:
            return -1
        plt.figure(ifig)
        if color==None: plt.plot(self.x,self.R,marker)
        if color!=None: plt.plot(self.x,self.R,marker,color=color)
    def write(self,filename):
        """ given the file name writes angle and R value """
        fout(filename,self.x,self.R)

def ex1(ext='cmb', dang=3.,filename=None, interaction=3):
    """
    Reads all texture/*.cmb files and writes the results 
    into those files whose names are R_* using 'class r_text'

    >>> vpsc_param.ex1(ext='cmb', dang=3.) 
                     or
    >>> vpsc_param.ex1(filename='texture/00500.cmb', dang=5.)

    >>> vpsc_param.plt_R (ifig=3) ---> should be performed 
                                       where matplotlib is accessible
            OR, YOU MAY USE THE BELOW
    >>> vpsc_param.plt_R2(ifig=3) ---> unfolded lankford profile
                * It is possible to do plt_R indepently,
                in which all the resulting files are located

       ex1 and plt_R can be run independently since the output
      of ex1 module is written in files. By invoking plt_R you'll
      read those files and plot them properly provided that you have
      properly installed matplotlib module.
    """
    if filename==None:
        FILES = glob.glob("%s%s*.%s"%('texture',os.sep, ext))  
    else:
        if type(filename).__name__=='str': FILES = [filename]
        elif type(filename).__name__=='list': FILES = filename
        
    #which together with ext='cmb' makes "texture/*.cmb"

    print 'files',FILES #;raw_input()
    jobs = []
    #FILES=["%s%s%s"%('texture',os.sep,'00500.cmb')]
    for i in FILES:
        jobs.append(r_text(ftex=i, dang=dang,interaction=interaction))
    for i in range(len(jobs)):
        jobs[i].run()
        jobs[i].pp()
        #jobs[i].write(filename='R_%s'%str(i).zfill(2))
        if len(glob.glob('R_*'))==0: 
            fout('R_%s'%str(i).zfill(2),jobs[i].x,jobs[i].R)
        else:
            ib = 0
            while True:
                if os.path.exists('R_%s'%str(i+ib).zfill(2)):
                    ib = ib + 1
                else:
                    fout('R_%s'%str(i+ib).zfill(2),
                         jobs[i].x,jobs[i].R)
                    break
            pass
    # This is a reminder in order to encourage the end-user
    # to plot the results without having to look at the manual once more.
    print ex1.__doc__
    print "\n"
    print "You just have performed module ex1"
    print "There should be resulting R_* files"
    print "in it, you will find the Lankford coefficients"
    print "for the first increments of tension"


    return jobs  #returns the vpsc class for further application


## ----------------------------------------------------------------- ##
## Here the below is the plotting procedures for resulting R values  ##
## ----------------------------------------------------------------- ##
def plt_R(ifig,marker=None,*args):
    """ 
    Plots all "R_*" files if no *args is given.
    Given *args will be plotted.

    --Arguments:
        ifig
        marker
        *args    : Supposedly, a list of string of file names
    """
    FILES=[]
    if len(args)==0: FILES = glob.glob("R_*")
    else: 
        for i in range(len(args)): FILES.append(args[i])
        
    print "\n Below is the files that you'd like to plot"
    print FILES
    raw_input(" Enter to proceed >>")
    for i in range(len(FILES)):
        if marker==None: fplot(filename=FILES[i], ix=0, iy=1, 
                               nhead=1, ifig=ifig)
        else: 
            x,y = fplot(filename=FILES[i], ix=0, iy=1, 
                        nhead=1, ifig=ifig)
            x=np.array(x[0]); y=np.array(y[0])
            figure=plt.figure(ifig)
            try: figure.Axes[0]
            except: figure.gca().plot(x,y,marker=marker)
            else:
                myax=figure.add_subplot(111)
                myax.plot(x,y,marker=marker)
    fig = plt.gcf()
    ax=fig.gca()
    ax.set_xticks(np.arange(0.,90.001,45.))
    ax.grid('on')
    ax.set_xlabel('angle from RD',dict(fontsize=25))
        
def plt_R2(ifig, marker=None, *args):
    """ 
    Plots all 'R_*' without folding them in half

    Note that the class r_text produces the results in the way that 
    the spanning angular axis is as folded to overlapped the result.
    Since this is sometimes uncesseary or uncomfortable, plt_R2 unfold
    the axis when plotting.
    """
    FILES = []
    if len(args)==0:
        FILES = glob.glob("R_*")
    else: 
        for i in range(len(args)): FILES.append(args[i])
    
    print "\n Below is the files that you'd like to plot"
    print FILES
    raw_input(" Enter to proceed >>")
    fig = plt.figure(ifig)
    ax=fig.add_subplot(111)
    for i in range(len(FILES)):
        x, y = fplot(filename=FILES[i], ix=0, iy=1, nhead=1, ifig=None)
        x = x[0]; y = y[0]
        nst = len(x)/2
        inc = 90/(nst-1)
        x1 = np.linspace( 0., 90., 90./inc + 1)
        x2 = np.linspace(90.,180., 90./inc + 1)
        newx = np.array([x1,x2]).reshape((len(x),))
        if marker==None: ax.plot(newx,y)
        else : ax.plot(newx,y, marker)
        fout('_%s'%FILES[i],newx,y)
        
        #182 
        # This is the number of element if one uses 1 
        # as the resolution angle increment for lnkf probing
        
        """
        x[0]/2     :# of steps 
        inc=90/(x[0]/2 - 1) : #degree increment
        x: 0~90 , x: 90~ 180
        x1=np.linspace(0.,90., 90./inc+1)
        x2=np.linspace(90.,180.,90./inc+1)
        newx=np.array([x1,x2]).reshape((len(x),))
        """
    fig = plt.gcf()
    ax=fig.gca()
    ax.set_xticks(np.arange(0.,180.001,45.))
    ax.grid('on')
    ax.set_xlabel('angle from RD',dict(fontsize=25))

def plt_R1(ifig):
    """ Plot only selected files """
    FILES = glob.glob("R_*")
    files = []
    for i in range(len(FILES)):
        print "Is '%s'  included?(y,n) "%FILES[i]
        if raw_input()=='y': files.append(FILES[i])
    for i in range(len(files)):
        fplot(filename=files[i], ix=0, iy=1, nhead=1, ifig=ifig)

def rotate_all(rang=0.):
    """ 
    Rotate all '*.cmb' files in the folder 
    
    End-users must know how this would be like. The files
    will have their own names without affecting anything else
    than the orientation of grains, which means one may not be able  
    to figure out whether or not the files has been rotated by just
    looking at the files

    Argument : rang=0. -rotating angle 
    """
    #ftx = "%s%s%s"%(self.myjob.cwd,os.sep,self.myjob.ftx[0])
    #print "texture file that you are rotating is '%s'"%ftx
    FILES = glob.glob('*.cmb')
    for i in FILES:
        _rot_(filename=i,ang=rang)

    print "\n You just have rotated %i"%len(FILES)
    print "texture file with an angle of %f"%rang

"""
500, 1000, 9000, 100
"""        

###
###  EXAMPLE #5 PROPROTIONAL LOADING and get the pcys evolution
###
###  Path dependency of yield surface can be found 
###      even between monotonic proportional loadings
###

class ev_pcys():
    """
    >>> myev = vpsc_param.ev_pcys(texture='texture/02000.cmb',
                                   stp=60, eqc=0.005, ict=7)
    >>> myev.run()
    >>> myev.pp()
    >>> myev.write()
    ---> You can plot the result in the same folder but in where you can
        access to matplotlib library as following
    >>> vpsc_param.plt_ev(ifig=6)
    """
    def __init__(self,texture='texture/00500.cmb',
                 stp=2, eqc=0.005, ict=7,interaction = 3,
                 fsx='sx/Neil.sx', u11=1.0, u22=-0.5):
        self.myjob = vpsc(mode='ev_pcys', 
                          stp=stp, eqc=eqc, ict=ict,
                          texture=texture, fsx=fsx,
                          u11=u11, u22=u22, 
                          interaction = interaction
                          )
    def run(self,):
        self.myjob.run()
        pass
    def pp(self,):
        self.myjob.pp()
    def write(self,):
        self.E11 = self.myjob.datamaster['E11']
        self.S11 = self.myjob.datamaster['S11']
        self.wrk = self.myjob.datamaster['wrk']
        sigx = self.myjob.datamaster['pcyspl'][0]
        sigy = self.myjob.datamaster['pcyspl'][1]
        #sigx = np.array(sigx); sigy = np.array(sigy)
        for i in range(len(sigx)):
            fout("YS_%s"%str(i).zfill(3), sigx[i], sigy[i])
            
def plt_ev(ifig=5):
    """
    plots... to be continued.
    """
    files = glob.glob("YS_*")
    fig = plt.figure(ifig)
    ax = fig.add_subplot(111)
    for i in range(len(files)):
        x,y = fplot(filename=files[i],ix=0,iy=1,nhead=1)
        x = x[0]; y = y[0]; x = np.array(x); y = np.array(y)
        ax.plot(x,y,'.-')
    ax.set_xlim(0.,); ax.set_ylim(0.,); ax.set_aspect('equal')
    pass

#############################################################################
##  example #6
##  R-VALUE linear fitting of E_x and E_y after a certain amout of elongation
##
##  * R-value is estimated just like in the experimental way.
##  For each direction which is incrementally increasing from +RD to -RD (0~pi)
##  uniaxial tension is carried out. The width strain and thickness strain 
##  which together provides the R-value out of the linear fitting of the two

##  Added feature(s)
##  1. Now it records Yield stress of the first increments.
##  2. Writes it down to YSprob files

"""
Here uniaxial tensions with given increments and steps are carried out. The
resulting longitudinal and transverse strain will produce a linearly fitted
slope that amounts to R-value. By so doing, the unstability of using only the
first increment is solved, which is done in VPSC code's SUBROUTINE LANKFORD.


class R_prob can be a single running upon a texture file.
Below, rpro_ex1 runs class R_probs on many different texture file
one by one, which in the end makes a R-value dependency upon different 
variables, for example, number of grains.
"""

class R_prob:
    """
    myR = vpsc_param.R_prob(dang=45, texture='texture/00500.cmb',
                            fsx=None, interaction=3,stp=3, eqc=0.005)
    myR.__run__(); myR.__pp__()
    myR.__calcR__(); myR.__yieldstress__()
    myR.__writes__(R_filename='Rprb', YS_filename='YSprob')
    """
    def __init__(self, dang=3, texture=None, fsx=None,
                 ## information how many probings
                 th_0=0., th_1=180.000001,
                 ## information on how each uniaxial 
                 ## loading is performed
                 eqc=0.005, stp=15, ict=7, interaction=3,
                 mode='ten_ang'):


        if interaction ==0: 
            print "Cannot impose FC"
            raise IOError
        self.jobs = []
        th_1 = th_1+0.00001  #for gaurantee the next line
        angs = np.arange(th_0,th_1,dang)
        """
        for i in range(len(angs)):
            angs[i] = angs[i]*np.pi/180.  #to radians
        """
        for i in range(len(angs)):
            print "*****************"
            print "R_prob's job #%i"%i
            print "*****************"
            self.jobs.append(vpsc(texture=texture, fsx=fsx,
                                  eqc=eqc, stp=stp, ict=ict,
                                  ang=angs[i],
                                  mode='ten_ang',
                                  interaction=interaction
                                  )
                             )
        self.angs = angs
        pass

    def __run__(self):
        for i in range(len(self.jobs)): self.jobs[i].run()

    def __pp__(self):
        for i in range(len(self.jobs)): self.jobs[i].pp()

    def __calcR__(self):
        """
        Calculates the R-value
        """
        E33 = []; E22 = [] ;R = []
        for i in range(len(self.jobs)):
            et = self.jobs[i].datamaster['E33'][0]
            ew = self.jobs[i].datamaster['E22'][0]
            try:
                r = __slope__(x=et,y=ew)
                R.append(r)
            except: 
                print "Error occured during calculation of slope"
                print "et and ew are repsectively as below"
                raw_input('Press Enter to proceed  >>')
                print "**et**\n",et
                print "**ew**\n",ew
                R.append('NA')
            else:
                print "\n***********************************"
                print "Ang %f loading's R: %f "%(self.angs[i],R[i])
                print "***********************************"
        self.R = R

    def __yieldstress__(self,inorm = 0):
        self.ys = []
        for i in range(len(self.jobs)):
            ##flattening the multidimentional S11 numpy.array
            ys = self.jobs[i].datamaster['S11'][0].flatten(1)
            ys = ys[0]
            self.ys.append(ys)
        self.ys = np.array(self.ys)
        self.ys = self.ys / self.ys[inorm]
    
    def __writes__(self,R_filename="Rprob", YS_filename="YSprob"):
        fout(R_filename,self.angs,self.R)
        fout(YS_filename,self.angs,self.ys)
        pass

"""
    def __init__(self, dang=3, texture=None,
                 th_0=0., th_1=180.000001,  ## information how many probings
                 eqc=0.005, stp=15, ict=7,  ## information on how each uniaxial loading is performed
                 mode='ten_ang'):
"""

def rpro_ex1(ext='cmb', filename=None,
             #For total increment
             th_0=0., th_1=180.00001, dang=3.,  

             #for each subset unixial tension 
             eqc=0.005, stp=15, ict=7,    
             message=None
             ):
    """
    Run the class R_prob for a number of texture.
    The texture file under investigation is fixed to be
    the return of the wild card file collector "glob", 
    glob.glob('texture/*.cmb'). The only argument is
    incremental in forming

    It is designed to return myjobs to further investigate directly.

    ** Arguments:
       dang :incremental angle
       ext='cmb'  :file extension in 'texture' folder
       filename=None : If filename is given then ext is not taken
       th_0=0., th_1=180.0001 : start and end theta from RD
       dang=3. :angle increment
       eqc=0.005 :strain increment
       stp=15 :nstep of strain increment
       ict=7 : strain control mode
       
    ** Exemplary use
    >>> mypro = vpsc_param.rpro_ex1(ext='cmb', dang=5., stp=30)
                   or
        mypro = vpsc_param.rpro_ex1(filename='texture/00500.cmb',
                                     dang=5., stp=30)


    >>> mypro.myjobs ---> has a list of the class vpsc

    """
    dang = float(dang)
    if filename==None: 
        filename = 'texture/*.%s'%(ext)
        files = glob.glob(filename)
    elif filename!=None: 
        if type(filename).__name__=='list':  files = filename
        elif type(filename).__name__=='str': files = [filename]

    myjobs = []
    for i in range(len(files)):
        myjobs.append(R_prob(texture=files[i],
                             dang=dang, th_0=th_0, th_1=th_1,
                             eqc=eqc, stp=stp, ict=ict
                             )
                      )
    for i in range(len(files)):
        #exectes each class vpsc
        myjobs[i].__run__()
    for i in range(len(files)):
        #post-execution process for each class vpsc
        myjobs[i].__pp__()
    for i in range(len(files)):
        #Calculates the slopes of Et vs. Ew
        try:
            myjobs[i].__calcR__()
        except:
            print "__calcR__ raised an Error for the job #%i"%i

        try:
            myjobs[i].__yieldstress__()
        except:
            print "__yieldstress__ raised an Error for the job#%i"%i
    for i in range(len(files)):
        #writes the individual results
        myjobs[i].__writes__(R_filename='Rprob_%s'%str(i).zfill(3),
                             YS_filename='YSprob_%s'%str(i).zfill(3)
                             )

    print rpro_ex1. __doc__  #Show the __doc__ recursively!!!
    print "\n***************************************************"
    print " Information of the last run"
    print " You ran vpsc_param.rpro_ex1 module"
    print " Below is the information of the last run\n"
    print " angle increment was %f"%dang,
    print " and extension of the texture files were %s"%ext
    print " You has total %i jobs executed\n"%len(files)
    print " Now you have the resulting files of Rprob_*"
    print "***************************************************"
    
    if message!=None: 
        print "\n You left the following message"
        print "**"
        print message
        print "**"

    return myjobs

def current_time():
    pass
    
def r_prob_plot():
    FILE = open('')
    pass


#############################################################################
##  example #7
##  Loading RD, TD, and normalize the stress based on their work level

class RDTD:
    """ I should come up with a better name for this """
    def __init__(self,texture=None,fsx=None,
                 stp=2, eqc=0.005, ict=7 ): 
        self.jobs = []
        #TD, RD, 
        self.modes = ['TD','RD']
        for mode in self.modes:
            cjob = vpsc(mode=mode, texture=texture, 
                        stp=stp, eqc=eqc, ict=ict, fsx=fsx)
            self.jobs.append(cjob)
        pass

    def __run__(self,):
        """ Run the jobs! """
        for job in self.jobs:
            job.run()
        pass
    def __pp__(self,):
        """ post-execution process """
        ## Global stress
        self.S11 = []; self.S22 = []
        self.S33 = []; self.S12 = []
        self.S23 = []; self.S13 = []
        ## Global Strain
        self.E11 = []; self.E22 = []
        self.E33 = []; self.E12 = []
        self.E23 = []; self.E13 = []

        self.wrk = [];
        for i in range(len(self.jobs)):
            self.jobs[i].pp()
            cs11 = self.jobs[i].datamaster['S11'].flatten(1)
            self.S11.append(cs11)
            cs22 = self.jobs[i].datamaster['S22'].flatten(1)
            self.S22.append(cs22)
            cs33 = self.jobs[i].datamaster['S33'].flatten(1)
            self.S33.append(cs33)
            cs12 = self.jobs[i].datamaster['S12'].flatten(1)
            self.S12.append(cs12)
            cs13 = self.jobs[i].datamaster['S13'].flatten(1)
            self.S13.append(cs13)

            ce11 = self.jobs[i].datamaster['E11'].flatten(1)
            self.E11.append(ce11)
            ce22 = self.jobs[i].datamaster['E22'].flatten(1)
            self.E22.append(ce22)
            ce33 = self.jobs[i].datamaster['E33'].flatten(1)
            self.E33.append(ce33)
            ce12 = self.jobs[i].datamaster['E12'].flatten(1)
            self.E12.append(ce12)
            ce13 = self.jobs[i].datamaster['E13'].flatten(1)
            self.E13.append(ce13)


            wrk = self.jobs[i].datamaster['wrk'].flatten(1)
            self.wrk.append(wrk)
        pass
    def __norm__(self,):
        ref = np.copy(self.S11[-1]) # RD
        for i in range(len(self.jobs)):
            self.S11[i]=self.S11[i]/ref
            self.S22[i]=self.S22[i]/ref
        pass
    def __plwrk__(self,*args):
        """
        Given the *args, the stress is interpolated
        """
        self.plwrk = []
        for i in range(len(args)): pass
            
    def __flow__(self,):
        pass
    

#---------------------------------------------------
# What professor Barlat asked me to do (2011-01-31)
# For checking the equivalent work contour for AZ31.
# He was suspicious about the results!!! But eventu-
# ally better understood and accepted it. This 
# taught me that explanations upon complex results
# should be accompanied with a systematic and easy 
# presentation.
#--------------------------------------------------

# The first task: Prob the stresses along certain paths
# to see how the full stress tensor is evoloving.
#    1. Select a region to be probed. (here 130~160 degree) of 
#      strain rate vector is probed. This region is discretized
#      into 4 paths. 
#    2. For each loading, stress and strain components are 
#      written down to 'loading*' and 'straining*' files 
#      together with work levels to make it ready to be plotted.
#      Making use of matplotlib library, simultaneous plotting
#      is performed. 
"""
>>> s, s11, s22, ys = vpsc_param.aa()

>>> plot(ys[0],ys[1])  # useful when whole strain vector plane is spanned

>>> vpsc_param.aaa(ifig=88)
"""
####
def aa(nloading=4,interaction = 3, ## I am always suffered from naming.
                                   ## If difficult I tend to put better-naming
                                   ## off. I don't know why.
       nincr=50,eqc=0.001,ini_ang=130,fil_ang=160.,
       ifig=1,mode='inplanebiaxial'):
       ## mode='ev_pcys', mode='inplanebiaxial', mode='inplanebiaxialfc'
       ## I'm woring on the case that mode is given as 'ev_pcys' 

    if interaction==0 and mode=='inplanebiaxial':
        print "Interaction and mode is not matched."
        return -1

    ## vpsc_param.mywrk class
    mywrk = wrk_cnt(ini_ang=ini_ang, fil_ang=fil_ang, nstep=nloading, stp=nincr,
                    ict=7, eqc= eqc, interaction=interaction,#fc
                    texture='texture/00500.cmb',fsx='sx/Neil.sx',
                    mode = mode, norm=False)

    fig = plt.figure(ifig); ysfig = plt.figure(ifig+1)
    ax = fig.gca() ; ysax = ysfig.gca()
    mywrk.run()
    mywrk.pp()
    s11 = []; s22 = []
    filenames = [];filenames_ = []; filenames__ = [];
    delf = glob.glob('loadings*')
    delfs = glob.glob('straining*')
    delfw = glob.glob('workings*')
    for i in delf:os.remove(i)
    for i in delfs:os.remove(i)
    for i in delfw:os.remove(i)

    for i in range(len(mywrk.jobs)):
        filenames.append("%s%s"%("loadings",str(i).zfill(2)))
        filenames_.append("%s%s"%("straining",str(i).zfill(2)))
        filenames__.append("%s%s"%("workings",str(i).zfill(2)))
        
    for i in range(len(mywrk.jobs)):
        ## stress
        S11=mywrk.jobs[i].datamaster['S11'].flatten()
        S22=mywrk.jobs[i].datamaster['S22'].flatten()
        S33=mywrk.jobs[i].datamaster['S33'].flatten()
        S12=mywrk.jobs[i].datamaster['S12'].flatten()
        S23=mywrk.jobs[i].datamaster['S23'].flatten()
        S13=mywrk.jobs[i].datamaster['S13'].flatten()

        ## strain
        E11=mywrk.jobs[i].datamaster['E11'].flatten()
        E22=mywrk.jobs[i].datamaster['E22'].flatten()
        E33=mywrk.jobs[i].datamaster['E33'].flatten()
        E12=mywrk.jobs[i].datamaster['E12'].flatten()
        E23=mywrk.jobs[i].datamaster['E23'].flatten()
        E13=mywrk.jobs[i].datamaster['E13'].flatten()

        ## work
        wrk=mywrk.jobs[i].datamaster['wrk'].flatten()

        ###
        S11=S11-S33
        S22=S22-S33
        S33=S33-S33
        ###
        s11.append(S11);s22.append(S22)
        #### The columns are consisted of wrk level
        ## and 6 components either strain or stress tensors
        fout(filenames[i],wrk, S11,S22,S33,S12,S23,S13)
        fout(filenames_[i],wrk, E11,E22,E33,E12,E23,E13)
        fout(filenames__[i],wrk)
    ax.plot([-10,10],[0,0],'--',alpha=0.2,color='k')
    ax.plot([0,0],[-10,10],'--',alpha=0.2,color='k')
    ysax.plot([-10,10],[0,0],'--',alpha=0.2,color='k')
    ysax.plot([0,0],[-10,10],'--',alpha=0.2,color='k')
    ##initial yield surface drawing
    ys1 = []
    ys2 = []
    for i in range(len(s11)):
        ys1.append(s11[i][0]);ys2.append(s22[i][0])
    ysax.plot(ys1,ys2,'--',label='YS',alpha=0.3,color='k')
    for i in range(len(s11)):
        ax.plot(s11[i],s22[i],'o',label=str(i))
    ax.set_aspect('equal')
    ##return
    return mywrk,s11,s22,[ys1,ys2]

def ppp(x,y,ifig,title=None):
    plt.figure(ifig).gca().plot(x,y)
    plt.figure(ifig).gca().set_title(title)

def aaa(ifig=88):
    """ 
    This is for plotting the results without having to do
    origin stuffs. But for systematic organization should be
    done, I guess, in origion for databasing.
    """
    files = glob.glob('loadings*')
    wfiles = glob.glob('wrk_cntr_*')
    efiles = glob.glob('straining*')
    figs = []
    for i in range(len(files)):
        s = fin(files[i], 1,   1,2,3,4,5,6)
        w = fin(wfiles[i],3,3)[0]
        e = fin(efiles[i],1,   1,2,3,4,5,6)
        for j in range(6): 
            fig = plt.figure(ifig+i)
            ax  = fig.add_subplot(1,2,1) ##Stress axis
            ax1 = fig.add_subplot(1,2,2) ##Strain axis
            ax.plot(w,s[j])
            ax1.plot(w,e[j])
            ax.set_xlabel('work')
            ax.set_ylabel(r'$\sigma$',dict(fontsize=20))
            ax.legend(('S11','S22','S33','S12','S23','S13'),
                      loc='best')
            ax.set_ylim(-15,20)
            ax.set_xlim(0,0.5)
            
            #fig.savefig("sig_%s"%str(i).zfill(2))
            ax1.set_xlabel('work')
            ax1.set_ylabel(r'$\varepsilon$',dict(fontsize=20))
            ax1.legend(('E11','E22','E33','E12','E23','E13'),
                       loc='best')
            ax1.set_ylim(-0.05,0.2)
            ax1.set_xlim(0,0.5)
            #fig.savefig("eps_%s"%str(i).zfill(2))        
        fig.savefig("se_%s"%str(i).zfill(2))


# The second task: Prob the stresses along certain paths and see how 
# yield surface is evoloving.
#    1. Select a path to be probed.
#    2. For each loading, stress and strain components are 
#      written down to 'loading*' and 'straining*' files 
#      together with work levels to make it ready to be plotted.
#      Addition to these, I also need to plot yield loci.
#      Making use of matplotlib library, simultaneous plotting
#      is performed. 


"""
I probably need to make use of one of the exisiting examples, class ev_pcys():
"""
def task2(ifig=1,interaction=3, u11=1.0, u22=-0.5):
    """
    
    """
    myev = ev_pcys(texture='texture/00500.cmb',
                   stp=20, eqc=0.0025, ict=7, fsx='sx/Neil.sx',
                   u11=u11, u22=u22,interaction = interaction
                   )
    #run, pp, and write'em down.
    myev.run();myev.pp();myev.write()
    # and even plot'em up!!!
    files = glob.glob('YS_*')
    fig = plt.figure(ifig)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    for i in range(len(files)):
        sx, sy = fplot(filename=files[i],ix=0,iy=1,nhead=1)
        sx = sx[0]; sy = sy[0];
        sx.append(sx[0]); sy.append(sy[0])
        ax.plot(sx,sy,'.-')

    s11 = myev.myjob.datamaster['S11']
    s22 = myev.myjob.datamaster['S22']
    s33 = myev.myjob.datamaster['S33']

    s11 = np.array(s11).flatten()
    s22 = np.array(s22).flatten()
    s33 = np.array(s33).flatten()
    s11 = s11 - s33; s22 = s22 - s33

    ## Writes the trajectory down
    fout('YS_trajectory',s11,s22)

    ax.plot(s11,s22,'o',color='k', alpha='0.8',label='path')
    ax.set_xlabel(r'$\sigma_{RD}$',dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{TD}$',dict(fontsize=20))
    ax.set_title("Yield surface evolution u11=%2.1f u22=%2.1f"%(u11,u22))
    fig.savefig('YS.png')

    return s11,s22,s33,myev
    

## Example #8
## Faster YS probing 
"""
>>> myys = vpsc_param.ysprob(interaction=3, fsx=None, texture=None, dang=5.0)
>>> myys.__plot__(ifig=3)
>>> myys.__write__()  --> writes down to the 'YSprob'
"""
class ysprob:
    def __init__(self,interaction=3, fsx=None, texture=None, dang=5.0):
        self.__makeprcs__(dang=dang,nstp=1)
        self.job = vpsc(interaction=interaction, 
                        fsx=fsx,
                        texture=texture,
                        ### process is directly given
                        prcs=self.prcs, iupdate=[0,0,0,0])
        ##-------------------
        ##running and pp
        ##-------------------
        self.job.run()
        self.__pp__()

    def __makeprcs__(self,dang,nstp=1):
        """   """
        histmaker = vpsc_in.histmaker
        rotfile = vpsc_in.rotfile
        self.prcs = []
        ang = np.arange(0.,180.001,dang)
        self.ang = ang
        for i in range(len(ang)): 
            self.prcs.append('0')
            self.prcs.append('hist/temp%s'%str(dang).zfill(3))
            self.prcs.append('4')
            self.prcs.append('rot/r%s'%str(dang).zfill(5))
            ##--------------------
            ##rotation file making
            ##--------------------
            rotfile(filename='r%s'%str(dang).zfill(5), ang=dang)
            histmaker(filename='hist/temp%s'%str(dang).zfill(3), 
                      nstep = nstp, ictrl = 7, eqincr=0.005)

    def __pp__(self):
        """
        """
        self.job.pp()
        s11 = self.job.datamaster['S11'][0]
        self.ys = []
        for i in range(len(s11)):
            if i%2==0: self.ys.append(s11[i])
        """
        for i in range(len(self.prcs)):
            if self.prcs[i*2]=='0':
                self.ys.append(s11[i])
            else: pass
        """
        self.ys = np.array(self.ys)

    def __plot__(self,ifig=1):
        fig = plt.figure(ifig)
        ax = fig.add_subplot(111)
        ax.plot(self.ang, self.ys)
        ax.set_xlabel(r'$\theta$'+' from RD', dict(fontsize=20))
        ax.set_ylabel('YS',dict(fontsize=20))
        ax.set_xticks(np.arange(0.,180.001,45.))
        ax.grid('on')
        pass

    def __write__(self,append=None):
        """
        append must be string
        """
        if append ==None:
            fout('YSprob', self.ang, self.ys)
        else:
            fout('%s%s'%('YSprob',append), self.ang, self.ys)
        
## ex
def myys():
    files = glob.glob('texture/*.cmb')
    myjobs = []

    for i in range(len(files)):
        myjobs.append(ysprob(texture=files[i], interaction=3))
        myjobs[i].__write__(append=files[i].split(os.sep)[1].split('.')[0])

    for i in range(len(files)):
        try: myjobs[i].__plot__(ifig=1)
        except: pass

    try: plt.gca().legend(files)
    except: pass

    return myjobs
    

#------------------
#Below is the plan
#------------------

###
### FLD test ... 
###
###

# how to? which model to start with?
# monotonic loadings 
# Which approach among others will be taken? -- MK theory???? 
# I need to read and study...

# path change..

"""
MANUAL upon vpsc_param.vpsc class 
myvpsc_class = vpsc_param.vpsc()
myvpsc_class.pp()        call post-execution process.
                         __hist__ , __texture__, __str_wrk_R__ .. and so on


 After calling the post-processing sripts you'll have a data class, datamaster.
 datamaster is actually a dictionary variable.
 datamater['E11-S11']
 datamaster['R']
 datamaster['SIG-EPS']
 datamaster['wrk']
"""
