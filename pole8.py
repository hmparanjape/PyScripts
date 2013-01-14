import os, glob, shutil
import numpy as np
import matplotlib.pyplot as plt
from sx_maker import cubic
from pf import poleinmaker

def exec_compile(fsrce = 'pole8.for', 
                 outputfile = 'pole_gfortran_linux',
                 compiler = 'gfortran',
                 flags ='-fdefault-real-8 -fdefault-double-8'):
    # fsrce ='pole8.for'
    # outputfile = 'pole_gfortran_linux'
    # compiler = 'gfortran'
    # flags ='-fdefault-real-8 -fdefault-double-8 '
    os.system('%s %s %s -o %s'%(
            compiler,   #compiler name
            fsrce,      #source code
            flags,      #compiler options
            outputfile  #name of the output file
            ))
    
def exec_compile_g95(fsrce='pole8.for',
                     outputfile='pole_g95_linux',
                     compiler='g95',
                     flags='-d8'
                     ):
    os.system('%s %s %s -o %s'%(
            compiler,
            fsrce,
            flags,
            outputfile))

def f2py_compile(compiler='g95'): 
    #os.system('%s %s %s')
    modulename='pf8'
    fsrce ='pole8.for'
    try:
        os.remove("%s.%s"%(modulename,'so'))
        os.remove("%s.%s"%(modulename,'pyf'))
    except: pass
    if compiler=='g95':
        modulename='pf8_g95'
        iflag = os.system(
            "sudo f2py -c --f90exec='%s' --f77exec='%s' --f90flags='-d8' --f77flags='-d8' -m %s %s "%(
                '/usr/bin/g95',  #--f90exec
                '/usr/bin/g95',  #--f77exec
                modulename, fsrce))
    elif compiler=='gfortran':
        modulename='pf8_gfortran'
        iflag = os.system(
            "sudo f2py -c --f90flags='%s' --f77flags='%s' -m %s %s"%(
                '-fdefault-double-8 -fdefault-real-8 -fno-align-commons', #f90flags
                '-fdefault-double-8 -fdefault-real-8 -fno-align-commons', #f77flags
                modulename,
                fsrce))
    
class pfplot:
    """
    Provided that the resulting pole figure files are created,
    plots the pole figure under a linux system. 

    >>> import pole8
    >>> mypf = pole8.plot(pfname='aaaaa')

    """
    osname = os.name
    def __init__(self, pfname='aaaaa'):
        """ Plotting options  """
        files = glob.glob('*%s*'%pfname)

        ## dat & cnt file creation ##
        self.cntfiles = []
        self.datfiles = []
        self.cnt_levels = []
        for i in range(len(files)):
            flg = files[i].split('.')[0][-3::1]
            try: int(flg)
            except: self.datfiles.append(files[i])
            else: 
                self.cntfiles.append(files[i])
                self.cnt_levels.append(flg)
        ##                         ##

        ## cnt file plot
        self.__plot_cnt__()

        ## dat file plot and Label
        self.__plot_dat__()
        
    def __plot_dat__(self,):
        """
        cc and lb
        """
        segblock = self.__line_seg_plot_file__(
            filename=self.datfiles[0]) ## '_CC'

        for each_loop in segblock:
            x, y = np.array(each_loop).transpose()
            self.ax.plot(x,y, color='grey')
            
        ##label
        lbf = open(self.datfiles[1],'r') ## '_LB' -- label
        lines = lbf.readlines()
        lbf.close()
        for i in range(len(lines)):
            cline = lines[i].split()
            self.ax.text(float(cline[0]), float(cline[1]), cline[2])

        ## xx
        lbf = np.genfromtxt(self.datfiles[3])
        lbf = lbf[0:-1:1]
        x,y = np.array(lbf).transpose()
        self.ax.plot(x,y,',',color='k', alpha=0.7)
            
    def __plot_cnt__(self, mode='color'):
        """
        mode = 'color' or 'bw'    
        """
        self.fig = plt.figure(figsize=(20,5))
        self.ax = self.fig.add_subplot(111,
                                       frameon=False,
                                       aspect='equal')
        self.ax.set_axis_off()
        self.idrawn=[]
        for i in range(len(self.cntfiles)):
            self.idrawn.append(False)

        ## color setting, either colorful or black&white ##
        if mode=='color':
            self.colors=['k','g','r','b','y','m','b','g',
                         'r','b','y','m','b','g','r','b',
                         'y','m','b','g','r','b','y','m']
        elif mode=='bw':
            self.colors=['1','0.9','0.8','0.7','0.6',
                         '0.5','0.4','0.3','0.2','0.1']        
        ##                                               ##

        for i in range(len(self.cntfiles)):
            segblock = self.__line_seg_plot_file__(
                filename=self.cntfiles[i])

            for j in range(len(segblock)): #each_loop in segblock:
                x, y = np.array(segblock[j]).transpose()
                p1 = self.ax.plot(x,y,color=self.colors[i])
                p1[0].set_clip_on(False)

                ## legend drawning
                if self.idrawn[i]==False:
                    l1 = self.ax.legend([p1], 
                                        [self.cnt_levels[i]], loc=1,
                                        bbox_to_anchor=(1.06, 1-i*0.1),
                                        #frameon=False
                                        )
                    self.ax.add_artist(l1)
                    self.idrawn[i]==True

    def __line_seg_plot_file__(self, filename):
        block = np.genfromtxt(filename, dtype='string')
        segments = []
        c_seg = []
        for line in block:
            try: xy = map(float, line)
            except:
                segments.append(c_seg)
                c_seg = []
            else: c_seg.append(xy)
        return segments

class pf_exec:
    """
    Execute the pole8_g95_linux,
    and constructs the pole figure plottings
    """
    def __init__(self, texture=None, fsx=None, miller=None):
        dfiles = glob.glob('*.dat') + glob.glob('*.DAT')
        for ifile in dfiles: os.remove(ifile)
        os.system("./%s"%'pole_g95_linux')
        pfs = pfplot('aaaaa') #plot pole figures...
        self.fig = pfs.fig
        self.ax = pfs.ax
        self.fig.savefig('pf.pdf')

    def __polein__(self, texture,
                   miller, fsx):
        #texture file
        if os.path.exist(texture):
            pass
        else:
            print "Given texture file is not existed"
            raise IOError
            pass

        #single crystal file
        if os.path.exist(fsx):
            pass
        else:  #default miller indicies : (100),(110),(111)
            fsx = 'cubic.sx' # SX file maker
            cubic(ss=[[[1,1,1],]  ,[[1,1,0],]],filename='cubic')
            pass
        poleinmaker(ftex = texture, ind=miller, fsx=fsx)
        pass
        
    
