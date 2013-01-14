"""
update version of pf.py scripts

author : Youngung Jeong
"""
import os, glob
import pf
import numpy as np
import matplotlib.pyplot as plt
pmaker = pf.poleinmaker

class pf:
    """
    """
    def __init__(self, pfname='aaaaa', 
                 ftex=['TEX_PH1.OUT'], fsx='sx/hijhiihb.sx',
                 pole8=None,
                 ind=[[1,0,0],[1,1,0],[1,1,1]],
                 ipfig=0,
                 idot=False,
                 ifig=3):
        """
        pfname = 'aaaaa'               /pole figure lable, typically 'aaaaa'
        ftex = [ 'TEX_PH1.OUT']        /texture files
        fsx = 'sx/hijhiihb.sx'         /single crystal file
        pole8 = None                     /name of pole8 program executable
        ind = [[1,0,0],[1,1,0],[1,1,1]]  /Miller indices of poles
        ipfig = 0                      /pf(0), ipf(-1) ... and so on
        idot = False                   /contour or dot type polefigure flag
        """
        ## initiation of global variables
        self.contourfiles = []; self.dfiles = []
        self.contour_levels = []; self.ftex = ftex
        self.fsx = fsx; self.ind = ind
        self.ipfig = ipfig; self.FILES = []
        self.pole8 = pole8; self.pfname = pfname
        ##-------------------------------
        
        self.run()

        for i in self.FILES:
            try: float(i.split(self.pfname)[1][0:3])
            except: self.dfiles.append(i)
            else: self.contourfiles.append(i)


        ## iso-level contour lines are fetched from self.xy
        self.x,self.y = self.xy(FILES=self.contourfiles)


        ## Sample axis(axes) and circle(s) are fetched from '_CC.DAT' file
        self.xc,self.yc = self.xy(FILES=["%s%s"%(pfname,'_CC.DAT')])
        self.xc = self.xc[0]; self.yc = self.yc[0]
        if len(self.xc)==len(self.yc):
            self.npf = len(self.xc)
        else: print "\nSomething Wrong!"; raise IOError
        self.xax=[]; self.yax=[]
        for i in range(self.npf):
            self.xax.append(self.xc[i][0:4])
            self.yax.append(self.yc[i][0:4])
        ## Decompose circular rims from the circle and axes mixed string
        for i in range(self.npf):
            self.xc[i]=self.xc[i][4:len(self.xc[i])-1]
            self.yc[i]=self.yc[i][4:len(self.yc[i])-1]

        ## plot circle(s)
        self.fig = plt.figure(ifig)
        self.myax = self.fig.add_subplot(1,1,1)
        self.myax.set_aspect('equal')
        for i in range(self.npf):
            self.myax.plot(self.xc[i], self.yc[i], 'k-')

        ## plot points indicating empty space
        self.xdots, self.ydots = self.xy(FILES=["%s%s"%(pfname,'_XX.DAT')])
        self.xdots = self.xdots[0]; self.ydots = self.ydots[0]
        for ipf in range(self.npf):
            self.myax.plot(self.xdots[ipf],self.ydots[ipf],',',
                           color='k', alpha=0.6)


        ## Material Axis(axes)
        for ipf in range(self.npf):
            self.myax.plot(self.xax[ipf],self.yax[ipf],
                           alpha=0.5,color='k')

        ## Iso-contour levels
        colors = ['gray','g','r','b','y','m','b','g',
                  'r','b','y','m','b','g','r','b',
                  'y','m','b','g','r','b','y','m']
        for icnt in range(len(self.contourfiles)):
            for ib in range(len(self.x[icnt])):
                self.myax.plot(self.x[icnt][ib],self.y[icnt][ib],color=colors[icnt])
        

        
                                    
        

            
        """
        self.xc,self.yc = self.cc()
        """

    def run(self):
        """
        Executes Tome's pole8 program
        If it was successful, returns 0, otherwise -1
        """
        if os.name=='nt':
            self.FILES = glob.glob('*.dat')
            if self.pole8==None:
                self.pole8 = 'p.exe'
        elif os.name=='posix': 
            self.FILES = glob.glob('*.dat') + glob.glob('*.DAT')
            if self.pole8==None:
                self.pole8 = './p'
                
        # deletes exisiting '*.DAT' files        
        #for i in FILES: os.remove(i)
        
        # pole8.IN file maker
        pmaker(ftex=self.ftex, fsx=self.fsx,
               ind=self.ind, ipfig=self.ipfig)
        #irst_flg = os.system(pole8)
        irst_flg=0
        if irst_flg == 0:
            print '\n\n*************************************'
            print '   Calling of pole8 was successful!  '
            print '*************************************\n\n'
            return 0
        elif irst_flg !=0 :
            #raise IOError
            print '\n\n*************************************'
            print '   Calling pole8 was not successful '
            print 'pf function ends with raising IOError'
            print '*************************************\n\n'
            if irst_flg == 512:
                print '\n\n************************************************'
            	print 'Error possibly related to input file has occured'
                print '     Check out the input files to POLE8 code    '
                print '************************************************\n\n'
                print
                raw_input("Press Enter  >>>")
            return -1


    def cc(self):
        FILENAME = '%s%s'%(self.pfname,'_CC.DAT')
        try: FILE = open(FILENAME, 'r')
        except:
            print "\nCould not find circle data file '*_CC.DAT' "
            raise IOError
        lines = FILE.readlines()
        xmaster=[];ymaster=[]
        x,y=[],[]
        for i in range(len(lines)):
            x,y=[],[]
            if i== 0: 
                x.append([])
                y.append([])
                ib = 0
            try: temp = map(float, lines[i].split())
            except: 
                """
                ib = ib + 1
                x.append([])
                x.append([])
                """
                pass
            else:
                try:
                    x[ib].append(temp[0]);y[ib].append(temp[1])
                except IndexError:
                    pass
                    #print 'temp = ', temp
        return x,y

    def lb(self):
        pass

    def xx(self):
        pass

    def xy(self,FILES=None):
        """
        Finds contour files and data files
        contourlevel - block level
        """

        
        #_cc.dat : circle and axes 
        #_lb.dat : axes label position
        #_xx.dat : spots on which the pole figure contour does not cover
        xmast=[];ymast=[]
        for fn in FILES:#self.contourfiles:
            #Loop over iso-contour levels(files)
            icon = fn.split(self.pfname)[1]
            """
            current_lev = float(icon[0:1])
            current_lev = current_lev + float(icon[1:3])*0.01
            self.contour_levels.append(current_lev*10)
            self.isolines = []
            """
            FILE = open(fn,'r')
            temp = FILE.read().split('     -         -    ')
            x=[];y=[]
            for ib in range(len(temp)):
                #contour over blocks
                ctemp = temp[ib].split('\n')
                xb = []; yb = []
                for ilines in range(len(ctemp)):
                    # if ctemp is empty then just pass it
                    if len(ctemp[ilines])==0: pass
                    else:
                        try: temp_xy = map(float, ctemp[ilines].split())
                        except: pass
                        else:
                            try: xb.append(temp_xy[0]);yb.append(temp_xy[1])
                            except: print 'temp_xy:', temp_xy
                if len(xb)*len(yb)!=0: x.append(xb);y.append(yb)
            if len(x)*len(y)!=0: xmast.append(x);ymast.append(y)
            
        #Reports how files are configured.
        #Number of contour lines as well as 
        """
        print "# of iso-contour levels: %i"%len(self.contour_levels)
        print " iso-pf intensity levels are as below"
        for i in range(len(self.contour_levels)):
            print "%6.2f  "%self.contour_levels[i],
        """


        return xmast,ymast

            

        
            
                
        
        
            
            
        
            
