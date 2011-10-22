# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 15:09:51 2010
@author: youngung

pf __ class pf: plotting pole figures from pole8 pf files
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os,glob
try:
    import sx_maker
    cubic = sx_maker.cubic
except:
    print "Could not find sx_maker.py module"
    print "Making the sx file in case there's no proper"
    print "one is given is impossible."
plt.ion()


def poleinmaker(ftex=['TEX_PH1.OUT'],
                ind=[[1,0,0],[1,1,0],[1,1,1]],
                fsx=None, ipfig=0):
    """
    pole8.IN file maker

    ftex = ['TEX_PH1.OUT']  / texture files
    ind = [[1,0,0],[1,1,0],[1,1,1]]  /  miller indices of poles
    fsx = None   / single crystal file
    ipfig = 0    / 0:polefigure , 1: inverse pole figure, -1 for *.SOD -2 for *.EPF
    """
    if fsx==None:
        print " You should give a proper name for single crystal file"
        return -1
    else:
        isfile = glob.glob(fsx)
        isfile = len(isfile)
        if isfile ==0:
            print "\n**********************************"
            print " Could not find '%s' file "%fsx
            print "   Please check this again"
            print "**********************************"
            return -1
    
    FILENAME='pole8.IN'
    FILE=open(FILENAME,'w')
    FILE.writelines('** # of ftex **\n')
    FILE.writelines('%i\n'%len(ftex))
    FILE.writelines('** texture files **\n')
    for fname in ftex:
        FILE.writelines('%s\n'%fname)
    FILE.writelines('** ipfig : 0 (pf),')
    FILE.writelines('1 (IPF), -1 (*.SOD), -2 (*.EPF) **\n')
    FILE.writelines('%i\n'%ipfig)
    FILE.writelines('** Single crystal unit cell **\n')
    FILE.writelines('%s\n'%fsx)
    FILE.writelines('** # pf poles (PF) or sample axes')
    FILE.writelines(' (IPF) **\n')
    FILE.writelines('%i\n'%len(ind))
    FILE.writelines('** Miller indices **\n')
    for i in range(len(ind)):
        for j in range(len(ind[i])):
            FILE.writelines(' %i  '%ind[i][j])
        FILE.writelines('\n')

class pf:
    """
    **pole figure plotting**
    arguments :
       pf_name (5digits)
       path = directory path where pf files are located
       pole8_exe_path = 'p.exe'
       idot = False
       ftex = 'tex_ph1.out'
       sx_file = 'sx/hijhiihb.sx'
    """
    def __init__ (self, pf_name='aaaaa', idot=False,
                  path=None,
                  ftex = 'TEX_PH1.out',
                  sx_file = 'sx/hijhiihb.sx',
                  pole8_exe_path =None,
                  ishow = True):
        """ documentation """
        """ POLE8.in preparation """
        sxfile = glob.glob(sx_file)
        if len(sx_file)==0:
            print "You could not find single crystal file"
            print "This automatically generates fcc crystal file"
            cubic(filename=sx_file)

        #os.system('rm *.DAT')
        #os.system('rm *.dat')
        if path == None: path == os.getcwd()
        if pole8_exe_path == None:
            if os.name=='nt': pole8_exe_path=os.getcwd()+'/p.exe'
            elif os.name=='posix': pole8_exe_path='./p'

        #Remove existing *.DAT file
        if os.name=='nt': #windows
            try: os.system('del *.dat')
            except: pass
        elif os.name=='posix': #Linux: Ubuntu installed using VM
            try:
                os.system('rm *.dat')
                os.system('rm *.DAT')
            except: pass
        else:
            print 'Error: Not allowed os system on which pf.py is called'
            raise IOError


        try:
            fpole = file('pole8.in', 'r')  #reads
        except: 
            print "You don't have proper pole8.in file in your folder"
            print "pf.py generates a pole.IN file under assumption that"
            print "your aggregate is cubic and the Miller indices of your"
            print "pole are (1,1,1), (1,0,0), (1,1,0)"
            if len(raw_input("press Enter if yes >>> "))==0:
                poleinmaker(ftex=[ftex],fsx=sx_file)
            else:
                print 'You should have pole8.in file in the folder'
                raise IOError

        lines = fpole.read()
        lines = lines.split('\n')
        fpole.close()
        fpole = file('pole8.in', 'w')  #writes
        for i in range(3):
            fpole.writelines(lines[i] + '\n')
        fpole.writelines(ftex+'\n')

        iline = 4
        while True:
            if iline==7:
                fpole.writelines(sx_file + '\n');
                iline = iline + 1
            else:
                try: lines[iline]
                except IndexError: break
                fpole.writelines(lines[iline] + '\n')
                iline = iline + 1
        fpole.close()

        """ POLE FIGURE DATA FILE GENERATION """
        irst_flg = os.system(pole8_exe_path)
        if irst_flg == 0:
            print '\n\n*************************************'
            print '   Calling of pole8 was successful!  '
            print '*************************************\n\n'
        elif irst_flg !=0 :
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
            raise IOError
            
        self.pf_name = pf_name
        #print 'path =', path; raw_input()
        self.datafile(path, pf_name)
        #print self.contour_level,
        
        self.myfigure = plt.figure(figsize=(15,5))#,frameon=False)
        self.myfigure.suptitle('Pole Figure___ '+ pf_name, 
                               fontsize=12, 
                               fontweight='bold')
        self.ax = self.myfigure.add_subplot(111)#,frame_on=False)
        self.ax.set_aspect('equal')
        self.ax.set_axis_off()
        
        colors=['k','g','r','b','y','m','b','g',
                'r','b','y','m','b','g','r','b',
                'y','m','b','g','r','b','y','m']
        #colors=['1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1']        
        #colors=[]
        #for i in range(10):
        #     colors.append(str(1-(i+1)*0.2))
        print 'path == ', path
        if path == None: path =os.getcwd()
        try:
            self.plot ("%s%s%s%s"%(path,os.sep,pf_name,'_CC.DAT'), 
                       color='black')
        except IOError:
            print '%s%s%s%s is missing'%(path,os.sep,pf_name,'_CC.DAT')
            return None
        try:
            self.splot("%s%s%s%s"%(path,os.sep,pf_name,'_XX.DAT'), 
                       color='white',ms=1)
        except IOError:
            print "%s%s%s%s is missing"%(path,os.sep,pf_name,'_XX.DAT') 
            return None
        for i in range(len(self.contour_level)):
            mylabel = self.contour_level[i]
            cnt = i + 1
            print 'LABEL = ', mylabel
            if idot==False:
                self.plot("%s%s%s%s%s"%(path, os.sep, pf_name, 
                                        self.contour_level[i], '.DAT'),
                          color=colors[i], label=mylabel, cnt=i)
            elif idot==True:
                self.plot("%s%s%s%s%s"%(path,os.sep,pf_name,
                                        self.contour_level[i],'.DAT'),
                          color=colors[i], label=mylabel, cnt=i, opt='spot')
            
        self.plotlabel("%s%s%s%s"%( path, os.sep, pf_name, '_LB'+'.DAT'), 
                       color='black')
        print 'Plotting finished!! '
        if ishow==True: self.show()

        


        ### legend making
        self.contourlevels = np.array( map(float,self.contour_level))*0.1
        """
        for icnt in range(len(self.contourlevels)):
            l1=self.ax.legend(([str(self.contourlevels[icnt])]),
                              loc=1,
                              bbox_to_anchor=(1.06,1-(cnt+1)*0.1),
                              frameon=False)
            self.ax.add_artist(l1)
        """
        """
        l1=self.ax.legend([p1],[label],loc=1,
                          bbox_to_anchor=(1.06,1-cnt*0.1),frameon=False)
        
        self.ax.add_artist(l1)
        """
            
            


    def plot(self, filename, color='r', ms=1, 
             opt='line', label='label-1',cnt=0):
        f = file(filename)
        kount = 0
        self.pfline=[[]]
        count = 0
        while True:
            temp=[]
            s = f.readline()
            if len(s)<5:
                #print 'an EOF is reached'
                break
            if s.split()[0]=='-':
                #print s
                count = count + 1
                pass
            else:
                if self.is_number(s.split()[0]):
                    if self.is_number(s.split()[1]):
                        self.pfline[count].append([float(s.split()[0]),
                                                   float(s.split()[1])])
                        kount = kount + 1
                        self.pfline.append([])
                else:
                    pass
        #print 'count=', count
        #print 'kount=', kount
        
        #line
        #print 'self.pfline= ',self.pfline
        if opt=='line':
            #label = float(label)/10.
            if self.is_number(label):
                label = float(label)/10.
                label = str(label)
            for ipf in range(count):
                x, y = [],[]
                for i in range(len(self.pfline[ipf])-1) :
                    #x = [self.pfline[ipf][i][0],self.pfline[ipf][i+1][0]]
                    x.append(self.pfline[ipf][i][0])
                    #y = [self.pfline[ipf][i][1],self.pfline[ipf][i+1][1]]
                    y.append(self.pfline[ipf][i][1])
                p1 = self.ax.plot(x,y,color=color)
                p1[0].set_clip_on(False)  



                #print 'label=',label
                try: float( label )
                except: pass
                else:
                    l1 = self.ax.legend([p1],[label],loc=1,
                                        bbox_to_anchor=(1.06,1-cnt*0.1),
                                        frameon=False)
                    self.ax.add_artist(l1)

                """
                if cnt==0: pass
                else:
                    l1 = self.ax.legend([p1],[label],loc=1,
                                      bbox_to_anchor=(1.06,1-cnt*0.1),
                                      frameon=False)
                    self.ax.add_artist(l1)
                """


                
                           
        #spot
        elif opt=='spot':
            for ipf in range(count):
                x , y = [], []
                for i in range(len(self.pfline[ipf])):
                    x.append(self.pfline[ipf][i][0])
                    y.append(self.pfline[ipf][i][1])
                plt.plot(x,y,',',color='k', alpha=0.7)
                plt.gray()
        f.close()

    def plotlabel(self,pfname,color='black'):
        f=file(pfname,'r')
        while True:
            s=f.readline()
            if len(s)<5:
                #print 'an EOF is reached'
                break
            plt.text(float(s.split()[0]),float(s.split()[1]),s.split()[2])
        f.close()
    
    def splot(self,filename,color='black',ms=3,opt='spot'):
        self.plot(filename,color,ms,opt)
    
    def is_number(self,st):
        try:
            float(st)
            return True
        except ValueError:
            return False
                
    def datafile(self,path,pfname):
        """
        
        """
        if path ==None: path =os.getcwd()
        files=os.listdir(path)
        filenames=[]
        if os.name=='nt':
            for i in range(len(files)):
                if files[i].split('.')[-1]=='DAT':
                    filenames.append(files[i])
        elif os.name=='posix':
            for i in range(len(files)):
                if files[i].split('.')[-1]=='dat':
                    filenames.append(files[i])
                elif files[i].aplit('.')[-1]=='DAT':
                    filenames.append(files[i])
        else: print 'ERR: Unexpected os system'; raise IOError
        idit=0
        tempname=[]
        tempfullname=[]
        self.contour_level=[]
        for i in range(len(filenames)):
            prefix=filenames[i].split('.')[0][0:len(pfname)]
            if prefix==pfname:
                #print 'filename=',filenames[i]
                #print 'prefix=',prefix
                if self.is_number(filenames[i].split('.')[0][-1]):
                    self.contour_level.append(filenames[i].split('.')[0][5:8])
                else:
                    #if they are labels or cc or some other meaningful things
                    pass
            else:
                pass

    def show(self):
        plt.show()


def plot(filename, color='r', ms=1, 
             opt='line', label='label-1',cnt=0):
    f = file(filename)
    kount = 0
    pfline=[[]]
    count = 0
    while True:
        temp=[]
        s = f.readline()
        if len(s)<5:
            #print 'an EOF is reached'
            break
        if s.split()[0]=='-':
            #print s
            count = count + 1
            pass
        else:
            if self.is_number(s.split()[0]):
                pfline[count].append([float(s.split()[0]),float(s.split()[1])])
                kount = kount + 1
                pfline.append([])
            else:
                pass
        #print 'count=', count
        #print 'kount=', kount
        
        #line
        #print 'self.pfline= ',self.pfline
    if opt=='line':
        #label = float(label)/10.
        if self.is_number(label):
            label = float(label)/10.
            label = str(label)
        for ipf in range(count):
            x, y = [],[]
            for i in range(len(self.pfline[ipf])-1) :
                #x = [self.pfline[ipf][i][0],self.pfline[ipf][i+1][0]]
                x.append(self.pfline[ipf][i][0])
                #y = [self.pfline[ipf][i][1],self.pfline[ipf][i+1][1]]
                y.append(self.pfline[ipf][i][1])
            p1 = plt.plot(x,y,color=color)
            p1[0].set_clip_on(False)  
            if cnt==0 :
                pass
            else:

                #loc='upper right',bbox_to_anchor=(1.,1-cnt*0.5)
                l1=plt.legend([p1],[label],loc=1,bbox_to_anchor=(1.06,1-cnt*0.1))
                    
                plt.gca().add_artist(l1)
                           
        #spot
    elif opt=='spot':
        for ipf in range(count):
            x , y = [], []
            for i in range(len(pfline[ipf])):
                x.append(pfline[ipf][i][0])
                y.append(pfline[ipf][i][1])
            plt.plot(x,y,',',color='k', alpha=0.7)
            plt.gray()
    f.close()
