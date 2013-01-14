"""
In-plane biaxial tester (kokusai) post processing python script

class biaxial (referring some of its methods to uniaxial's)
class uniaixal (some fundamental features and methods + R_value)

Started with uniaxial class thus biaxial class makes use of some
legacy defs from the class uniaxial.

author: Youngung Jeong at MML, GIFT, POSTECH, Korea
date of development : 2010-July-29 ~ November-20
  Ver. 2 (Oct-21)
     Detects the axiality :
       either uniaxial or biaxial based on sig or eps ratio
     def sr :
       Instantaneous strain rate
     Actively works with bpp_tot script
  Ver. 3 (Nov-11)
     Enhanced stress, strain, and direction of plastic strain estimation
       Linear Interpolation
"""

import os
import matplotlib.pyplot as plt
import matplotlib.font_manager
import scipy.integrate as integrate
from scipy.interpolate import splprep, splev
from scipy.interpolate import UnivariateSpline
from numpy import linspace
import numpy as np
import math
os.sys.path.append('c:\\python26\\myworkplace\\myscripts')
from plot_xy import pl

if plt.isinteractive()==False:
   plt.ion()


print __doc__


################################################################################
#   Class biaxial:
#         Class biaxial includes calling of class uniaxial
#         Class biaxial includes several features
#            - Calculation of work (either plastic, total)
#            - Estimate the stresses at the given work
#            - Estimate the direction of strain (either total or plastic)
#            - Plotting data using "matplotlib" library
#                - Deviation of two load cells
#                - Strain and stress raw data with respect to time flow
#                - Stress-strain curves
#                - Plotting pseudo-yield surface with direction of strain rate
#            - Extraction of data in certain formats
#                - Stress-strain curves
#                - Work curve w.r.t time flow
#                - Instantaneous strain rate w.r.t time flow
#                - etc. to be added
################################################################################

class biaxial:
    """
    Post-processes in-plane biaxial results.
    Plots Flow stress-strain curves and changes in strain or stress 
    (To be done but not yet) A flag to choose load_control or strain_control
    The class biaxial is under consideration.

    Arguments :

      path        (=None)
      filename    (=None)
      iplot       (=False)
      modulusx    (=None)
      modulusy    (=None)
      lsx         (=20.)
      usx         (=50.)
      lsy         (=20.)
      usy         (=50.)

    ** ls and us represents respectively lower and upper stress level
       from which the modulus is calculated.
    
    """

    def __init__(self, path=None, filename=None, iplot=False,
                 lsx=20., usx=50., lsy=20., usy=50.,
                 modulusx = None, modulusy=None, iflip=False,
                 str_channel = None, maxtime=None):
        """
        Aurgments:
           path
           filename
           iplot
           lsx = 20.
           usx = 50.
           lsy = 20.
           usy = 50.
           modulusx = None
           modulusy = None
           iflip = False
           str_channel = None

        Notes
           * str_channel is going to be used only for uniaxial axaility case.
            This variable will be a clue to do right analysis with lack of
            full channels. For instance, if you have one strain channel for
            uniaxial tension along x direction, x1 and x2 will be equalized
            by the active channel during the test. Thus, the active channel
            must be indicated so that the end-user can rightly input.

           * The str_channel variable will be just passed to class uniaxial from
            class biaxial, in the first place.

           * str_channel is basically a list variable possibly containing indices
            of strain, e.g. it can be [1,2] or [2].
        """
        if path == None:
            print 
            print '**************************'
            print 'No path is given '
            print 'Current path is hardwired!'
            print '**************************'
            print
            path = os.getcwd()
        if filename == None :
            print
            print '********************'
            print 'No filename is given'
            print 
            return 1
        print '** Upper and lower stress levels'
        print 'from which the moduli are calculated'
        print 'x-axis'
        print 'Upper : ', usx
        print 'Lower : ', lsx
        print 'y-axis'
        print 'Upper : ', usy
        print 'Lower : ', lsy
        self.filename = filename
        self.path = path
        self.mxt = maxtime

        ##---------------------------------
        f = open(filename, 'r')
        # if path[-1]=='\\':
        #    f = file(path+filename)
        # else:
        #    f = file(path+'\\'+filename)
        source = f.read()
        f.close()
        ##---------------------------------

        
        lines = source.split('\n')
        self.uni = uniaxial(path=self.path, filename=self.filename, iplot=iplot,
                            modulusx = modulusx, modulusy=modulusy,
                            low_sigx=lsx, up_sigx=usx, low_sigy=lsy, up_sigy=usy,
                            iflip=iflip, str_channel=str_channel, maxtime=maxtime)
        #uniaxial(iplot=True) plottings are in the figure (10)

        print 
        print 'Sig_ratio ' , self.uni.stressratio[0], self.uni.stressratio[1]
        print 'Eps_ratio ' , self.uni.strainratio[0], self.uni.strainratio[1]
        self.__default_plots__(iplot=iplot, fig_id=20, ind_plot=331)

        
    def total_work(self):
        """
        Returns the plastic total work, plastic strains, elastic strains
        As a bunus, total elastic work is printed. (but not returned)
        No Arguments are necessary since only default arguments are allowed for
        inner calling of self.work definition.
        In that calling you have several arugments set with default values:
           sig_label = 'avg'
           eps_label = 'avg'
           filename  = 'temp.out'
           iplot = False
           sigxl = 'avg'
           sigyl = 'avg'
           epsxl = 'plavg'
           epsyl = 'plavg'
        """
        wkx, wky, wktp = self.work(sig_label='avg', eps_label='avg', filename='temp.out',
                                  iplot=False, sigxl='avg',
                                  sigyl='avg', epsxl='plavg', epsyl='plavg')
        try:
            print '** Max plastic work:', wktp[self.uni.break_indx]
            ind =self.uni.break_indx
        except IndexError:
            print 'The index raising the Error is : ', self.uni.break_indx
            ind = int(raw_input('please enter your own index(trivial)'))
            raw_input('please press enter to proceed')
        
        wkx, wky, wkt = self.work(sig_label='avg', eps_label='avg', filename='temp.out',
                                  iplot=False, sigxl='avg',
                                  sigyl='avg', epsxl='avg', epsyl='avg')
        try:
            print '** Max total work:', wkt[self.uni.break_indx]
            ind = self.uni.break_indx
        except IndexError:
            print 'The index raising the Error is : ', self.uni.break_indx
            ind = int(raw_input('please enter your own index(trivial)'))
            raw_input('please enter to proceed')

        #self.Epsy_pl
        return wktp[ind], self.uni.epsx[ind], self.uni.epsy[ind],self.uni.Epsx_pl[ind], self.uni.Epsy_pl[ind], self.uni.sigx[ind], self.uni.sigy[ind],

    def sumup(self, x, y):
        """
        Returns summedup of x and y
        Arguments :
           x, y in the list type
        """
        temp = []
        for i in range(len(x)):
            temp.append(x[i]+y[i])
        return temp
        
    
    def out(self, filename='pp.out', path = 'c:\\python26', uni=False):
       """
       Writes some selected post-processed outputs into the file
          Arguments-
            filename
            path
            uni = False : this is for uniaxial tests from biaxial machine.
                          It's still got problems.

       Note that only average terms are concerned here.
       Results with different load and strain is not of concern in this method
       """
       
       if uni==False: ### biaxial case
           workx, worky, workt = self.work(sig_label='avg', eps_label='avg',
                                           filename='temp.out',
                                           iplot=False, sigxl='avg',
                                           sigyl='avg', epsxl='avg', epsyl='avg')
           wx = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.epsx)
           wy = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.epsy)
           wtot = self.sumup(x=wx, y=wy)
           wxpl = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.Epsx_pl)
           wypl = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.Epsy_pl)
           wpltot = self.sumup(x=wxpl, y=wypl)
           
       elif uni[0]=='x':  ### uniaxial along x-axis
           workx, worky, workt = self.work(sig_label='avg', eps_label='avg',
                                           filename='temp.out',
                                           iplot=False, sigxl='avg',
                                           sigyl='avg',epsxl='1', epsyl='avg')
           workt = workx
           if uni[1]=='1':
               wx = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.epsx1)
               wxpl = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.Epsx_pl11)
           elif uni[1]=='2':
               wx = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.epsx2)
               wxpl = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.Epsx_pl22)
           #wx = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.epsx2)
           wy = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.epsy)
           wtot = wx
           #wxpl = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.Epsx_pl11)
           wypl = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.Epsy_pl)
           wpltot = wxpl
           
       elif uni[0]=='y':  ### uniaxial along y-axis
           workx, worky, workt = self.work(sig_label='avg', eps_label='avg',
                                           filename='temp.out',
                                           iplot=False, sigxl='avg',
                                           sigyl='avg',epsxl='avg',epsyl='1')
           workt = worky
           if uni[1]=='1':
               wy = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.epsy1)
               wypl = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.Epsy_pl11)
           elif uni[1]=='2':
               wy = integrate.cumtrapz(y=self.uni.sigy,  x=self.uni.epsy2)
               wypl = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.Epsy_pl22)
           wx = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.epsx)               
           #wy = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.epsy1)
           wtot = wy
           wxpl = integrate.cumtrapz(y=self.uni.sigx, x=self.uni.Epsx_pl)
           #wypl = integrate.cumtrapz(y=self.uni.sigy, x=self.uni.Epsx_pl11)
           wpltot = wypl
       
       if path[-1]=='\\': pass
       else : filename = path + '\\' + filename


       ### WRITES HEADER TO THE OUTPUT FILE ###
       f = open(filename, 'w')
       f.writelines('** post-processed data **\n')
       f.writelines(' time        Ex         Ex_pl      Sx')
       f.writelines('         Ey         Ey_pl      Sy')
       f.writelines('         Srx        Sry')
       f.writelines('        Work_tot   Work_x     Work_y')
       f.writelines('     Wtot_pl    Wx_pl      Wy_pl\n')
       
       i = 0
       delT = self.uni.acq_rate/10**3  # time increment

       while True:
         try:
           time = self.uni.time_flow[i]

           if uni==False:  ## biaxial
               expl = self.uni.Epsx_pl[i]
               eypl = self.uni.Epsy_pl[i]
               ex = self.uni.epsx[i]
               ey = self.uni.epsy[i]
               
           elif uni[0] =='x':  ## uniaxial along x
               eypl = self.uni.Epsy_pl[i]
               ey = self.uni.epsy[i]
               if uni[1]=='1':
                   expl = self.uni.Epsx_pl11[i]
                   ex = self.uni.epsx1[i]
               elif uni[1]=='2':
                   expl = self.uni.Epsx_pl22[i]
                   ex = self.uni.epsx2[i]
                   
           elif uni[0] =='y':  ## uniaxial along y
               expl = self.uni.Epsx_pl[i]
               ex = self.uni.epsx[i]
               if uni[1]=='1':
                   eypl = self.uni.Epsy_pl11[i]
                   ey = self.uni.epsy1[i]
               elif uni[1]=='2':
                   eypl = self.uni.Epsy_pl22[i]
                   ey = self.uni.epsy2[i]

           sigmax = self.uni.sigx[i]
           sigmay = self.uni.sigy[i]

           
           if i==0: srx = 0
           else: srx = (expl - self.uni.Epsx_pl[i-1]) / delT
           if i==0: sry = 0
           else: sry = (eypl - self.uni.Epsy_pl[i-1]) / delT
           #f.write('%10.4e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n'% (time, ex, expl, sigmax, ey, eypl, sigmay, srx, sry, wtot[i], wx[i], wy[i], wpltot[i], wxpl[i], wypl[i]))
           f.write('%10.4e '%(time))
           f.write(' %10.3e'%(ex))
           f.write(' %10.3e'%(expl))
           f.write(' %10.3e'%(sigmax))
           f.write(' %10.3e'%(ey))
           f.write(' %10.3e'%(eypl))
           f.write(' %10.3e'%(sigmay))
           f.write(' %10.3e'%(srx))
           f.write(' %10.3e'%(sry))
           f.write(' %10.3e'%(wtot[i]))
           f.write(' %10.3e'%(wx[i]))
           f.write(' %10.3e'%(wy[i]))
           f.write(' %10.3e'%(wpltot[i]))
           f.write(' %10.3e'%(wxpl[i]))
           f.write(' %10.3e\n'%(wypl[i]))
         except IndexError:
             #print 'IndexError has been raised'
             #raw_input('Press Enter')
             break
         else:
           i = i + 1
           pass

       f.close()
       
    def sr(self, epsxl=None, epsyl=None, iplot=False, ifig=82, plot_ind=111, delt= 20):
      """
      Calculates instantaneous strain rate
      Arguments:
        epsxl = None
        epsyl = None
        iplot = False, yet to be made (2010-Sep 2nd)
        ifig, plot_ind = 82, 111
        delt = 20
            Note that delt is size of the bin in which data are selected
             to be used for linearization sort of a thing.
      """
      #Epsilon X
      if epsxl=='1': Ex = self.uni.epsx1
      elif epsxl=='2': Ex = self.uni.epsx2
      elif epsxl=='avg': Ex = self.uni.epsx
      elif epsxl=='pl11': Ex = self.uni.Epsx_pl11
      elif epsxl=='pl12': Ex = self.uni.Epsx_pl12
      elif epsxl=='pl21': Ex = self.uni.Epsx_pl21
      elif epsxl=='pl22': Ex = self.uni.Epsx_pl22
      elif epsxl=='plavg': Ex = self.uni.Epsx_pl
      else:
         print 'Wrong label in def sr'
         print 'The EPSX is corrected as epsx (avg) '
         Ex = self.uni.epsx
         raw_input('Please press Enter')

      #Epsilon Y
      if epsyl =='1': Ey = self.uni.epsy1
      elif epsyl =='2': Ey = self.uni.epsy2
      elif epsyl =='avg': Ey = self.uni.epsy
      elif epsyl=='pl11': Ey = self.uni.Epsy_pl11
      elif epsyl=='pl12': Ey = self.uni.Epsy_pl12
      elif epsyl=='pl21': Ey = self.uni.Epsy_pl21
      elif epsyl=='pl22': Ey = self.uni.Epsy_pl22
      elif epsyl=='plavg': Ey = self.uni.Epsy_pl          
      else:
         print 'Wrong label in def sr'
         print 'The EPSY is corrected as epsy (avg)'
         Ey = self.uni.epsy
         raw_input('Please press Enter')

      delt1 = delt
      delt2 = delt1
      fname = self.filename.split('.')[0]+'.sr'
      if self.path[-1]=='\\': f = open(self.path + self.filename.split('.')[0]+'.sr', 'w')
      else : f = open(self.path + '\\' + self.filename.split('.')[0]+'.sr', 'w')
      f.writelines('     time       norm     xslope     yslope\n')

      if iplot==True:
         fig = plt.figure(ifig)
         ax = fig.add_subplot(plot_ind)
         ax2 = ax.twinx()

      y1, y2, y3, x = [], [], [], []
      """
      for i in range(len(Ex)- delt1 - delt2 - 10):
         xslope = self.uni.__slope__(x=self.uni.time_flow[i+delt1:i+delt1+delt2],
                                     y=Ex[i+delt1:i+delt1+delt2])
         yslope = self.uni.__slope__(x=self.uni.time_flow[i+delt1:i+delt1+delt2],
                                     y=Ey[i+delt1:i+delt1+delt2])
         t = self.uni.time_flow[i+int((delt1+delt2)/2.)]
         x.append(t)
         y1.append(xslope)
         y2.append(yslope)
         y3.append(math.sqrt(xslope**2 + yslope**2))
      """
      for i in range(len(self.uni.time_flow)):
         if i-delt > 1:
            if i+delt < int(len(self.uni.time_flow)-2):
               xslope = self.uni.__slope__(
                  x=self.uni.time_flow[i-delt:i+delt], y=Ex[i-delt:i+delt])
               yslope = self.uni.__slope__(
                  x=self.uni.time_flow[i-delt:i+delt], y=Ey[i-delt:i+delt])
            else: xslope, yslope='NA', 'NA'
         else: xslope, yslope='NA', 'NA'
         y1.append(xslope)
         y2.append(yslope)
         if xslope=='NA': y3.append('NA')
         else: y3.append(math.sqrt(xslope**2 + yslope**2))
      for i in range(len(self.uni.time_flow)):
         f.write('%9.3e'% (self.uni.time_flow[i]))
         if any([y1[i],y2[i],y3[i]][j]=='NA' for j in range(3)):
            f.write('         Na         Na         Na\n')
         else: f.write('  %9.3e  %9.3e  %9.3e\n'%(y3[i],y1[i],y2[i]))

      if iplot==True:
         #ax.plot(x, y1, marker='o', color='red', alpha=0.2)
         #ax2.plot(x, y2, marker='o', color='blue', alpha=0.2)
         ax.plot(x, y1, color='red', alpha=0.1)
         ax.plot(x, y2, color='blue', alpha=0.1)
         ax2.plot(x, y3, marker='o',  color='gray', alpha=0.2)

      f.close()
      print
      print '***************************************'
      print 'Instantaneous strain rate is written in'
      print fname.rjust(20)
      print '***************************************'

      return x, y1, y2, y3

    def __multyieldpoints__(self, w, ifig=1, plot_ind=331):
        """
        Runs yield_vector with all possible combination of 
        two load signals and two strain signals
        for each axis.
        Arguments :
           w : worklevel
           ifig = 1
           plot_ind = 331
        """
        labels=[['1','1','1','1'],
                ['1','1','1','2'],
                ['1','1','2','1'],
                ['1','2','1','1'],
                ['2','1','1','1'],
                ['2','2','1','1'],
                ['2','1','2','1'],
                ['2','1','1','2'],
                ['1','2','2','1'],
                ['1','2','1','2'],
                ['1','1','2','2'],
                ['1','2','2','2'],
                ['2','1','2','2'],
                ['2','2','1','2'],
                ['2','2','2','1'],]

        for i in range(len(labels)):
            self.yield_vector(sigxl=labels[i][0], sigyl=labels[i][1],
                              epsxl=labels[i][2], epsyl=labels[i][3],
                              ifig=ifig, w=w, plot_ind=plot_ind)
            """
            sx,sy,ix = self.plot_yield(w=w, sigxl=labels[i][0], 
                       sigyl=labels[i][1], epsxl=labels[i][2], 
                       epsyl=labels[i][3], ifig=ifig)
            temp, sr= self.strain_rate(epsxl=labels[i][2], epsyl=labels[i][3],
                                       iplot=False)
            vect_norm=20.
            x0 = sx
            y0 = sy
            x1 = x0 + vect_norm*math.cos(sr[ix]*math.pi/180.)
            y1 = y0 + vect_norm*math.sin(sr[ix]*math.pi/180.)
            fig=plt.figure(ifig)
            ax=fig.add_subplot(111)
            ax.plot([x0,x1],[y0,y1],ls='-',color='black')
            """

    def yield_vector(self, sigxl, sigyl, epsxl, epsyl, w, ifig=1,
                     instant=True, plot_ind=111, marker='o', color='black',
                     label='', markersize=20, alpha=0.7, linewidth=2.0,
                     ls='-', iplot=True):
        """
        Plots yield surface + strain vector
        Yield surface is plotted using another def. in the class
        which is plot_yield. Thus, some of the arguments are just
        passed to that method.
        
        Arguments:
             sixgl, sigyl, epsxl, epsyl, w, ifig
             instant = True : True means that instantaneous strain rate is used
                   for determining the direction of strain
             ifig, plot_ind
             iplot = True
             marker = None
             color  = 'k'
             label =' '
             markersize = 20
             alpha = None
             linewidth = 2.0
             ls = '-'
        """
        ## checking the labels // debugging
        labels = [sigxl,sigyl,epsxl,epsyl]
        print labels
        if any(labels[i]==None for i in range(len(labels))):
           print labels
           raw_input
        ##
           
        sx, sy, ix, ibreak = self.plot_yield(sigxl=sigxl, sigyl=sigyl, epsxl=epsxl,
                                             epsyl=epsyl, w=w, ifig=ifig,
                                             plot_ind=plot_ind, marker=marker,
                                             color=color, label=label,
                                             markersize=markersize, alpha=alpha,
                                             ls=ls, iplot=iplot)
        
        vect_norm = 20.
        if instant==True:
            temp, theta = self.instantaneous_strain_rate(epsxl=epsxl,
                          epsyl=epsyl, ind=ix, delt=10, iplot=iplot, k=1)
            x0 = sx
            y0 = sy
            x1 = x0 + vect_norm * math.cos(theta * math.pi/180.)
            y1 = y0 + vect_norm * math.sin(theta * math.pi/180.)
            
        elif instant==False:
            temp, sr = self.acc_strain_rate(epsxl=epsxl, epsyl=epsyl,
                                        iplot=iplot, ifig=30)
            x0 = sx
            y0 = sy
            x1 = x0 + vect_norm * math.cos(sr[ix]*math.pi/180.)
            y1 = y0 + vect_norm * math.sin(sr[ix]*math.pi/180.)
            
        if iplot==True:
            pl(x=[x0,x1],y=[y0,y1], color='black', ls=ls,
               title='Yield Surface', xlabel='Sigx (MPa)', ylabel='Sigy (MPa)',
               ifig=ifig)
        else: pass
        print 'sx, sy, theta  =', sx, sy, theta



        maxtime = self.mxt
        if maxtime != None:
           if self.uni.time_flow[ix-1] >= maxtime:
              return -1, -1, -1, ix, self.uni.time_flow[ix]
           else:
              return sx, sy, theta, ix, self.uni.time_flow[ix]
        elif maxtime == None:
           return sx, sy, theta, ix, self.uni.time_flow[ix]
        else: print 'Something wrong with mxtime given in def yield_vector'
        


        #if mxtime > self.uni.time_flow[ix]:


        
        #return sx, sy, theta, ix, self.uni.time_flow[ix]


      
        #   elif ibreak == -1: return sx, sy, theta, ix, self.uni.time_flow[ix]
        #else :
        #   print
        #   print ' *************************************'
        #   print ' Beyond the assigned maxtime'
        #   print ' The assignmed maxtime is', str(mxtime).rjust(4),'seconds'
        #   print ' *************************************'
        #   print
        #   return -sx, -sy, 720, ix, self.uni.time_flow[ix]
        #   return -10, -10, 720, ix, self.uni.time_flow[ix]
        """
        if self.uni.axiality[0:3]=='uni':
            return sx,sy,theta, ix, self.uni.time_flow[ix]
        else:
            return 0.0001,0.0001,0.0001, 0, 0.0001
        """
         
            

    def instantaneous_strain_rate(self, ind, epsxl='plavg', epsyl='plavg',
                                  delt=10, iplot=False, ifig=30,
                                  ilegend=False, plot_ind=332, k=1):
        """
        Calculate instantaneous strain rate at the certain level of work (index)
        Arguments:
            epsxl, epsyl : strain label must be among '1', '2', and 'avg'
            ind: the index in which the strain is sought
            delt: strain rate will be calculated in [ind-delt:ind+delt] range
            iplot = False
            ifig = 30
            ilegend = False
            plot_ind = 332
            k = 1
        """
        x = self.uni.time_flow
        if epsxl=='1': Ex = self.uni.epsx1
        elif epsxl=='2': Ex = self.uni.epsx2
        elif epsxl=='avg': Ex = self.uni.epsx
        elif epsxl=='pl11': Ex = self.uni.Epsx_pl11
        elif epsxl=='pl12': Ex = self.uni.Epsx_pl12
        elif epsxl=='pl21': Ex = self.uni.Epsx_pl21
        elif epsxl=='pl22': Ex = self.uni.Epsx_pl22
        elif epsxl=='plavg': Ex = self.uni.Epsx_pl
        else:
           print 'Wrong label given[1]'
           raw_input('Press Enter...')
           return -1
        if epsyl=='1': Ey = self.uni.epsy1
        elif epsyl=='2': Ey = self.uni.epsy2
        elif epsyl=='avg': Ey = self.uni.epsy
        elif epsyl=='pl11': Ey = self.uni.Epsy_pl11
        elif epsyl=='pl12': Ey = self.uni.Epsy_pl12
        elif epsyl=='pl21': Ey = self.uni.Epsy_pl21
        elif epsyl=='pl22': Ey = self.uni.Epsy_pl22
        elif epsyl=='plavg': Ey = self.uni.Epsy_pl        
        else:
           print 'Wrong label given[2]'
           raw_input('Press Enter...')
           return -1
        if ind<delt:
           print 'The index give is too small that strain rate cannot be estimated'
           raw_input('Press Enter...')
           return -1
        else:
           x = x[ind-delt:ind+delt]
           y1 = Ex[ind-delt:ind+delt]
           y2 = Ey[ind-delt:ind+delt]
        xs, y1s = self.uni.__UnivariateSpline__(x=x, y=y1,
                          xdown=x[0], xup=x[-1], k=k)
        xs, y2s = self.uni.__UnivariateSpline__(x=x, y=y2,
                          xdown=x[0], xup=x[-1], k=k)
        if iplot==True:
            fig = plt.figure(ifig)
            ax = fig.add_subplot(plot_ind)
            ax.plot(x, y1, marker='o', label='y1')
            ax.plot(x, y2, marker='o', label='y2')
            ax.plot(xs, y1s, ls='--', label='y1s')
            ax.plot(xs, y2s, ls='--', label='y2s')
            ax.set_xlabel('time')
            ax.set_ylabel('Strain')
            ax.set_title(label='Instantaneous strain rate fitting')
            if ilegend==True:
                ax.legend(loc='upper left')
        dExx = y1s[-1] - y1s[0]
        dEyy = y2s[-1] - y2s[0]
        """
        dt = xs[-1]-xs[0]
        Exdot = dExx/dt
        Eydot = dEyy/dt
        """
        try :
             #E = Exdot/Eydot
             E = dEyy / dExx
        except ZeroDivisionError:
             print 'ZeroDivisionError occured in instantaneous_strain_rate'
             raw_input('Press Enter...')
             E = - 1
        theta = math.atan(E)*180. / math.pi
        if dExx < 0: theta = theta + 180
        else : pass
        #print 'Theta from instantaneous strain rate :' , theta
        """
        print 'y1s[-1] = ', y1s[-1]
        print 'y1s[0]  = ', y1s[0]
        print 'y2s[-1] = ', y2s[-1]
        print 'y2s[0]  = ', y2s[0]
        print 'E       = ', E
        """
        if theta == np.nan:
           print theta
           theta = 0.
        
        xs = np.mean(xs)

        
        ##  If theta is numpy.nan because of singularity (possibly for uniaxial cases)
        ##  theta is checked and converted to 0.
        ##  To detect this situation following is done.
        try:
           int(theta)
        except ValueError:
           theta=0.

           
        return xs, theta

    def acc_strain_rate(self, epsxl='1', epsyl='1', iplot=True,
                    k=3, s=1., ifig=30, plot_ind=334):
        """
        Calculates accumulated strain response with respect to the time flow.
        Plotting is given when 'iplot' = True

        Arguments :
             epsxl = '1'
             epsyl = '1'
             iplot = True
             k = 3
             s = 1.
             ifig = 30
             plot_ind = 334
        """
        time_flow = self.uni.time_flow
        if epsxl=='1': epsx = self.uni.epsx1
        elif epsxl=='2': epsx = self.uni.epsx2
        elif epsxl=='avg': epsx = self.uni.epsx
        elif epsxl=='pl11': epsx = self.uni.Epsx_pl11
        elif epsxl=='pl12': epsx = self.uni.Epsx_pl12
        elif epsxl=='pl21': epsx = self.uni.Epsx_pl21
        elif epsxl=='pl22': epsx = self.uni.Epsx_pl22
        elif epsxl=='plavg': epsx = self.uni.Epsx_pl
        else: 
           print 'wrong label given[3]'
           raw_input('Press Enter...')
           return -1

        if epsyl=='1': epsy = self.uni.epsy1
        elif epsyl=='2': epsy = self.uni.epsy2
        elif epsyl=='avg': epsy = self.uni.epsy
        elif epsyl=='pl11': epsy = self.uni.Epsy_pl11
        elif epsyl=='pl12': epsy = self.uni.Epsy_pl12
        elif epsyl=='pl21': epsy = self.uni.Epsy_pl21
        elif epsyl=='pl22': epsy = self.uni.Epsy_pl22
        elif epsyl=='plavg': epsy = self.uni.Epsy_pl
        else: 
           print 'wrong label given[4]'
           raw_input('Press Enter...')
           return -1

        txs = time_flow
        tys = time_flow
        epsxs, epsys = epsx, epsy

        theta=[]
        for i in range(len(txs)):
           try:
             E = epsys[i]/epsxs[i]
           except ZeroDivisionError:
             print ' Zerodivision Error occured in def strain_rate. -1 returns to E '
             print ' epsys[i], epsxs[i] = ', epsys[i], epsxs[i]
             print ' index at this point =', i
             raw_input('Press Enter...')
             E = -1
           temp_theta = math.atan(E)*180. / math.pi
           if epsxs<=0. :  theta.append(temp_theta + 180)#theta.append(math.atan(E)*180./math.pi)
           else : theta.append(temp_theta)

        if iplot==True:
           fig = plt.figure(ifig)
           ax = fig.add_subplot(plot_ind)
           ax.plot(txs, epsxs, label='Ex', marker='o')
           ax.plot(tys, epsys, label='Ey', marker='o')
           ax.legend('upper center')
           ax2 = ax.twinx()
           plt.ylim([0, 90])
           ax2.plot(txs, theta, label='Thet')

           #plotting
           pass
        return txs, theta
        
    def plot_yield(self, sigxl=None, sigyl=None, epsxl=None, epsyl=None,
                   w=None, ifig=None, plot_ind=None, marker=None,
                   color='gray', label=None, markersize=None, alpha = 1.,
                   ls=None, iplot=True ):
        """
        Plots the stress state corresponding to the work
        Arguments:
               sigxl, sigyl, epsxl, epsyl : labels
               w                           : work amount (w can be 0.4, 1.0 .. and so on)
               ifig
               plot_ind
               marker
               color = 'gray'
               label
               markersize
               alpha = 1.
               ls
               iplot
        """
        ## Calculate work by integrate stress strain curve of the given labels.
        wx, wy, worktot = self.work(sigxl=sigxl, sigyl=sigyl, epsxl=epsxl, epsyl=epsyl, iplot=iplot)

        ## Assign local work variables to global ones
        self.WX = wx
        self.WY = wy
        self.WT = worktot

        ## index where the work is estimated
        ix, ibreak, wrk  = self.__ind_at_workof__(work=worktot, time_flow=self.uni.time_flow, w=w)
        print 'w = ', w
        print 'wrk = ', wrk
        if iplot==True:
            fig = plt.figure(ifig)
            ax = fig.add_subplot(plot_ind ) #, aspect='equal')
        else: pass

        ## Assining requested label stress and stress
        if sigxl=='1': sigx = self.uni.sigx1
        if sigxl=='2': sigx = self.uni.sigx2
        if sigxl=='avg': sigx = self.uni.sigx
        if epsxl=='1': epsx = self.uni.epsx1
        if epsxl=='2': epsx = self.uni.epsx2
        if epsxl=='avg': epsx = self.uni.epsx
        if epsxl=='plavg': epsx = self.uni.Epsx_pl
        if epsxl=='pl11': epsx = self.uni.Epsx_pl11
        if epsxl=='pl12': epsx = self.uni.Epsx_pl12
        if epsxl=='pl21': epsx = self.uni.Epsx_pl21
        if epsxl=='pl22': epsx = self.uni.Epsx_pl22
        
        if sigyl=='1': sigy = self.uni.sigy1
        if sigyl=='2': sigy = self.uni.sigy2
        if sigyl=='avg': sigy = self.uni.sigy
        if epsyl=='1': epsy = self.uni.epsy1
        if epsyl=='2': epsy = self.uni.epsy2
        if epsyl=='avg': epsy = self.uni.epsy
        if epsyl=='plavg': epsy = self.uni.Epsy_pl
        if epsyl=='pl11': epsy = self.uni.Epsy_pl11
        if epsyl=='pl12': epsy = self.uni.Epsy_pl12
        if epsyl=='pl21': epsy = self.uni.Epsy_pl21
        if epsyl=='pl22': epsy = self.uni.Epsy_pl22

        if iplot==True:
            ax.plot(sigx[ix], sigy[ix], marker, markerfacecolor=color,
                    alpha=alpha, label=label, markersize=markersize)
            #plt.ylim(-100,400)
            #plt.xlim(-100,400)
            ax.legend()
            #plt.xticks(np.arange(0,400,100))
            #plt.yticks(np.arange(0,500,100))
        else: pass

        ## interpolation of the stress
        w0, w1 = worktot[ix-1], worktot[ix]
        sx0, sx1 = sigx[ix-1], sigx[ix]
        sy0, sy1 = sigy[ix-1], sigy[ix]

        xslp = (sx1 - sx0) / (w1 - w0)
        yslp = (sy1 - sy0) / (w1 - w0)

        XS = xslp * (w - w0) + sx0
        YS = yslp * (w - w0) + sy0
        
        """
        print 'XS = ', XS
        print 'YS = ', YS
        print 'work, wrk, w', wrk, worktot[ix], w
        """
        if ibreak ==1 :
            #return sigx[ix], sigy[ix], ix,  1
            return XS, YS, ix, 1
        elif ibreak == -1 :
            #return sigx[ix], sigy[ix], ix, -1
            return XS, YS, ix, 1

    def stat(self):
        """
        Shows statistical information for the raw data obtained.
        Currently under construction.
        As of now (2010-Nov 02 no further progesse was made)
        """
        pass

    def __ind_at_workof__(self, work, time_flow, w=0, iplot=True,
                          ifig=1, plot_ind=336):
        """
        Returns index in the reponse to the work level
           and if work level request is possibly over the break_point
           returns ibreak as -1 . (other wise ibreak = 1)
        
        Arguments:
           work : work column
           time_flow : time column
           w : level of the work
           iplot = True : flag whether or not to plot
           ifig = 1
           plot_ind = 336
        """
        time, wrk = self.__time_at_workof__(work, time_flow, w=w,
                                       iplot=iplot, ifig=ifig,
                                       plot_ind=plot_ind)
        ibreak = 1
        for i in range(len(time_flow)):
           if time_flow[i]==time:
              if i > self.uni.break_indx:
                  print
                  print 'Request on seeking index for corresponding work'
                  print 'level is in risk since the index is possibly'
                  print 'in the post break point'
                  print '** FYI, the given work level is ',w
                  #raw_input('Press Enter...')
                  ibreak = -1
              return i, ibreak, wrk
              break
        
    def __time_at_workof__(self, work, time_flow, w, plot_ind=337,
                           iplot=True, ifig=1, ismooth=False):
        """
        Returns the time at the accumulated work that is supposed to be
        calculated from preceding methods embedded in the prior post-processings

        Once 'ismooth' flag is True, 
        the given work column undergoes UnivariateSpline method to interpolate.
        Note that this flag is false as a default.

        Arguments:
           work : work
           time_flow
           w : level of the work
           iplot = True
           ifig = 1
           ismooth = False
        """
        if abs(len(work)-len(time_flow))==1:
           print "Input work level: ", w #, "\n"
           if ismooth==True:
                 xnew, ynew = self.uni.__UnivariateSpline__(x=time_flow[0:-1],
                              y=work, xdown=time_flow[0], xup=time_flow[-2],
                              spacing=len(time_flow))
                 for i in range(len(ynew)):
                   if ynew[i]>=w:
                      break
                 if iplot==True:
                   fig = plt.figure(ifig)
                   ax = fig.add_subplot(plot_ind)
                   ax.plot(xnew, ynew, ls='-', marker='o')
                   ax.plot(time_flow[0:-1], work, ls='-', marker='o')                    
                 return time_flow[i]

           else:
                 for i in range(len(work)):
                      if work[i]>=w: break
                 return time_flow[i], work[i]
        else:
            print "The dimension of work is not properly given"
            print "or is not the one the author conceived"
            print "Thus, __time_at_workof__ method only returns 0"
            return -1, -1

    def __show__(self):
        """
        Matplotlib.pyplot method plot.show becomes active.
        """
        plt.show()

         
    def __default_plots__(self, iplot=False, breakdetect=True,
                          fig_id=20, ind_plot=331):
        """
        Plots various curves for the current data file.
        1. Stress-strain curves
        sigx1-epsx1  sigy1-epsy1 as a default set
        """      
        if breakdetect==True:
            sigx1 = self.__upto__(col=self.uni.sigx1, index=self.uni.break_indx)
            sigx2 = self.__upto__(col=self.uni.sigx2, index=self.uni.break_indx)
            sigy1 = self.__upto__(col=self.uni.sigy1, index=self.uni.break_indx)
            sigy2 = self.__upto__(col=self.uni.sigy2, index=self.uni.break_indx)
            sigy = self.__upto__(col=self.uni.sigy, index=self.uni.break_indx)
            sigx = self.__upto__(col=self.uni.sigx, index=self.uni.break_indx)
            epsx1 = self.__upto__(col=self.uni.epsx1, index=self.uni.break_indx)
            epsx2 = self.__upto__(col=self.uni.epsx2, index=self.uni.break_indx)
            epsy1 = self.__upto__(col=self.uni.epsy1, index=self.uni.break_indx)
            epsy2 = self.__upto__(col=self.uni.epsy2, index=self.uni.break_indx)
            epsx = self.__upto__(col=self.uni.epsx, index=self.uni.break_indx)
            epsy = self.__upto__(col=self.uni.epsy, index=self.uni.break_indx)
        else:
            sigx1, sigx2 = self.uni.sigx1, self.uni.sigx2
            sigy1, sigy2 = self.uni.sigy1, self.uni.sigy2
            sigx , sigy  = self.uni.sigx , self.uni.sigy
            epsx1, epsx2 = self.uni.epsx1, self.uni.epsx2
            epsy1, epsy2 = self.uni.epsy1, self.uni.epsy2
            epsx , epsy  = self.uni.epsx , self.uni.epsy
        #print '** Max time: ',self.uni.break_time,'** Max indx: ',self.uni.break_indx
        plot=self.uni.plot_xy_range
        if iplot==True:
            plot(x=epsx1, y=sigx1, fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            xlabel='True strain', ylabel='True stress (MPa)', label='Sx1-Ex1', 
            title='Stress-strain curves')
            plot(x=epsx2, y=sigx2, fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sx2-Ex2')
            plot(x=epsx1, y=sigx2, fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sx2-Ex1')
            plot(x=epsx2, y=sigx1, fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sx1-Ex2')
            plot(x=epsx , y=sigx , fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sx -Ex ')
            plot(x=epsy1, y=sigy1, fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sy1-Ey1',ls='--')
            plot(x=epsy2, y=sigy2, fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sy2-Ey2',ls='--')
            plot(x=epsy1, y=sigy2, fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sy2-Ey1',ls='--')
            plot(x=epsy2, y=sigy1, fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sy1-Ey2',ls='--')
            plot(x=epsy , y=sigy , fig_id=fig_id+1, ind_plot=111, loc='lower right', 
            label='Sy -Ey ',ls='--')
        elif iplot==False:
            pass
        
    def __smooth__(self, x,y,s=3.0, k=3, nest=-1, mult=10,
                   iplot=True, ind_plot=339, fig_id=1):
        """
        Returns a smoothed line
        Arguments:
            x,y       : x and y columns
            s=3.0     : smoothness parameter (1<= s <= 5)
            k=3       : cube (order)
            nest=-1   :
            mult=10   :
            iplot=True:  Flag whether or not to plot the resulting xnew,ynew
        """
        xnew,ynew = self.uni.__spline__(x=x, y=y, s=s, k=k, nest=nest, mult=mult*len(x))
        if iplot==True:
            if xnew!=0:
                self.uni.plot_xy_range(x=xnew, y=ynew, ind_plot=ind_plot,
                                       fig_id=35, title='smoothed curve')
            else: 
                print 'Your request on curve smoothing failed.'
                print 'It is supposed to be due to fail in spline method you choose'
        return xnew, ynew

    def __upto__(self, col, index=10):
        temp=[]
        for i in range(index):
            temp.append(col[i])
        return temp

    def work(self, sig_label=None, eps_label=None, filename='work.out',
             iplot=False, breakdetect=False, marker='-',
             ifig=20, plot_ind=335,
             sigxl=None, sigyl=None, epsxl=None, epsyl=None):
         """
         Calculates the work (defined as sig_i x eps_i)
         Returns the total work summing x and y axis works up.
         And it also saves the work as a list variable.
         (accessible through calling self.potential)

         Arguments are sig_label, eps_label, filename, iplot, and breakdetect.
         The default arguments for sig_label and eps_label are '1' and '1'.
         sig_label and eps_label are to be given as among followings.
           -'1','2' and 'avg' for sig and eps, respectively.
           -'avg' stands for the average ('1'+'2')/2.

         sigxl,sigyl,epsxl, and epsy1 are optional arguments.
         They are None as default, but once the end-user activates
         More specific sig and eps label set can be given.

         Integration scheme here uses the trapezoidal method.
         Trapezoidal method was implemented from one the Scipy package subsets.
         The package is readily available in the python society.

         Arguments:
             sig_label, eps_label = '1', '1'
             filename='work.out'
             iplot=True
             breakdetect = False,
             marker='-'
             ifig=20, plot_ind=335
             sigxl=None, sigyl=None, epsxl=None, epsyl=None
         """
         if sigxl==None:
              if   sig_label=='1':
                 Sx = self.uni.sigx1
                 Sy = self.uni.sigy1
              elif sig_label=='2':
                 Sx = self.uni.sigx2
                 Sy = self.uni.sigy2
              elif sig_label=='avg':
                 Sx = self.uni.sigx
                 Sy = self.uni.sigy
              else: 
                 print 'Wrong label was given for sig[5]'
                 print 'Returns -1 (1)'
                 raw_input('Press Enter...')
                 return -1
              
              if   eps_label=='1':
                 Ex = self.uni.epsx1
                 Ey = self.uni.epsy1
              elif eps_label=='2':
                 Ex = self.uni.epsx2
                 Ey = self.uni.epsy2
              elif eps_label=='avg':
                 Ex = self.uni.epsx
                 Ey = self.uni.epsy
              else:
                 print 'Wrong label was given for sig[6]'
                 print 'Returns -1 (2)'
                 raw_input('Press Enter...')
                 return -1
                
         elif sigxl!=None:
              if sigxl=='1': Sx = self.uni.sigx1
              elif sigxl=='2': Sx = self.uni.sigx2
              elif sigxl=='avg': Sx = self.uni.sigx
              else:
                   print 'Wrong label was given for sigxl[7]'
                   print 'label?'
                   print 'Returns -1 (3)'
                   raw_input('Press Enter...')
                   return -1
                
              if sigyl=='1': Sy = self.uni.sigy1
              elif sigyl=='2': Sy = self.uni.sigy2
              elif sigyl=='avg': Sy = self.uni.sigy
              else:
                   print 'Wrong label was give for sigyl[8]'
                   print 'Returns -1 (4)'
                   raw_input('Press Enter...')
                   return -1
              
              if epsxl=='1': Ex = self.uni.epsx1
              elif epsxl=='2': Ex = self.uni.epsx2
              elif epsxl=='avg': Ex = self.uni.epsx
              elif epsxl=='plavg' : Ex = self.uni.Epsx_pl
              elif epsxl=='pl11': Ex = self.uni.Epsx_pl11
              elif epsxl=='pl12': Ex = self.uni.Epsx_pl12
              elif epsxl=='pl21': Ex = self.uni.Epsx_pl21
              elif epsxl=='pl22': Ex = self.uni.Epsx_pl22
              else:
                   print 'Wrong label was given for epsxl[9]'
                   print 'Returns -1 (5)'
                   raw_input('Press Enter...')
                   return -1
                
              if epsyl=='1': Ey = self.uni.epsy1
              elif epsyl=='2': Ey = self.uni.epsy2
              elif epsyl=='avg': Ey = self.uni.epsy
              elif epsyl=='plavg' : Ey = self.uni.Epsy_pl
              elif epsyl=='pl11': Ey = self.uni.Epsy_pl11
              elif epsyl=='pl12': Ey = self.uni.Epsy_pl12
              elif epsyl=='pl21': Ey = self.uni.Epsy_pl21
              elif epsyl=='pl22': Ey = self.uni.Epsy_pl22
              else:
                   print 'Wrong label was given for epsyl[10]'
                   print 'Returns -1 (6)'
                   raw_input('Press Enter...')
                   return -1              

         self.SXX=Sx
         self.EXX=Ex
         self.SYY=Sy
         self.EYY=Ey
         print epsyl
         workx = integrate.cumtrapz(y=Sx, x=Ex)
         worky = integrate.cumtrapz(y=Sy, x=Ey)
         workTotal = []
         for i in range(len(workx)):
             workTotal.append(abs(workx[i])+abs(worky[i]))
         if self.path[-1]=='\\': name = self.path+filename
         else : name = self.path+'\\'+filename
         f = file(name, 'w')
         
         f.write('Work per unit volume : integration of the stress-strain curve \n' )
         f.write('time flow(sec)  workx   worky   total work \n')

         time = []
         for i in range(len(self.uni.time_flow)-1):
              time.append(self.uni.time_flow[i]+self.uni.acq_rate*0.001/2.)

         if breakdetect==True: mxind = self.uni.break_indx
         else: mxind = len(workTotal)

         for i in range(mxind-1):
              f.write('%8.3f  %8.3e  %8.3e  %8.3e \n' %(time[i],
                                                        workx[i],
                                                        worky[i],
                                                        workTotal[i]))
         f.close()

         plot = self.uni.plot_xy_range
         if iplot==True:
              #ifig=33, plot_ind=111,
              #xlabel='X(sig:'+sig_label+' ,'+'eps:'+eps_label+')'
              #ylabel='Y(sig:'+sig_label+' ,'+'eps:'+eps_label+')'
              #totlab='total: sig'+sig_label+' eps'+eps_label
              plot(x=time, y=workx, ind_plot=plot_ind,
                   fig_id=ifig, loc='upper left',
                   #xlabel='time flow (sec)',
                   #ylabel='Work per unit vol.(10^6 J/m^3) ',
                   #label=xlabel,
                   breakdetect=breakdetect, ls=marker)
              plot(x=time, y=worky, ind_plot=plot_ind, fig_id=ifig,
                   loc='upper left',
                   #label=ylabel,
                   breakdetect=breakdetect, ls=marker)
              plot(x=time, y=workTotal, ind_plot=plot_ind, fig_id=ifig,
                   loc='upper left',
                   #label=totlab,
                   breakdetect=breakdetect, ls=marker)
         plt.draw()
         return workx, worky, workTotal

    def write_output(self, path=os.getcwd(), filename='biaxial_output.out',
                     label=['xload1','xload2','yload1','yload2',
                            'xstr1','xstr2','ystr1','ystr2'],
                     nind= [0,1,2,3,4,5,6,7,8,-1]):
        
        """
        Writes an output file listing load and str which come from the test as raw data.
        Arguments:
            path
            filename
            label
            nind
        """
        self.uni.write_output(path=path, filename=filename,
                              label=label, nind=nind)

    def monitor_plot(self):
        """
        Plot raw data for monitoring the signals
        """
        fig = plt.figure(20)
        ax1 = fig.add_subplot(111)
        x = self.uni.time_flow
        y1 = self.uni.sigx1
        y2 = self.uni.sigx2
        y3 = self.uni.sigy1
        y4 = self.uni.sigy2

        y5 = self.uni.epsx1
        y6 = self.uni.epsx2
        y7 = self.uni.epsy1
        y8 = self.uni.epsy2

        ax1.plot(x,y1,label='sigx1')
        ax1.plot(x,y2,label='sigx2')
        ax1.plot(x,y3,label='sigy1')
        ax1.plot(x,y4,label='sigy2')
        plt.legend(loc='upper left')
        ax1.set_xlabel('Time flow(sec)')
        ax1.set_ylabel('True Stress (MPa)'  )

        ax2 = ax1.twinx()
        ax2.plot(x,y5,label='epsx1', ls='--')
        ax2.plot(x,y6,label='epsx2', ls='--')
        ax2.plot(x,y7,label='epsy1', ls='--')
        ax2.plot(x,y8,label='epsy2', ls='--')
        ax2.set_ylabel('True Strain')
        plt.legend(loc='upper center')
        plt.draw()




















##################################################################################
#
#          class uniaxial
#               developed initially for the uniaxial test post-processing from
#              Kokusai in-plane biaxial process.
#               Now it is extended to span general biaxial testings.
#
##################################################################################

class uniaxial:
    """
    post-process uniaxial results (flow stress and R-values)
    optional but important arguments : path, filename
    defaults for the arguments are below:
    """
    default_filename = "2nd.csv"

    def __init__(self, low_sigx, up_sigx, low_sigy, up_sigy,
                 path=os.getcwd(), filename=default_filename,
                 iplot=False, modulusx = None, modulusy=None,
                 iflip = False, str_channel = None, maxtime=None):
        """
        Class uniaixal __init__ def.
        This def includes opening file, read the source data, 
        handling header of it, and extracting columns of data.
        One may ask for intrinsic plottings if necessary.

        time_flow : The time flow to the data calculated from acqusition rate

        Arguments:
           low_sigx
           up_sigx
           low_sigy
           up_sigy
           path
           filename
           iplot=False
           modulusx
           modulusy
           iflip
           str_channel
        """
        print 
        print "*  The given file: '", filename,"'",
        print "at the path of ...", path[-15:len(path)]

        self.filename = filename
        self.path = path
        self.ipyplot = iplot
        self.str_channel = str_channel   # str_channel is inherited as a global value for class uniaxial

        # if path[-1]=='\\': f = file(path + filename)
        # else: f = file(path + os.sep + filename)
        f = open(filename, 'r')
        source = f.read()
        f.close()
        self.source = source
        self.lines = source.split('\n')

        # Header part is returned
        self.header, self.h_line = self.__get_header__(self.lines)
        col_labels = self.lines[self.h_line + 1]
        self.acq_rate = float(self.__get_param_header__
                          (label='Sampling Rate(ms)',
                           source=self.header))
        self.stressratio = []
        self.strainratio = []
        for i in self.header:
            if i.split(',')[0]=='Stress Ratio':
                self.stressratio.append(int(i.split(',')[1]))
                self.stressratio.append(int(i.split(',')[2]))
                self.strainratio.append(0)
                self.strainratio.append(0)
            if i.split(',')[0]=='All Strain Ratio':
                self.strainratio.append(int(i.split(',')[1]))
                self.strainratio.append(int(i.split(',')[2]))
                self.stressratio.append(0)
                self.stressratio.append(0)
        if iflip==True:
            temp = self.stressratio[0]
            self.stressratio[0] = self.stressratio[1]
            self.stressratio[1] = temp
            temp = self.strainratio[0]
            self.strainratio[0] = self.strainratio[1]
            self.strainratio[1] = temp
        ctrl, axiality = self.__ratio__(sigma=self.stressratio, epsilon=self.strainratio)
        self.ctrl = ctrl
        self.axiality = axiality
        print 
        print '*******'
        print axiality
        print '*******'
        print 
             
        # ---------------------------------------------
        # Return column values with the indices of'em.
        self.__init_columns__(modulusx=modulusx, modulusy=modulusy,
                              low_sigx=low_sigx, up_sigx=up_sigx,
                              low_sigy=low_sigy, up_sigy=up_sigy,
                              axiality=axiality, ctrl=ctrl, iflip=iflip,
                              maxtime=maxtime)

        if iplot==True:
            print 
            print '******************************'
            print ' Intrinsic plot flag is True  '
            print ' Intrinsic plots are selected '
            print '******************************'
            print 
            self.__intrinsic_plots()
        """
        print '____________________'
        print 'Plot def explanation'
        print 'plot_xy_range()'
        print ' --arguments'
        print '   x,y            : x (horizontal), y (vertical) axis data'
        print '   ind_plot       : the index of matplotlib.pyplot model'
        print '   fig_id         : id # of the figure'
        print '   xrange,yrange  : list of two elements range of x and y'
        print "   loc            : location of label (e.g. 'upper right') "
        print "   xlabel, ylabel : labels of two axes (e.g. 'Time flow')  "
        print "   label          : Legend of the curve (e.g. '316L')      "
        print "   title          : title of the axes                      "
        print ''
        print "write_output()"
        print ' --arguments'
        print '   filename, path'
        """
    def __ratio__(self, sigma, epsilon):
        """
        Returns either 'sig_ctrl' or 'eps_ctrl'
        as well as the axiaility: 'uniaxial' or 'biaxial'
        """
        if any(sigma[i]!=0 for i in range(2)) :
            if any(sigma[i]==0 for i in range(2)):
                print
                print '*********************'
                print 'Uniaxial load control'
                print '*********************'
                if sigma[0]==0: return 'sig_ctrl', 'uniy'
                elif sigma[1]==0: return 'sig_ctrl', 'unix'
                else: return 'sig_ctrl', -1
            else:
                print
                print '******************'
                print 'Load ratio control'
                print '******************'
                return 'sig_ctrl', 'biaxial'

        elif any(epsilon[i]!=0 for i in range(2)) :
            if epsilon[0] == 0 : #any(epsilon[i]==0 for i in range(2)):
                print
                print '*******************************'
                print 'Plane strain control (eps[0]=0)'
                print '*******************************'
                return 'eps_ctrl', 'uniy'
            elif epsilon[1] == 0:
                print
                print '*******************************'
                print 'Plane strain control (eps[1]=0)'
                print '*******************************'
                return 'eps_ctrl', 'unix'
            else:
                print
                print '********************'
                print 'Strain ratio control'
                print '********************'
                return 'eps_ctrl', 'biaxial'
        else: return -1, -1

    def __intrinsic_plots(self, fig_id=10):
        """
        Plots some common data
        """
        # ------------------
        #Intrinsic plots
        self.plot_str_str(stress=self.x1stress_engi, strain=self.xstr1, 
                       label='x1 engineering', ind_plot=321, fig_id=fig_id)
        self.plot_str_str(stress=self.x1stress_true, strain=self.x1str_true, 
                       label='x1 true', ind_plot=322,fig_id=fig_id)
        self.plot_str_str(stress=self.x2stress_engi, strain=self.xstr1,
                       label='x2 engineering', ind_plot=321, fig_id=fig_id)
        self.plot_str_str(stress=self.x2stress_true, strain=self.x1str_true,
                       label='x2 true', ind_plot=322, fig_id=fig_id)
        self.plot_R(l_str=self.x1str_true, w_str=self.y1str_true, 
                 label='x1y1 R', ind_plot=323, fig_id=fig_id)
        self.plot_R(l_str=self.x1str_true, w_str=self.y2str_true, 
                 label='x1y2 R', ind_plot=323, fig_id=fig_id)

        self.plot_xy(x=self.time_flow, y=self.xload1, ind_plot=324,
                  fig_id=fig_id, label='xload1')
        self.plot_xy(x=self.time_flow, y=self.xload2, ind_plot=324,
                  fig_id=fig_id, label='xload2')
        devi=[]
        for i in range(len(self.xload1)):
         devi.append( self.xload1[i]-self.xload2[i]  )

        rel_devi=[]
        for i in range(len(devi)):
         if  self.xload1[i] == 0 : rel_devi.append(0)
         else:rel_devi.append(devi[i]/self.xload1[i])

        self.plot_xy_range (x=self.time_flow, y=devi, ind_plot=325, 
                        fig_id=fig_id, label='deviation btwn 1&2',
                        xrange=[10,160], loc='upper center')
        self.plot_xy_range (x=self.time_flow, y=rel_devi , ind_plot=326, 
                        fig_id=fig_id, label='rel_devi', xrange=[10,160], 
                        loc='upper left', ylabel='Relative deviation',
                        xlabel='time (s)' )


        ###
        # raw data curves like the ones simultaneously plotted
        # in the course of experiments:
        #      Load   x1, x2, y1, y2
        #      Strain x1, x2, y1, y2
        # with respect to the testing time flow.
        self.plot_xy(x=self.time_flow, y=self.xload1, ind_plot=111, fig_id = 63, label ='xload1')
        self.plot_xy(x=self.time_flow, y=self.xload2, ind_plot=111, fig_id = 63, label ='xload2')
        self.plot_xy(x=self.time_flow, y=self.yload1, ind_plot=111, fig_id = 63, label ='yload1')
        self.plot_xy(x=self.time_flow, y=self.yload2, ind_plot=111, fig_id = 63, label ='yload2')
        self.plot_xy(x=self.time_flow, y=self.xstr1,  ind_plot=111, fig_id = 64, label ='xstr1')
        self.plot_xy(x=self.time_flow, y=self.xstr2,  ind_plot=111, fig_id = 64, label ='xstr2')
        self.plot_xy(x=self.time_flow, y=self.ystr1,  ind_plot=111, fig_id = 64, label ='ystr1')
        self.plot_xy(x=self.time_flow, y=self.ystr2,  ind_plot=111, fig_id = 64, label ='ystr2')

        """
        self.xload1 = self.__get_col__(self.lines, self.h_line, 0)
        self.xload2 = self.__get_col__(self.lines, self.h_line, 1)
        self.yload1 = self.__get_col__(self.lines, self.h_line, 2)
        self.yload2 = self.__get_col__(self.lines, self.h_line, 3)
        self.xstr1  = self.__get_col__(self.lines, self.h_line, 4)
        self.xstr2  = self.__get_col__(self.lines, self.h_line, 5)
        self.ystr1  = self.__get_col__(self.lines, self.h_line, 6)
        self.ystr2  = self.__get_col__(self.lines, self.h_line, 7)
        """
                                         

    def __spline__(self, x, y, s=3.0, k=3, nest=-1, mult=1):
         """
         Returns interpolated spline calculated
         Arguments:
            s=3.0   : smoothness parameter
            k=3     : spline order (default: cube)
            nest=-1 : estimate of number of knots needed (-1=maximal)
            mult=10 : Multiplication factor for the points #
         """
         try :
            tckp, u = splprep([x,y], s=s, k=k, nest=nest)
         except SystemError:
            print "************************* Exception occurs! *****************************"
            print "* SystemError Occured in '__spline__' method under the class 'uniaxial' *"
            print "* Accordingly, it only returns zeroes                                   *" 
            print "*************************************************************************"
            return 0,0
         else:
            xnew, ynew = splev(linspace(0, 1, mult), tckp)
            return xnew, ynew  

    def __UnivariateSpline__(self, x, y, xdown, xup, s=1, spacing=1000, k=3):
        """
        One-dimensional smoothing spline fit to a given set of data points.
        *implemented from scipy.interpolate package

        Returns xs, ys which are sampled

        Arguments:
        x, y : sequence data for x and y. 
               x must be increasing
        s=1  : positive smoothing factor
        xdown, xup : sampled x limit
        spacing : spacing to be used for sample
        """
        s = UnivariateSpline(x, y, s=s, k=k)
        xs = linspace(xdown, xup, spacing)
        ys = s(xs)
        return xs, ys
    
    def __Eps_pl__(self, strain, stress,
                   lowstress = 20., upstress = 80.,
                   iplot=False, ifig=20, plot_ind=339,
                   modulus=None):
        """
        Decomposes the strain into elastic and plastic
        Returns plastic strain, and elastic strain(E_tot - E_pl)
        Arguments:
            strain
            stress
            lowstress = 20.
            upstress = 50.

        if modulus is given, it is to be used.
        """
        """
        Sig/Modulus = E_tot - E_pl
        -> E_pl = E_tot - Sig/Modulus
        if E_pl < 0 then just regard this as 0 (pure elastic regime)
        """
        if modulus == None:
            # getting modulus firstly.
            i_low = self.__indSeek__(value = lowstress, col = stress)
            i_up = self.__indSeek__(value = upstress, col = stress)
            modulus = self.__slope__(x = strain[i_low:i_up], y = stress[i_low:i_up])
            #print 'Modulus calculated = ', modulus
            if modulus < 0 :
                inegloading = True
            else: inegloading = False
        else :
            #print 'Modulus given = ' , modulus
            if modulus < 0:
                inegloading = True
            else: inegloading = False

            
        E_pl = []
        E_el = []
        for i in range(len(strain)):
            temp = strain[i] - stress[i] / modulus

            
            ## below is the consideration of imposing different treatment for different loadings.
            ## but it is decided not to make things complex here thus it is just suspended for the moment.
            ## Nevertheless, it is expected that final results may not vary much.
            """
            if inegloading == True:
                if temp <= 0.:
                   E_pl.append(temp)
                else:
                   E_pl.append(0.)
            else:  # Thus positive loading
                if temp <= 0.:
                   E_pl.append(0.)
                else:
                   E_pl.append(temp)
            """
            E_pl.append(temp)            
            E_el.append(strain[i] - E_pl[i])
            
        if iplot==True:
            fig = plt.figure(ifig)
            ax = fig.add_subplot(plot_ind)
            ax.plot(E_pl, stress)
            ax.set_xlabel('Plastic strain')
            ax.set_ylabel('Stress')
            ax.set_title('Str-Str curve')
            
        return E_pl, E_el, modulus
    
    def __slope__(self, x, y, iplot=False, ifig=30, plot_ind=339):
        """
        Obtain the slope of y w.r.t x
        """
        z = np.polyfit(x,y,1)        
        if iplot==True:

            fig = plt.figure(ifig)
            ax = fig.add_subplot(plot_ind)
            ax.plot(x,y,marker='o')
            ax.set_title('Slope fitting requested representation')
        else: pass
        y1 = []
        for i in range(len(x)):
           y1.append(z[0]*x[i]+z[1])
        if iplot==True: ax.plot(x,y1)
        else: pass
        return z[0]
    
    def __indSeek__(self,value,col):
        """
        returns the corresponding index close to the stress given
        """
        i = 0
        while True:
           
           try:
              if col[i]>=value: pass
           except IndexError:break
           else:
              if col[i]>=value:
                 return i
                 break
              else: i = i + 1
            
           """
           if col[i] >= value:
               return i
           elif i > 100000:
               print "Too many interation in __ind_exceeding_value_of__"
               break
           else: i = i + 1
           """

    def __make_zero__(self,column=None):
        temp = []
        for i in range(len(column)):
           temp.append(0)
        return temp
    
    def __flip__(self, col1, col2):
        return col2, col1
    
    def __init_columns__(self, modulusx=None, modulusy=None, axiality=None, ctrl=None,
                         low_sigx=20., up_sigx=50., low_sigy=20., up_sigy=50., iflip=None,
                         maxtime=None):
        """
        Sets default column values for some post-processed data.
        """
        # Header info for post process, (e.g. cross section area, width etc.)
        Ay=self.__cross_sectionY__(self.header)
        Ax=self.__cross_sectionX__(self.header)
        self.Ay=Ay
        self.Ax=Ax
        if iflip ==True:
            temp = self.Ay
            self.Ay = self.Ax
            self.Ax = temp

            
        tmp_col = self.__get_col__(self.lines, self.h_line, 0)
        self.time_flow = self.__time_xcol__(rate=self.acq_rate, 
                                         size=(len(tmp_col)) )
        # ---------------------------------------------
        # Return column values with index of it.
        # if maxtime is given, the columns size needs be constrained accordingly.
        
        if maxtime==None:
           self.xload1 = self.__get_col__(self.lines, self.h_line, 0)
           self.xload2 = self.__get_col__(self.lines, self.h_line, 1)
           self.yload1 = self.__get_col__(self.lines, self.h_line, 2)
           self.yload2 = self.__get_col__(self.lines, self.h_line, 3)
           self.xstr1  = self.__get_col__(self.lines, self.h_line, 4)
           self.xstr2  = self.__get_col__(self.lines, self.h_line, 5)
           self.ystr1  = self.__get_col__(self.lines, self.h_line, 6)
           self.ystr2  = self.__get_col__(self.lines, self.h_line, 7)
        else:
           i = 0
           while True:
              if self.time_flow[i] > maxtime:
                 idx = i - 1
                 if idx < 0: idx = 0
                 break
              else:
                 i = i + 1
                 pass
           self.xload1 = self.__get_col__(self.lines, self.h_line, 0)[0:idx]
           self.xload2 = self.__get_col__(self.lines, self.h_line, 1)[0:idx]
           self.yload1 = self.__get_col__(self.lines, self.h_line, 2)[0:idx]
           self.yload2 = self.__get_col__(self.lines, self.h_line, 3)[0:idx]
           self.xstr1  = self.__get_col__(self.lines, self.h_line, 4)[0:idx]
           self.xstr2  = self.__get_col__(self.lines, self.h_line, 5)[0:idx]
           self.ystr1  = self.__get_col__(self.lines, self.h_line, 6)[0:idx]
           self.ystr2  = self.__get_col__(self.lines, self.h_line, 7)[0:idx]
           
        self.time_flow = self.__time_xcol__(rate=self.acq_rate, 
                                         size=(len(self.xload1)))
        print 'maximum time of this file is ', self.time_flow[-1]





        if iflip==True:
            self.xload1, self.yload1 = self.__flip__(self.xload1, self.yload1)
            self.xload2, self.yload2 = self.__flip__(self.xload2, self.yload2)
            self.xstr1 , self.ystr1   = self.__flip__(self.xstr1, self.ystr1)
            self.xstr2 , self.ystr2   = self.__flip__(self.xstr2, self.ystr2)

        # Axiality tells about whether or not the test was performed under biaxial or uniaxial condition
        if axiality[0:3]=='uni':
           if ctrl=='eps_ctrl':
              pass
              # str_channel
              if axiality[-1]=='x':    # in case of 'unix'
                 self.ystr2 = self.__get_col__(self.lines, self.h_line, 5)   # ystr2 <- xstr2 
                 self.xstr2 = self.__get_col__(self.lines, self.h_line, 7)   # xstr2 <- ystr2
                 self.xstr1 = self.__get_col__(self.lines, self.h_line, 4)   # xstr1 <- xstr1
              elif axiality[-1]=='y':  # in case of 'uniy'
                 self.xstr2 = self.__get_col__(self.lines, self.h_line, 7)   # xstr2 <- ystr2
                 self.ystr2 = self.__get_col__(self.lines, self.h_line, 5)   # ystr2 <- xstr2
                 self.ystr1 = self.__get_col__(self.lines, self.h_line, 6)   # ystr1 <- ystr1
              else: pass
        

        # ------------------
        # True strains
        # strain must be in the unit of ust (strain x 10^-6)
        self.x1str_true = self.__strain_true__(str=self.xstr1)  
        self.x2str_true = self.__strain_true__(str=self.xstr2) 
        self.y1str_true = self.__strain_true__(str=self.ystr1) 
        self.y2str_true = self.__strain_true__(str=self.ystr2)
        
        self.epsx1 = self.x1str_true
        self.epsx2 = self.x2str_true
        self.epsy1 = self.y1str_true
        self.epsy2 = self.y2str_true

        # ------------------
        # Engineering Stress
        self.x1stress_engi = self.__stress_engi__(load=self.xload1, 
                                               cross_section=self.Ax)
        self.x2stress_engi = self.__stress_engi__(load=self.xload2, 
                                               cross_section=self.Ax)
        self.y1stress_engi = self.__stress_engi__(load=self.yload1, 
                                               cross_section=self.Ay)
        self.y2stress_engi = self.__stress_engi__(load=self.yload2, 
                                               cross_section=self.Ay)
        self.SIGx1=self.x1stress_engi
        self.SIGx2=self.x2stress_engi
        self.SIGy1=self.y1stress_engi
        self.SIGy2=self.y2stress_engi

        # ------------------
        # True Stress
        self.x1stress_true = self.__stress_true__(self.x1stress_engi, 
                                               self.xstr1)
        self.x2stress_true = self.__stress_true__(self.x2stress_engi, 
                                               self.xstr2)
        self.y1stress_true = self.__stress_true__(self.y1stress_engi,
                                               self.ystr1)
        self.y2stress_true = self.__stress_true__(self.y2stress_engi,
                                               self.ystr2)
        self.sigx1=self.x1stress_true
        self.sigx2=self.x2stress_true
        self.sigy1=self.y1stress_true
        self.sigy2=self.y2stress_true


        if axiality=='biaxial':
            pass
        elif axiality[0:3]=='uni':
            if axiality[-1]=='x':
                if ctrl=='sig_ctrl':
                    self.yload1 = self.__make_zero__(column = self.yload1)
                    self.yload2 = self.__make_zero__(column = self.yload2)
                    self.SIGy1 = self.__make_zero__(column = self.SIGy1)
                    self.SIGy2 = self.__make_zero__(column = self.SIGy2)
                    self.sigy1 = self.__make_zero__(column = self.sigy1)
                    self.sigy2 = self.__make_zero__(column = self.sigy2)
                elif ctrl=='eps_ctrl':
                    pass
                    #self.ystr1 = self.__make_zero__(column = self.ystr1)
                    #self.ystr2 = self.__make_zero__(column = self.ystr2)
                    #self.epsy1 = self.__make_zero__(column = self.epsy1)
                    #self.epsy2 = self.__make_zero__(column = self.epsy2)                   
                else: print 'Err: Wrong ctrl'
            elif axiality[-1]=='y':
                if ctrl=='sig_ctrl':
                    self.xload1 = self.__make_zero__(column = self.xload1)
                    self.xload2 = self.__make_zero__(column = self.xload2)
                    self.SIGx1 = self.__make_zero__(column = self.SIGx1)
                    self.SIGx2 = self.__make_zero__(column = self.SIGx2)
                    self.sigx1 = self.__make_zero__(column = self.sigx1)
                    self.sigx2 = self.__make_zero__(column = self.sigx2)
                elif ctrl=='eps_ctrl':
                    pass
                    #self.xstr1 = self.__make_zero__(column = self.xstr1)
                    #self.xstr2 = self.__make_zero__(column = self.xstr2)
                    #self.epsx1 = self.__make_zero__(column = self.epsx1)
                    #self.epsx2 = self.__make_zero__(column = self.epsx2)
                else: print 'Err: Wrong ctrl'
            else: print 'Err: Wrong axiality string error occured in __init__columns__'
        else:
            print 'Err: Wrong axiality string passed to class uniaxial'
        
        # Average terms of x and y loads and strain.
        self.xload_avg = self.__get_data_avg__(self.xload1, self.xload2)
        self.yload_avg = self.__get_data_avg__(self.yload1, self.yload2)
        self.xstr_avg  = self.__get_data_avg__(self.xstr1, self.xstr2)
        self.ystr_avg  = self.__get_data_avg__(self.ystr1, self.ystr2)
        self.sigx = self.__get_data_avg__(self.sigx1, self.sigx2)
        self.sigy = self.__get_data_avg__(self.sigy1, self.sigy2)
        self.epsx = self.__get_data_avg__(self.epsx1, self.epsx2)
        self.epsy = self.__get_data_avg__(self.epsy1, self.epsy2)
        self.SIGx = self.__get_data_avg__(self.SIGx1, self.SIGx2)
        self.SIGy = self.__get_data_avg__(self.SIGy1, self.SIGy2)


        # Some plastic elastic decompositions
        if ctrl=='sig_ctrl':
            if axiality=='unix':
                #print 'axiality ==', axiality, 'unix'               
                ### x strains

                #average
                self.Epsx_pl, self.Epsx_el, mod = self.__Eps_pl__(strain=self.epsx,
                                                          stress=self.sigx,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                          modulus = modulusx)
                if modulusx==None: print 'x Modulus (pl) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl) =',round(mod/10**3, 2) ,'GPa'

                #pl11
                self.Epsx_pl11, self.Epsx_el11, mod = self.__Eps_pl__(strain=self.epsx1,
                                                              stress=self.sigx1,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                              modulus = modulusx)
                if modulusx==None: print 'x Modulus (pl11) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl11) =', round(mod/10**3, 2) ,'GPa'

                #pl12
                self.Epsx_pl12, self.Epsx_el12, mod = self.__Eps_pl__(strain=self.epsx1,
                                                              stress=self.sigx2,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                              modulus = modulusx)
                if modulusx==None: print 'x Modulus (pl12) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl12) =', round(mod/10**3, 2) ,'GPa'

                #pl21
                self.Epsx_pl21, self.Epsx_el21, mod = self.__Eps_pl__(strain=self.epsx2,
                                                              stress=self.sigx1,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                              modulus = modulusx)
                if modulusx==None: print 'x Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl21) =', round(mod/10**3, 2) ,'GPa'

                #pl22
                self.Epsx_pl22, self.Epsx_el22, mod = self.__Eps_pl__(strain=self.epsx2,
                                                              stress=self.sigx2,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                              modulus=modulusx)
                if modulusx==None: print 'x Modulus (pl22) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl22) =', round(mod/10**3, 2) ,'GPa'
                
                tmp_e = self.__make_zero__(column=self.Epsx_pl)
                self.Epsy_pl, Epsy_el = tmp_e, tmp_e
                self.Epsy_pl11, self.Epsy_pl12, self.Epsy_pl21, self.Epsy_pl22 = tmp_e,tmp_e,tmp_e,tmp_e
                modulusy = 1.
                
            elif axiality=='uniy':
                ### y strains
                
                #average
                print 'axiality ==', axiality, 'uniy'
                self.Epsy_pl, self.Epsy_el, mod = self.__Eps_pl__(strain=self.epsy,
                                                          stress=self.sigy,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                          modulus=modulusy)

                if modulusy==None: print 'y Modulus (pl) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl) =', round(mod/10**3, 2) ,'GPa'

                #pl11
                self.Epsy_pl11, self.Epsy_el11, mod = self.__Eps_pl__(strain=self.epsy1,
                                                              stress=self.sigy1,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                              modulus=modulusy)
                if modulusy==None: print 'y Modulus (pl11) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl11) =', round(mod/10**3, 2) ,'GPa'

                #pl12
                self.Epsy_pl12, self.Epsy_el12, mod = self.__Eps_pl__(strain=self.epsy1,
                                                              stress=self.sigy2,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                              modulus=modulusy)
                if modulusy==None: print 'y Modulus (pl12) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl12) =', round(mod/10**3, 2) ,'GPa'

                #pl21
                self.Epsy_pl21, self.Epsy_el21, mod = self.__Eps_pl__(strain=self.epsy2,
                                                              stress=self.sigy1,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                              modulus=modulusy)
                if modulusy==None: print 'y Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl21) =', round(mod/10**3, 2) ,'GPa'

                #pl22
                self.Epsy_pl22, self.Epsy_el22, mod = self.__Eps_pl__(strain=self.epsy2,
                                                              stress=self.sigy2,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                              modulus=modulusy)
                if modulusy==None: print 'y Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl21) =',round(mod/10**3, 2) ,'GPa'

                tmp_e = self.__make_zero__(column=self.Epsy_pl)
                self.Epsx_pl, Epsx_el = tmp_e, tmp_e
                self.Epsx_pl11, self.Epsx_pl12, self.Epsx_pl21, self.Epsx_pl22 = tmp_e,tmp_e,tmp_e,tmp_e
                modulusy = 1.

                
            elif axiality=='biaxial':
                ### x strains


                # Names of Epsx_pl sometimes have a trailing 
                
                #average
                self.Epsx_pl, self.Epsx_el, mod = self.__Eps_pl__(strain=self.epsx,
                                                          stress=self.sigx,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                          modulus = modulusx)
                if modulusx==None: print 'x Modulus (pl) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl) =',round(mod/10**3, 2) ,'GPa'

                #pl1x  (using the averaged x-stress)
                self.Epsx_pl1x, self.Epsx_el1x, mod = self.__Eps_pl__(strain=self.epsx1,
                                                                  stress=self.sigx,
                                                                  lowstress=low_sigx,
                                                                  upstress=up_sigx,
                                                                  modulus=modulusx)
                if modulusx==None: print 'x Modulus (pl1x) = ', round(mod/10**3, 2), 'GPa'
                else: print 'Given x Modulus (pl1x) =', round(mod/10**3, 2), 'GPa'

                #pl2x (using the averaged x-stress)
                self.Epsx_pl2x, self.Epsx_el2x, mod = self.__Eps_pl__(strain=self.epsx2,
                                                                  stress=self.sigx,
                                                                  lowstress=low_sigx,
                                                                  upstress=up_sigx,
                                                                  modulus=modulusx)
                if modulusx==None: print'x Modulus (pl2x) = ', round(mod/10**3, 2), 'GPa'
                else: print 'Given x Modulus (pl2x) =', round(mod/10**3, 2), 'GPa'

                #pl11
                self.Epsx_pl11, self.Epsx_el11, mod = self.__Eps_pl__(strain=self.epsx1,
                                                              stress=self.sigx1,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                              modulus = modulusx)
                if modulusx==None: print 'x Modulus (pl11) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl11) =', round(mod/10**3, 2) ,'GPa'

                #pl12
                self.Epsx_pl12, self.Epsx_el12, mod = self.__Eps_pl__(strain=self.epsx1,
                                                              stress=self.sigx2,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                              modulus = modulusx)
                if modulusx==None: print 'x Modulus (pl12) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl12) =', round(mod/10**3, 2) ,'GPa'

                #pl21
                self.Epsx_pl21, self.Epsx_el21, mod = self.__Eps_pl__(strain=self.epsx2,
                                                              stress=self.sigx1,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                              modulus = modulusx)
                if modulusx==None: print 'x Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl21) =', round(mod/10**3, 2) ,'GPa'

                #pl22
                self.Epsx_pl22, self.Epsx_el22, mod = self.__Eps_pl__(strain=self.epsx2,
                                                              stress=self.sigx2,
                                                          lowstress=low_sigx,
                                                          upstress=up_sigx,
                                                              modulus=modulusx)
                if modulusx==None: print 'x Modulus (pl22) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given x Modulus (pl22) =', round(mod/10**3, 2) ,'GPa'
                
                ### y strains
                #average
                self.Epsy_pl, self.Epsy_el, mod = self.__Eps_pl__(strain=self.epsy,
                                                          stress=self.sigy,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                          modulus=modulusy)        
                if modulusy==None: print 'y Modulus (pl) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl) =', round(mod/10**3, 2) ,'GPa'

                #pl1y  (using the averaged y-stress)
                self.Epsy_pl1y, self.Epsy_el1y, mod = self.__Eps_pl__(strain=self.epsy1,
                                                                      stress=self.sigy,
                                                                      lowstress=low_sigy,
                                                                      upstress=up_sigy,
                                                                      modulus=modulusy)
                if modulusy==None: print 'y modulus (pl1y) = ', round(mod/10**3, 2) ,'GPa'
                else: print 'Given y modulus (pl1y) =', round(mod/10**3, 2) ,'GPa'

                #pl2y (using the averaged y-stress)
                self.Epsy_pl2y, self.Epsy_el2y, mod = self.__Eps_pl__(strain=self.epsy2,
                                                                      stress=self.sigy,
                                                                      lowstress=low_sigy,
                                                                      upstress=up_sigy,
                                                                      modulus=modulusy)
                if modulusy==None: print 'y modulus (pl2y) = ', round(mod/10**3, 2) ,'GPa'
                else: print 'Given y modulus (pl2y) =', round(mod/10**3, 2) ,'GPa'

                #pl11
                self.Epsy_pl11, self.Epsy_el11, mod = self.__Eps_pl__(strain=self.epsy1,
                                                              stress=self.sigy1,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                              modulus=modulusy)
                if modulusy==None: print 'y Modulus (pl11) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl11) =', round(mod/10**3, 2) ,'GPa'

                #pl12
                self.Epsy_pl12, self.Epsy_el12, mod = self.__Eps_pl__(strain=self.epsy1,
                                                              stress=self.sigy2,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                              modulus=modulusy)
                if modulusy==None: print 'y Modulus (pl12) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl12) =', round(mod/10**3, 2) ,'GPa'

                #pl21
                self.Epsy_pl21, self.Epsy_el21, mod = self.__Eps_pl__(strain=self.epsy2,
                                                              stress=self.sigy1,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                              modulus=modulusy)
                if modulusy==None: print 'y Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl21) =', round(mod/10**3, 2) ,'GPa'

                #pl22
                self.Epsy_pl22, self.Epsy_el22, mod = self.__Eps_pl__(strain=self.epsy2,
                                                              stress=self.sigy2,
                                                               lowstress=low_sigy,
                                                               upstress=up_sigy,
                                                              modulus=modulusy)
                if modulusy==None: print 'y Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
                else : print 'Given y Modulus (pl21) =',round(mod/10**3, 2) ,'GPa'


        elif ctrl=='eps_ctrl':
            ### x strains
            #average
            self.Epsx_pl, self.Epsx_el, mod = self.__Eps_pl__(strain=self.epsx,
                                                      stress=self.sigx,
                                                      lowstress=low_sigx,
                                                      upstress=up_sigx,
                                                      modulus = modulusx)
            if modulusx==None: print 'x Modulus (pl) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given x Modulus (pl) =',round(mod/10**3, 2) ,'GPa'

            #pl11
            self.Epsx_pl11, self.Epsx_el11, mod = self.__Eps_pl__(strain=self.epsx1,
                                                          stress=self.sigx1,
                                                      lowstress=low_sigx,
                                                      upstress=up_sigx,
                                                          modulus = modulusx)
            if modulusx==None: print 'x Modulus (pl11) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given x Modulus (pl11) =', round(mod/10**3, 2) ,'GPa'

            #pl12
            self.Epsx_pl12, self.Epsx_el12, mod = self.__Eps_pl__(strain=self.epsx1,
                                                          stress=self.sigx2,
                                                      lowstress=low_sigx,
                                                      upstress=up_sigx,
                                                          modulus = modulusx)
            if modulusx==None: print 'x Modulus (pl12) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given x Modulus (pl12) =', round(mod/10**3, 2) ,'GPa'

            #pl21
            self.Epsx_pl21, self.Epsx_el21, mod = self.__Eps_pl__(strain=self.epsx2,
                                                          stress=self.sigx1,
                                                      lowstress=low_sigx,
                                                      upstress=up_sigx,
                                                          modulus = modulusx)
            if modulusx==None: print 'x Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given x Modulus (pl21) =', round(mod/10**3, 2) ,'GPa'

            #pl22
            self.Epsx_pl22, self.Epsx_el22, mod = self.__Eps_pl__(strain=self.epsx2,
                                                          stress=self.sigx2,
                                                      lowstress=low_sigx,
                                                      upstress=up_sigx,
                                                          modulus=modulusx)
            if modulusx==None: print 'x Modulus (pl22) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given x Modulus (pl22) =', round(mod/10**3, 2) ,'GPa'

            ### y strains
            #average
            self.Epsy_pl, self.Epsy_el, mod = self.__Eps_pl__(strain=self.epsy,
                                                      stress=self.sigy,
                                                           lowstress=low_sigy,
                                                           upstress=up_sigy,
                                                      modulus=modulusy)        
            if modulusy==None: print 'y Modulus (pl) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given y Modulus (pl) =', round(mod/10**3, 2) ,'GPa'

            #pl11
            self.Epsy_pl11, self.Epsy_el11, mod = self.__Eps_pl__(strain=self.epsy1,
                                                          stress=self.sigy1,
                                                           lowstress=low_sigy,
                                                           upstress=up_sigy,
                                                          modulus=modulusy)
            if modulusy==None: print 'y Modulus (pl11) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given y Modulus (pl11) =', round(mod/10**3, 2) ,'GPa'

            #pl12
            self.Epsy_pl12, self.Epsy_el12, mod = self.__Eps_pl__(strain=self.epsy1,
                                                          stress=self.sigy2,
                                                           lowstress=low_sigy,
                                                           upstress=up_sigy,
                                                          modulus=modulusy)
            if modulusy==None: print 'y Modulus (pl12) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given y Modulus (pl12) =', round(mod/10**3, 2) ,'GPa'

            #pl21
            self.Epsy_pl21, self.Epsy_el21, mod = self.__Eps_pl__(strain=self.epsy2,
                                                          stress=self.sigy1,
                                                           lowstress=low_sigy,
                                                           upstress=up_sigy,
                                                          modulus=modulusy)
            if modulusy==None: print 'y Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given y Modulus (pl21) =', round(mod/10**3, 2) ,'GPa'

            #pl22
            self.Epsy_pl22, self.Epsy_el22, mod = self.__Eps_pl__(strain=self.epsy2,
                                                          stress=self.sigy2,
                                                           lowstress=low_sigy,
                                                           upstress=up_sigy,
                                                          modulus=modulusy)
            if modulusy==None: print 'y Modulus (pl21) = ', round(mod/10**3, 2) ,'GPa'
            else : print 'Given y Modulus (pl21) =',round(mod/10**3, 2) ,'GPa'
        else: print 'Wrong ctrl string pass to __init_columns__'




        #print self.Epsy_pl[0:4]
        #print max(self.Epsy_pl)





        # Set time flow column calculated from data aquisition rate
        # Note that the below is only avalid after calling __inti_columns def.
        self.time_flow = self.__time_xcol__(rate=self.acq_rate, 
                                         size=(len(self.xload1)) ) 

        # ------------------
        # Get the flow time when value is maximum
        # This can be a clue for detecting the break point

        self.break_time, self.break_indx = self.__mxtimepoint__()
        print 
        print '** Max time: ',self.break_time,'** Max indx: ',self.break_indx

    def __mxtimepoint__(self):
        times = []
        times.append(self.__mxtime__(self.epsx1))
        times.append(self.__mxtime__(self.epsx2))
        times.append(self.__mxtime__(self.epsy1))
        times.append(self.__mxtime__(self.epsy2))
        times.append(self.__mxtime__(self.sigx1))
        times.append(self.__mxtime__(self.sigx2))
        times.append(self.__mxtime__(self.sigy1))
        times.append(self.__mxtime__(self.sigy2))
        indices = []
        indices.append(self.__mxpoint__(self.epsx1))
        indices.append(self.__mxpoint__(self.epsx2))
        indices.append(self.__mxpoint__(self.epsy1))
        indices.append(self.__mxpoint__(self.epsy2))
        indices.append(self.__mxpoint__(self.sigx1))
        indices.append(self.__mxpoint__(self.sigx2))
        indices.append(self.__mxpoint__(self.sigy1))
        indices.append(self.__mxpoint__(self.sigy2))

        return min(times), min(indices)

    def __mxpoint__(self,col1):
        # Get absolute value to find the max value
        val=[]
        for i in range(len(col1)):
           val.append(abs(col1[i]))
        mx = max(val)

        iflag = 0
        for i in range(len(val)):
            if mx==val[i]:
                iflag = i
            else: pass
        return iflag
    
    def __mxtime__(self,col1):
        indx = self.__mxpoint__(col1)
        return self.time_flow[indx]

    def __return_columns(self,lines,h_line):
        temp=[]
        for i in range(8):
            temp.append(self.__get_col__(lines,h_line,i))
        return temp

    def plot_xy_range(self,x,y,fig_id=1,ind_plot=111,label='', 
                      xrange=[0,0], yrange=[0,0], loc='upper right', 
                      xlabel='', ylabel='', title='',ls='-', breakdetect=False):
        """
        Basically, it plots x-y curve using matplotlib.pyplot package
        
        ____________________'
        plot_xy_range
        --arguments'
           x,y            : x (horizontal), y (vertical) axis data'
           ind_plot       : the index of matplotlib.pyplot model'
           fig_id         : id # of the figure'
           xrange,yrange  : list of two elements range of x and y'
           loc            : location of label (e.g. 'upper right') "
           xlabel, ylabel : labels of two axes (e.g. 'Time flow')  "
           label          : Legend of the curve (e.g. '316L')      "
           title          : title of the axes                      "
        """
        prop = matplotlib.font_manager.FontProperties(size=5) #font size!!!          
        if breakdetect==True:
            tempx, tempy = [], []
            for i in range(self.break_indx):
                 tempx.append(x[i])
                 tempy.append(y[i])
            x=tempx
            y=tempy
        else: pass
        if self.__is_same__(xrange[0],xrange[1]):
             if self.__is_same__(yrange[0],yrange[1]):
                 fig=plt.figure(fig_id)
                 ax=fig.add_subplot(ind_plot)
                 ax.plot(x,y,label=label,ls=ls)                     
                 if len(xlabel)!= 0: ax.set_xlabel(xlabel)
                 if len(ylabel)!= 0: ax.set_ylabel(ylabel)
                 if len(label)>0: plt.legend(loc=loc, prop=prop)
             else: 
                 for i in range(len(y)):
                     if y[i]>yrange[0]:
                         down = i
                     if y[i]<yrange[1]:
                         up   = i
                 fig=plt.figure(fig_id)
                 ax=fig.add_subplot(ind_plot)
                 ax.plot(x[down:up],y[down:up],label=label,ls=ls)
                 if len(label)>0: plt.legend(loc=loc, prop=prop)
                 if len(xlabel)!=0: ax.set_xlabel(xlabel)
                 if len(ylabel)!=0: ax.set_ylabel(ylabel)
                 # y range case
        elif self.__is_same__(yrange[0],yrange[1])==True:
             for i in range(len(x)):
                 if x[i] <= xrange[0]:
                    down= i
                 if x[i] < xrange[1]:
                    up  = i
             print 'down and up= ', down, up
             fig=plt.figure(fig_id)
             ax=fig.add_subplot(ind_plot)
             ax.plot(x[down:up], y[down:up], label=label,ls=ls)
             if len(label)>0: plt.legend(loc=loc, prop=prop)
             if len(xlabel)!=0: ax.set_xlabel(xlabel)
             if len(ylabel)!=0: ax.set_ylabel(ylabel)
             pass #x range case
        else : 
             print 'Wrong range setting, will return -1'
             raw_input('Press Enter...')
             return -1
        if len(title)>0 :
             ax.set_title(label=title)

    def __is_same__(self,x,y):
        if x==y: return True
        else   : return False

    def plot_xy(self, x, y, fig_id=1, ind_plot=111, label='', cutx_off=0, 
                cuty_off=0, loc='upper right' ):
        """
        Plots simple x vs y plot

        Arguments
        x, y : data
        fig_id : id for the figure
        ind_plot : indicies for the plot, can be 111, 221, 332, etc. 
                   Refer to the matplotlib.pyplot manual
        cutx_off : cut the x below this number (default = -1)
        cuty_off : cut the y below this number (default = -1)
        """
        prop = matplotlib.font_manager.FontProperties(size=10) #font size!!!
        fig1 = plt.figure(fig_id)
        ax = fig1.add_subplot(ind_plot)
        x_temp = []
        y_temp = []
        cut_ind = 0

        if x !=0  :
            if y !=0: pass
               #print
               #print ' ***'
               #print 'Problem x & y cannot be cut-off simutaneously.\n',
               #'X is prefered in this case'
               #print ' ***'
               #print
            for i in range(len(x)) :
               if x[i] < cutx_off : cut_ind=i
               else: pass

        elif y!=0 :
             if x!=0: pass
             if x==0:
                 for i in range(len(y)):
                     if y[i] < cuty_off: cut_ind=i
        else: pass

        ax.plot(x[cut_ind:-1],y[cut_ind:-1],label=label)

        if len(label) >0 :
           plt.legend(loc=loc, prop=prop)
        else: pass

    def __time_xcol__(self, rate= 50, size=3 ):
         temp=[]
         for i in range(size):
            temp.append(rate*i*10**(-3))
         return temp

    def __stress_engi__(self,load,cross_section):
         """
         Changes the given load to engineering (nominal) stress
         """
         value = []
         for i in range(len(load)):
             value.append ( load[i] /cross_section * (10.**3))
         return value

    def __stress_true__(self,stress,strain):
         """
         Returns true stress
         """
         value = []
         for i in range(len(stress)):
             value.append(stress[i]*(1+strain[i]))
         return value

    def plot_str_str(self,stress,strain, label='',ind_plot=221,fig_id=1):
         prop = matplotlib.font_manager.FontProperties(size=5) #font size!!!
         fig1=plt.figure(fig_id)
         ax = fig1.add_subplot(ind_plot)
         ax.plot(strain,stress, label=label)
         plt.legend(prop=prop, loc='lower right')

    def plot_R(self, w_str, l_str,label='',ind_plot=222, fig_id=1):
         prop = matplotlib.font_manager.FontProperties(size=5) #font size!!!        
         fig1=plt.figure(fig_id)
         ax = fig1.add_subplot(ind_plot)

         t_str=[]
         r    =[]
         for i in range(len(w_str)):
            t_str.append( -w_str[i]-l_str[i] )
            try:
             r.append(  w_str[i]/t_str[i] )
            except ZeroDivisionError:
             r.append(0)

         ax.plot(l_str[1000:],r[1000:],label=label)
         plt.legend(prop=prop, loc='lower right')

    def plot_strain_rate(self, strain):
        """
        This is not completed, not even started yet.
        """
        self.acq_rate
        pass

    def __strain_true__(self, str):
        """
        Changes the given engineering strain in the unit of ust 
        into true strain in the unist of st
        """
        for i in range(len(str)):
            value = str[i]*(10**-6)+1.
            if value <= 0:
                 print "Error! strain given is below the value of -1. "
                 str[i] = 0.
            else:
                 str[i]=math.log(value)
        return str

    def __get_param_header__(self, source, label = 'Sampling Rate(ms)'):
        for i in range(len(source)):
            if source[i].split(',')[0] == label:
                return source[i].split(',')[1]
            else: pass
        pass

    def __get_data_avg__(self, list1, list2):
        """
        Returns an averaged columns from two lists
        Argument: list1, list2
        """
        temp = []
        for i in range(len(list1) ):
           temp.append( (list1[i]+list2[i]) / 2.)
        return temp

    def __get_col__(self,lines,h_line,indx):
         """
         Gets column values with given index of column
         Returns a column as in a list varaible.
          
         Arguments are as below
         lines  : Line format of source data
         h_line : The number of header lines
         """
         x=[]
         for i in range(len(lines[(h_line+2):])):
              current_line = self.lines[i+h_line+2].split(',')
              if len(current_line) < 3 :
                      pass
              else:
                   x.append(float(current_line[indx]))

         del current_line
         return x

    def write_output(self,path,filename='biaxial_output.out', 
                     label=['xload1','xload2','yload1','yload2',
                            'xstr1','xstr2','ystr1','ystr2'], 
                     nind=[0,1,2,3,4,5,6,7,8,-1]):
         """
         Write the post-processed data

         Arguments as below:

         Optional arguments:
         label = ['xload1','xload2','yload1','yload2','xstr1','xstr2','ystr1','ystr2']
         nind  = [       0,       1,       2,       3,      4,      5,      6,      7]
         filename = 'biaxial_output.out'

         Requisite argument(s):
         path  =  Requisite

         Note!
         if nind[-1] < 0 (any negative value):
             ----> not disregard nind but label
         """

         if path[-1] =='\\':
               f=file(path+filename,'w')
         else: f=file(path+'\\'+filename, 'w')

         if nind[-1]< 0:
            # label given by text
            nind=[]
            for i in range(len(label)):
                nind.append( self.__getlabel_return_inNumber__(label[i]) )
         else:
            pass
            # label given by nind

         f.write('** post-processed biaxial test output \n')

         #label writing
         f.write('units \n')
         f.write('time: second \n load: kN \n strain: ust(10^-6 mm/mm)\n')
         f.write( 'time'.rjust(8))
         for i in nind:
             f.write( self.__getnumber_return_inLabel__(i).rjust(12) )
         f.write('\n')

         temp = []
         for i in nind:
              temp.append(self.__get_col__(self.lines, self.h_line, i ))

         for i in range(len(temp[0]) ):
            #for k in nind :
            f.write( '%8.4f'%(self.acq_rate*10**(-3)*i))
            for k in range(len(nind)):
                f.write( '%12.4e'%(temp[k][i]))
            f.write('\n')

         del temp
         f.close()

    def __cross_sectionX__(self, header, index='x'):
         """
         Returns the cross section 
         """
         for i in range(len(header)):
             if header[i].split(',')[0] == 'CrossSectionX':
                 value_x = float( header[i].split(',')[1] )
             if header[i].split(',')[0] == 'CrossSectionY':
                 value_y = float( header[i].split(',')[1] )

         if index == 'x': return value_x
         if index == 'y': return value_y

    def __cross_sectionY__(self, header):
         """
         Returns the cross section Y
         """
         return self.__cross_sectionX__(header, index='y')

    def __get_header__(self,lines,indicator='[Data]'):
         """
         Returns header written until the indicator is shown up.
         """
         i=0
         while True:
           if len(lines[i])==0:
             pass
           elif lines[i].split()[0]==indicator:
             #print indicator, 'block is reached'
             break
           i=i+1

         header = []
         for k in range(i):
            header.append(lines[k])

         return header,i

    def __getlabel_return_inNumber__(self,label):
         if label =='xload1':
            return 0
         if label =='xload2':
            return 1
         if label =='yload1':
            return 2
         if label =='yload2':
            return 3
         if label =='xstr1' :
            return 4
         if label =='xstr2' :
            return 5
         if label =='ystr1' :
            return 6
         if label =='ystr2' :
            return 7

    def __getnumber_return_inLabel__(self, i):
         if i ==0:
            return 'xload1'
         if i ==1:
            return 'xload2'
         if i ==2:
            return 'yload1'
         if i ==3:
            return 'yload2'
         if i ==4:
            return 'xstr1'
         if i ==5:
            return 'xstr2'
         if i ==6:
            return 'ystr1'
         if i ==7:
            return 'ystr2'
         else: print 'wrong argument range or type in __getnumber_inLabel__'





#### prints __doc__ for biaxial and uniaxial
print '######################'
print 'class biaxial'
print '######################'
print biaxial.__doc__
print '\n\n'
print '######################'
print 'class uniaxial'
print '######################'
print uniaxial.__doc__
