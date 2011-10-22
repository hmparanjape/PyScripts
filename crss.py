"""
Created on Tue Oct 12 2010

CRSS analysis:
  Analyzes criticial resolved shear stresses
  of grains resulting from VPSC run.
  CRSS is a representative critical value of
  a slip system block. One a set of slip systems
  is in a block, just one criticial resolved
  shear stress is evoloving.

Some statistical analysis is to be done
"""

## libraries
import glob
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os

class crss:
    def __init__(self, path=None):
        """
        --------
        Argument
        --------
        path=None   ; if path==None: The current path is input.
        """
        if path == None:
            tglob = 'crss*.out'
        else: tglob = os.getcwd()+'/'+path +'/crss*.out'
        filenames = glob.glob(tglob)

        ## sorting the list of filenames!
        filenames.sort()

        filecrss = []
        self.master_source = []
        self.nstep = len(filenames)
        for ifile in range(len(filenames)):
            #Appending file object to file crss
            filecrss.append(file(filenames[ifile], 'r'))
            self.master_source.append([]) #master_source[ifile]
            igr = 0
            iline = 0
            while True:
                if iline==0:
                    #print 'iline = 0'
                    #print 'iline --? '
                    #print iline
                    #raw_input()
                    tmp_line=filecrss[ifile].readline()
                    tmp_line=filecrss[ifile].readline()

                if iline==2:
                    self.nmode=len(tmp_line.split())

                tmp_line = filecrss[ifile].readline()   #line input
                 
                if len(tmp_line)<10:                     #check e.o.f
                    #print 'end of line'
                    #print tmp_line
                    #raw_input()
                    #print 'tmp_line'
                    #print tmp_line
                    #raw_input()
                    break

                self.master_source[ifile].append([])  #self.master_source[ifile][igr]
                crss_values = map(float,tmp_line.split())
                self.master_source[ifile][igr] = crss_values
                #print crss_values
                #print 'igr=',igr
                #print tmp_line
                igr = igr + 1
                iline = iline + 1
                
            filecrss[ifile].close()

        """
        for i in range(self.nstep):
            filecrss.append(open(filenames[i], 'r'))
            self.master_source.append(filecrss[i].readlines())
            if i == 0:
                self.ngr = len(self.master_source[0]) - 2 # Excludes the header lines
                self.nmode = len(self.master_source[0][2].split())
        """
    def data(self, igrain=0, imode=0, istep=0):
        """
        Returns the crss under the given conditions
        ---------
        Arguments
        ---------        
        igrain=0
        imode =0
        istep =0
        """
        data = self.master_source
        """
        temp = self.master_source[istep][igrain]
        temp = temp.split()
        return float(temp[imode])
        """
        return self.master_source[istep][igrain][imode]

    def step(self, igrain=0, istep=0):
        """
        with fixed step (deformation)
        Returns the crss of slip systems of a grain

        ---------
        Arguments
        ---------        
        igrain = 0
        istep = 0
        """
        data = self.master_source
        r2 = []
        for i in range(len(data[istep][igrain])):
            r2.append(data[istep][igrain][i])
        return r2
            
    def grain(self, igrain=1):
        """
        Returns crss of the slip systems for each deformation step
        --------
        Argument
        --------
        igrain  = 1

        usage:
           crss = grain(igrain=i)
           crss[j][k]
           j:mode, k:nstep
        """
        data=[]
        for j in range(self.nmode):        
            data.append([])
            for i in range(self.nstep):
                data[j].append(self.data(igrain=igrain, imode=j, istep=i))
        return data

    def plot(self, igrain, nmode, marker=True):
        """
        Plots the crss evolution of the given grain for given slip system

        Arguments:
            igrain
            nmode
            marker = True
        """
        data = self.grain(igrain=igrain)
        data = data[nmode]
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        x = []
        for i in range(len(data)):  # Loop over number of steps
            x.append(i)
        if marker==True: ax.plot(x, data, 'o',alpha=0.4)
        else : ax.plot(x,data)
        
    def show(self):
        try: plt.show()
        except AttributeError:
            print 'Please do use plot method before command show method'
        else: pass


def info_crss(filename='crss_00001.out'):
    """
    Returns the information of the given crss file

    ngr, nstep, nmode  (number of grain, step and slip system)
    """
    a=crss()
    ngr = len(a.master_source[0])
    nstep = a.nstep
    nmode = a.nmode
    return ngr, nstep, nmode
    
def print_crss(filename='pp.out', increment=None):
    """
    Print crss data into a sinlge file of which file name is given

    ---------
    Arugments
    ---------
    filename  = 'pp.out'
    increment = 0.005
    """
    if increment==None: raise IOError
    a=crss()
    fout=file(filename,'w')
    ngr = len(a.master_source[0])
    print '*******************************'
    print 'number of grains =', ngr
    print 'number of steps =', a.nstep
    print 'number of slip system =', a.nmode
    print 'output file name =', filename
    print '*******************************'
    #raw_input()
    for i in range(a.nstep):
        fout.writelines(' %12.8e  '%(increment*(i+1)))
    fout.writelines('\n')
    for k in range(a.nmode):           #nmode
        for i in range(ngr):           #ngrain
            for j in range(a.nstep):   #nstep
                fout.writelines(' %12.8e  '%(
                        a.data(igrain=i, imode=k, istep=j)))
                pass
            fout.writelines('\n')
            pass
        pass
    fout.close()
