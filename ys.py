"""
A VPSC code wrapper for equivalent plastic work contour.
This code wraps the VPSC code, written in Fortran, to run
VPSC in the in-plane strain space.
Resulting values are returned, and as an option are plotted.

Project starts from 08-23-2010

Author: Youngung Jeong
        Materials Mechanics Laboratory,
        Graduate Institute of Ferrous Technology,
        Pohang University of Science and Technology.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import scipy.integrate as integrate
import smtplib
if plt.isinteractive()==False:
    plt.ion()
path = os.getcwd()

class ys:
    """
    yield surface class
    Arguments:
    ialert = False    :inform the process via email
    quad = 'all'      :choice on the stress space ('quad' or '1st')
    """
    def __init__(self, inp=path+'\\vpsc_pcpc.in',
                 vpscin=path+'\\vpsc7.in', ialert=False, quad='all'):
        #os.system('del *.out')
        if ialert==True:
            fromaddr = raw_input('From address>> ')
            toaddrs = raw_input('To address>> ')
            msg = raw_input('msg>> ')
            username = raw_input('user-name>> ')
            password = raw_input('password>> ')
        fig = plt.figure(1)
        ax = fig.add_subplot(111, aspect='equal')
        fcondi = file(inp)
        fvpscin = file(vpscin)
        lines = fcondi.read()
        lines = lines.split('\n')
        ntheta = int(lines[3])
        eqincr = float(lines[5].split(',')[0])
        nsteps = int(lines[5].split(',')[1])
        self.nsteps = nsteps
        w = []
        k = 7
        """
        while True:
            if len(lines[k]) < 2:
                break
            else:
                for j in range(len(lines[k].split(','))):
                    w.append(float(lines[k].split(',')[j]))
                k = k + 1
        """
        for j in range(len(lines[7].split(','))):
            w.append(float(lines[7].split(',')[j]))
        print 'Work levels:', w
        #---- VPSC7.in 
        self.vpscin(fin = fvpscin)

        fout = file(path + '\\ys.rst', 'w')
        fout.writelines('sigx sigy epsx epsy\n')
        fwork = file(path +'\\work.rst', 'w')
        fwork.writelines('Wx   Wy   Wtot')
        #---- Main Loop over ntheta
        workfiles = []
        for j in range(len(w)):
            temp_file = open(path+'\\YS_'+str(w[j]).split('.')[0]+str(w[j]).split('.')[1]+'.rst','w')
            workfiles.append(temp_file)
            workfiles[j].writelines('Work contour file for w = ' + str(w[j]) + '\n')
            workfiles[j].writelines(' Sigx     Sigy       Epsx      Epsy\n')
            
        for i in range(ntheta):
            if quad=='all':
                theta = 360. / ntheta * i
            elif quad=='1st':
                theta = -30. + 180./ntheta * i
            ex, ey = self.__ratio__(theta = theta)
            print 'theta =', theta

            ### ---- HISOTRY FILE
            self.__fdeform__(ictrl=0, eqincr=eqincr, sratio = [ex,ey]) 

            ### ---- VPSC executable run
            print '\n\n\n***** VPSC RUN'
            print i+1,'/',ntheta + 1
            self.vpscrun(path=os.getcwd())
            
            fstr_str = file(path + '\\str_str.out')
            S11, S22, E11, E22, work, workx, worky = self.str_str(f=fstr_str)
            
            for m in range(len(work)):
                fwork.write('%.3e %.3e %.3e \n'% (workx[m], worky[m], work[m])  )
            
            print 'max work = ', work[-1]
            print 'min work = ', work[0]
            color = ['k','gray','b','r','yellow','green']
            for j in range(len(w)):
                if work[-1] > w[j] :
                    if work[0] > w[j]:
                        print 'work level is too low'
                        #raw_input('Press Enter...')
                    #else:
                    for k in range(len(work)):
                            #print 'len(work) = ', len(work)
                            if work[k] > w[j]:
                                ind1 = k - 1
                                ind2 = k
                                print 'ind1, ind2=', ind1, ind2
                                break
                            else: pass
                    Sx = self.interpolate(y1=S11[ind1], y2=S11[ind2],
                                              x1=work[ind1], x2=work[ind2],
                                              x=w[j])
                    Sy = self.interpolate(y1=S22[ind1], y2=S22[ind2],
                                              x1=work[ind1], x2=work[ind2],
                                              x=w[j])
                    Ex = self.interpolate(y1=E11[ind1], y2=E11[ind2],
                                              x1=work[ind1], x2=work[ind2],
                                              x=w[j])
                    Ey = self.interpolate(y1=E22[ind1], y2=E22[ind2],
                                              x1=work[ind1], x2=work[ind2],
                                              x=w[j])
                        
                    print 'Sx1, Sx2, Sx = ', S11[ind1], S11[ind2], Sx
                    print 'Sy1, Sy2, Sy = ', S22[ind1], S22[ind2], Sy

                else:
                    print
                    print '*** Alert! ***'
                    print ' the level of work is not achived '
                    print
                    raw_input('Press Enter...')
                    
                fout.write('%8.3e %8.3e %8.3e %8.3e \n' % (Sx, Sy, Ex, Ey))
                workfiles[j].writelines('%8.3e %8.3e %8.3e %8.3e \n'% (Sx, Sy, Ex, Ey))
                ax.plot(Sx,Sy, marker='o', color= color[j], alpha = 0.5)
                plt.draw()
            fwork.write('\n')
            #plt.xlim([0,20])
            #plt.ylim([0,20])
        fout.close()
        for i in range(len(workfiles)):
            workfiles[i].close()
        if ialert == True:
            server = smtplib.SMTP('smtp.gmail.com:587')
            server.starttls()
            server.login(username, password)
            server.sendmail(fromaddr, toaddrs, msg)
            server.quit()









    def interpolate(self, x1, x2, y1, y2, x):
        return (y2-y1)/(x2-x1)*(x-x1)+y1
            
    def str_str(self, f):
        """
        Gets str_str.out file and analyze it
        """
        lines = f.read()
        lines = lines.split('\n')
        EVM, SVM, E11, E22, E33, S11, S22 = [], [], [], [], [], [], []
        nline = len(lines) - 1
        for i in range(nline + 5):
            try:
                EVM.append(float(lines[i+1].split()[0]))
                SVM.append(float(lines[i+1].split()[1]))
                E11.append(float(lines[i+1].split()[2]))
                E22.append(float(lines[i+1].split()[3]))
                E33.append(float(lines[i+1].split()[4]))
                S11.append(float(lines[i+1].split()[8]))
                S22.append(float(lines[i+1].split()[9]))
            except IndexError:
                break
        workx = integrate.cumtrapz(y=S11, x=E11)
        worky = integrate.cumtrapz(y=S22, x=E22)
        workTotal = []
        for k in range(len(workx)):
            workTotal.append(workx[k]+worky[k])
        return S11, S22, E11, E22, workTotal, workx, worky

    def vpscin(self, fin):
        """
        MAKE 'vpsc7.in' to be dedicated only to HISTORY deformation
        """
        source = fin.read()
        source = source.split('\n')
        #modify vpsc7.in accordingly
        fname = fin.name
        fin.close()
        
        #vfile = file(fname.split('.')[0]+'_ex.'+fname.split('.')[1],'w')
        vfile = file(fname, 'w')
        for i in range(len(source[0:31])):
            vfile.writelines(source[i]+'\n')
        vfile.writelines('1\n*IVGVAR AND PATH\NAME OF FILE FOR EACH PROCESS (dummy if ivgvar=2,3)')
        vfile.writelines('\n0\nhistory')
        vfile.close()
            
    def __ratio__(self, theta):
        """
        Returns the epsx and epsy ratio for the give theta
        """
        return math.cos(theta * math.pi / 180.), math.sin(theta * math.pi / 180.)

    def vpscrun(self, path, filename='vpsc7.exe'):
        """
        Run the vpsc!
        """
        os.chdir(path)
        os.system(filename)
        
    def __fdeform__(self, ictrl, eqincr, sratio):
        """
        Writes history file!
        """
        hist = file(os.getcwd()+'\\history', 'w')         #VPSC DEFORMATION HISTORY FILE
        hist.writelines(str(self.nsteps) + '   ' + str(ictrl) + '  ' + str(eqincr) +'  298')
        hist.writelines('\n* boundary conditions')
        hist.writelines('\n    1       0       0           iudot    |    ')
        hist.writelines('flag for vel.grad.')
        hist.writelines('\n    1       1       0                    |    ')
        hist.writelines('(0:unknown-1:known)')
        hist.writelines('\n    1       1       0                    |    ')
        hist.writelines('\n                                         |')
        hist.writelines('\n  '+str(sratio[0]))
        hist.writelines('      0.      0.          udot     |    vel.grad')   
        hist.writelines('\n    0.     '+str(sratio[1]))
        hist.writelines('     0.                   |')
        hist.writelines('\n    0.      0.    '+str(-(sratio[0]+sratio[1]))+'                 |')
        hist.writelines('\n                                         |')
        hist.writelines('\n    0       1       1           ')
        hist.writelines('iscau    |    flag for Cauchy')
        hist.writelines('\n            0       1                    |')
        hist.writelines('\n                    1                    |')
        hist.writelines('\n                                         |')
        hist.writelines('\n    0.      0.      0.          scauchy  |    Cauchy stress')
        hist.writelines('\n            0.      0.                   |')
        hist.writelines('\n                    0.                   @')
        #iudot (flag for vle. gradient 0:unknown 1:known)
        #velocity gradient
        #flag for Cauchy stress
        #Cauchy stress: (0,0,0\n0,0,0\n0,0,0)
        hist.close()        
