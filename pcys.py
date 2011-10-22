"""
A VPSC code wrapper for equivalent plastic work contour.
This code wraps the VPSC code, written in Fortran, to run
VPSC in the in-plane strain space.
Resulting values are returned, and as an option are plotted.


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
try: import scipy.integrate as integrate
except:  print'You need scipy installed'; raise IOError

if plt.isinteractive()==False:
    plt.ion()
try:
    fcondi = file('vpsc_pcpc.in')  #WRAPPER INPUT FILE
except IOError:
    print "Could not find 'vpsc_pcpc.in' file"
    fcondi = file(path+'/vpsc_pcpc.in','w')
    fcondi.writelines('** VPSC polycrystal plastic work contour wrapper input file\n')
    fcondi.writelines('* HEADER \n')
    ntheta = int(raw_input(' number of theta? >>'))
    fcondi.writelines('* ntheta \n %i\n'%(ntheta))
    print "strain increment, and # of increment (delimiter: ',') "
    inc = map(float,raw_input( ' >> ').split(','))
    fcondi.writelines('* Plastic strain increment')
    fcondi.writelines(", # of increment\n %f, %i\n "%(inc[0],inc[1]))
    print "Type the work levels that you want "
    print "with ',' as the delimiter"
    plw = map(float,raw_input(' >> ').split(','))
    fcondi.writelines("* work contour \n")
    for i in range(len(plw)):
        fcondi.writelines('%f, '%(plw[i]))
    fcondi.close()
    fcondi = file('vpsc_pcpc.in')


lines = fcondi.readlines(); fcondi.close()
print lines ; raw_input()

fvpscin = file('vpsc7.in')      #VPSC INPUT FILE
#os.system('vpsc7.exe')       VPSC EXECUTABLE
#hist = file(path+'\\history', 'w')         #VPSC DEFORMATION HISTORY FILE

if os.name=='nt':      del_file = glob.glob('*.out')
elif os.name=='posix': del_file=glob.glob('*.out') + glob.glob('*.OUT')

for i in del_file:
    os.remove(i)


sratio = []
i = 3
while True:
    r1 = lines[i].split(',')[0]
    r2 = lines[i].split(',')[1]
    try:
        r1 = float(r1)
        r2 = float(r2)
    except ValueError:
        #print ' An EOF is reached'
        print 
        print '** Strain ratios are ... '
        for i in range(len(sratio)):
            print sratio[i][0], sratio[i][1]
        print 
        break
    else:
        sratio.append([])
        sratio[i-3].append(r1)
        sratio[i-3].append(r2)
        i = i + 1

eqincr = lines[i+5].split(',')[0]      #Plastic strain increment (eg. 0.01)
print 'eqincr=',eqincr
nsteps = lines[i+5].split(',')[1]          #Number of plastic strain increments
                                           #in string type as-received
print 'nsteps=',nsteps

def vpscin(fin):
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
    
def fdeform(ictrl, eqincr, sratio):#, sratio, nsteps, ictrl, eqincr):
    """
    Make history file in accordance to strain ratio.
    Arguments
        sratio: strain ratio sequence, e.g. [-1,3]
        nsteps
        ictrl
        eqincr
          The above three are the default variable or flag in
          original history file of VPSC code.
    """
    hist = file(os.getcwd()+'\\history', 'w')         #VPSC DEFORMATION HISTORY FILE
    hist.writelines(str(nsteps) + '   ' + str(ictrl) + '  ' + eqincr + '  298')
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
    

def vpscrun():
    os.system('vpsc7.exe')



#Brings the vpsc.in file in then,
vpscin(fin=fvpscin)

"""
for i in range(len(sratio)):
    #Make Deformation history file
    #inctrl, ratio = historycondition(sratio = sratio[0],

    fdeform(ictrl=0, eqincr=eqincr, sratio=sratio[i])

    #Run VPSC!!!!
    vpscrun()

    #Rename file store
    print '*** Rename str_str.out file ***'
    a = os.system('ren str_str.out str_str'+str(i)+'.out')
    #b = os.system('cls')
    print a #, b
"""

ntheta = 18
for i in range(ntheta):
    theta = -45 + 180./ntheta * i
    sratio = [math.cos(theta*math.pi/180.), math.sin(theta*math.pi/180.)]
    fdeform(ictrl=0, eqincr=eqincr, sratio = sratio)
    #Run VPSC!!
    vpscrun()
    a = os.system('ren str_str.out str_str'+str(i)+'.out')
    print a


f = glob.glob('str_str*')
print 'files..', f
files = []
for i in range(len(f)):
    files.append(file(os.getcwd()+'\\'+f[i]))

#----------------
#plotting the results using matplotlib.pyplot
fig = plt.figure(1)
ax = fig.add_subplot(111)


evm, svm, e11, e22, e33, s11, s22 = [],[],[],[],[],[],[]
for i in range(len(files)):
    source = files[i].read()
    lines = source.split('\n')
    EVM, SVM = [], []
    E11, E22, E33 = [], [], []
    S11, S22 = [], []
    #print lines
    count = 0
    for j in range(len(lines)):
        if len(lines[j]) > 3:
         try:
            float(lines[j].split()[0])
         except ValueError:
            pass
         else:
             #print 'j=',j
             EVM.append(float(lines[j].split()[0]))
             SVM.append(float(lines[j].split()[1]))
             E11.append(float(lines[j].split()[2]))
             E22.append(float(lines[j].split()[3]))
             E33.append(float(lines[j].split()[4]))  #just in case
             S11.append(float(lines[j].split()[8]))
             S22.append(float(lines[j].split()[9]))
             count = count + 1
    evm.append(EVM)
    svm.append(SVM)
    e11.append(E11)
    e22.append(E22)
    e33.append(E33)
    s11.append(S11)
    s22.append(S22)
    workx = integrate.cumtrapz(y=S11, x=E11)
    worky = integrate.cumtrapz(y=S22, x=E22)
    workTotal = []
    for k in range(len(workx)):
        workTotal.append(workx[k]+worky[k])
    fout = file('pcpc_'+str(i)+'.out','w')
    fout.writelines('** Plastic work **')
    fout.writelines('\n workx       worky      workTotal\n')
    for k in range(len(workx)):
        fout.write('%9.3e  %9.3e  %9.3e \n' % (workx[k], worky[k], workTotal[k]))
    fout.close()

            


files = glob.glob('pcpc_*')
print files






frslt = file('YS.rst','w')
frslt.writelines('** Results files **')
frslt.writelines('\n id  sigx        sigy        epsx        epsy ')

os.system('cls')
w = [0.005] #[0.002, 0.003, 0.004, 0.006] #, 0.008, 0.01]   #0.002, 0.004, 0.006
color = ['k','gray','b','r','yellow','green']
for iw in range(len(w)):
    for i in range(len(files)):
        f = file(files[i])
        source = f.read()
        f.close()
        lines = source.split('\n')
        temp = []
        kount = 0
        #for j in range(len(lines)):
        j = 0
        while True:
            try:
                wtot = float(lines[j+2].split()[2])
                
            except IndexError:
                break
            if wtot > w[iw]:
                wtot_1 = float(lines[j+1].split()[2])
                print 'wtot, wtot_1 = ', wtot, wtot_1
                ind1 = j - 1
                ind2 = j
                print 'ind1, ind2=', ind1, ind2
                break
            j = j + 1

        
        #Interpolate the sigx, sigy, epsx, and epsy (Simpler, easier)
        print s11[i][ind2], s11[i][ind1], w[iw]
        print s22[i][ind2], s22[i][ind1], w[iw]
        sigx = (s11[i][ind2] - s11[i][ind1]) / (wtot - wtot_1) * (w[iw] - wtot_1) + s11[i][ind1]
        sigy = (s22[i][ind2] - s22[i][ind1]) / (wtot - wtot_1) * (w[iw] - wtot_1) + s22[i][ind1]
        epsx = (e11[i][ind2] - e11[i][ind1]) / (wtot - wtot_1) * (w[iw] - wtot_1) + e11[i][ind1]
        epsy = (e22[i][ind2] - e22[i][ind1]) / (wtot - wtot_1) * (w[iw] - wtot_1) + e22[i][ind1]
        #sigx = (s11[i][ind1] + s11[i][ind2]) / 2.
        #sigy = (s22[i][ind1] + s22[i][ind2]) / 2.
        #epsx = (e11[i][ind1] + e11[i][ind2]) / 2.
        #epsy = (e22[i][ind1] + e22[i][ind2]) / 2.
        print 'sigx, sigy = ', sigx, sigy
        print 'epsx, epsy = ', epsx, epsy
        print
        frslt.write('\n %i %11.3e %11.3e %11.3e %11.3e' %(i, sigx, sigy, epsx, epsy))
        ax.plot([sigx],[sigy],marker='o', color=color[iw])
        plt.xlim([0,20])
        plt.ylim([0,20])
        """
        fig = plt.figure(i)
        ax = fig.add_subplot(111)
        ax.plot(sigx,sigy)
        """
    frslt.write('\n')

frslt.close()
#plt.show()

