from scipy import *
#from pylab import *
import os
import matplotlib.pyplot as plt
path = os.getcwd()
plt.ion()
fig = plt.figure(2)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
plt.draw()
def voce(parameter=[1.0,0.1,3.0,2.], step=0.02, path=None, sxfile=None):
     p = parameter
     if p[1]==0. :
        p[1]==0.00001
     t = arange(0.,1.,step)
     stress = p[0]+(p[1]+p[3]*t)*(1-exp(-t*abs(p[2]/p[1])))
     ax2.plot(t,stress)
     #plot(t,stress)
     #---------------------------------------
     # VPSC singly crystal file modification
     #     & run VPSC
     #---------------------------------------
     sx_mod(param=parameter, path='\\sx', fname = sxfile)
     #path='c:\python26\myworkplace\\VPSC\\VPSC_COMPILE\\'
     run(path=path, fname='VPSC7.exe')
     print 'run finished'



     f = open(path+'\\str_str.out')

     line= f.readline()

     x,y=[],[]
     while True:
         temp=[]
         s=f.readline()
         if len(s)<10:
             print 'An EOF is reached'
             break
         else:
             x.append(float(s.split()[0]))
             y.append(float(s.split()[1]))
     f.close()
     #plot(x,y)
     ax1.plot(x,y)
     """
     source=f.read()
     lines = source.split('\n')
     x,y=[],[]
     for i in range(len(lines)):
         x.append(lines[i+1].split()[0])
         y.append(lines[i+1].split()[1])
     """

     #print 'x=', x
     #print 'y=', y


     plt.draw()
     f.close()
     return stress

def run(path=None, fname='VPSC7.exe' ):
     os.chdir(path)
     os.system('VPSC7.exe')

def sx_mod(path=None,fname=None,param=[1.0,0.1,1.0,0.]):
     f      = open(path+'\\'+fname)
     source = f.read()
     lines  = source.split('\n')
     p      = map(float,lines[19].split())
     for i in range(4):
        p[i] = param[i]
     f.close()
     f      = open(path+'\\'+fname, 'w')

     for i in range(len(lines)):
         if i == 19:
             for j in range(6):
                f.write(str(p[j])+'   ')
             f.write('\n')
         else:
             f.write(lines[i])
             f.write('\n')

     f.close()

















