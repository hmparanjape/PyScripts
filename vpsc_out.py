# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 10:01:03 2010 /n

@author: youngung /n

This file gives back default plots of the three basic output files of a VPSC run.
The files are Lankford, str_str, and PCYS.

"""

ncol=6
inc_col=2
nrow=6

import matplotlib.pyplot as plt
plt.ion()
fig= plt.figure(figsize=(12,6), frameon=False)
fig.suptitle('VPSC output plotting', fontsize=14, fontweight='bold')
plt.draw()

def is_number(st):
    try:
        float(st)
        return True
    except ValueError:
        return False

"""
LANKFORD PLOTTING
"""
icol=1
fname='c:\Python26\myWorkplace\VPSC\LANKFORD.OUT'
f=open(fname)
s=f.readline()
x=[]
y=[]
plt.ion()
while True:
    s=f.readline()
    if len(s)<5:
        print 'an EOF is reached'
        break
    else:
        x.append(float(s.split()[0]))
        y.append(float(s.split()[2]))
    #print 'x=',x
    #print 'y=',y

f.close()
ax_lnkf=fig.add_subplot(nrow,ncol,icol)
#ax_lnkf=fig.add_axes([0.1,0.1,1/4.,1/4.])
ax_lnkf.set_title('Lankford probing')
ax_lnkf.set_xlabel('R-value')
ax_lnkf.set_ylabel('Angle from RD')
plt.axis([0,90.,0,3.])
plt.xticks([0.,30.,60.,90.])
ax_lnkf.set_aspect('auto','datalim')
line, = ax_lnkf.plot(x,y,'ro')
#del ax_lnkf.lines[0]
plt.draw()


"""
STR_STR PLOTTING
"""
icol=icol+inc_col
fname='c:\Python26\myWorkplace\VPSC\STR_STR.OUT'
f=open(fname)
s=f.readline()
x=[]
y=[]
while True:
    s=f.readline()
    if len(s)==0:
        print 'an EOF is reached'
        break
    else:
        x.append(s.split()[2] )   #-E11
        y.append(s.split()[12])   #-S11
    #print 'x=',x
    #print 'y=',y


f.close()
ax_str=fig.add_subplot(nrow,ncol,icol)
ax_str.set_title('Stress Strain Curve')
ax_str.set_ylabel(r'$\sigma_1$'  +'$_1$')
ax_str.set_xlabel(r'$\epsilon_1$'+'$_1$')
ax_str.set_aspect('auto')
ax_str.plot(x,y,'bo')

"""
PCYS plotting
"""
icol=icol+inc_col
fname='c:\Python26\myWorkplace\VPSC\pcys_pl.OUT'
f=open(fname)
s=f.readline()
#plt.subplot(333)
ax_pcys=fig.add_subplot(nrow,ncol,icol)
x=[[]]
y=[[]]
iset = 0
while True:
    s=f.readline()
    if len(s)==0:
        print 'an EOF is reached'
        break
    else:
        if is_number(s.split()[0]) :
            x1=float(s.split()[0])
            y1=float(s.split()[1])
            x[iset].append(x1)   #-E11
            y[iset].append(y1)   #-S11
        else:
            x.append([])
            y.append([])
            iset=iset+1
            pass


p=[]
for i in range(iset+1):
    p.append(ax_pcys.plot(x[i],y[i],label='test'+str(iset)))

plt.legend(('undeformed','deformed'),bbox_to_anchor=(1.05,1),loc=2,fancybox=True)
ax_pcys.set_aspect('equal')
ax_pcys.set_title('PCYS')
ax_pcys.set_xlabel('S11')
ax_pcys.set_ylabel('S22')
#plt.tick_params(direction='both',length=6,width=2,color='blue')
#plt.xticks([0,30,60,90])
f.close()
#plt.plot(x,y)



