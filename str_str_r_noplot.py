import os
os.sys.path.append('c:\\python26\\myworkplace\\myscripts')
import math
#import matplotlib.pyplot as plt
import glob

print
print '*****************************************'
print 'Simple *.TRA plotting and post processing'
print 'But probably, maybe, possibly, presumably'
print 'can be applied to any of format of data  '
print '*****************************************'
print

"""
path = raw_input('Path or Eneter(current working directory)  >>')
if len(path) == 0: path = os.getcwd()
else :
    try:
        os.chdir(path)
    except:
        print 'Exception raised possibly due to wrong path input'
"""

files = glob.glob('*.tra')
if len(files) == 0:
    print
    print '*********************************'
    print 'No files in the current directory'
    print 'Current directory is ', os.getcwd()
    print '*********************************'
    print
for i in range(len(files)):
    print files[i], '\n'
"""
iext = 0
iforce = 1
iwid = 2
"""
def param():
    t = raw_input('Thickness of the sample, (default = 0.42) >>')
    L0 = raw_input('Gauge length in which the extension is obtained (default = 50) >>')
    W0 = raw_input('Transverse width of the sample (default = 12.5) >>')
    if len(t) == 0: t = 0.42
    else: t = float(t)
    if len(L0) == 0: L0 = 50.
    else: L0 = float(L0)
    if len(W0) == 0: W0 = 12.5
    else: W0 = float(W0)
    area = t * W0
    print 't, L0, W0, area = ', t, L0, W0, area
    return t, L0, W0, area
    
def column():
    iext = raw_input('Column # of extension (default = 0) >>')
    ifor = raw_input('Column # of standard force (default = 1) >>')
    iwid = raw_input('Column # of change in width (default = 2) >>')
    if len(iext) == 0: iext = 0
    else: iext = int(iext)
    if len(ifor) == 0: ifor = 1
    else: ifor = int(iext)
    if len(iwid) == 0: iwid = 2
    else: iwid = int(iwid)
    return iext, ifor, iwid
    
def delimiter():
    delimt = raw_input("Delimiter separating the data columns(default=';')")
    if len(delimt) == 0:
        delimt =';'
    print 'The given delimiter is ', delimt
    return delimt
"""
t = 0.42 # unit: mm
L0 = 50.  # unit: mm
W0 = 12.5 # unit: mm
area = t * L0 # unit: mm^2
"""
if len(files) > 0:
    print
    print '*******************'
    print 'Intrinsic functions'
    print '*******************'
    print
    print '*******************'
    t, L0, W0, area = param()
    iext, ifor, iwid = column()
    delimt = delimiter()

for i in range(len(files)):
    print 'filename = ', files[i]
    #print path+'\\'+files[i]
    f = open(path+'\\'+files[i], 'r')
    source = f.read()
    f.close()
    fout = open(path+'\\'+files[i].split('.')[0]+'.str', 'w') # write file
    fout.writelines('   Strain            Stress                R-value\n')
    fout.writelines('                    (MPa) \n')
    
    lines = source.split('\n')
    header = lines[0]
    E_l = []
    S = []
    E_w = []
    E_t = []
    R = []
    for j in range(len(lines)-1):
        try:
            cline = map(float,lines[j+1].split(delimt))
        except ValueError:
            #print lines[j+1].split(',')
            print 'EOF at the', j,'th line.'
            print
        else:
            ext = cline[iext]
            force = cline[ifor]
            chg_wid = cline[iwid]
            el = math.log(ext / L0 + 1)
            E_l.append(el)
            ew = math.log(chg_wid / W0 + 1)
            ew = -ew
            E_w.append(ew)
            stress = force * (ext / L0 + 1) / 10.
            S.append(stress)
            E_t.append(-(E_l[j]+E_w[j]))
            R_ = E_w[j]/E_t[j]
            R.append(R_)
            fout.write('%13.3e  %13.3e  %17.3e  \n'%(el, stress,  R_ ))

        #E_l.append(math.log(lines[j+1].split(',')[0]/L0 + 1))
        """
        S.append(lines[j+1].split(',')[1]/area *(1 + lines[j+1].split(',')[0]/L0))
        E_w.append(math.log(liens[j+1].split(',')[2]/W0+1))
        E_w[k] = -E_w[k]
        E_t.append(-(E_l[k]+E_w[k]))
        R.append(E_w[k]/E_t[k])
        """
    """
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(E_l,S)
    plt.xlim(0,0.8)
    
    fig = plt.figure(2)
    ax2 = fig.add_subplot(111)
    ax2.plot(E_l, R)
    ax2.set
    plt.ylim(0,1.5)
    plt.xlim(0,0.8)
    """


    fout.close()
"""
if len(files)==0: pass
else:  plt.show()
"""
