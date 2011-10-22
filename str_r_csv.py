import os
os.sys.path.append('c:\\python26\\myworkplace\\myscripts')
import math
import np
import matplotlib.pyplot as plt
import glob

print
print '*****************************************'
print 'Simple CSV plotting and post processing'
print 'But probably, maybe, possibly, presumably'
print 'can be applied to any of format of data  '
print '*****************************************'
print

path = raw_input('Path or Enter(returns os.getcwd() )  >>  ')
extension = raw_input('Extension of files ... (e.g. txt)')
if len(path) == 0: path = os.getcwd()
else:
    try:
        os.chdir(path)
    except:
        print '** Exception raised possibly due to wrong path input **'
files = glob.glob('*.' + extension)
if len(files) == 0:
    print
    print '*********************************'
    print 'No files in the current directory'
    print 'whose extension is ' + extension
    print 'Current directory is ', os.getcwd()
    print '*********************************'
    print
for i in range(len(files)):
    print files[i], '\n'
    
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

def get_info():
    t, L0, W0, area = param()
    iext, ifor, iwid = column()
    delimt = delimter()

for i in range(len(files)):
    print 'filename = ', files[i]
    f = open(path + '\\'


