"""
Estimates the maximum data in terms of stress and strain
until the achievable values to the maximum
"""
import os
import scipy.integrate as integrate
import math
import glob

#os.sys.path.append('c:\\python26\\myworkplace\\myscripts')

import biaxial_pp

modulusx = raw_input('Modulus along RD >> ') 
modulusy = raw_input('Modulus along TD >> ') 
filename = raw_input('Filename >> ')

lsx = raw_input('Lower sigma X bound >>')
lsy = raw_input('Lower sigma Y bound >>')
usx = raw_input('Upper sigma X bound >>')
usy = raw_input('Upper sigma Y bound >>')

if len(modulusx)==0: modulusx=None
else: modulusx = float(modulusx)
if len(modulusy)==0: modulusy=None
else: modulusy = float(modulusy)

if len(lsx)==0: lsx = 20.
else: lsx = float(lsx)
if len(lsy)==0: lsy = 20.
else: lsy = float(lsy)
if len(usx)==0: usx = 50.
else: usx = float(usx)
if len(usy)==0: usy = 50.
else: usy = float(usy)

def sumup(x, y):
    temp = []
    for i in range(len(x)):
        temp.append(x[i]+y[i])
    return temp

#files = glob.glob('*.csv')
files = [filename]
unia = []
for i in files:
    unia.append(biaxial_pp.uniaxial(
        path=os.getcwd(), low_sigx=lsx, up_sigx=usx,
        low_sigy=lsy, up_sigy=usy, filename=i,
        iplot=False, modulusx=modulusx, modulusy=modulusy,
        #isuni=False
        ))

for i in range(len(files)):
    fout = open(os.getcwd() + os.sep + files[i].split('.')[0] + '.bpp', 'w')
    fout.writelines('** post-processed data **\n')
    fout.writelines(' time        Ex         Ex_pl      Sx')
    fout.writelines('         Ey         Ey_pl      Sy')
    fout.writelines('         Srx        Sry')
    fout.writelines('        Work_tot   Work_x     Work_y')
    fout.writelines('     Wtot_pl    Wx_pl      Wy_pl\n')
    j = 0
    delT = unia[i].acq_rate/10**3
    wx = integrate.cumtrapz(y=unia[i].sigx, x=unia[i].epsx)
    wy = integrate.cumtrapz(y=unia[i].sigy, x=unia[i].epsy)
    wtot = sumup(wx, wy)
    wxpl = integrate.cumtrapz(y=unia[i].sigx, x=unia[i].Epsx_pl)
    wypl = integrate.cumtrapz(y=unia[i].sigy, x=unia[i].Epsy_pl)
    wpltot = sumup(wxpl, wypl)
    while True:
        try:
            time = unia[i].time_flow[j]
            expl = unia[i].Epsx_pl[j]
            eypl = unia[i].Epsy_pl[j]
            ex = unia[i].epsx[j]
            ey = unia[i].epsy[j]
            sigmax = unia[i].sigx[j]
            sigmay = unia[i].sigy[j]
            if i==0: srx = 0
            else: srx = (ex - unia[i].epsx[j-1]) / delT
            if i==0: sry = 0
            else: sry = (ey - unia[i].epsy[j-1]) / delT
            fout.write('%10.4e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n'% (time, ex, expl, sigmax, ey, eypl, sigmay, srx, sry, wtot[j], wx[j], wy[j], wpltot[j], wxpl[j], wypl[j]) )
        except IndexError:
            break
        else:
            j = j + 1
    fout.close()
