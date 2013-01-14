"""
Replacement of old biaxial_pp and its family
developed by youngung.jeong@MML in GIFT, POSTECH

a simpler in plane biaxial tester post-processor!!

input of the file : Due to its trailing tail for control values
attached to the file, np.loadtxt method is not applicable.
"""
bpp_manual="""

## How to complete the pseudo yield surface with the error bar
## Just follow steps below

1. Convert the MTS result file into '*.str' files
    using pp() under tr_str_r.py module.
2. Convert the resulting '*.str' files into '*.bpp' files
    by UniMultiLoader() under biaxial.py module.
3. Convert the in-plane biaxial tester's '*.csv' files
    into '*.bpp' files using multiloader(l=3,h=10)
    under biaxial.py module
4. Once all '*.bpp' files are collected correctly,
    run eqwcplot(wlv = [ ], confidence=0.95, ifig=3)
    under biaxial.py module

    
under the folder where '*.txt' MTS files are located
>>> import str_str_r
>>> str_str_r.pp()


>>> import biaxial
>>> biaxial.UniMultiLoader()
*.bpp files are created.


Move to the folder where '*.csv' in-plane-biaxial test files are located
>>> biaxial.multiloader()


## follow the instructions popping up over the operation.
## move all *.bpp files into the same folder


>>> data = biaxial.eqwcplot(wlv=[0.2, 1.0, .... ], confidence=0.95, ifig=2)
"""
uni_manual="""
Given the bpp files of uniaxial tension that may have been obtained accoring to
the "bpp_manual", mean value and its error bar featured stress-strain& R-strain curves
are saved and plotted.

All that need be done is as below

>>> import biaxial (as you already did)
>>> rst = biaxial.TotUniTensioin(wildcard = '*RD*.bpp', iplot=True,
                                 confidence=0.95, stp=10,
                                 sfout='eps-sig_meancurve',
                                 rfout='eps-rvl_meancurve')

>>> rst[0] --> x_points (String)
>>> rst[1] --> sig_points with error of each point
>>> rst[2] --> rvl_points with error of each point
"""

## list of modules

## intrinsic modules
import os, glob

# Scientific Python modules
import numpy as np
import scipy as sp
import scipy.integrate as integrate

## the plotting has to be optional
## (despite of its superior advantages...)
try:
    import matplotlib.pyplot as plt
except: iplt = False
else: iplt = True


def menual():
    print bpp_manual
    print uni_manual
    pass

def loader(filename=None, delimiter=',', thick=None,
           StressAtLow=5, StressAtHigh=30,bin=100):
    """
    file loader

    Argument:
    filename
    delimiter
    thick
    StressAtLow, StressAtHigh
    
    returns
    avg force, avg ystrains, avg xstrains
    """
    if filename==None:
        print "Give filename "
        return -1
    
    dlines = open(filename, 'r').readlines()  ## data input
    head = dlines[:26]
    
    #--------------------------------------------------------------------#
    ## dimension and control environments
    area, dt, sratio = header(
        head, delimiter) #Extract info from header
    if thick!=None: area = thick * 60 #dimention 60 is fixed
    #--------------------------------------------------------------------#

    #--------------------------------------------------------------------#
    ## data trim
    dlines = dlines[26::] #trim the header
    dd = trim(dlines, delimiter).transpose()
    xl1,xl2,yl1,yl2 = dd[0][::],dd[1][::],dd[2][::],dd[3][::]
    xe1,xe2,ye1,ye2 = dd[4][::],dd[5][::],dd[6][::],dd[7][::]
    time = np.linspace(start=0., stop=dt*(len(xe1)-1), num=len(xe1))
    #print len(time), len(xe1);raw_input(); return time
    #--------------------------------------------------------------------#

    #--------------------------------------------------------------------#
    ## post-process
    # Plastic decomposition, Youngung's modulus and so on.
    # print 'xl1, xl2, yl1, yl2', xl1, xl2, yl1, yl2
    # print 'xe1,xe2,ye1,ye2',xe1,xe2,ye1,ye2
    # raw_input()
    

    
    xs1,xs2,ys1,ys2,xep1,xep2,yep1,yep2,xrep1,xrep2,yrep1,yrep2,xrs1,xrs2,yrs1,yrs2 = pp(
        xl1,xl2,yl1,yl2,xe1,xe2,ye1,ye2,
        area,StressAtLow,StressAtHigh,time,iplt,bin
        )
    #--------------------------------------------------------------------#

    #--------------------------------------------------------------------#
    ## data average
    xs = (xs1+xs2)/2.  ; ys = (ys1+ys2)/2.
    xep = (xep1+xep2)/2.; yep= (yep1+yep2)/2.
    xrep = (xrep1+xrep2)/2.; yrep = (yrep1+yrep2)/2.
    yrs = (yrs1+yrs2)/2.; xrs = (xrs1+xrs2)/2.
    #--------------------------------------------------------------------#

    #--------------------------------------------------------------------#
    ## strain direction!
    ed = np.arctan2(yrep, xrep)
    ## stress direction!
    sd = np.arctan2(yrs, xrs)
    #--------------------------------------------------------------------# 
    
    ## plot the average value if iplt
    if False:
        figure = plt.figure(98)
        ax = figure.add_subplot(111)
        ax.plot(xep,xs,label='X-axis');ax.plot(yep,ys,label='Y-axis')
        ax.legend(loc='best')
        pass

    # print '\n\n______________________________'
    # print '       Data trim option       '
    # print '______________________________\n'

    istop = True
    while True:
        emt = float(raw_input('type in the maximum time(s)  >>  '))
        #99 and 98 are reserved on showing purpose
        try:
            plt.close(99)
            plt.close(98)
        except: pass

        iid=-1
        for i in range(len(time)):
            if time[i]>emt:
                iid = i
                break
            pass
        
        #--------------------------------------------------------------------#
        ## trims the data
        Xep = xep[:iid:]; Yep = yep[:iid:]
        Xs = xs[:iid:]; Ys = ys[:iid:]
        Xrep = xrep[:iid:]; Yrep = yrep[:iid:]
        Ed = ed[:iid:]; Time=time[:iid:]

        #--------------------------------------------------------------------#
        ## Plastic work
        workx = integrate.cumtrapz(x=Xep, y=Xs)
        worky = integrate.cumtrapz(x=Yep, y=Ys)
        workx = np.append(np.array([0]),workx)
        worky = np.append(np.array([0]),worky)
        work = workx + worky
        #--------------------------------------------------------------------#

        fig = plt.figure(99)
        ax1 = fig.add_subplot(231, label='sig-eps')
        ax2 = fig.add_subplot(232, label='work - sig');
        ax3 = fig.add_subplot(233, label='time - eps'); 
        ax4 = fig.add_subplot(234, label='time - eps_d')
        ax5 = fig.add_subplot(235, label='eps-load')

        ax1.plot(Xep, Xs);ax1.plot(Yep,Ys)
        ax1.set_xlabel(r'$\varepsilon_{avg}^{pl}$',dict(fontsize=20))
        ax1.set_ylabel(r'$\sigma_{avg}$',dict(fontsize=20))
        ax1.legend()
        
        ax2.plot(work,Xs);ax2.plot(work,Ys)
        ax2.set_xlabel(r'$work^{pl}$', dict(fontsize=20))
        ax2.set_ylabel(r'$\sigma_{avg}$', dict(fontsize=20))
        ax2.legend()
        
        ax3.plot(Time,Xs);ax3.plot(Time,Ys)
        ax3.set_xlabel('time [s]', dict(fontsize=20))
        ax3.set_ylabel(r'$\varepsilon_{avg}^{pl}$', dict(fontsize=20))
        ax3.legend()
        
        ax4.plot(Time,Ed*180./np.pi)
        ax4.set_xlabel('time',dict(fontsize=20))
        ax4.set_ylabel(r'$\theta_{avg}$', dict(fontsize=20))
        ax4.legend()        

        ax5.plot(xe1, xl1); ax5.plot(ye1,yl1); ax5.plot(xe2, xl2); ax5.plot(ye2,yl2)
        ax5.set_xlabel(r'$\epsilon$',dict(fontsize=20))
        ax5.set_ylabel('Force [KN]',dict(fontsize=20))
        ax5.legend()

        ## temporal save to pdf file ##
        try: fig.savefig('temp.pdf')
        except: print 'you are not authorized to write '
        else: pass
            
            
        #-----------------------------#
        
        iflag = raw_input('Type (y) close the figure or (n) iterate again  (default = y) >>>')
        if iflag=='y' or len(iflag)==0:
            plt.close(99)
            break
        elif iflag=='n': plt.close(99)
        else:
            print "It accepts the stroke as 'n' "
            plt.close(99)
        pass
    
    return Xep, Yep, Xs, Ys, workx, worky, work, Xrep, Yrep, Ed, sratio

def pp(xl1,xl2,yl1,yl2,xe1,xe2,ye1,ye2,
       area,StressAtLow,StressAtHigh,time,iplt,bin):
    """
    What a tedius post-processing!
    
    Arguements:
    xl1,xl2,yl1,yl2,xe1,xe2,ye1,ye2,
    area,StressAtLow,StressAtHigh,time
    """
        ## engineering stress
    xs1 = xl1 / area * 10**3 ##area[mm] , load [KN]  
    xs2 = xl2 / area * 10**3
    ys1 = yl1 / area * 10**3
    ys2 = yl2 / area * 10**3  
    
    ## engineering strain
    xe1 = xe1/10**6
    xe2 = xe2/10**6
    ye1 = ye1/10**6
    ye2 = ye2/10**6    

    ## True stress
    xs1 = xs1 * ( 1 + xe1)
    xs2 = xs2 * ( 1 + xe2)
    ys1 = ys1 * ( 1 + ye1)
    ys2 = ys2 * ( 1 + ye2)    

    ## true strain
    xe1 = np.log(xe1+1)
    xe2 = np.log(xe2+1)
    ye1 = np.log(ye1+1)
    ye2 = np.log(ye2+1)

    ## Elastic moduls calculation
    # StressAtLow, StressAtHigh  = 10, 80 [Mpa]
    # print 'xe1, xs1'; print xe1, xs1; raw_input()
    #print len(xe1) ; raw_input()

    while True:
        try: Ex1 = slope(xe1, xs1, StressAtLow, StressAtHigh)
        except:
            print "Young's modulus along x1 has a problem"
            print " type mannual!"
            Ex1 = float(raw_input('Ex1 in [MPa] level >> '))
            pass
        try: Ex2 = slope(xe2, xs2, StressAtLow, StressAtHigh)
        except:
            print "Young's modulus along x2 has a problem"
            print " type mannual!"
            Ex2 = float(raw_input('Ex2 in [MPa] level >> '))
            pass
        try: Ey1 = slope(ye1, ys1, StressAtLow, StressAtHigh)
        except:
            print "Young's modulus along y1 has a problem"
            print " type mannual!"
            Ey1 = float(raw_input('Ey1 in [MPa] level >> '))
            pass
        try: Ey2 = slope(ye2, ys2, StressAtLow, StressAtHigh)
        except:
            print "Young's modulus along y2 has a problem"
            print " type mannual!"
            Ey2 = float(raw_input('Ey2 in [MPa] level >> '))
            pass        
        
        print "Young's modulus as follow"
        print "Ex1: %6.3f Ex2: %6.3f Ey1: %6.3f Ey2: %6.3f [GPa]"%(
            Ex1/1000., Ex2/1000., Ey1/1000., Ey2/1000.)
        young = np.array([Ex1/1000., Ex2/1000., Ey1/1000., Ey2/1000.])
        if any(abs(young[i]) < 100. for i in range(4)):
            print "Too low slope estimation"
            print "Do you want to proceed or not"
            iflag = raw_input("n:again  y:pass  >>>")
            if iflag=='n':
                print "\ntwo option,\n1. chage stress at low and high"
                print "2. Manuually input the slope"
                if int(raw_input("Type the option (1 or 2 (any else number) >>"))==1:
                    print "Current value as below"
                    print "Low: %6.3f High: %6.3f "%(StressAtLow, StressAtHigh)
                    StressAtLow = float(raw_input("StressAtLow [MPa] >> "))
                    StressAtHigh = float(raw_input("StressAtHigh [MPa]>> "))
                else:
                    Ex1 = float(raw_input("type Ex1[Gpa]"))* 1000
                    Ex2 = float(raw_input("type Ex2[Gpa]"))* 1000
                    Ey1 = float(raw_input("type Ey1[Gpa]"))* 1000
                    Ey2 = float(raw_input("type Ey2[Gpa]"))* 1000
                    break
                pass
            elif iflag=='y': break
        else: break
        pass

    xep1 = xe1 - xs1/Ex1
    xep2 = xe2 - xs2/Ex2    
    yep1 = ye1 - ys1/Ey1
    yep2 = ye2 - ys2/Ey2

    ## Instantaneous strain rate calculation
    xrep1 = sr(xep1, time, bin)
    xrep2 = sr(xep2, time, bin)    
    yrep1 = sr(yep1, time, bin)
    yrep2 = sr(yep2, time, bin)


    xrs1 = sr(xs1, time, bin)
    xrs2 = sr(xs2, time, bin)    
    yrs1 = sr(ys1, time, bin)
    yrs2 = sr(ys2, time, bin)    

    
    if iplt==True:
        stp = 100
        fig = plt.figure(99); ax = fig.add_subplot(111)
        ax.plot(time[::stp], xs1[::stp], 'd', mfc='None',
                label=r'$\varepsilon_{x1}$')
        ax.plot(time[::stp], xs2[::stp], 'o', mfc='None',
                label=r'$\varepsilon_{x2}$')
        ax.plot(time[::stp], ys1[::stp], 'd', mfc='None',
                label=r'$\varepsilon_{y1}$')
        ax.plot(time[::stp], ys2[::stp], 'o', mfc='None',
                label=r'$\varepsilon_{y2}$')
        ax.legend(loc='best')
        
        # ax.plot(xep1[::stp],  xs1[::stp], 'd', mfc='None',
        #         label=r'$\sigma_{x1}\varepsilon_{x1}$')
        # ax.plot(xep2[::stp], xs2[::stp], 'o', mfc='None',
        #         label=r'$\sigma_{x2}\varepsilon_{x2}$')
        
        # ax.plot(yep1[::stp],  ys1[::stp], 'd', mfc='None',
        #         label=r'$\sigma_{y1}\varepsilon_{y1}$')
        # ax.plot(yep2[::stp], ys2[::stp], 'o', mfc='None',
        #         label=r'$\sigma_{y2}\varepsilon_{y2}$')
        # ax.legend(loc='best')
        pass
    
    return xs1,xs2,ys1,ys2,xep1,xep2,yep1,yep2,xrep1,xrep2,yrep1,yrep2,xrs1,xrs2,yrs1,yrs2

def slope(x,y,yl,yh):
    """
    Calculates the slope dy/dx in the prefined y bound by yl and yh

    argument:
    x,y,   yl, yh
    returns slope
    """
    idxl, idxh = 0, 0
    for i in range(len(y)):
        if y[i]<yl: idxl = i
        if y[i]<yh: idxh = i
        if y[i]>yh : break
        pass
    #print 'idxl, idxh =', idxl, idxh; raw_input()
    z = np.polyfit(x[idxl:idxh:], y[idxl:idxh:], 1)
    return z[0]

def trim(dlines,delimiter):
    dd = []
    for i in range(len(dlines)):
        cline = dlines[i].split(delimiter) # delimiter
        try:
            map(float, cline)
        except:
            return np.array(dd)
        else: dd.append(map(float,cline))
        pass
    """
    the raw csv file does not follow the convention that I thought.
    The trailing control vavlue block does not exist? Yes! 
    """
    #raise IOError
    return np.array(dd)
    

def header(head, delimiter=','):
    """
    returns trivial information from the header block
    """
    date = head[1].split(delimiter)[1].split()[0] #string
    time = head[1].split(delimiter)[1].split()[1] #string
    exptime = int(head[2].split(delimiter)[1]) #integer
    name = head[4].split(delimiter)[1]
    operator = head[5].split(delimiter)[1]
    sample = head[8].split(delimiter)[1]
    matl   = head[9].split(delimiter)[1]
    area   = float(head[14].split(delimiter)[1]) #float
    dt = float(head[16].split(delimiter)[1]) #float [micro second]
    dt = dt * (10**-3) # now it is in the second unit.
    ratio = map(float, [head[20].split(',')[1],head[20].split(',')[2]])  #ratio

    ## only some selected information is retrieved.
    return area, dt, ratio

def sr(e,t,b):
    """
    Calculates instantanouse strain rate vector
    Arguments
    e: strain
    t: time
    b: bin size
    """
    insr = []
    for i in range(len(e)):
        if i-b > 0:
            if i+b < int(len(e)):
                temp = np.polyfit( #linear fitting
                    t[i-b:i+b:],
                    e[i-b:i+b:],
                    1)[0]
                insr.append(temp)
            else: insr.append(0.)
        else: insr.append(0.)
        pass
    return np.array(insr)

def multiloader(l=5,h=15):
    """
    Load *.csv file in the current working directory

    ## to-be-returned variables from a single loading
    xep, yep, xs, ys, workx, worky, work, xrep, yrep, ed

    l=5,
    h=15
    """
    files = glob.glob('*.csv')
    files.sort()
    for f in files:
        ## post-process
        print 'Current file: %s'%f
        rst = loader(filename=f, delimiter=',',
                     StressAtLow=l, StressAtHigh=h)
        sratio  = rst[-1]
        Rst = rst[:-1:]
        print '\n'
        
        #----------------------------------------------
        #post-processed data saver.
        fn = f.split('.csv')[0] + '.bpp'
        header = 'sratio, %i, %i'%(sratio[0], sratio[1])
        fwrite(filename=fn, dat=np.array(Rst).transpose(), header=header) #Nothing but writes!
        pass
    pass

def fwrite(filename, dat, header):
    """
    Writes file with header!!
    """
    f = open(filename, 'w')
    f.write(' %s \n'%header)
    for c in range(len(dat)):
        for d in range(len(dat[c])):
            f.write(' %e '%dat[c][d])
            pass
        f.write('\n')
        pass
    f.close()
    print "filename: %s was saved\n\n"%filename
    pass
    
def analysis():
    """
    Provides the information of the generated '.bpp' files
    """
    files = glob.glob('*.bpp')
    maxwork=0.; minwork = 100000.
    indx = 6
    for f in files:
        dat = np.loadtxt(f, dtype='float',skiprows=1).transpose()
        if maxwork < max(dat[indx]):
            maxworkf = f
            maxwork = max(dat[indx])
            pass
        if minwork > max(dat[indx]):
            minwork = max(dat[indx])
            minworkf = f
            pass
        pass

    print 'maximum work: %s'%maxwork
    print 'maximum work filename: %s'%maxworkf
    print 'minimum work: %s'%minwork
    print 'minimum work filename: %s'%minworkf
    return maxwork, minwork, maxworkf, minworkf
    pass

    indx = rst
    for f in files:
        dat = np.loadtxt(f, dtype='float').transpose()
        if maxwork < max(dat[indx]):
            maxworkf = f
            maxwork = max(dat[indx])
            pass
        if minwork > max(dat[indx]):
            minwork = max(dat[indx])
            minworkf = f
            pass
        pass

def weq(wlv=None, stp=5, dump=-1000000, confidence=0.90):
    """
    Provides the work equivalent segments

    Arguments:
    wlv=None: work level
    stp=5   : number of work level (used only when wlv!=None)
    """
    if wlv==None:
        mx = analysis()[0] #maxwork
        wlv = np.linspace(0.00001, mx*0.999, stp)
        pass

    print "work levels are as below"
    print wlv #;raw_input('Press enter >>>>')

    ## ----------------------------------------------------------------------
    weq = dict()
    for iwlv in range(len(wlv)):
        weq['%6.3f'%wlv[iwlv]] = dict()
        #weq['%f'%wlv[i]]['data_0'%?]
        pass
    dat = data_dict() ## call def.data_dict()

    ## loop over iso work segment
    for iwlv in range(len(wlv)):
        ## loop over stress ratio
        for irat in range(len(dat.keys())):
            rat = 'ratio%i'%irat
            current_data_dict = dat['ratio%i'%irat]
            sratio = dat[rat]['sratio']
            nkeys = len(current_data_dict.keys())
            if nkeys<=1: pass
            else:
                for idata in range(nkeys-1):
                    if current_data_dict[
                        'data_%i'%idata]['mxw'] < wlv[iwlv]:
                        pass
                    else:
                        work = current_data_dict[
                            'data_%i'%idata]['work']
                        ys = current_data_dict[
                            'data_%i'%idata]['ys']
                        xs = current_data_dict[
                            'data_%i'%idata]['xs']
                        xe = current_data_dict[
                            'data_%i'%idata]['xe']
                        ye = current_data_dict[
                            'data_%i'%idata]['ye']
                        ed = current_data_dict[
                            'data_%i'%idata]['ed']
                        for index in range(len(work)):
                            if work[index]>wlv[iwlv]:break
                            pass

                        ## interpolate the values:
                        ## yield stress, Strain direction(ED)
                        x0 = work[index-1]; x1=work[index]
                        y0 = xs[index-1]; y1=xs[index-1]
                        ys_x = interpolate(
                            x0,x1,y0,y1, wlv[iwlv])
                        
                        y0 = ys[index-1]; y1=ys[index-1]
                        ys_y = interpolate(
                            x0,x1,y0,y1, wlv[iwlv])
                        
                        y0 = ed[index-1]; y1=ed[index-1]
                        ED = interpolate(
                            x0,x1,y0,y1, wlv[iwlv])
                        
                        #if the current ratio is not in the key-list
                        if not(weq['%6.3f'%wlv[iwlv]].has_key(rat)): 
                            weq['%6.3f'%wlv[iwlv]][rat] = {}
                            weq['%6.3f'%wlv[iwlv]][
                                rat]['ratio'] = sratio
                            weq['%6.3f'%wlv[iwlv]][
                                rat]['raw_data']={}
                            weq['%6.3f'%wlv[iwlv]][
                                rat]['raw_data']['ys']=[]
                            weq['%6.3f'%wlv[iwlv]][
                                rat]['raw_data']['ed']=[]
                            pass
                        else: pass
                        weq['%6.3f'%wlv[iwlv]][
                            rat]['raw_data'][
                            'ys'].append([ys_x, ys_y])
                        weq['%6.3f'%wlv[iwlv]][
                            rat]['raw_data']['ed'].append(ED)
                        pass
                    pass
                pass
            pass
        pass
    ## ----------------------------------------------------------------------
    # eliminate the empty work brance to the main dictionary
    
    keys = weq.keys()
    for i in keys:
        if len(weq[i])==0: weq.pop(i)
        pass
    for iwl in weq.keys():
        for irat in weq[iwl].keys():
            ysx = []; ysy = []; thet = []
            for ix in range(len(weq[iwl][irat]['raw_data']['ys'])):
                ysx.append(weq[iwl][irat]['raw_data']['ys'][ix][0])
                ysy.append(weq[iwl][irat]['raw_data']['ys'][ix][1])
                thet.append(weq[iwl][irat]['raw_data']['ed'][ix])
                pass
            # print ysx, ysy, thet
            # raw_input()
            weq[iwl][irat]['ysx_avg'] = np.mean(ysx)
            weq[iwl][irat]['ysy_avg'] = np.mean(ysy)
            weq[iwl][irat]['ed_avg'] = np.mean(thet)
            weq[iwl][irat]['ysx_conf'] = mean_confidence_interval(ysx,confidence)
            weq[iwl][irat]['ysy_conf'] = mean_confidence_interval(ysy,confidence)
            weq[iwl][irat]['ed_conf'] = mean_confidence_interval(thet,  confidence)
            pass
    
    ## ----------------------------------------------------------------------    
    return weq
    pass

def mean_confidence_interval(data, confidence=0.90):
    """
    T-distribution
    """
    from scipy import stats
    import numpy
    data = np.array(data)
    n = len(data)
    m, se = np.mean(data), stats.stderr(data)
    h = se * sp.stats.t._ppf((1+confidence)/2. , n-1)
    return h*2

def data_dict():
    """
    Returns the data dictionary.

    data are input from '*.bpp'
    and sorted based on stress ratio

    Returns the YSS_dict.
    YSS_dict dictionary has stress ratio attributes
    under which multiple number of data resides
    """

    files = glob.glob('*.bpp')
    files.sort()

    ssr = [] # stress ratio list
    for f in range(len(files)):
        tmpfile = open(files[f], 'r')
        header = tmpfile.readline(); tmpfile.close()
        sratio = [int(header.split(',')[1]),
                  int(header.split(',')[2])]
        #Only unique set of stress ratio is found
        # and put into ssr among '*.bpp' files
        if any(ssr[i]==sratio for i in range(len(ssr))):pass
        else: ssr.append(sratio) 
        pass
    print "total sratio as below"
    print ssr
    ssr.sort() #  [[1,2],[1,4],... ]
    #raw_input()
    #print 'sratio: %i %i'%(sratio[0], sratio[1])

    
    #recieve data and sorting based on their stress ratio
    ## dictionary generator
    ## yield stresses are saved
    ## according to the stress ratios
    YSS_dict = {}
    for i in range(len(ssr)):
        YSS_dict['ratio%i'%i]= dict()
        YSS_dict['ratio%i'%i]['sratio'] = ssr[i]
        pass

    ## loop over each of '*.bpp' files in the folder
    for f in range(len(files)):
        dat = dict()
        
        raw_data = np.loadtxt(
            files[f], dtype='float',
            skiprows=1).transpose()
        
        tmpfile = open(files[f], 'r')
        header = tmpfile.readline(); tmpfile.close()
        header = header.split(',')[1::]
        sratio = map(int, [header[0],header[1]])

        ## add keys to the dictionary
        dat['sratio'] = sratio
        dat['filename'] = files[f]
        dat['mxw'] = max(raw_data[6])
        dat['work'] = raw_data[6].copy()
        
        dat['xe'] = raw_data[0]
        dat['ye'] = raw_data[1]        
        dat['xs'] = raw_data[2]
        dat['ys'] = raw_data[3]
        dat['ed'] = raw_data[9]

        ## loop over each of the yield stress dictionary
        for i in range(len(YSS_dict)):
            ## if encouters the same stress ratio
            if YSS_dict[
                'ratio%i'%i]['sratio']==dat['sratio']:
                j = 0
                while True:
                    if YSS_dict[
                        'ratio%i'%i].has_key('data_%i'%j):
                        ## if yss alread has the data
                        j = j + 1
                        pass
                    else:
                        ## or it is first inciden to encouter
                        ## the data, add dat dictionary
                        ## to YSS_dict dictionary
                        YSS_dict[
                            'ratio%i'%i]['data_%i'%j] = dat
                        break
                    pass
                pass
            pass
        pass
    ## in YSS_dict dictionary, there is a 'ratio' attribute
    ## under which the corresponding data reside.
    return YSS_dict

def interpolate(x0, x1, y0, y1, x):
    """
    Interpolates the stress, strain, whatever with given points
    """
    slope = (y1-y0)/(x1-x0)
    value = slope * (x-x0) + y0
    return value


def ysplot(wlv=[0.1,10, 20, 30],
           confidence=0.95, ifig=None):
    """
    Provides the work level
    based on wlv
    Provides the errorbar based on confidence
    level.
    plots them in the matplotlib figure
    whose id is indicated as ifig
    unless ifig==None
    """
    if ifig==None: pass
    else:
        fig = plt.figure(ifig)
        ax = fig.add_subplot(111)
        pass

    ## iso work segment dictionary
    dat = weq(wlv, confidence)
    al = dat.keys() 
    al.sort() # All levels of work is sorted

    # clist = np.log(np.linspace(0, 3, len(wlv)))
    # colors= []
    # for i in range(len(clist)):
    #     colors.append(plt.cm.cmap_d['winter_r'](clist[i]))
    #     pass

    colors = ['b','g','r','c','m','y','k']

    eqwc = [] #eq-work contour
    for wlv in range(len(al)):
        rl = dat[al[wlv]].keys(); rl.sort()
        eqwc.append([])
        for ir in rl: #loop over stress ratio
            xm = dat[al[wlv]][ir]['ysx_avg']
            ym = dat[al[wlv]][ir]['ysy_avg']
            hx = dat[al[wlv]][ir]['ysx_conf']
            hy = dat[al[wlv]][ir]['ysy_conf']

            # assigns the mean value and error to a
            # iso-work segment
            eqwc[wlv].append([[xm,ym],[hx,hy]])

            if ifig==None: pass
            else: ax.errorbar(
                xm, ym, xerr=hx,
                yerr=hy, fmt='o',#mfc='None',
                color=colors[wlv],
                mec=colors[wlv],
                mfc='None')
            pass
        pass
    if ifig==None:pass
    else:
        ax.set_xlim(-10,)
        ax.set_ylim(-10,)
        ax.set_aspect('equal')
        pass
    return eqwc

def eqwcplot(wlv=[0.1, 2, 3, 8], confidence=0.95, ifig=1):
    """
    plot the equivalnet work contour based on given work level list.
    Basically, it read all '*.bpp' files in the folder
    and sorts them by their stress ratio in the header (1st line of the file)
    Then, calculates the mean value and the error based on the Student's
    t-distribution with a confidence level passed as 'confidence' argument.

    Then plots the data in a figure whose id is given as 'ifig'
    using the matplotlib library for the Python interpreter.
    The normalized work contour will be plotted in the ifig+1 figure.
    
    Arguments:
    wlv = [0.1, 2, 3,8 .. ]
    confidence = 0.95
    ifig = 1
    """
    colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
    markers = ['o','^','d','p','h','*','+','o','^','d','p','h','*','+']
    # markers = ['+','*',',','.','1','2','3',
    #            '4','<','>','D','H','^','_','h']

    ## eqwc is consisted of iso-work segment
    ## for all possible stress ratio
    eqwc = ysplot(wlv, confidence, ifig=None) 
    fig = plt.figure(ifig,figsize=(8,6))
    fignorm = plt.figure(ifig+1,figsize=(8,6))

    ## ax = fig.add_subplot(111)
    ax = fig.add_axes((0.12,0.12,0.65,0.8))
    axnorm = fignorm.add_axes((0.12,0.12,0.65,0.8))

    ## xnorm
    for i in range(len(eqwc)):
        cwc = np.array(eqwc[i]).transpose()
        x = cwc[0][0]
        y = cwc[1][0]
        xe = cwc[1][1]
        ye = cwc[0][1]
        ## finds the loading along rd from which the
        ## normalization will take place.

        for j in range(len(y)):
            if y[j]==0:
                xnorm = x[j]
                break
            pass
        ax.errorbar(x, y, xe, ye, marker=markers[i],
                    color='Black', #color = colors[i],
                    #mec = colors[i],
                    label='%5.2f'%wlv[i],
                    mfc='None',
                    ls=' ', #no line
                    markersize=10.,
                    markeredgewidth=1.,
                    #label=r'work$^{pl}$=%6.2f'%wlv[i]
                    )

        axnorm.errorbar(x/xnorm, y/xnorm, xe/xnorm, ye/xnorm, marker=markers[i],
                        color='Black', #color = colors[i],
                        #mec = colors[i],
                        label='%5.2f'%wlv[i],
                        mfc='None',
                        ls=' ', #no line
                        markersize=10.,
                        markeredgewidth=1.,
                        #label=r'work$^{pl}$=%6.2f'%wlv[i]
                        )            
        pass
    ax.legend(loc=2, bbox_to_anchor=(1.05, 1.0))
    ax.set_xlabel(r'$\sigma_{RD}$ [MPa]',dict(fontsize=20))
    ax.set_ylabel(r'$\sigma_{TD}$ [MPa]',dict(fontsize=20))
    ax.set_aspect('equal')
    ax.set_xlim(-30,);ax.set_ylim(-30,)
    ax.text(x=5, y=25, s='confidence: %3.1f'%(confidence*100)+'%')
    axnorm.legend(loc=2, bbox_to_anchor=(1.05, 1.0))
    axnorm.set_xlabel(r'$\sigma_{RD}/\bar{\sigma}^{YS}_{RD}$',dict(fontsize=20))
    axnorm.set_ylabel(r'$\sigma_{TD}/\bar{\sigma}^{YS}_{RD}$',dict(fontsize=20))
    axnorm.set_aspect('equal')
    axnorm.set_xlim(-0.05,);axnorm.set_ylim(-0.05,)
    axnorm.text(x=5, y=25, s='confidence: %3.1f'%(confidence*100)+'%')
    pass
        
def uniaxial(filename):
    """
    Provides the dat file in the same format
    that is input to the biaxial case.
    '*.bpp' with proper heading (the stress ratio,
    in uniaxial case, it is either 1:0 or 0:1)
    """
    print 'file: %s'%filename
    while True:
        iaxis = raw_input("axis?  (x or RD: 0,  y or TD: 1) >>")
        if iaxis=='0' or iaxis=='1':
            break
        else: print "Type correctly again"
        pass

    if iaxis=='0': sratio = [1,0]
    else: sratio = [0,1]
    
    dat = np.loadtxt(filename, skiprows=2, dtype='string')
    dat = dat.transpose()
    time = dat[4] ; 'time'    
    if iaxis=='0':
        #xe = dat[0] ; 'True axial strain'
        xe = dat[1].copy() ; 'plasitic strain'
        ye = dat[9].copy() ; 'plasitic strain'
        xs = dat[2].copy() ; 'True stress'
        ys = np.zeros((len(xs),)) ; 'True stress along y'
        workx = dat[3].copy() ; 'plastic work along x'
        worky = np.zeros((len(workx))) ; 'plastic work along y'
        work = workx.copy()
        Xrep = dat[5].copy() ; 'plastic strain rate'
        Yrep = dat[10].copy() ; 'plastic strain rate'
        pass
    
    elif iaxis=='1':
        ye = dat[1].copy() ; 'plasitic strain along y'
        xe = dat[9].copy() ; 'plasitic strain along x'
        ys = dat[2].copy() ; 'True stress'
        xs = np.zeros((len(ys),)) ; 'True stress along y'
        worky = dat[3].copy() ; 'plastic work along x'
        workx = np.zeros((len(worky),)) ; 'plastic work along y'
        work = worky.copy()
        Yrep = dat[5].copy() ; 'plastic strain rate'
        Xrep = dat[10].copy() ; 'plastic strain rate'
        pass

    for i in range(len(Xrep)):
        if Xrep[i]=='NA': Xrep[i]='0'
        if Yrep[i]=='NA': Yrep[i]='0'
        if ye[i]=='NA': ye[i] ='0'
        if xe[i]=='NA': xe[i] ='0'
        if ys[i]=='NA': ys[i] ='0'
        if xs[i]=='NA': xs[i] ='0'
        pass
    
    header = 'sratio, %i, %i'%(sratio[0], sratio[1])
    
    Yrep = np.array(Yrep, dtype='float')
    Xrep = np.array(Xrep, dtype='float')
    
    xe = np.array(xe, dtype='float')
    ye = np.array(ye, dtype='float')
    
    xs = np.array(xs, dtype='float')
    ys = np.array(ys, dtype='float')
    
    workx = np.array(workx, dtype='float')
    worky = np.array(worky, dtype='float')    
    work = np.array(work, dtype='float')

    ed = np.arctan2(Yrep, Xrep)
    
    data = np.array((xe, ye, xs, ys, workx, worky, work, Xrep, Yrep, ed), dtype='float')
    fn = '%s.bpp'%filename.split('.str')[0]
    fwrite(fn, data.transpose(), header)
    
    return data


def UniMultiLoader():
    """
    Uniaxial result multi loader
    Loads *.str files
    """
    files = glob.glob('*.str')
    for i in files:
        uniaxial(filename=i)
        pass
    pass


def TotUniTension(wildcard, stp=10, iplot=False, confidence=0.95,
                  sfout='eps-sig_meancurve',
                  rfout='eps-rvl_meancurve',
                  ):
    """
    Average the uni file ".bpp"
    
    arguments
    wildcard = '*RD*.bpp'
    stp = 100
    iplot = False
    confidence = 0.95
    sfout ='eps-sig_meancurve'
    rfout ='eps-rvl_meancurve'
    """
    files = glob.glob(wildcard)
    eps_f = 0.2;eps_i = 0.001
    eps = []; sig = []; rvl = []
    
    for i in files:
        e, s, r = UniTension(i) #eps, sig, R-value
        if iplot==True:
            ifig = 96
            fig = plt.figure(ifig,figsize=(6.0,8.575))
            ax0 = fig.add_subplot(211)
            ax1 = fig.add_subplot(212)        
            ax0.plot(e,s)
            ax1.plot(e,r)
            ax1.set_ylim(0.,2.)
            pass
        for ii in range(len(e)):
            indf = -1
            if e[ii]>eps_f: indf = ii; break
            pass
        for ii in range(len(e)):
            indi = 0
            if e[ii]>eps_i: indi = ii; break
            pass

        eps.append(e[indi:indf:stp])
        sig.append(s[indi:indf:stp])
        rvl.append(r[indi:indf:stp])

        if iplot==True:
            # ifig = 95
            # fig = plt.figure(ifig,figsize=(6.0,8.575))
            # ax0 = fig.add_subplot(211)
            # ax1 = fig.add_subplot(212)        
            ax0.plot(e[indi:indf:stp],s[indi:indf:stp],'o')
            ax1.plot(e[indi:indf:stp],r[indi:indf:stp],'o')
            ax1.set_ylim(0., 2.)
            pass
        pass

    ## close the plots
    flag=raw_input('press enter to close the figures>>>')
    if flag=='y' or flag=='Y' or len(flag)==0: plt.close(ifig)

    ## Average curve calculatation
    npoints= 50
    x_points = np.linspace(eps_i, eps_f, npoints)
    sig_points = np.zeros((len(x_points),2))
    rvl_points = np.zeros((len(x_points),2))

    for ix in range(len(x_points)):
        sig_mean = 0
        rvl_mean = 0
        sig_tmp = []; rvl_tmp = []
        for seg_id in range(len(eps)):
            x = x_points[ix]
            #y = unknown: aim
            current_segment_strain = eps[seg_id]
            current_segment_stress = sig[seg_id]
            current_segment_rvalue = rvl[seg_id]
            low_index=0
            for each_point in range(len(current_segment_strain)):
                if current_segment_strain[each_point] > x:
                    low_index=each_point ; break
                    pass
                pass
            # x linear space
            x0 = current_segment_strain[low_index-1]
            x1 = current_segment_strain[low_index]

            ## sigma
            y0 = current_segment_stress[low_index-1]
            y1 = current_segment_stress[low_index]
            slope = (y1 - y0)/(x1 - x0)
            sig_tmp.append(slope * (x-x0) + y0)

            ## R-value        
            y0 = current_segment_rvalue[low_index-1]
            y1 = current_segment_rvalue[low_index]
            slope = (y1 - y0)/(x1 - x0)
            rvl_tmp.append(slope * (x-x0) + y0)
            pass

        sig_mean = np.array(sig_tmp).mean()
        rvl_mean = np.array(rvl_tmp).mean()

        sig_err = mean_confidence_interval(sig_tmp, confidence)
        rvl_err = mean_confidence_interval(rvl_tmp, confidence)
        
        sig_points[ix][0] = sig_mean
        sig_points[ix][1] = sig_err
        rvl_points[ix][0] = rvl_mean
        rvl_points[ix][1] = rvl_err
        pass

    if iplt==True:
        fig = plt.figure(ifig,figsize=(6.0,8.575))
        ax0 = fig.add_subplot(211)
        ax1 = fig.add_subplot(212)

        ax0.errorbar(x_points,sig_points.transpose()[0], yerr=sig_points.transpose()[1])
        ax1.errorbar(x_points,rvl_points.transpose()[0], yerr=rvl_points.transpose()[1])
        ax0.set_ylim(0.,)
        ax1.set_ylim(0.,)
        flag==raw_input("Press enter to close >>> ")
        if flag=='y' or len(flag)==0: plt.close(ifig)
        else: pass
        pass
    
    ## sig_points and rvl points include the confidence interval

    mean_sig  = sig_points.transpose()[0].copy()
    err_sig   = sig_points.transpose()[1].copy()
    mean_rvl  = rvl_points.transpose()[0].copy()
    err_rvl   = rvl_points.transpose()[1].copy()
    mean_eps = x_points.copy()
    sdata_to_file = np.array([mean_eps, mean_sig, err_sig])
    rdata_to_file = np.array([mean_eps, mean_rvl, err_rvl])
    print "The mean stress-strain Written : '%s'"%sfout
    print "The mean Rvalue-strain Written : '%s'"%rfout
    np.savetxt(sfout, sdata_to_file.transpose())
    np.savetxt(rfout, rdata_to_file.transpose())
    return x_points, sig_points, rvl_points # x _points: eps (strain)

def UniTension(filename):
    """
    Averages the uni Tension curves.
    stress-strain curve
    R-value - strain curve
    """
    F = open(filename)
    header = F.readline();F.close()
    tmp = header.split(',')
    sratio = map(int,[tmp[1],tmp[2]])
    dat = np.loadtxt(filename, skiprows=1).transpose()

    if sratio==[1,0]:
        'RD'
        eps = dat[0]
        sig = dat[2]
        eps_rl = dat[7]
        eps_rw = dat[8] #transverse strain rate
        eps_rt = -(eps_rl + eps_rw)
        pass
    elif sratio==[0,1]:
        'TD'
        eps = dat[1]
        sig = dat[3]
        eps_rl = dat[8]
        eps_rw = dat[7] #transverse strain rate
        eps_rt = -(eps_rl + eps_rw)
        pass

    r = eps_rw/ eps_rt
    #r = 1./r
    return eps, sig, r

    
    
    
    
    
def eps_work(filename=None,mode='RD'):
    data = np.loadtxt(filename, skiprows=1).transpose()
    if mode=='RD':
        return data[0], data[6]
    elif mode=='TD':
        return data[1], data[6]
    pass



def epsatwork(wildcard='RD?.bpp', mode='RD',
              worklevels = [],confidence=0.95):
    """
    Calculates the strain at certain work levels given through worklevels=[]
    For the prepartion of MSMSE
    """
    
    files = glob.glob('*%s*'%wildcard)
    worklevels = np.array(worklevels)

    w_data = []
    yerr_data = []
    #mean value + yerr (t-distribution)
        
    for worklevel in worklevels:
        interpolated_works = []
        for filename in files:
            eps, work = eps_work(filename, mode)
            ## Interpolates the work---------
            index = 0
            for iy in range(len(work)):
                if max(work)< worklevel:
                    raise IOError, "Work level exceeding..."
                if work[iy]>worklevel:
                    index = iy
                    break
                pass
            w0 = work[index-1];  w1 = work[index]
            eps0 = eps[index-1]; eps1 = eps[index]
            slope = (eps1-eps0)/(w1-w0)
            y = slope * (worklevel - w0) + eps0
            ##------------------------------
            interpolated_works.append(y)
            pass
        
        mean = np.array(interpolated_works).mean()
        error = mean_confidence_interval(
            interpolated_works,
            confidence)
        w_data.append(mean); yerr_data.append(error)
        pass

    fig = plt.figure(33)
    ax =fig.add_subplot(111)
    # print 'worklevels', worklevels
    # print 'w_data', w_data
    # print 'yerr_data', yerr_data
    # raw_input()
    ax.errorbar(x=worklevels, y=w_data, yerr=yerr_data)
    
    return worklevels, w_data, yerr_data
        
    
    
    
