"""
Automatizes the biaxial test results post-processing
Written by Youngung
"""

import os
import glob
import math
import biaxial_pp; reload(biaxial_pp)
#os.sys.path.append('c:\\python26\\scripts')
os.sys.path.append('/home/youngung/Dropbox/py')

try:
    biaxial_pp
except NameError:
    import biaxial_pp
else:
    reload(biaxial_pp)

def bpp_tot():
    lsx = raw_input('lsx >>')
    usx = raw_input('usx >>')
    lsy = raw_input('lsy >>')
    usy = raw_input('usy >>')
    print
    print '*******************************************************************'
    print '              Strain and load cell selection to be used '
    print " Defaults are 'plavg' and 'avg' respectively for strain and stress "
    print '*******************************************************************'
    print
    epsxl = raw_input('epsxl >>')
    epsyl = raw_input('epsyl >>')
    sigxl = raw_input('sigxl >>')
    sigyl = raw_input('sigyl >>')
    n   = raw_input('vector size? (MPa) >>')
    if len(lsx) == 0: lsx = 20.
    else: lsx = float(lsx)
    if len(usx) == 0: usx = 50.
    else: usx = float(usx)
    if len(lsy) == 0: lsy = 20.
    else: lsy = float(lsy)
    if len(usy) == 0: usy = 50.
    else: usy = float(usy)
    if len(n)   == 0: n = 10.
    else: n = float(n)

    if len(epsxl) == 0 : epsxl ='plavg'
    if len(epsyl) == 0 : epsyl ='plavg'
    if len(sigxl) == 0 : sigxl ='avg'
    if len(sigyl) == 0 : sigyl ='avg'

    #files = glob.glob('*.csv')
    bpp = []
    #w=[1.0, 2.0, 4.0, 8.0, 10.0, 14.0]
    files, w = [],[]
    times = []
    delt = raw_input('Bin size in terms of number of data points (default:20) >>')
    if len(delt)==0: delt = 20
    else: delt = int(delt)

    while True:
        if raw_input('Input file (bpp_tot.in) ? (y or n) >>') =='n':
            ifile = False
            tmp = raw_input('type the files >>')
            if len(tmp) == 0:
                break
            else:
                files.append(tmp)
        else:
            ifile = True
            f = open('bpp_tot.in', 'r')
            f.readline()
            while True:
                source = f.readline()
                if len(source) < 3:
                    break
                else:
                    files.append(source.split()[0])
                    try: source.split()[1]
                    except IndexError:
                        times.append(None)
                    else:
                        times.append(float(source.split()[1]))
            break
            
                
    if len(files) == 0: files = glob.glob('*.csv')

    while True:
        tmp = raw_input('type the work level >> ')
        if len(tmp)==0:
            break
        else:
            tmp = float(tmp)
            w.append(tmp)

    if len(w) == 0:
        w = [1.0, 2.0, 4.0, 8.0, 10.0, 14.0, 20, 36, 46, 56, 66]

    for i in range(len(files)):
        """
        if files[i] == 'No10_uniy.csv':
            raw_input('flipping happened.. press the enter')
            bpp.append(biaxial_pp.biaxial(filename = files[i], path=os.getcwd(), iflip=True, lsx= lsx, usx=usx, lsy=lsy, usy=usy))
            bpp[i].out(filename=bpp[i].filename.split('.')[0]+'__pp.out',path=os.getcwd(), uni='x1')       
        elif files[i] =='No14_RD_unix.csv':
            raw_input('flipping happened.. press the enter')
            bpp.append(biaxial_pp.biaxial(filename = files[i], path=os.getcwd(), iflip=True, lsx= lsx, usx=usx, lsy=lsy, usy=usy))
            bpp[i].out(filename=bpp[i].filename.split('.')[0]+'__pp.out',path=os.getcwd())
        elif files[i] =='10-1_316L_cruciform_23':
            bpp.append(biaxial_pp.biaxial(filename = files[i], path=os.getcwd(), iflip=True, lsx= lsx, usx=usx, lsy=lsy, usy=usy))
            bpp[i].out(filename=bpp[i].filename.split('.')[0]+'__pp.out',path=os.getcwd())        
        else: bpp.append(biaxial_pp.biaxial(filename = files[i], path=os.getcwd(), iflip=False, lsx= lsx, usx=usx, lsy=lsy, usy=usy))
        """

        ### Declare a looping bpp object considering the maxtime given as augmented possibly through 'bpp_pp.in' file.
        if ifile==True:
            bpp.append(biaxial_pp.biaxial(filename = files[i], path=os.getcwd(), iflip=False,
                                          lsx= lsx, usx=usx, lsy=lsy, usy=usy, maxtime=times[i]))
        elif ifile == False:
            bpp.append(biaxial_pp.biaxial(filename = files[i], path=os.getcwd(), iflip=False,
                                          lsx= lsx, usx=usx, lsy=lsy, usy=usy))
        else:
            print 'Wrong ifile input'

        trash = bpp[i].sr(epsxl='plavg', epsyl='plavg', delt=delt)
        bpp[i].out(filename=bpp[i].filename.split('.')[0]+'__pp.out',
                   path=os.getcwd(), uni=False)
        del trash

    for i in range(len(w)):
        f = open('w_'+str(w[i])+'.weq', 'w')
        f.writelines('work contour at the work of '+str(w[i]) + '\n')
        f.write('filename                   Sx   Sy  Ex   Ey        sig_x        sig_y')
        f.write('      sig_x`      sig_y`        theta   indx         time\n')
        
        for j in range(len(files)):
            sx, sy, theta, tx, t = bpp[j].yield_vector(w=w[i], sigxl=sigxl,
                                                       sigyl=sigyl, epsyl=epsyl,
                                                       epsxl=epsxl, iplot=False)
            print 'assigned maxtime of this file is ', times[j]
            #if t >= times[j]:  pass
            if times[j]*0.99 < t:   #bpp[j].uni.time_flow[-1]:
                print
                print '****************************************'
                print 'It is beyond the maxtime of the raw data'
                print '****************************************'
                print
                pass
            else:
                f.write(bpp[j].filename.split('.')[0][0:25].ljust(25))
                f.write('%4i %4i'%(bpp[j].uni.stressratio[0] , bpp[j].uni.stressratio[1]))
                f.write('%4i %4i'%(bpp[j].uni.strainratio[0] , bpp[j].uni.strainratio[1])) 
                f.write('   %10.4e   %10.4e  %10.4e  %10.4e %12.4e'%
                        (sx,sy,sx+n*math.cos(theta*math.pi/180.),
                         sy+n*math.sin(theta*math.pi/180.), theta))
                f.write('   %4i   %10.2f \n'%(tx, t))
            
        f.close()
        
    del bpp
    del files
    del w        

def __ratio__(extension='weq', headern=2, filename=None):
    """
    Returns possible ratioes from those files whose extension is given.
    For example, once 'weq' is given as the extension argument,
    glob.glob() brings all the files into a list of files, and
    'def ratio()' scan through all the files listed to gather all the possible ratioes
    in the active folder.
    """
    # files
    weq_files = []
    if filename==None:
        weq_filelist = glob.glob('*.'+extension)
    else:
        weq_filelist = [filename]
    
    # Group of the sigma or epsilon ratios
    ratio = []
    
    for i in range(len(weq_filelist)):
        weq_files.append(open(weq_filelist[i], 'r'))

        for n in range(headern):
            weq_files[i].readline()   #two liens of header are treated to be dummy lines
        
        while True:
            line = weq_files[i].readline()
            if len(line) < 3:
                break
            else:
                tmp = line.split()
                rat = tmp[1] + tmp[2] + tmp[3] + tmp[4]
            if any(ratio[i] == rat for i in range(len(ratio))): pass
            else: ratio.append(rat)
        weq_files[i].close()

    """
    # Based on the ratio list, declare variable dynamically
    # rat_var includes variables like rat_1100, rat_1200, .. and so on
    # which would have been dynamically made
    rat_var = []
    for i in range(len(ratio)):
        globals()['rat_' + ratio[i]] = []
        
        rat_var.append(globals()['rat_' + ratio[i]])
    # total kinds of ratio is, then, len(rat_var)
    """
    #print 'total number of ratio', 'and rat_var ', " and 'ratio' "
    return ratio
            
def __weq_eqw_combine__(uni_ext = 'eqw', bi_ext='weq'):
    """
    Finds results from tension_r written down to files with weq extension
    and append them into weq files.
    This is to be done before average the data with avg method.

    ** Arguments:
       uni_ext = 'eqw'
       bi_ext = 'weq'

       These are all for files with different extensions
    """
    tension_files = glob.glob('*.'+uni_ext)
    biaxial_files = glob.glob('*.'+bi_ext)
    if len(tension_files) == 0:
        print
        print '*******************************************'
        print 'NOTICE!!'
        print 'no uniaxial results in the working directory'
        print '*******************************************'
        print
        return -1
    else:
        if len(tension_files) == len(biaxial_files):
            for i in range(len(tension_files)):
                temp_files = file(tension_files[i].split(uni_ext)[0] + 'tmp', 'w')
                bi_f = file(biaxial_files[i], 'r')
                bi_source = bi_f.read()  #bi_source is a bulky block of data
                uni_f = file(tension_files[i], 'r')
                uni_source = uni_f.readlines()  #uni_source is a data-line list.
                
                temp_files.write(bi_source)

                #while True:
                for k in range(len(uni_source)):
                    if k == 0: pass
                    else:
                        if len(uni_source[k]) < 3: break
                        else:
                            temp_files.writelines(uni_source[k].split()[0].ljust(25))
                            temp_files.writelines(
                                '%4i %4i'%
                                (int(uni_source[k].split()[1]),
                                 int(uni_source[k].split()[2])))
                            temp_files.writelines(
                                '%4i %4i'%
                                (int(uni_source[k].split()[3]),
                                 int(uni_source[k].split()[4])))
                            temp_files.writelines(
                                '   %10.4e   %10.4e\n'%
                                (float(uni_source[k].split()[5]),
                                 float(uni_source[k].split()[6])))
                            #temp_files.writelines(uni_source[k])
                            """
                            f.write(bpp[j].filename.split('.')[0][0:25].ljust(25))
                            f.write('%4i %4i'%(bpp[j].uni.stressratio[0] , bpp[j].uni.stressratio[1]))
                            f.write('%4i %4i'%(bpp[j].uni.strainratio[0] , bpp[j].uni.strainratio[1])) 
                            f.write('   %10.4e   %10.4e  %10.4e  %10.4e %12.4e'%(sx,sy,sx+n*math.cos(theta*math.pi/180.), sy+n*math.sin(theta*math.pi/180.), theta))
                            f.write('   %4i   %10.2f \n'%(tx, t))
                            """
                # closing the current working file
                
                tmp_name = temp_files.name
                tmp_name_bi = bi_f.name
                temp_files.close()
                bi_f.close()
                uni_f.close()
                """
                print tmp_name
                print tmp_name_bi
                raw_input()
                os.system('del '+tmp_name_bi)
                os.system('ren '+tmp_name+' '+tmp_name_bi)
                os.system('del '+tmp_name)
                """
        else:
            print
            print '****************************************************************************'
            print ' The number of work levels are different between uni- and bi- results'
            print ' Please make corrections for this prior to further post-processing the data'
            print '****************************************************************************'
            print

def __avg_out__(extension='weq', headern=2, filename = None):
    """
    Finds files whose extension is given and average the values
    whose ratio, consisting of stress and strain ones,
    is the same.

    Finally the averaged values for respective ratio is again written down
    to the files of the same name but different extension, '.wvg'

    Arguments:
       extension = 'weq'
       headern = 2   : the number of header lines in the files of the given extension
       filename = None
          One can assign a filename if a particular file is of interest
    """

    ## files
    weq_files = []    
    if filename==None:
        
        weq_filelist = glob.glob('*.'+extension)

        for i in range(len(weq_filelist)):
            weq_files.append(open(weq_filelist[i],'r'))

    else:
        weq_files.append(open(filename, 'r'))

    ## available ratio combinations
    rat = __ratio__(extension, headern, filename=filename)   # rat is now a list includes all the possible ratios

            
    for i in range(len(rat)):
        globals()['rat_'+rat[i]] = []
        
    for i in range(len(weq_files)):
        fout = open(weq_files[i].name.split(extension)[0][0:-1]+'.wvg','w')
        #print weq_files[i].readline()
        #print weq_files[i].readline()
        for n in range(headern):
            weq_files[i].readline()
            
        fout.writelines('Avg contour\n ')
        fout.writelines(' Sx  Sy  Ex  Ey      sigx      sigy\n')        
        print
        print '******************'
        print weq_files[i].name
        print '******************'
        while True:
            line = weq_files[i].readline()
            if len(line)<3:
                break
            else:
                tmp = line.split()
                r = tmp[1]+tmp[2]+tmp[3]+tmp[4]
                sx, sy = float(line.split()[5]), float(line.split()[6])
                globals()['rat_'+r].append([sx,sy])

        sigma = []
        epsilon = []
        
        for j in range(len(rat)):
            sum_sx, sum_sy = 0, 0
            for k in range(len(globals()['rat_'+rat[j]])):
                #print globals()['rat_'+rat[j]]
                #print globals()['rat_'+rat[j]]
                sum_sx = sum_sx + globals()['rat_'+rat[j]][k][0]
                sum_sy = sum_sy + globals()['rat_'+rat[j]][k][1]
                pass
            try:  avg_sx = sum_sx/len(globals()['rat_'+rat[j]])
            except ZeroDivisionError: avg_sx = 'NA'
            try: avg_sy = sum_sy/len(globals()['rat_'+rat[j]])
            except ZeroDivisionError: avg_sy = 'NA'
            
            #print rat[j]
            print ' Sx  Sy  Ex  Ey'
            print rat[j][0].rjust(3), rat[j][1].rjust(3), rat[j][2].rjust(3), rat[j][3].rjust(3)
            #print rat[j]
            print avg_sx, avg_sy
            fout.writelines(
                rat[j][0].rjust(4) + rat[j][1].rjust(4) + rat[j][2].rjust(4) + rat[j][3].rjust(4))
            try: fout.writelines('%10.2f%10.2f \n'%(avg_sx, avg_sy))
            except TypeError:
                print 'avg_sx, avg_sy =', avg_sx, avg_sy
                break
            #raw_input()
        weq_files[i].close()
        fout.close()
    
    pass

def avg(extension='weq', uni_ext = 'eqw', headern=2):
    """
    Returns the averaged out values for each of files resulting from the testings
    with a multiple loadings for each possibility of ratio.
    Detects files through def avg_out().
    
    Arguments:
       extension = 'weq'
       uni_ext = 'eqw'      :extension of uniaxial tension post-processed result files
       headern = 2
    """

    # detects if uniaxial results are in the current working directory.
    # if so, append the results to the files having the given extension.
    files = glob.glob('*.'+uni_ext)
    if len(files) == 0:
        pass
    elif len(files) > 0:
        print
        print 'Appending uniaxial results to biaxial reult files'
        print
        raw_input()
        __weq_eqw_combine__(uni_ext = uni_ext, bi_ext=extension)
        extension = 'tmp'
    
    # average the data having the same control ratio
    files = glob.glob('*.'+extension)
    for i in range(len(files)):
        __avg_out__(filename=files[i], extension=extension, headern=headern)

def norm(extension ='wvg', headern = 2):
    """
    normalize the resulting wvg files
    """
    filelist = glob.glob('*.'+extension)
    files = []
    for i in range(len(filelist)):
        isgood=False
        files.append(open(filelist[i], 'r'))
        
        #file info
        print 'file name: ', files[i].name
        tmp_line = []
        iline = 0
        for n in range(headern):
            files[i].readline()
        while True:
            line = files[i].readline()
            #print ' line =', line
            r = []
            if len(line.split())==0:
                if isgood:
                    break
                else:
                    print '*************************'
                    print 'No uniaxial loading found'
                    print 'Returns -1'
                    print '*************************'
                    return -1
            else:
                iline = iline + 1
                tmp_line.append(line)
                r.append(line.split()[0])
                r.append(line.split()[1])
                r.append(line.split()[2])
                r.append(line.split()[3])
                #print 'r=', r
                if r == ['0','0','1','0']:
                    nrm = float(line.split()[4])
                    nline = iline
                    print 'nrm: ', nrm
                    isgood=True
                else:
                    pass
        fout = open(filelist[i].split(extension)[0][0:-1]+'.neq', 'w')
        fout.writelines(' Normalized EQ plastic work\n')
        fout.writelines(' Sx Sy Ex Ey     sigx     sigy\n')
        for j in range(len(tmp_line)):
            #print 'line = '
            #print line
            #raw_input()
            try:
                tmp_line[j].split()[0].rjust(3)
                tmp_line[j].split()[1].rjust(3)
                tmp_line[j].split()[2].rjust(3)
                tmp_line[j].split()[3].rjust(3)
            except IndexError: pass
            else:
                fout.writelines(tmp_line[j].split()[0].rjust(3) + tmp_line[j].split()[1].rjust(3))
                fout.writelines(tmp_line[j].split()[2].rjust(3) + tmp_line[j].split()[3].rjust(3))                
                stressx = float(tmp_line[j].split()[4])/nrm
                stressy = float(tmp_line[j].split()[5])/nrm
                fout.writelines(' %8.3f %8.3f \n' %(stressx, stressy))
