import biaxial_pp as bpp
import os
import glob

try:
    reload(bpp)
except NameError:
    pass


class worklevel:
    def __init__(self, workmax=0.1, sigxl='avg', sigyl='avg', epsxl='plavg', epsyl='plavg',
                 inc=0.1, filename=None, path=None, modulusx=None, modulusy=None, iplot=True):
        step = int(workmax/inc)
        current_work = 0
        if path == None: path = os.getcwd()
        pp = bpp.biaxial(filenae=filename, path=path, modulusx=modulusx, modulusy=modulusy)
        fout = open(path + '\\' + filename.split('.')[0] + '.wrk', 'w')
        fout.writelines('** Work vs Sig and theta contour\n')
        fout.writelines('**\n')
        fout.writelines('  six     sigy      work       theta\n')
        while True:
            current_work = current_work + inc
            if current_work > workmax: break
            sigx, sigy, theta = pp.yield_vector(sigxl=sigxl, sigyl=sigyl, epsxl=epsxl,
                                                epsyl=epsyl, w = current_work, color='k',
                                                alpha=0.5, marker='o', iplot=iplot)
            fout.write('%8.4f  %8.4f  %8.4f  %8.4f  \n' %(sigx, sigy, current_work, theta))
        fout.close()

class bpp_exe:
    """
    Returns and write maximum work, Eps, Sig available from the files
    """
    def __init__(self, files = None):
        bpps = []
        fout = open(os.getcwd()+'\\'+ 'worklevels.out', 'w')
        fout.writelines('** plastic maximum work and plastic strain contour **\n')
        fout.writelines('filename                 Sig_ratio  Eps_ratio    ')
        fout.writelines('Epsx        Epsy        Epsxpl        Epsypl        ')
        fout.writelines('Sigx          Sigy          max pl_work   \n')
        if files == None:
            files = glob.glob('*.csv')
        else : pass

        for i in range(len(files)):
            temp = bpp.biaxial(filename=files[i],path=os.getcwd(),
                               modulusx = 337585., modulusy = 297914.)
            maxplwork, maxEx, maxEy, maxExpl, maxEypl, maxSx, maxSy = temp.total_work()
            stressratio = []
            strainratio = []
            try:
                stressratio.append(temp.uni.stressratio[0])
            except IndexError:
                stressratio.append('0')
                stressratio.append('0')
                try:
                    strainratio.append(temp.uni.strainratio[0])
                except IndexError:
                    print 'no proper ratio was found'
                    raw_input('press enter ')
                else:
                    strainratio.append(temp.uni.strainratio[1])
            else:
                stressratio.append(temp.uni.stressratio[1])        
                strainratio.append('0')
                strainratio.append('0')
                
            fout.write(files[i][0:29].ljust(30))
            fout.write(str(stressratio[0]).ljust(3))
            fout.write(str(stressratio[1]).ljust(5))
            fout.write(str(strainratio[0]).ljust(3))
            fout.write(str(strainratio[1]).ljust(5))
            fout.write('{0:.5f}'.format(maxEx).rjust(10))
            fout.write('{0:.5f}'.format(maxEy).rjust(12))
            fout.write('{0:.5f}'.format(maxExpl).rjust(12))
            fout.write('{0:.5f}'.format(maxEypl).rjust(14))
            fout.write('{0:.5f}'.format(maxSx).rjust(14))
            fout.write('{0:.5f}'.format(maxSy).rjust(15))
            fout.write('{0:.4f}'.format(maxplwork).rjust(12)+'\n')

        del temp, maxplwork, maxEx, maxEy, maxExpl, maxEypl
        fout.close()
