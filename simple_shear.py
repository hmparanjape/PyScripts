"""
Simple shear post-process

Cyclic test is positively reproduced.
"""
import numpy as np
import glob
import matplotlib.pyplot as plt
from vpsc_param import fout

class ss():
    """
    Post-processor to Simple shear test from MTS
    """
    def __init__(self, wildcard=None, filename=None,
                 t=1.0, l=50., g=3., mode=0):
        """
        Arguments:
        wildcard = 'txt'
        t, l, g are thickness, length, gauge width of the sample. [mm]
        mode (0: mirror-flip)
        mode (1: mirror-flip + translate to the origin point)
        """
        self.mode = mode
        
        if wildcard==None and filename==None:
            print "One of wildcard and filename must be given"
            raise IOError
        elif wildcard!=None and filename!=None:
            print "wildcard and filename is complementary to each other!"
            print "Select only one!"
            raise IOError

        if wildcard!=None:
            files = glob.glob('*%s*'%wildcard); files.sort()
            pass

        else:files = [filename]
        self.files = files
        
        self.t = t #[mm]
        self.l = l #[mm]
        self.area = t * l /(10.**6) #[m^2]
        self.g = g  #[mm]
        
        print 'height:', l, '[mm]\n thickness:', t
        print '[mm]\n gauge: ', g,'[mm]'
        pass
    
    def __file__(self,filename=None):
        """
        File is loaded and returns the data

        Argument
           filename=None
        """
        if filename==None:
            print "A filename should be given!"
            raise IOError
        data = np.loadtxt(filename, skiprows=8, delimiter=',', dtype='float')
        data = data.transpose()
        data[0] = data[0] - data[0][0] #time trim
        data[3] = data[3] - data[3][0] #displacement trim
        return data
    
    def __flip__(self, load, ext):
        """
        If any load of the data is found to be negative,
        mirror-flip the flow curve post to that negative point (self.mode==0)
        or translate it to the original point (self.mode=1)
        
        This will be called several times on a cyclic test result

        Arguments:
          load  (load array)
          ext   (extension array)
        """
        # self.mode= 0 or 1 (0 continuous, 1 translate to the origin)


        ## Dict list data (deformation segment) ------------- ##
        ext_dict  = dict()
        load_dict = dict()
        rvs_points = list()
        ## -------------------------------------------------- ##
        if len(load)!=len(ext): raise IOError

        iseg = 0
        rvs_points.append(0)
        for i in range(len(load)):
            if load[i]<0:
                iseg = iseg + 1
                rvs_points.append(i)
                ## change of the load sign (crossed the border)
                try:
                    x0 =  ext[i-1]; x1 =  ext[i]
                    y0 = load[i-1]; y1 = load[i]
                    pass
                ## Probably due to wrong index for ext or load
                except:
                    raise IOError, 'unexpected index problem'
                else:
                    slope = (y1 - y0) / (x1 - x0)
                    ## mirror-flip
                    if self.mode==0:
                        point = - y0 / slope + x0
                        ext[i:len(ext)] = point + point - ext[i:len(ext)]
                        load[i:len(load)] = - load[i:len(load)]
                        pass
                    ## mirror-flip and translates
                    elif self.mode==1:
                        ## perform the flip! ------------------------- ##
                        ext[i:len(ext)] = - (ext[i:len(ext)] - ext[i])
                        load[i:len(load)] = - load[i:len(load)]
                        ## ------------------------------------------- ##
                        
                        ## dictionary type data for each segment ----- ##
                        seg_start_ind = rvs_points[len(rvs_points)-2]
                        seg_finish_ind = rvs_points[-1]
                        # print seg_start_ind
                        # print seg_finish_ind
                        # print 'id: %i %i'%(seg_start_ind, seg_finish_ind)
                        
                        dict_name = '%s_seg'%str(iseg).zfill(2)
                        ext_dict[dict_name] = ext[seg_start_ind:seg_finish_ind]
                        load_dict[dict_name] = load[seg_start_ind:seg_finish_ind]
                        ## ------------------------------------------- ##
                        pass
                    pass
                pass
            pass
        if self.mode==0: return load, ext
        elif self.mode==1: return load_dict, ext_dict
        else: raise IOError
            

    def __pp__(self, filename):
        """
        Flips the load upon extension or displacement
        and returns the flipped load and extension
        
        Argument:
          filename
        """
        data = self. __file__(filename=filename)
        ext = data[2].copy()
        load = data[1].copy()
        disp = data[3].copy()
    
        #Flip the reverse loading upon extension
        load, ext = self.__flip__(load=load, ext=ext)
        return load, ext

    def __plotnwrites__(self, ifig=None):
        """
        Plot the stress-strain curves and writes 
        the post-processed data to files.
        
        Argument
          ifig = None
        """

        fig = plt.figure(ifig,figsize=(9,5))
        ax = fig.add_axes((0.10, 0.15, 0.60, 0.8),
                          label='str_str')
        
        for i in self.files:
            #load and either extension or displacement
            try: load, ext = self.__pp__(filename=i)
            except: print 'Error'; return 0
            if len(load)!=len(ext): raise IOError

            if self.mode==0:
                ## shear stress strain calculation ---------------- ##
                stress = (load.copy()*(10.**3))/self.area #[Pa]
                stress = stress/(10.**6) #[Mpa]
                strain = ext.copy()/self.g
                ## ------------------------------------------------ ##
                try:
                    # Plots the results.
                    ax.plot(strain, stress, label=i,
                            marker='.', color='black', ls='None'
                            )
                except:
                    print "Error occured during plotting!"
                    return -1
                # Writes the results.
                else: fout('%s.pp'%i.split('.')[0],
                           strain, stress)
                
                pass
            elif self.mode==1:
                keys = load.keys()
                keys.sort()
                for j in range(len(keys)):
                    load_seg = np.array(load[keys[j]])
                    ext_seg = np.array(ext[keys[j]])
                    ## stress and strain calculation ##
                    stress = (load_seg * (10.**3))/self.area
                    stress = stress/(10.**6)
                    strain = ext_seg/self.g
                    ##
                    try:
                        if len(self.files)==1: label = keys[j]
                        else: label = '%s_%s'%(keys[j], self.files[i])
                        ax.plot(strain, stress, label=label)
                        pass
                    except: raise IOError, "Error occured during plotting"
                    pass
                pass
            else: raise IOError, 'Unexpected mode'
            pass

        ax.set_xlabel(r'$\gamma$',
                      dict(fontsize=20))
        ax.set_ylabel(r'$\tau$'+' [MPa]',
                      dict(fontsize=20))
        ax.legend(bbox_to_anchor=(1.05,1), loc=2)
        ax.set_xlim(0.,);ax.set_ylim(0.,)
        pass
        
    pass # end of the class ss

def main(mode, filename, wildcard, t, l, g, ishow):
    """
    command script

    Argument
       mode
       filename
       wildcard
       t, l, g
       ishow
    """
    import matplotlib.pyplot as plt
        
    if ishow==True: ifig=1
    elif ishow==False:ifig=False
    else: raise IOError, "Unexpected ifig"

    ## 
    MyClass = ss(
        wildcard = wildcard, filename = filename,
        t = t, l = l, g=g, mode = mode)

    ## plots and writes down
    MyClass.__plotnwrites__(ifig = ifig)
    if ishow==True:
        try: plt.show()
        except:
            print '\n** Warning) Plot is not visible **\n'
            pass
        pass
            
    pass

if __name__=='__main__':
    #import matplotlib.pyplo as plt
    import getopt, sys

    ## arguments ------------------------- ##
    try: opts, args = getopt.getopt(
        sys.argv[1:], 'm:i:w:l:t:g:s')
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
        pass
    ## ----------------------------------- ##

    ## default options -------------------------------------------- ##
    mode = 0          # mode (0: mirror, 1: mirror+translation)
    filename = None   # single filename to be post-processed
    wildcard = None   # wildcard for multiple number of files
    t = 1.0           # thickness of the sample
    l = 50           # sample length along the loading direction 
    g = 3.            # gauge length
    ishow = False     # plotting    
    ## ------------------------------------------------------------ ##
    

    for o, a in opts:
        if o in ('-m'): mode = int(a)
        elif o in ('-i'): filename = a
        elif o in ('-w'): wildcard = a
        elif o in ('-t'): t = float(a)
        elif o in ('-l'): l = float(l)
        elif o in ('-g'): g = float(g)
        elif o in ('-s'): ishow = True
        pass

    
    ## excute the main function ------------------- ##
    main(mode, filename, wildcard, t, l, g, ishow)
    ## -------------------------------------------- ##
    pass
