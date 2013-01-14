"""
The COD section view module for RVE
Youngung Jeong, 2011-May 27~ Jun 10

(log)
  1. (2011-May-28)
  The projection at the boundary is not completed. 
  The phi1 and phi2 seem to work correcly.
  For phi projection, the normalization has to be reconsidered...
  I wonder if phi projection plane can be sample-symmetry-applied.

  2. (2011-May-29)
  The above conflict is now completed.

  3. (2011-Jun-10)
  Multiple section view (def cod)
"""

import math
import euler
euler = euler.euler
import numpy as np
import matplotlib.pyplot as plt

def spread(eul=None, w=None, nseg=0):
    """
    Given the the Euler angle set (including the volume),
    spread the current angle into nseg number of poles
    based on a distribution with the scatter angle is given as w

    eul = [ph1, ph, ph2, vf]
    """
    rst = []
    
    vf = eul[3] / nseg
    a = euler(ph = eul[0], th = eul[1], tm = eul[2], echo=False)
    for i in range(nseg):

        ## Generate a random rotation axis, about which the given grain rotates
        delta, phi = rot_axis()

        ## Calculates the spread (rotation) matrix b
        b = spread_matrix(
            delta, phi,
            gaussian(w))
        ## Rotates the rotation matrix a

        c = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    c[i,j] = c[i,j] + b[m,i] * a[m,j]
                    pass
                pass
            pass
        #c = np.dot(b.transpose(), a)

        p1, p, p2 = euler(a=c, echo=False)

        rst.append([p1, p, p2, vf])
        pass

    ## Returns the spread angle sets
    return rst

def rot_axis():
    """
    Random rotation axis generator
    """
    acos = math.acos
    delta = random.uniform(0.,1.)
    delta = acos(delta)
    phi = random.uniform(0.,2.*math.pi)
    return delta, phi   #radians...

def spread_matrix(delta, phi, w):
    """
    delta, phi, w = delta*math.pi/180. , phi*math.pi/180. , w*math.pi/180.
    Generates the rotation matrix
    """
    sin = math.sin
    cos = math.cos
    w = w * math.pi/180.
    d1 = sin(delta) * cos(phi)
    d2 = sin(delta) * sin(phi)
    d3 = cos(delta)
    cw = cos(w)
    sw = sin(w)
    p=[[None,None,None],[None,None,None],[None,None,None]]
    p[0][0] = ( 1 - d1**2 ) * cw + d1**2
    p[0][1] = d1 * d2 * ( 1 - cw ) + d3 * sw
    p[0][2] = d1 * d3 * ( 1 - cw ) - d2 * sw
    p[1][0] = d1 * d2 * ( 1 - cw ) - d3 * sw
    p[1][1] = ( 1 - d2**2 ) * cw + d2**2
    p[1][2] = d2 * d3 * ( 1 - cw ) + d1 * sw
    p[2][0] = d1 * d3 * ( 1 - cw ) + d2 * sw
    p[2][1] = d2 * d3 * ( 1 - cw ) - d1 * sw
    p[2][2] = ( 1 - d3**2 ) * cw + d3**2
    return np.array(p)

def gaussian(w0):
    """
    Returns the spread angle
    """
    gauss = random.gauss
    dp = gauss(mu=0., sigma=w0)
    return dp
                             

def main(phi1=None, phi2=None, phi=None,
         resolution=7.5, gr=None,
         nseg=1, ismooth=False, ireverse=False,
         ixlabel=True, iylabel=True, iclose=True,
         ifig=1, echo=True, x=0., y=0., levels=None):
    """
    A single COD section view
    (For multiple section view, refer to the def cod)
    
    Arguments:
    phi1, phi2, Phi (Bunge convention)
    resolution
    gr
    spread = None (an integer)
    nseg = 3
    """
    pi = math.pi
    rad = pi / 180.
    cos = math.cos
    sin = math.sin
    dang = resolution
    ngrain = len(gr)
    gr0 = gr.copy()
    gr = np.array(gr).transpose()
    wgt = gr[3].sum()
    if echo==True:
        print '***  number of grain: %i  ***'%ngrain
        print 'Total volume: %3.1f'%wgt
        pass
    
    # self.gr  # grains
    # gr = self.gr.transpose() #phi1, phi, phi2, volume fraction

    ## Boundary of the Euler space in which the COD is sampled. ##
    ## 90 is assumed to be common factor for the Euler space
    factor = 90 #common factor
    mx = np.array([max(gr[0]), max(gr[1]), max(gr[2])])
    mn = np.array([min(gr[0]), min(gr[1]), min(gr[2])])
    mx = np.array(map(round, mx/factor), dtype='int') * factor
    mn = np.array(map(round, mn/factor), dtype='int') * factor
    if echo==True:
        print 'max:  %3.1f %3.1f %3.1f'%(mx[0],mx[1],mx[2])
        print 'min:  %3.1f %3.1f %3.1f'%(mn[0],mn[1],mn[2])
        pass
    ## -------------------------------------------------------- ##

    ## spread poles ------------------------------------------- ##
    w = resolution/10. #spreading angle
    if nseg==1:pass
    else:
        ## Multiplication of grains ------------------- ##
        temp = []
        for i in range(ngrain):
            rst = spread(eul=gr0[i], w=w, nseg=nseg)
            count = 0
            for n in range(len(rst)):
                if rst[n][0]>=mn[0] and rst[n][0]<=mx[0] and rst[n][1]>=mn[1] and rst[n][1]<=mx[1] and rst[n][2]>=mn[2] and rst[n][2]<=mx[2]:
                    temp.append(rst[n])
                    pass
                else: count = count + 1
                pass
            if count!=0:
                # print 'count: %i'%count
                # print 'volume %5.3f'%(gr0[i][3]/nseg * count)
                temp.append([gr0[i][0], gr0[i][1], gr0[i][2], gr0[i][3]/nseg*count])
                pass
            pass
        gr = np.array(temp).transpose()
        ngrain = len(gr)
        pass
        ## -------------------------------------------- ##
    ## -------------------------------------------------------- ##

    ## -------------------------------------------------------- ##
    del_ph1 = (mx[0] - mn[0]) * rad
    del_ph  = cos(mx[1]* rad) - cos(mn[1] * rad)
    del_ph2 = (mx[2] - mn[2]) * rad

    FUL = abs(del_ph1 * del_ph * del_ph2)
    if echo==True : print 'FULL INTEGRATION VALUE: %3.1f * pi * pi'%(
        FUL/pi/pi)
    ## -------------------------------------------------------- ##

    ## x, y, z axis insertion & rdi calculation --------------- ##
    if None in (phi1, phi2, phi):
        if phi1!=None:
            if echo==True:print 'section at phi1 = %4.1f'%phi1
                
            x_axis = gr[2] #phi2
            y_axis = gr[1] #phi
            z_axis = gr[0]
            
            xmax = mx[2]; xmin = mn[2]
            ymax = mx[1]; ymin = mn[1]
            
            z = phi1
            rdi = del_ph2 * del_ph * (resolution * rad)
            pass

        elif phi!=None:
            if echo==True:print 'section at phi = %4.1f'%phi
            x_axis = gr[0] #phi1
            y_axis = gr[2] #phi2
            z_axis = gr[1]

            xmax = mx[0]; xmin = mn[0]
            ymax = mx[2]; ymin = mn[2]
            
            z = phi
            rdi = del_ph1 * del_ph2 * cos(
                (phi - resolution / 2.) * rad) - cos(
                (phi + resolution / 2.) * rad)
            pass
        elif phi2!=None:
            if echo==True:print 'section at phi2 = %4.1f'%phi2
            x_axis = gr[0] #phi1
            y_axis = gr[1] #phi
            z_axis = gr[2]
            
            xmax = mx[0]; xmin = mn[0]
            ymax = mx[1]; ymin = mn[1]
            
            z = phi2
            rdi = del_ph1 * del_ph * (resolution * rad)
            pass
        rdi = abs(rdi)
        rdi = rdi #/ FUL #= abs(del_ph1 * del_ph * del_ph2)
        if echo==True:print 'rdi: %f'%rdi
        pass
    else: raise IOError
    ## -------------------------------------------------------- ##


    ## The cubic cell maker (f) ------------------------------- ##
    dx = dang # grid angle along the x-axis
    dy = dang # grid angle along the y-axis
    
    mx = round(max(x_axis) / factor) * factor
    nx = round(min(x_axis) / factor) * factor
    my = round(max(y_axis) / factor) * factor
    ny = round(min(y_axis) / factor) * factor
    mz = round(max(z_axis) / factor) * factor
    nz = round(min(z_axis) / factor) * factor
    
    if echo==True:
        print 'x: %3.1f ~ %3.1f'%(nx, mx)
        print 'y: %3.1f ~ %3.1f'%(ny, my)
        print 'z: %3.1f ~ %3.1f'%(nz, mz)
        pass

    # if (z + resolution/2.) > mz: raise IOError
    # if (z - resolution/2.) < nz: raise IOError
    
    nnx = abs(mx - nx) / dang #nnx grid
    nny = abs(my - ny) / dang #nny grid
        
    # The cubic cell in the Euler space
    f = np.zeros((nnx, nny))
    ## -------------------------------------------------------- ##

    ## Assigns the volume fraction to each of the cubic celll - ##
    for n in range(ngrain):
        xvalue = x_axis[n]
        yvalue = y_axis[n]
        zvalue = z_axis[n]
        vf = gr0[n][3]

        # cell indicies
        ix = int(xvalue / resolution)
        iy = int(yvalue / resolution)

        ## if the projection plane is close to
        ## the Given Euler space boundary
        if z + resolution/2. > mz:
            zdiff = z + resolution/2. - mz
            if zvalue < nz + zdiff or zvalue > z - resolution / 2.:
                f[ix, iy] = f[ix, iy] + vf
                pass
            pass
        elif z - resolution/2. < nz:
            zdiff = nz - (z - resolution/2.)
            if zvalue > mz - zdiff or zvalue < z + resolution / 2.:
                f[ix, iy] = f[ix, iy] + vf
                pass                
            pass
        ##
        ## ---------------------------------- ##

        ## ---------------------------------- ##        
        else:
            if abs(zvalue - z) < resolution/2. :
                f[ix, iy] = f[ix, iy] + vf
                pass
            else: pass
            pass
        ## ----------------------------------##

        
        pass
    ## -------------------------------------------------------- ##


    ### Calculation of the cubic weight if were random ###
    ### Normalization ###
    
    
    if phi1!=None:
        for m in range(len(f)): #phi2
            for n in range(len(f[m])): #phi
                ## cube volume ----------------------------- ##
                dp1 = resolution * rad
                dp2 = resolution * rad
                current_ph = (n * resolution  + resolution / 2. )*rad                
                ph_up   = current_ph + rad * resolution * 0.5
                ph_down = current_ph - rad * resolution * 0.5
                dp = cos(ph_up) - cos(ph_down)
                cube_vol = abs(dp1 * dp * dp2)
                ## ----------------------------------------- ##
                
                if cube_vol==0: raise IOError

                f[m,n] = f[m,n] / cube_vol * FUL / wgt
                pass
            pass
        pass
    
    elif phi!=None:
        ## delta phi calculation
        
        ## cube volume ----------------------------- ##
        dp1 = resolution * rad
        dp2 = resolution * rad

        ##### phi is closer to the boudnary ----- #####
        if phi + resolution/2. > mz:
            diff = mz - phi
            
            ph_down = (phi - resolution/2.) * rad
            ph_up = mz * rad
            dp = cos(ph_up) - cos(ph_down)

            ph_down = nz * rad
            ph_up = (nz + (resolution/2. - diff)) * rad
            dp = dp + cos(ph_up) - cos(ph_down)
            pass
        elif phi - resolution/2. < nz:
            diff = phi - nz

            ph_down = (mz - (resolution/2. - diff)) * rad
            ph_up = mz * rad
            dp = cos(ph_up) - cos(ph_down)

            ph_down = nz * rad
            ph_up = (phi + resolution/2.) * rad
            dp = dp + cos(ph_up) - cos(ph_down)
            pass
        else:
            ph_up = (phi + resolution/2.) * rad
            ph_down = (phi - resolution/2.) * rad
            dp = cos(ph_up) - cos(ph_down)
            pass
        ##### ----------------------------------- ####
        cube_vol = abs(dp1 * dp2 * dp)
        ## ----------------------------------------- ##
        for m in range(len(f)):
            for n in range(len(f[m])):
                f[m,n] = f[m,n] / cube_vol * FUL / wgt
                pass
            pass        
        pass
    
    elif phi2!=None:
        for m in range(len(f)):
            for n in range(len(f[m])):
                ## cube volume ----------------------------- ##
                dp1 = resolution * rad
                dp2 = resolution * rad
                current_ph = (n + 0.5) * rad * resolution
                ph_up = current_ph + rad * resolution * 0.5
                ph_down = current_ph - rad * resolution * 0.5
                dp = cos(ph_up) - cos(ph_down)
                cube_vol = abs(dp1 * dp2 * dp)
                ## ----------------------------------------- ##

                f[m,n] = f[m,n] / cube_vol * FUL / wgt
                pass
            pass        
        pass
    ### ------------------------------------------ ###

    # ## mesh grid for the cell f ---------------- ##
    # x_lin = np.arange(xmin + resolution/2., xmax, dx)
    # y_lin = np.arange(ymin + resolution/2., ymax, dy)
    # print 'x_axis', x_lin
    # print 'y_axis', y_lin

    # x, y = np.meshgrid(y=x_lin, x=y_lin)
    # print x.shape
    # print y.shape
    # print f.shape
    # ## ----------------------------------------- ##
    
    # fig = plt.figure(1)
    # ax = fig.add_subplot(211)
    # cnt = ax.contourf(x, y, f ) #, cmap=plt.cm.cmap_d['gray_r'])
    # level = cnt._levels
    # print 'levels: ', level

    # ## making nodes from cell ------------------ ##
    nodes = np.zeros((nnx+1, nny+1))
    # ## ----------------------------------------- ##

    ## assuming x-axis and y-axis boundaries are the mirror-axes
    for ix in range(int(nnx+1)):
        for iy in range(int(nny+1)):
            if ix==0:
                if iy==0:
                    #at the first corner
                    nodes[ix, iy] = f[ix, iy] * 4. / 4.
                    pass
                elif iy==nny:
                    nodes[ix, iy] = f[ix, iy-1] * 4. / 4.
                    pass
                else:
                    #at along the x=0 edge
                    temp = f[ix, iy-1] * 2.
                    temp = temp + f[ix, iy] * 2.
                    nodes[ix, iy] = temp / 4.
                    pass
                pass
            elif ix > 0 and ix < nnx:
                ## lower (iy=0) boundary
                if iy==0:
                    temp = f[ix-1 ,iy]*2.
                    temp = temp + f[ix, iy] * 2.
                    nodes[ix, iy] = temp/4.
                    pass
                ## upper (iy=ny) boundary
                elif iy==nny:
                    temp = f[ix-1, iy-1] * 2.
                    temp = temp + f[ix, iy-1] * 2.
                    nodes[ix, iy] = temp /4.
                    pass
                ## the rest
                else:
                    temp = f[ix-1, iy-1] + f[ix, iy-1]
                    temp = temp + f[ix-1, iy] + f[ix, iy]
                    nodes[ix,iy] = temp / 4.
                    pass
                pass
            elif ix==nnx:
                if iy==0:
                    nodes[ix,iy] = f[ix-1, iy] * 4. / 4.
                    pass
                elif iy==nny:
                    nodes[ix,iy] = f[ix-1, iy-1] * 4. / 4.
                    pass
                else:
                    temp = (f[ix-1,iy-1] + f[ix-1,iy]) * 2.
                    nodes[ix, iy] = temp / 4.
                    pass
                pass
            pass
        pass
    ## ------------------------------------------------------ ##

    ## average along the phi axis -------------- ##
    if phi==None:
        xlabel = r'$\Phi$'
        if ismooth==True:
            temp = nodes.transpose()[0].mean()
            nodes = nodes.transpose()
            nodes[0] = temp
            nodes = nodes.transpose()
            pass
        if phi1!=None: ylabel=r'$\phi_2$'
        elif phi2!=None: ylabel=r'$\phi_1$'
        pass
    if phi!=None:
        xlabel=r'$\phi_2$'
        ylabel=r'$\phi_1$'
        pass
    ## ----------------------------------------- ##


    ## plotting the section based on nodes ----- ##
    nplot=ifig
    fig = plt.figure(nplot,figsize=(6,5))
    tiny = 0.00001
    

    x_lin = np.arange(xmin, xmax + tiny, dx)
    y_lin = np.arange(ymin, ymax + tiny, dy)
    x, y = np.meshgrid(y = x_lin, x = y_lin)
    axt = fig.add_axes(
        (0.15, 0.05, 0.75, 0.80),
        aspect='equal'
        )#add_subplot(111,)
    ax = axt.twiny()

    if levels==None: cnt = ax.contour(
        x,y, nodes, cmap = plt.cm.cmap_d['gray_r'])
    elif levels!=None: cnt = ax.contour(
        x,y, nodes, levels, cmap = plt.cm.cmap_d['gray_r'])
    
    if ixlabel:ax.set_xlabel(xlabel, dict(fontsize=20))
    if iylabel:axt.set_ylabel(ylabel, dict(fontsize=20))
    level = cnt._levels
    
    ax.legend()

    ## x and y limit, and level plotting ------- ##
    ax.set_xlim(ymin, ymax); 
    if ireverse:
        if ixlabel==False:
            ax.set_yticks(np.arange(xmax, xmin-0.001,-15.))
            ax.set_xticks(np.arange(ymin, ymax+0.001, 15.))
            for tick in ax.xaxis.get_major_ticks():
                tick.label1On=False;tick.label2On=False
                pass
            pass
        if iylabel==False:
            for tick in axt.yaxis.get_major_ticks():
                tick.label1On=False;tick.label2On=False
                pass
            pass
        ax.set_ylim(xmax, xmin)
        for tick in axt.xaxis.get_major_ticks():
            tick.label1On=False
            tick.label2On=False
            pass
        if echo==True: print 'xmax: %f xmin: %f'%(xmax,xmin)
        pass
    else:
        ax.set_ylim(xmin, xmax)
        ax.set_yticks(np.arange(xmin, xmax+0.001, 15.))
        ax.set_xticks(np.arange(ymin, ymax+0.001, 15.))
        pass
    
    
    plt.clabel(cnt, inline=False, fontsize=9., fmt = '%4.2f', inline_spacing=1. )

    # tcolors = cnt.tcolors
    # for i in range(len(tcolors)):
    #     cc = tcolors[i][0][0:3]
    
    ## ----------------------------------------- ##
    if echo==True: print 'levels: ', level
    ## ------------------------------------------ ##
    
    return x, y, nodes, ax, axt, level


def cod(phi1=None, phi2=None, phi=None,
        resolution=7.5, gr=None,
        nseg=1,ismooth=False, nr=None, ifig=10, niso=7,
        cmap = None, figsize=8,
        ):
    """
    Plotting the sections along a particular axis among the
    three Euler axes.

    Arguments
     phi1=None, phi2=None, phi=None, 
     resolution = 7.5 (cell resolution in the Euler angle)
     gr = None  (polycrystal aggregate)
     nseg = 1 (spreading the pole: As of Jun 10, still under development)
     ismooth= False (Enforce the same intensity along Phi=0 axis)
     nr=None (incremental angle for the cross section, in degree)
     ifig=10   
     niso = 7 (number of iso-level contour)
     cmap = None
    """
    nplot = ifig
    pi = math.pi
    rad = pi / 180.
    cos = math.cos
    sin = math.sin
    ##grid (l,m,n)
    dang = resolution
    ngrain = len(gr)
    gr0 = gr.copy()
    gr = np.array(gr).transpose()
    wgt = gr[3].sum()
    print 'Total volume: %3.1f'%wgt
    
    ## Boundary of the Euler space in which the COD is sampled. ##
    ## 90 is assumed to be common factor for the Euler space
    factor = 90 #common factor
    mx = np.array([max(gr[0]), max(gr[1]), max(gr[2])])
    mn = np.array([min(gr[0]), min(gr[1]), min(gr[2])])
    mx = np.array(map(round, mx/factor), dtype='int') * factor
    mn = np.array(map(round, mn/factor), dtype='int') * factor
    print 'max:  %3.1f %3.1f %3.1f'%(mx[0],mx[1],mx[2])
    print 'min:  %3.1f %3.1f %3.1f'%(mn[0],mn[1],mn[2])

    ## -------------------------------------------------------- ##
    del_ph1 = (mx[0] - mn[0]) * rad
    del_ph  = cos(mx[1]* rad) - cos(mn[1] * rad)
    del_ph2 = (mx[2] - mn[2]) * rad

    if nr==None: nr = resolution
    tiny = 0.00001    
    ## x, y, z axis insertion & rdi calculation --------------- ##
    if None in (phi1, phi2, phi):
        if phi1!=None:
            print 'section at phi1 = %4.1f'%phi1
            x_axis = gr[2] #phi2
            y_axis = gr[1] #phi
            z_axis = gr[0]
            
            xmax = mx[2]; xmin = mn[2]
            ymax = mx[1]; ymin = mn[1]
            zmax = mx[0]; zmin = mn[0]

            zlin = np.arange(zmin, zmax + tiny, nr)
            phi1 = zlin
            z = phi1
            rdi = del_ph2 * del_ph * (resolution * rad)
            pass

        elif phi!=None:
            print 'section at phi = %4.1f'%phi
            x_axis = gr[0] #phi1
            y_axis = gr[2] #phi2
            z_axis = gr[1]

            xmax = mx[0]; xmin = mn[0]
            ymax = mx[2]; ymin = mn[2]
            zmax = mx[1]; zmin = mn[1]
            zlin = np.arange(zmin, zmax + tiny, nr)
            
            z = phi
            rdi = del_ph1 * del_ph2 * cos(
                (phi - resolution / 2.) * rad) - cos(
                (phi + resolution / 2.) * rad)
            pass
        
        elif phi2!=None:
            print 'section at phi2 = %4.1f'%phi2
            x_axis = gr[0] #phi1
            y_axis = gr[1] #phi
            z_axis = gr[2]
            
            xmax = mx[0]; xmin = mn[0]
            ymax = mx[1]; ymin = mn[1]
            zmax = mx[2]; zmin = mn[2]
            zlin = np.arange(zmin, zmax + tiny, nr)
            
            z = phi2
            rdi = del_ph1 * del_ph * (resolution * rad)
            pass
        
        else:
            print phi1, phi2, phi
            raise IOError
        
        rdi = abs(rdi)
        rdi = rdi #/ FUL #= abs(del_ph1 * del_ph * del_ph2)
        print 'rdi: %f'%rdi
        pass
    ## -------------------------------------------------------- ##
    

    ## -------------------------------------------------------- ##
    ## find relevant levels
    mxlev = 0.
    mnlev = 1. # level 1 amounts to random texture's intensity
    for v in zlin:
        print 'value = ' , v
        if phi!=None: x, y, nodes, ax, axt, level = main(
            phi1=phi1, phi2=phi2, phi=v,
            resolution=resolution, gr=gr0,
            nseg=nseg, ismooth=ismooth, ifig=nplot,
            ireverse=True,echo=False)
        elif phi1!=None: x, y, nodes, ax, axt, level = main(
            phi1=v, phi2=phi2, phi=phi,
            resolution=resolution, gr=gr0,
            nseg=nseg, ismooth=ismooth, ifig=nplot,
            ireverse=True,echo=False)
            
        elif phi2!=None:x, y, nodes, ax, axt, level = main(
            phi1=phi1, phi2=v, phi=phi,
            resolution=resolution, gr=gr0,
            nseg=nseg, ismooth=ismooth, ifig=nplot,
            ireverse=True,echo=False)
            
        else: raise IOError
        if mxlev < max(level): mxlev = max(level)
        if mnlev > min(level): mnlev = min(level)
        pass
    plt.close(nplot) # close the plot
    ## -------------------------------------------------------- ##

    ## -------------------------------------------------------- ##
    axes,axest = [],[]    
    tlevel = np.linspace(mnlev, mxlev + tiny, niso)

    xmaster = [ ]
    ymaster = [ ]
    nmaster = [ ]
    
    for v in zlin:
        print 'value = ', v
        if phi!=None:
            lang = r'\Phi'
            x, y, nodes, ax, axt, level = main(
                phi1=phi1, phi2=phi2, phi=v,
                resolution=resolution, gr=gr0,
                nseg=nseg, ismooth=ismooth, ifig=nplot,
                ireverse=True, echo=False, levels=tlevel)
            pass
        elif phi1!=None:
            lang = r'\phi_1'
            x, y, nodes, ax, axt, level = main(
                phi1=v, phi2=phi2, phi=phi,
                resolution=resolution, gr=gr0,
                nseg=nseg, ismooth=ismooth, ifig=nplot,
                ireverse=True,echo=False, levels=tlevel)
            pass
        elif phi2!=None:
            lang = r'\phi_2'
            x, y, nodes, ax, axt, level = main(
                phi1=phi1, phi2=v, phi=phi,
                resolution=resolution, gr=gr0,
                nseg=nseg, ismooth=ismooth, ifig=nplot,
                ireverse=True,echo=False, levels=tlevel)
            pass
        else: raise IOError
        if mxlev < max(level): mxleve = max(level)
        xmaster.append(x)
        ymaster.append(y)
        nmaster.append(nodes)
        pass

    ## -- xlabel and ylabel assignment
    if phi==None:
        xlabel = r'$\Phi$'
        if phi1!=None: ylabel=r'$\phi_2$'
        elif phi2!=None: ylabel=r'$\phi_1$'
        pass
    if phi!=None:
        xlabel=r'$\phi_2$'
        ylabel=r'$\phi_1$'
        pass
    ## -------------------------------- ##
    
    #plt.close(nplot) # close the plot
    fig = plt.figure(nplot)
    fig.clf()
    fig.set_figheight(figsize); fig.set_figwidth(figsize)

    ## -------------------------------------------------------- ##

    ## -------------------------------------------------------- ##
    ## allocate relevant section for all COD Section in the plt figure
    nax = len(zlin)
    nxax = int(round(math.sqrt(nax))) + 1

    w = 0.99; h = 0.99; l = 0.99
    xoff = 0.1; yoff = 0.1
    w = w - xoff; h = h - yoff
    lx, ly = w / nxax, h / nxax # local cell's size
    ox, oy = xoff, 1 - yoff# - ly #the first local cell's origin
    off = 0.005
    dl = off; db = off; dw = 1 - off; dh = 1 - off# fraction in a local cell
    dl, db, dw, dh = dl * lx, db * ly, dw * lx, dh * ly # the first local db, dl, dw, dh

    pos = []; oox = ox; ooy = oy
    for i in range(len(zlin)):
        if i%nxax==0:
            oox = ox; ooy = ooy - ly
            pass
        pos.append([oox + dl, ooy + db, dw, dh])
        oox = oox + lx
        pass
    
    #fig = plt.figure(nplot, figsize=(8, 8))
    axes = []; LEV = []; COL = [] # master level and col through out the section

    
    for i in range(len(zlin)):
        print '%i th axis'%(i + 1)
        # Initiate the axis and allocate a region in the plt.figure object canvas.
        if i==len(zlin)-1:
            ## --- an additional plt axis object for plotting iso-levels
            ax0 = fig.add_axes((
                    pos[i][0], pos[i][1],
                    pos[i][2], pos[i][3]*0.9))
                # (xoff, (1 - yoff + db + (ly - dh)) * (
                #         1. + 0.1), dw * 0.2, dh * 0.8 ))#(xoff,h+ly-ddw,dh))
            ax0.set_axis_off()
            ## -------------------------------------------------------- ##
            pass
        else:
            axes.append(fig.add_axes(pos[i]))
            if cmap==None:
                cnt = axes[i].contour(
                    xmaster[i], ymaster[i],
                    nmaster[i], tlevel)
                pass

            elif cmap!=None:
                cnt = axes[i].contour(
                    xmaster[i], ymaster[i],
                    nmaster[i], tlevel,
                    cmap=plt.cm.cmap_d[cmap])
                pass
            for s in range(len(cnt._levels)):
                if not(cnt._levels[s] in LEV):
                    # save levels and colors for latter use.
                    LEV.append(cnt._levels[s])   
                    COL.append(cnt.tcolors[s][0])
                    pass
                pass
            axes[i].set_ylim(xmax, xmin) # impose the tight margin to the plt axis canvas.
            ax = axes[i]
            ## -- hide labels and ticks of the cod sections to facilitate the visualization.
            for tick in ax.xaxis.get_major_ticks():
                tick.label1On=False
                tick.tick1On=False;tick.tick2On=False #tick removal
                pass
            for tick in ax.yaxis.get_major_ticks():
                tick.label1On=False
                tick.tick1On=False;tick.tick2On=False #tick removal
                pass
            ## ---

            if i==0:
                axes[0].set_ylabel(ylabel,dict(fontsize=15))
                at = axes[0].twiny() # upper axis
                at.set_xlabel(xlabel, dict(fontsize=15))
                at.set_xticks([ymin, ymax])
                for tick in axes[0].xaxis.get_major_ticks():
                    tick.label1On=False;tick.tick1On=False
                    pass
                for tick in axes[0].yaxis.get_major_ticks():
                    tick.tick1On=True;tick.label1On=True;
                    pass
                axes[0].set_yticks([xmax, xmin])
                pass
            pass
        pass
    
    lev_temp = LEV[::] ## copy to dummy lev_temp
    LEV.sort() ## sorting the level according to their intensities

    ## -- plotting the levels and write the relevant level string next to each of
    ## -- iso-level segment.
    for i in range(len(LEV)):
        for s in range(len(COL)):
            if LEV[i] == lev_temp[s]:
                c = COL[s]
                break
            pass
        ax0.plot([1.28, 1.35], [1.-i*0.2, 1.-i*0.2],color=c)
        ax0.text(x=1.37, y= 1. - i*0.2 - 0.05,
                 s='%3.2f'%(LEV[i]), fontsize=8.)
        pass
    ax0.text(x=1.60, y = 1. - i * 0.2 - 0.05,
             s=r'$\Delta%s:$ %3.2f'%(lang, nr),
             fontsize=15.)
    ax0.set_xlim(1.2, 1.8)
    ## -- 
    pass  # end of the current def cod


## command prompt
if __name__ == '__main__':
    import getopt, sys
    import matplotlib.pyplot as plt
    ## arguments ---------------------------------------- ##
    try: opts, args = getopt.getopt(
        sys.argv[1:],'p:v:i:r:')
    except getopt.GEtoptError, err:
        print str(err)
        sys.exit(2)
        pass
    ## -------------------------------------------------- ##
    
    ## default options
    phi1 = None; phi = None; phi2 = None
    inputfile = None
    resolution = 7.5
    which_section = 1

    for o, a in opts:
        if o in ('-p'): which_section = int(a)
        elif o in ('-v'): value = float(a)
        elif o in ('-i'): inputfile = a
        elif o in ('-r'): resolution = float(a)
        else: assert False, 'unhandled option'
        pass

    ## projection plane 
    if which_section==0: phi = value
    elif which_section==1: phi1 = value
    elif which_section==2: phi2 = value

    ## RVE input (Following the VPSC texture file format)
    gr = np.loadtxt(inputfile, skiprows=4)
    
    x,y, nodes, ax, axt, level = main(
        phi1=phi1, phi2=phi2, phi=phi,
        resolution=resolution, gr=gr,
        nseg=1, ismooth=False)
    
    plt.show()
    pass
