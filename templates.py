"""
Template modules in which a passed matplotlib axis
is decorated with predefined axis label


#1 Stress-strain curve 

"""

def str_str(axis):
    """
    """
    axis.set_xlabel(r'$\varepsilon$', dict(fontsize=15))
    axis.set_ylabel(r'$\sigma$', dict(fontsize=15))
    pass

def str_work(axis):
    """
    """
    axis.set_xlabel(r'$\work^{pl}$', dict(fontsize=15))
    axis.set_ylabel(r'$\sigma$', dict(fontsize=15))    

def ys(axis):
    """
    """
    axis.set_xlabel(r'$\sigma_{RD}$')
    axis.set_ylabel(r'$\sigma_{TD}$')
    pass

def ys3d(axis):
    """
    """
    axis.set_xlabel(r'$\sigma_{RD}$', dict(fontsize=15))
    axis.set_ylabel(r'$\sigma_{TD}$', dict(fontsize=15))
    axis.set_zlabel(r'$\tau$', dict(fontsize=15))
    pass

def r_str(axis):
    """    
    """
    axis.set_xlabel('$\varepsilon_{axial}$',dict(fontsize=15))    
    axis.set_ylabel('r value', dict(fontsize=15))
    
    pass


def str_str_r():
    l = 0.15  #left blank
    b = 0.15  #bottom blank
    w = 0.75  #width
    h = 0.80  #height 

    fig_width=7
    fig_height=5
    fig_rel_size= 0.9 
    
    fig_width, fig_height = np.array([fig_width,fig_height]) * fig_rel_size


    markers = ['o','^','d','p','h','*','+']
        
    if iplot==True:
        #fig = plt.figure(ifig, [fig_width *len(self.modes),fig_height])
        fig = plt.figure(ifig, [fig_width*2.,fig_height])  
        # each figure's length and width
        el = l/2.
        ew = w/2.
        fullscale = 1.0
        incr = fullscale / len(self.modes)

        ax0 = fig.add_axes((el,      b,ew,h))
        ax1 = fig.add_axes((el+ incr,b,ew,h))
                
        # for i in range(len(self.modes)):
        #     current_left = el + incr * i
        #     current_b    = b
        #     current_wid  = ew #+ incr * i                
        #     current_h    = h
        #     axes.append(
        #         fig.add_axes(
        #             (current_left, current_b,
        #              current_wid, current_h),label='ax #%i'%i
        #             )
        #         )
        #     pass
        pass
            
