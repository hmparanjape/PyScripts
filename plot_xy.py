try: reload(plt)
except: import matplotlib.pyplot as plt

class pl:
    def __init__(self, x, y, ifig, ind_plot=111, label='', loc='upper right',
                 xlabel='', ylabel='', title='', ls='-', marker=None,
                 ilegend=True, alpha=1., markersize=12, ylim=None, xlim=None,
                 color=None):
        figure = plt.figure(ifig)
        ax = figure.add_subplot(ind_plot)
        ax.plot(x, y, label=label, ls=ls, marker=marker,
                alpha=alpha, markersize=markersize,
                color=color)
        if len(xlabel)!=0: ax.set_xlabel(xlabel)
        if len(ylabel)!=0: ax.set_ylabel(ylabel)
        if ilegend==True: plt.legend(loc=loc)
        if len(title)>0: ax.set_title(label=title)
        plt.draw()
        if ylim!=None :
            plt.ylim(ylim[0],ylim[-1])
        if xlim!=None :
            plt.xlim(xlim[0],xlim[-1])
