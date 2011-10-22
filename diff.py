"""
Differentiate 
"""

from scipy.interpolate import UnivariateSpline
from numpy import linspace
import glob


def __UnivariateSpline__(x, y, xdown, xup, s=1, spacing=100, k=3):
    """
    One-dimensional smoothing spline fit to a given set of data points.
    *implemented from scipy.interpolate package

    Returns xs, ys which are sampled

    Arguments:
    x, y : sequence data for x and y. 
           x must be increasing
    s=1  : positive smoothing factor
    xdown, xup : sampled x limit
    spacing : spacing to be used for sample
    """
    s = UnivariateSpline(x, y, s=s, k=k)
    xs = linspace(xdown, xup, spacing)
    ys = s(xs)
    return xs, ys


def diff(x, y, ind, n_bin = 100):
    """
    Receives x and y, then calculate the differenciate of them
    by calculate instantaneous linear fitting to the curve.
    Over the linear fitting data of given binsize are sampled.
    """
    if ind <= n_bin:
        return -1
    elif ind + n_bin >= len(x):
        return -1
    else:
        xi = x[ind-n_bin:ind+n_bin+1]
        yi = y[ind-n_bin:ind+n_bin+1]
        #print 'xi **'
        #print xi
        #print 'yi **'
        #print yi
        #raw_input()
        xs, ys = __UnivariateSpline__(x=xi, y=yi, xdown=xi[0], xup=xi[-1], k=1)
        return (ys[-1] - ys[0]) / (xs[-1]-xs[0])


filename = raw_input('Type the full name of the file with extension included >>')

try: file(filename, 'r')
except IOError:
    while True:
        filename = raw_input('Type the filename again (Type exit to exit)>>')
        if filename == 'exit': break
        try : file(filename,'r')
        except IOError:
            pass
        else:
            f_in = file(filename, 'r')
            break
else:
    f_in = file(filename , 'r')


#n_head = float(raw_input('Type the number of header lines in your file >>'))


lines = f_in.readlines()
x = []
y = []

for i in range(len(lines)):
    try:
        float(lines[i].split()[0])
    except:
        pass
    else:
        x.append(float(lines[i].split()[0]))
        y.append(float(lines[i].split()[1]))



















