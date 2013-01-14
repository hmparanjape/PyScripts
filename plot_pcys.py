"""
Plotting PCYS 2D projection
"""
import matplotlib.pyplot as plt

def data(filename='pcys_pl.out'):
    f = file(filename)
    source = f.read()
    lines=source.split('\n')

    s11 = []
    s22 = []
    s12 = []
    iline = 0
    iz = -1
    while True:
        if len(lines[iline].split())<2: break
        try:
            float(lines[iline].split()[0])
        except:   #corresponding to header lines, i.e. starting of a block
                  #, starting of another height level along S12
            iz = iz + 1
            s11.append([])
            s22.append([])
            s12.append([])        
            pass
        else:
            col = map(float, lines[iline].split())
            #print col
            #raw_input()
            s11[iz].append(col[0])
            s22[iz].append(col[1])
            s12[iz].append(col[2])

        iline = iline + 1
    return s11, s22, s12

def plot_2D(x,y, marker=None, fig_id = None):
    if fig_id ==None: fig=plt.figure()
    else: fig=plt.figure(fig_id)
    ax=fig.add_subplot(111)
    for i in range(len(x)):
        if marker==None: ax.plot(x[i],y[i])
        else: ax.plot(x[i],y[i], marker)
    #plt.show()
