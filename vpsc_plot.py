"""
VPSC with instant plots
"""
import pylab
import os
class run:
    """
    Running the code
    """
    def __init__(self, filename='vpsc7.exe'):
        ios = os.system(filename)
        if ios !=0:
            print 'Running ',filename, 'failed'
            raw_input('press button to run a.exe')
            os.system('a.exe')
            
    def plot(self, plt='str_str', ix = 0, iy = 1):
        try:
            open(os.getcwd()+'\\'+plt + '.out')
        except IOError:
            print ' No such a file found'
            return -1
        else:
            f = open(os.getcwd()+'\\'+plt + '.out')
        source = f.readlines()
        self.source=source
        y, x = [], []
        for i in range(len(source)):
            try: map(float, source[i].split())
            except ValueError:
                pass
            else:
                current_line = map(float, source[i].split())
                y.append(current_line[iy])
                x.append(current_line[ix])
            pylab.plot(x, y, 'o')
            
    def show(self):
        try:
            pylab.show()
        except AttributeError:
            print 'Error: Probably any plot request was successful so far'
            print 'Try to see if any plot was successful'
            return -1
        pylab.show()
        
        
                    
                
        
