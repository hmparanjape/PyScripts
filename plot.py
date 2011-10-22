"""
This python script is for interactive plotting of the Lankford output file(s) of
the VPSC code. 

@author : youngung 

"""
import matplotlib.pyplot as plt
import os

def is_number(st):
    try:
        float(st)
        return True
    except ValueError:
        return False

class plot:
    """
    class plot:
    
    Explanation upon arguments
    
    path - The folder path including output file (default = c:\Python26\myWorkp
    lace\VPSC\\)
    plotid - output plot identification (default = ['LANK'] )
    suptitle - title of the figure.     (default = 'VPSC OUTPUT')
    """
    def __init__(self,
                 path='/mnt/hgfs/10_19_simple_shear/VPSC', 
                 plotid=['LANK'], 
                 suptitle='VPSC OUTPUT'):
        plt.ion()
        self.fig=[]
        self.fig.append(plt.figure(num=1,figsize=(8,8), frameon=False))
        ifig=0 # id for figure
        self.fig[ifig].suptitle('VPSC output', fontsize=14,fontweight='bold')
        plt.draw()
        #for i in range(len(plotid)):
        #    self.plot(plotid[i],path,ifig)
        self.plot(plotid,path,ifig)
        plt.draw()
        
    def plot(self, label,path,ifig):
        if   label[0]=='LANK'   :
            print 'yes'
            self.lankford(path,ifig)
        elif label=='STR_STR':
            self.str_str(path,ifig,xcol='E11',ycol='SCAU11')
            pass
        elif label=='PCYS'   :
            pass
        """ trial """
        self.str_str(path,ifig,xcol='E11',ycol='SCAU11')
        self.pcys   (path,ifig,xcol='S1' ,ycol='S2'    )

    def pcys(self,path,ifig,xcol='S1',ycol='S2'):
        #self.files.append(self.filedetection(path,'PCYS'))
        
        self.files=self.filedetection(path,'PCYS')
        #print self.files
        x=[]
        y=[]
        for i in range(len(self.files)):
            x.append([])
            y.append([])
            count=-1
            f=open(path+self.files[i])
            while True:
                temp=[]
                s=f.readline()
                col_label=s.split()
                print s.split()[0:4]
                templine=s.split()
                for j in range(len(templine)):
                    if templine[j]==xcol:
                        xlb=j
                    elif templine[j]==ycol:
                        ylb=j
                    else:
                        pass
                
                if len(s)<2:
                    print 'an EOF is reached'
                    break
                else:
                    if is_number(s.split()[0]):
                        #print i
                        x[i][count].append(float(s.split()[xlb]))
                        y[i][count].append(float(s.split()[ylb]))
                    else:
                        print s.split()
                        x[i].append([])
                        y[i].append([])
                        count=count+1
                        pass
            f.close()
        ax_pcys=self.fig[ifig].add_subplot(2,2,3)
        ax_pcys.set_title('Polycrystal Yield Surface')
        ax_pcys.set_xlabel(xcol)
        ax_pcys.set_ylabel(ycol)
        ax_pcys.set_aspect('auto','datalim')
        line=[]
        
        for i in range(len(x)):
            for j in range(len(x[i])):
                line.append(ax_pcys.plot(x[i][j],y[i][j],label='test-'+str(i)))
        

        for label in ax_pcys.xaxis.get_ticklabels():
            label.set_fontsize(8)
            label.set_rotation(90)

        plt.draw()
        
        
    def str_str(self,path,ifig,xcol='E11', ycol='SCAU11',y_col='E33'):
        self.files=self.filedetection(path,'STR_STR')
        print self.files
        x =[]
        y =[]
        y_=[]
        for i in range(len(self.files)):
            x.append([])
            y.append([])
            y_.append([])
            f=open(path+self.files[i])
            while True:
                temp=[]
                s=f.readline()
                col_label=s.split()     #assign column labels to col_labels for possible uses.
                print s.split()[0:2]
                templine=s.split()
                for j in range(len(templine)):
                    if templine[j]==xcol:
                        xlb=j
                    elif templine[j]==ycol:
                        ylb=j
                    elif templine[j]==y_col:
                        y_lb=j
                    else:
                        pass
                        
                if len(s)<2:
                    print 'an EOF is reached'
                    break
                else:
                    if is_number(s.split()[0]):
                        #print i
                        x[i].append(float(s.split()[xlb]))
                        y[i].append(float(s.split()[ylb]))
                        if float(s.split()[xlb])== 0 :
                             y_[i].append(0)
                        else:
                             y_[i].append( -(float(s.split()[xlb]) + float(s.split()[y_lb]))/float(s.split()[y_lb])  )
                    else:
                        pass
            f.close()
            
        ax_str_str=self.fig[ifig].add_subplot(2,2,2)
        ax_str_str.set_title('Stress Strain curve')
        ax_str_str.set_xlabel(xcol)
        ax_str_str.set_ylabel(ycol)
        #ax_str_str.set_aspect('auto','datalim')
        line=[]
        for i in range(len(x)):
            line.append(ax_str_str.plot(x[i],y[i],'o',label='test-'+str(i)))

        for label in ax_str_str.xaxis.get_ticklabels():
            label.set_fontsize(7)
            label.set_rotation(90)


	ax_R = self.fig[ifig].add_subplot(2,2,4)
	ax_R.set_title('R-value evolution ')
        ax_R.set_xlabel('Strain')
        ax_R.set_ylabel('R-value')
        ax_R.set_aspect('auto','datalim')
        line=[]
        for i in range(len(x)):
            line.append(ax_R.plot(x[i],y_[i],'o',label='test-'+str(i)))

        for label in ax_R.xaxis.get_ticklabels():
            label.set_fontsize(7)
            label.set_fontsize(7)

        plt.draw()
        
        
    def lankford(self, path, ifig):
        self.files=self.filedetection(path,'LANK')
        x=[]
        y=[]
        for i in range(len(self.files)):
            x.append([])
            y.append([])
            f=open(path+self.files[i])
            while True:
                temp=[]
                s=f.readline()
                print s.split()[0:4],' ... '
                if len(s)<2:
                    print 'an EOF is reached'
                    break
                else:
                    if is_number(s.split()[0]):
                        x[i].append(float(s.split()[0]))
                        y[i].append(float(s.split()[2]))
            f.close()

        ax_lnkf=self.fig[ifig].add_subplot(2,2,1)
        ax_lnkf.set_title('Lankford probing')
        ax_lnkf.set_xlabel('Angle from RD')
        ax_lnkf.set_ylabel('R')
        plt.axis([0,90.,0,3.])
        plt.xticks([0.,30.,60.,90.])
        ax_lnkf.set_aspect('auto','datalim')
        line=[]
        for i in range(len(x)):
            line.append(ax_lnkf.plot(x[i],y[i],'o',label='test-'+str(i)))

        for label in ax_lnkf.xaxis.get_ticklabels():
            label.set_fontsize(8)
            label.set_rotation(90)

        plt.draw()

    def filedetection(self, path, label):
        tempfiles=os.listdir(path)
        #print tempfiles
        files=[]
        for i in range(len(tempfiles)):
            #if tempfiles[i].split('.')[0:len(label)]==label:
            #print tempfiles[i].split('.')[0][0:5]
            if tempfiles[i].split('.')[0][0:len(label)]==label:
                files.append(tempfiles[i])
            else:
                pass
        print 'Being returned files=', files
        return files

    def close(self):
        plt.close()
