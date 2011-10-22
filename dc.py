"""
disc-compression test pp
"""
import os, glob
import numpy as np

def dc(filename):
    FILE = open(filename,'r')
    lines = FILE.read()
    lines = lines.split('\n')
    #return lines
    force =[]; lr = []; lt=[]; t = []

    for i in range(len(lines)):
        if len(lines[i]) < 3: pass
        else:
            try:
                cwl = map(float, lines[i].split()[0:4])
            except: pass
            else:
                force.append(cwl[0]); lr.append(cwl[1])
                lt.append(cwl[2]); t.append(cwl[3])


    force = np.array(force); lr=np.array(lr)
    lt=np.array(lt);t=np.array(t)
    return force, lr, lt, t

def l2eps(l):
    if type(l).__name__!='ndarray':
        print 'The length must be the type of numpy array object'
    elif type(l).__name__=='list':
        print 'The input length is converted into array'
        l = np.array(l)
    eng = (l - l[0])/l[0]
    eps = np.log(eng + 1)
    return eps
    
def l2eps_comp(l):
    if type(l).__name__!='ndarray':
        print 'The length must be the type of numpy array object'
    elif type(l).__name__=='list':
        print 'The input length is converted into array'
        l = np.array(l)

    eng = (l - l[0])/l[0]
    eng = -eng
    eps = np.log(eng+1)
    return -eps
    
def area(a,b):
    return a/2. * b/2. * np.pi
    
def volume(a,b,t):
    return area(a,b)*t

def stress(force, area):
    stress = force / area
    #       kN    / mm^2
    # 10^3 N / 10^-6  mm^2
    # 10^9   N/m^2
    # 10^3 M Pa
    stress = stress * 1000
    return stress

def estimate_thick(a,b):
    return -(a+b)


def dct(filename='8#1.txt'):
    force, lr, lt, t = dc(filename)
    e_rd = l2eps(lr)
    e_td = l2eps(lt)
    e_t = -(e_rd+e_td)
    e_thick = l2eps_comp(t)
    a = area(lr,lt)
    v = volume(lr,lt,t)
    sig = stress(force = force, area = a)
    dv = e_rd + e_td + e_thick
    
    return -e_t, -e_thick, sig,v, dv

def plot():
    x,dum,y,dum,dum = dct(filename='8#1.txt')
    x,dum,y,dum,dum = dct(filename='8#2.txt')
    x,dum,y,dum,dum = dct(filename='10#1.txt')
    x,dum,y,dum,dum = dct(filename='10#2.txt')
    x,dum,y,dum,dum = dct(filename='15#1.txt')
    x,dum,y,dum,dum = dct(filename='15#2.txt')

    


