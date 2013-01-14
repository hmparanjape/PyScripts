"""
Change basis subroutine from VPSC 
into python script version
"""
import math
import numpy as np
def chg_basis(ce2=None, c2=None, ce4=None, c4=None, iopt=0, kdim=5):
    """
    parameter(sqr2 = math.sqrt(2.))
    parameter(rsq2 = 1./math.sqrt(2.))
    parameter(rsq3 = 1./math.sqrt(3.))

    ce2 = array(kdim)
    c2 = array(3,3)
    ce4 = array(kdim,kdim)
    c4 = array(3,3,3,3)

    iopt = 1:
        ce2  -- >  c2
    iopt = 2:
        c2   -- >  ce2e
    iopt = 3:
        ce4  -- > c4
    iopt = 4:
        c4   -- > ce4
    """
    sqr2 = math.sqrt(2.)
    rsq2 = 1./math.sqrt(2.)
    rsq3 = 1./math.sqrt(3.)
    rsq6 = 1./math.sqrt(6.)

    b = np.resize(np.array((0.)), (3,3,6))
    if ce2 == None: ce2 = np.resize(np.array((0.)), (kdim))
    if c2  == None: c2  = np.resize(np.array((0.)), (3,3))
    if ce4 == None: ce4 = np.resize(np.array((0.)), (kdim,kdim))
    if c4  == None: c4  = np.resize(np.array((0.)), (3,3,3,3))

    b[0][0][1] = -rsq6
    b[1][1][1] = -rsq6
    b[2][2][1] = 2.*rsq6

    b[0][0][0] = -rsq2
    b[1][1][0] = rsq2

    b[1][2][2] = rsq2
    b[2][1][2] = rsq2

    b[0][2][3] = rsq2
    b[2][0][3] = rsq2

    b[0][1][4] = rsq2
    b[1][0][4] = rsq2

    b[0][0][5] = rsq3
    b[1][1][5] = rsq3
    b[2][2][5] = rsq3

    if iopt == 0:
        return b

    elif iopt == 1:
        #print 'c2 = \n', c2
        #print 'ce2 = \n', ce2
        #print 'b = \n' , b
        for i in range(3):
            for j in range(3):
                c2[i][j] = 0.0
                for n in range(kdim):
                    c2[i][j] = c2[i][j]+ce2[n]*b[i][j][n]
        return c2

    elif iopt == 2:
        for n in range(kdim):
            ce2[n] = 0.0
            for i in range(3):
                for j in range(3):
                    ce2[n] = ce2[n] + c2[i][j]*b[i][j][n]
        return ce2
    
    elif iopt == 3:
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        c4[i][j][k][l] = 0.
                        for n in range(kdim):
                            for m in range(kdim):
                                c4[i][j][k][l] = c4[
                                    i][j][k][l] + ce4[n][m] * b[
                                    i][j][n] * b[k][l][m]
        return c4

    elif iopt == 4:
        for n in range(kdim):
            for m in range(kdim):
                ce4[n][m] = 0.
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            for l in range(3):
                                ce4[n][m] = ce4[
                                    n][m]+c4[i][j][k][l]* b[
                                    i][j][n]*b[k][l][m]
        return ce4


    
                
    

