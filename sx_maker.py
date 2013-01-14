"""
VPSC SINGLE CRYSTAL PARAMETER FILE MANIPULATER
   UNIT CELL AXES, ANGLES
   ELASTIC STIFFNESS MATRIX
   THERMAL EXPANSION COEFFICIEINTS
   NMODESX
   NMODES
   MODE(i)   label of the modes to be used
   modex, nsmx, nrsx, iopsysx

   EITHER fcc, bcc

   tau0, tau1
   thet0, thet1
   hlat(nmodes)   : latent hardening coefficient
               - self and latent hardening coefficients
   
"""
import numpy as np

def sbsn(sb=[1,1,1], sn=[1,1,0]):
    """
    Provided slip directions as well as slip plane normals,
    unique slip system are returned
    """
    b, n = [], []
    #----- mutiplication of possible sets of slip system -----
    b.append([sb[0],sb[1],sb[2]])
    b.append([sb[0],sb[2],sb[1]])
    b.append([sb[1],sb[0],sb[2]])
    b.append([sb[1],sb[2],sb[0]])
    b.append([sb[2],sb[0],sb[1]])
    b.append([sb[2],sb[1],sb[0]])

    sb = __uniq__(b) #Filtering the uniqueness for slip direction

    n.append([sn[0],sn[1],sn[2]])
    n.append([sn[0],sn[2],sn[1]])
    n.append([sn[1],sn[0],sn[2]])
    n.append([sn[1],sn[2],sn[0]])
    n.append([sn[2],sn[0],sn[1]])
    n.append([sn[2],sn[1],sn[0]])
     
    sn = __uniq__(n) #Filtering the uniqueness for slip plane normal


    # Once again, multiplication is applied, however, this time by minus sign.
    # sign vector
    sign = [[ 1, 1, 1],[-1, 1, 1],[ 1,-1, 1],[ 1, 1,-1],
            [-1,-1, 1],[-1, 1,-1],[ 1,-1,-1],[-1,-1,-1]]

    # temp slip normal and direction
    tempsn , tempsb  = [], []
    tempsn_, tempsb_ = [], []

    # Each equivalent slip normal is multiplied by the sign vector defined above
    # There are still unique vector. (tempsn_)
    iss = 0    
    for i in range(len(sn)):
        tempsn.append([])
        for j in range(len(sign)):
            tempsn[i].append([])
            for k in range(3):  tempsn[i][j].append(sn[i][k]*sign[j][k])
            if any(tempsn[i][j]==tempsn_[k] for k in range(len(tempsn_))): pass
            else: tempsn_.append(tempsn[i][j])
            iss == iss + 1
            pass
        pass
    
    #As of here, equivalent slip directions are created.
    # The unique positive slip direction is multiplied by the sign vector
    iss = 0
    for i in range(len(sb)):
        tempsb.append([])
        for j in range(len(sign)):
            tempsb[i].append([])
            for k in range(3): tempsb[i][j].append(sb[i][k]*sign[j][k])
            if any(tempsb[i][j]==tempsb_[k] for k in range(len(tempsb_))): pass
            else: tempsb_.append(tempsb[i][j])
            iss == iss + 1
            pass
        pass
    sb = tempsb_

    ## the opposite normal to the exisiting slip plane is not accounted.
    sn_ind=[]
    sn = tempsn_
    for i in range(len(sn)):
        if any(sn[i] == [-sn_ind[j][0],
                         -sn_ind[j][1],
                         -sn_ind[j][2]] for j in range(len(sn_ind))): pass
        else: sn_ind.append(sn[i])
    sn = sn_ind

    ## check if the slip directions are perpendicular to the slip plane.
    ss = []
    for i in range(len(sn)):
        ss.append([sn[i]])
        for j in range(len(sb)):
            temp = sn[i][0]*sb[j][0]+sn[i][1]*sb[j][1]+sn[i][2]*sb[j][2]
            if abs(temp) < 0.00001 :  ss[i].append(sb[j]) #They must be normal!
            pass
        pass
    ##------------------------------------------------------------------
    
    tmp = []
    iss = 0
    for i in range(len(ss)):
        if len(ss[i])==1: pass #if there is only one slip direction!
        else:
            t = ss[i][1:len(ss[i])]
            tmp.append([])
            tmp[iss].append(ss[i][0])
            for i in range(len(t)):
                if any(t[i]==[-tmp[iss][k][0],-tmp[iss][k][1],
                               -tmp[iss][k][2]] for k in range(len(tmp[iss]))): pass
                else: tmp[iss].append(t[i])
            iss = iss + 1
            pass
        pass

    ##  tmp is a slip system in which the opposite (reverse) slip direction
    ##  has been intentionally is ommited.
    
    return tmp, ss

def __uniq__(s):
    """
    Provided a set of the vector, (expecting a list type of variable)
    only returns the unique vector.
    """
    # the expected type of s will be like
    # --   s = [  [1,2,3], [2,1,3], [3,2,1]]
    temp = []
    for i in range(len(s)):
        if i==0: temp.append(s[i])
        else: 
            if any(s[i]==temp[j] for j in range(len(temp))): pass
            else : temp.append(s[i])
            pass
        pass
    return temp

def __ntot__(ss):
    """
    Returns the total number of slip system depending upon iopsysx
    """
    nss = 0 
    for i in range(len(ss)):
        nss = nss + len(ss[i])-1
    return nss

def cubic(filename = None, ihardlaw='voce',
          hii = 1.0, hij = 1.4, hb = -0.4, hp = 0.9, # coplanar
          tau0 = 1.0, tau1= 0.2, thet0 = 1.0, thet1 = 0.05,
          hpfac=0., gndfac = 0., header ='** material info: ',
          c11=205., c12=138., c44 =126., ##elastic moduli matrix needs only three componen
                                         ## due to the crystallographic symmetry
          
          iopsysx = 1, #iopsysx(1): a slip mode includes two glide direction (+ and -)
          dependency = 'indv',
          nrsx = 20,
          twshx = 0., isectw = 0, thres1 =0., thres2 =0.,
          interaction ='Bauschinger Latent', b=[[1,1,0]], n=[[1,1,1]],
          copl = False,  #coplanar?

          ## dicloation hardeing single crystal file ('rauch' case)
          idislmode=1,  #idislmode for all is active(1) or active-inactive mode(0)
          ibau=1,       #ibau(1):sudden drop, ibau(0): no such drop
          pp=0.8,       #used to drop the sudden density level.
          ilat=1,       #ilat(1):latent hardening based on dislocation density

          ##
          
          grsize=3.0e-05,
          burgers=2.46e-10,
          fdisln=180,
          shearm=8.5e+04,
          ftau0=38.,
          portion=0.8,
          ftherm=2.8
          ):
    """
    ** INFO
    filename = 'mnt/hgfs/10_19_simple_shear/VPSC/sx/hiihijhb.sx',
    ihardlaw='voce' or 'rauch'
    header = '**materialinfo',
    slipsystem information should be given
    len(ss) = # of slip system set
    [  [ [b1,b2,b3],[n1,n2,n3]    ]   [    [b1,b2,b3],[n1,n2,n3]    ]   ] 
    
    ** HARDENING MATRIX
    hii = 1.0, hij = 1.0, hb = 1.0, hp=0.9 # colplanar coefficient is added.
    hpfac = 0., gndfac = 0.,

    ** VOCE HARDENING PARAMETERS, RATE-DEPENDENT EXPONENT
    tau0 = 1.045e2, tau1 = 70., thet0 = 2.6e2, thet1 = 0.95e2,
    nrsx = 20,

    ** SINGLE CUBIC CRYSTAL ELASTIC COEFFICIENTS
    c11 = 205, c12 = 138., c44 = 126.,

    ** INTERACTION AND LATENT HARDENING PARAMETERS
    iopsysx = 1, 0 (1): a slip mode includes two glide direction (+ and -)
                   (0): one slip mode for one glide direction [110] != [-1-10]
    dependency = 'iso', 'indv'
    interaction = 'Bauschinger', 'Bauschinger Latent', 'Latent'

    ** TWINNING RELATED PARAMTERS
    twshx = 0., isectw = 0, thres1 = 0., thre2 = 0.,

    
    ** dicloation hardeing single crystal file
    grsize=3.0e-05, burgers=2.46e-10, fdisln=180, shearm=8.5e+04,
    ftau0=38., portion = 0.8, ftherm = 2.8
    """
    f = open(filename ,'w')
    f.writelines(header + '\n')
    f.writelines('cubic\n') #           crysym \n')
    f.writelines('%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f'%(1.0, 1.0, 1.0, 90., 90., 90.,))
    f.writelines('  unit cell axes and angles \n')
    # SINGLE CRYSTAL ELASTIC STIFFNESS MATRIX :
    f.writelines('** Elastic stiffness (single crystal [GPa];')
    f.writelines('scaled = 0.85 x interpolated)\n')
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %
                 (c11, c12, c12, 0., 0., 0.))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %
                 (c12, c11, c12, 0., 0., 0.))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %
                 (c12, c12, c11, 0., 0., 0.))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %
                 (0., 0., 0., c44, 0., 0.,))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %
                 (0., 0., 0., 0., c44, 0. ))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %
                 (0., 0., 0., 0., 0., c44))
    # THERMAL EXPANSION COEFFICIENTS (single crystal in crystal axis)
    f.writelines('** Thermal expansion coefficients')
    f.writelines('(single crystal axis): \n')
    f.writelines(' %9.1e%9.1e%9.1e%9.1e%9.1e%9.1e \n' %
                 ( 10.0e-6,10.0e-6,10.0e-6, 0., 0., 0.))
    f.writelines('** Slip and Twinning modes (system) information: \n')
    ## BASIC SLIP SYSTEM DECLARATION OF cubic system
    ## NOW sbsn method is being used here.
    ## provided slip direction and slip plane normal,

    ss = []
    sss = []
    
    for i in range(len(n)): #loop over slip plane
        tmps1, tmps2 = sbsn(sb=b[i], sn=n[i])  #slip direction and slip plane normal
        ## slip system in which the opposite direction is not included.
        ss.append(tmps1) 
        ## slip system in which the opposite direction is independently included
        sss.append(tmps2)
        pass

    ss_ = []

    #if iopsysx == 0: ss = sss
    #elif iopsysx == 1: ss = ss
    #else: print 'Error: inappropriate iopsysx input'; raise IOError

    nss = 0
    iss = 0
    ntot = 0
    ss_header = []

    for n in range(len(ss)): ntot = ntot +  __ntot__(ss = ss[n])

    if iopsysx==0: ntot = ntot * 2 ## if yes, total number of slip system is doubled
                                   ## This differentiate the reverse slip direction from
                                   ## the forward direction.

    if dependency=='indv':
        ## number of total slip system
        for i in range(len(ss)): nss = nss + len(ss[i])

        ## slip system writing loop
        for i in range(len(ss)):  # slip system (either on (111) or (112).. )
            for j in range(len(ss[i])): # slip plane 
                if iopsysx==1:
                    for k in range(len(ss[i][j])-1): #slip direction
                        ### slip system info: -------------------------------
                        tmp = (
                            '(%2i%2i%2i) <%2i%2i%2i> system \n'%
                            (ss[i][j][0][0],
                             ss[i][j][0][1],
                             ss[i][j][0][2],
                             ss[i][j][k+1][0],
                             ss[i][j][k+1][1],
                             ss[i][j][k+1][2]))
                        tmp = tmp + (
                            '%6i%6i%6i%6i   %s \n'%
                            (iss+1, 1,nrsx, iopsysx,
                             'modex, nsmx, nrsx, iopsysx'))
                        tmp = tmp + (
                            '%6.2f%6i%6.2f%6.2f   %s\n'%
                            (twshx, isectw, thres1, thres2,
                             'twshx,isectw,thres1,thres2'))
                        if ihardlaw=='voce':
                            tmp = tmp + (
                                '%8.3e  %8.3e  %8.3e  %8.3e  %8.1e  %8.1e \n'%
                                (tau0, tau1, thet0, thet1, hpfac, gndfac))
                        elif ihardlaw=='rauch':
                            tmp = tmp + ('%6i%6i%6i%6.2f'%(idislmode,ilat, ibau, pp))
                            tmp = tmp + ('  idislmode, ilat, ibau, pp\n')
                            tmp = tmp + (
                                '%8.3e  %8.3e  %8.3e  %8.3e '%(
                                    grsize, burgers, fdisln, shearm))
                            tmp = tmp + ('%8.3e  %8.3e  %8.3e  %8.1e  %8.1e '%(
                                    ftau0, portion, ftherm, hpfac, gndfac))
                            tmp = tmp + ('%s'%
                                         'grsize, burgers, fdisln, shearm,')
                            tmp = tmp + ('%s\n'%
                                         'ftau0, portion, ftherm, hpfac, gndfac')

                        else:
                            raise IOError,'wrong ihardlaw input'
                        ## -----------------------------------------------------------------
                        #  Writes hardening matrix components
                        # iopsysx=1: Bidirectional

                        ssi = 0
                        for kk in range(len(ss)): # loop over kinds of slip system 
                            for ii in range(len(ss[i])): # loop over (111) (11-1), (1-11)...
                                for jj in range(len(ss[i][ii])-1):  # loop over direction
                                    spn = np.array(ss[i][ii][0])  ## slip plane normal
                                    sd = np.array(ss[i][ii][jj+1]) # slip direction
                                    
                                    # if it is the same slip mode (self-hardening)
                                    if ssi==iss: tmp = tmp + ('%5.2f  '%(hii)) #self-hardening
                                    #ss           --> slip system( (111), (112) ... )
                                    #ss[0][0]     --> slip directions
                                    #ss[0][0][0]  --> the fist slip direction:representative
                                    #ss[0][0][1:] --> the rest slip directions...

                                    # if on the same slip plane
                                    elif all(spn[q]==ss[i][j][0][q] for q in range(3)): #j:slip plane
                                        # if coplanar slip direction has a different hardening rate
                                        if copl==True: tmp = tmp + ('%5.2f  '%(hp))
                                        # if not, impose the latent hardening rate.
                                        else: tmp = tmp + ('%5.2f  '%(hij))
                                        pass
                                    #not even in the sample plane                                            
                                    else:
                                        tmp = tmp + ('%5.2f  '%(hij))
                                        pass
                                    ssi = ssi + 1 # counts the slip system
                                    pass
                                pass
                            pass

                        # the old way below:
                        # for l in range(ntot):
                        #     ## Either self-hardening or latent-hardening
                        #     if l==iss: tmp = tmp + ('%5.2f  '%(hii)) #iss: slip mode index
                        #     else:
                        #         tmp = tmp + ('%5.2f  '%(hij))
                        #         pass
                        #     pass
                        
                        ## -----------------------------------------------------------------
                        tes = ss[i][j][0] + ss[i][j][k+1] #combination of slip plane normal and direction
                        tmp = tmp + (
                            '\n %3i%3i%3i  %3i%3i%3i \n'%(
                                tes[0],tes[1],tes[2],
                                tes[3],tes[4],tes[5]))
                        ss_header.append(tmp)
                        iss = iss + 1
                        pass
                    pass
                        
                elif iopsysx==0:
                    for k in range(len(ss[i][j])-1):
                        #Slip direction in forward-wise
                        tmp = (
                            '(%2i%2i%2i) <%2i%2i%2i> system \n'%(
                                ss[i][j][0][0],
                                ss[i][j][0][1],
                                ss[i][j][0][2],
                                ss[i][j][k+1][0],
                                ss[i][j][k+1][1],
                                ss[i][j][k+1][2]))
                        tmp = tmp + (
                            '%6i%6i%6i%6i   %s \n'%(
                                iss+1, 1, nrsx, iopsysx, 
                                'modex, nsmx, nrsx, iopsysx'))
                        tmp = tmp + (
                            '%6.2f%6i%6.2f%6.2f   %s \n'%(
                                twshx, isectw, thres1, thres2,
                                'twshx,isectw,thres1,thres2'))
                        if ihardlaw=='voce':
                            tmp = tmp + (
                                '%8.3e  %8.3e  %8.3e  %8.3e  %8.1e  %8.1e \n'%
                                (tau0, tau1, thet0, thet1, hpfac, gndfac))
                        elif ihardlaw=='rauch':
                            tmp = tmp + ('%6i%6i%6i%6.2f'%(idislmode,ilat, ibau, pp))
                            tmp = tmp + ('  idislmode, ilat, ibau, pp\n')
                            tmp = tmp + (
                                '%8.3e  %8.3e  %8.3e  %8.3e '%(
                                    grsize, burgers, fdisln, shearm))
                            tmp = tmp + ('%8.3e  %8.3e  %8.3e  %8.1e  %8.1e '%(
                                    ftau0, portion, ftherm, hpfac, gndfac))
                            tmp = tmp + ('%s'%
                                         'grsize, burgers, fdisln, shearm,')
                            tmp = tmp + ('%s\n'%
                                         'ftau0, portion, ftherm, hpfac, gndfac')
                            
                        else:
                            raise IOError,'wrong ihardlaw input'
                        
                        """
                        j stands for the slip plane index.
                        """
                        ## -----------------------------------------------------------------
                        #  Writes hardening matrix components 
                        ssi = 0
                        for kk in range(len(ss)): # loop over kinds of slip system 
                            for ii in range(len(ss[i])): # loop over (111) (11-1), (1-11)...
                                for jj in range(len(ss[i][ii])-1):  # loop over direction
                                    for ll in range(2):
                                        spn = np.array(ss[i][ii][0])  ## slip plane normal
                                        sd = np.array(ss[i][ii][jj+1])
                                        if ll%2==0: pass
                                        else: sd = sd * -1
                                        if ssi==iss: tmp = tmp + ('%5.2f  '%(hii))
                                        elif all(spn[q]==ss[i][j][0][q] for q in range(3)): #
                                            if ssi==iss + 1: tmp = tmp + ('%5.2f  '%(hb))
                                            else:
                                                if copl==True: tmp = tmp + ('%5.2f  '%(hp))
                                                else: tmp = tmp + ('%5.2f  '%(hij))
                                        else: tmp = tmp + ('%5.2f  '%(hij))
                                        ssi = ssi + 1
                                        pass
                                    pass
                                pass
                            pass
                        ## Hardening matrix writing ends
                        ## -----------------------------------------------------------------
                        
                        
                        tes = ss[i][j][0] + ss[i][j][k+1]
                        tmp = tmp + (
                            '\n %3i%3i%3i  %3i%3i%3i \n'%(
                                tes[0],tes[1],tes[2],
                                tes[3],tes[4],tes[5]))
                        ss_header.append(tmp)
                        iss = iss + 1
                        
                        #Slip direction in backward-wise
                        tmp = (
                            '(%2i%2i%2i) <%2i%2i%2i> system \n'%(
                                ss[i][j][0][0],
                                ss[i][j][0][1],
                                ss[i][j][0][2],
                                ### Accordingly, negative direction is enforced
                                -ss[i][j][k+1][0],  
                                -ss[i][j][k+1][1],
                                -ss[i][j][k+1][2]))
                        tmp = tmp + (
                            '%6i%6i%6i%6i   %s \n'%(
                                iss+1, 1, nrsx, iopsysx,
                                'modex, nsmx, nrsx, iopsysx'))
                        tmp = tmp + (
                            '%6.2f%6i%6.2f%6.2f   %s \n'%(
                                twshx, isectw, thres1, thres2,
                                'twshx,isectw,thres1,thres2'))
                        if ihardlaw=='voce':
                            tmp = tmp + (
                                '%8.3e  %8.3e  %8.3e  %8.3e  %8.1e  %8.1e \n'%
                                (tau0, tau1, thet0, thet1, hpfac, gndfac))
                        elif ihardlaw=='rauch':
                            tmp = tmp + ('%6i%6i%6i%6.2f'%(idislmode, ilat, ibau, pp))
                            tmp = tmp + ('  idislmode, ilat, ibau, pp\n')
                            tmp = tmp + (
                                '%8.3e  %8.3e  %8.3e  %8.3e '%(
                                    grsize, burgers, fdisln, shearm))
                            tmp = tmp + ('%8.3e  %8.3e  %8.3e  %8.1e  %8.1e '%(
                                    ftau0, portion, ftherm, hpfac, gndfac))
                            tmp = tmp + ('%s'%
                                         'grsize, burgers, fdisln, shearm,')
                            tmp = tmp + ('%s\n'%
                                         'ftau0, portion, ftherm, hpfac, gndfac')
                            
                        else:
                            raise IOError,'wrong ihardlaw input'                        
                        # for l in range(ntot):
                        #     if l==iss: tmp = tmp + ('%5.2f  '%(hii))
                        #     elif l==iss-1: tmp = tmp + ('%5.2f  '%(hb))
                        #     else: tmp = tmp + ('%5.2f  '%(hij))A
                        """
                        j stands for the slip plane index.
                        """
                        ssi = 0
                        for kk in range(len(ss)): # loop over kinds of slip system 
                            for ii in range(len(ss[i])): # loop over (111) (11-1), (1-11)...
                                for jj in range(len(ss[i][ii])-1):  # loop over direction
                                    for ll in range(2):
                                        spn = np.array(ss[i][ii][0])  ## slip plane normal
                                        sd = np.array(ss[i][ii][jj+1])
                                        if ll%2==0: pass
                                        else: sd = sd * -1
                                        if ssi==iss: tmp = tmp + ('%5.2f  '%(hii))
                                        elif all(spn[q]==ss[i][j][0][q] for q in range(3)): #
                                            if ssi==iss - 1: tmp = tmp + ('%5.2f  '%(hb))
                                            else:
                                                if copl==True: tmp = tmp + ('%5.2f  '%(hp))
                                                else: tmp = tmp + ('%5.2f  '%(hij))
                                        else: tmp = tmp + ('%5.2f  '%(hij))
                                        ssi = ssi + 1
                                        

                        tes = ss[i][j][0] + [-ss[i][j][k+1][0],
                                             -ss[i][j][k+1][1],
                                             -ss[i][j][k+1][2]]
                        tmp = tmp + (
                            '\n %3i%3i%3i  %3i%3i%3i \n'%
                            (tes[0],tes[1],tes[2],tes[3],tes[4],tes[5]))
                        ss_header.append(tmp)
                        iss = iss + 1
                else: 
                    print 'Error: Inappropriate iopsysx'
                    raise IOError

    elif dependency=='iso':
        for i in range(len(ss)): #Unit of block 
        #(loop over block consisting of several equivalent slip systems)
            nss = 0
            for j in range(len(ss[i])):
                for k in range(len(ss[i][j])-1): 
                    nss = nss + 1
                
            if iopsysx==0: nss = nss * 2
            elif iopsysx==1 : pass
            else: 
                print 'Error: Inappropriate iopsysx'
                raise IOError
            
            tmpsn = ss[i][0][0]
            tmpsb = [abs(ss[i][0][1][0]),
                     abs(ss[i][0][1][1]), 
                     abs(ss[i][0][1][2])]

            tmp = ('(%2i%2i%2i) <%2i%2i%2i>  system\n'%
                   (tmpsn[0],tmpsn[1],tmpsn[2],
                    tmpsb[0],tmpsb[1],tmpsb[2]))
            tmp = tmp + ('%6i%6i%6i%6i   %s \n'%
                         (i+1,nss, nrsx, iopsysx,
                          ' modesx, nsmx, nrsx, iopsysx'))
            tmp = tmp + ('%6.2f%6i%6.2f%6.2f   %s \n'%
                         (twshx, isectw, thres1, thres2,
                          ' twshx, isectw, thres1, thres2'))

            if ihardlaw=='voce':
                tmp = tmp + (
                    '%8.3e  %8.3e  %8.3e  %8.3e  %8.1e  %8.1e \n'%
                    (tau0, tau1, thet0, thet1, hpfac, gndfac))
            elif ihardlaw=='rauch':
                tmp = tmp + ('%6i%6i%6i%6.2f'%(idislmode, ilat, ibau, pp))
                tmp = tmp + ('  idislmode, ilat, ibau, pp\n')                
                tmp = tmp + (
                    '%8.3e  %8.3e  %8.3e  %8.3e '%(
                        grsize, burgers, fdisln, shearm))
                tmp = tmp + ('%8.3e  %8.3e  %8.3e  %8.1e  %8.1e '%(
                        ftau0, portion, ftherm, hpfac, gndfac))
                tmp = tmp + ('%s'%
                             'grsize, burgers, fdisln, shearm,')
                tmp = tmp + ('%s\n'%
                             'ftau0, portion, ftherm, hpfac, gndfac')
                
            else:
                raise IOError,'wrong ihardlaw input'
            
            for j in range(len(ss)):
                if j == i: tmp = tmp + ('%5.2f  '%(hii))
                else: tmp = tmp + ('%5.2f  '%(hij))
                
            tmp = tmp + '\n'
            for j in range(len(ss[i])):
                for k in range(len(ss[i][j])-1):
                    tes = ss[i][j][0]+ss[i][j][k+1]
                    tmp = tmp + (' %3i%3i%3i  %3i%3i%3i \n'%
                                 (tes[0],tes[1],tes[2],
                                  tes[3],tes[4],tes[5]))
                    if iopsysx==0:
                        tmp = tmp + (' %3i%3i%3i  %3i%3i%3i \n'%
                                     (tes[0], tes[1], tes[2],
                                      -tes[3], -tes[4], -tes[5]))
            ss_header.append(tmp)
    else: 
        print 'Error: Inappropriate dependency input'
        print dependency
        raise IOError


    if dependency=='indv': pass
    elif dependency=='iso': ntot = len(ss)
    f.writelines('%4i    %s\n'%
                 (ntot, 'nmodesx (total # of modes listed in file)'))
    f.writelines('%4i    %s\n'%
                 (ntot, 
                  'nmodes  (# of modes to be used in the calculation)'))
    for i in range(ntot):
        f.writelines('%4i'%(i+1))
    f.writelines('  %s\n'%
                 ('mode(i)    (label of the modes to be used)'))
           
    for i in range(len(ss_header)):
        f.writelines(ss_header[i])
        pass
    pass
        
def invariant(filename='sx/invariants/alpha.sx',
              header = '*VARIANTS',
              crysym = 'cubic',
              b=[[1,1,0]], n=[[1,1,1]],
              iopsysx=0 #0 for 24, 1 for 12
              ):
    """
    Provided invariant system will be written down to the file
    whose name is given through filename argument.
    
    Mother and Child phases are distinguished in that
    mothers usually have multiple choices while childs have
    only one resulting crystallographic system.


    currently only one of invariant family set is possible.
    """
    FILE=open(filename,'w')
    if crysym =='cubic':
        a=[1. ,1. ,1. ]   #unit cell axes
        d=[90.,90.,90.]#unit cell angles
    else:
        print 'No else that cubic system is prepared yet'
        raise IOError

    
    ss = sbsn(sb=b[0], sn=n[0])   # Later this can be in a loop
    
# def sbsn(sb=[1,1,1], sn=[1,1,0]):
#     """
#     Provided slip directions as well as slip plane normals,
#     unique slip system are returned
#     """
    
    if   iopsysx ==  0:
        ss = ss[1]
    elif iopsysx ==  1:
        ss = ss[0]
    elif iopsysx == -1:
        ss = [ss[0][0][0:2]]
    else: print 'Error: Inappropriate iopsysx input'; raise IOError


    nm = 0
    for i in range(len(ss)):
        for j in range(len(ss[i])):
            if j==0: pass
            else: nm = nm + 1  #len(ss[i]) minus the plane

    nf = 1
    FILE.writelines(header+'\n')
    FILE.writelines('%s       crysym \n'%(crysym))
    FILE.writelines('%4.2f %4.2f %4.2f '%(a[0],a[1],a[2]))
    FILE.writelines('%4.2f %4.2f %4.2f '%(d[0],d[1],d[2]))
    FILE.writelines('    unit cell axes and angles \n')
    FILE.writelines('*Transformation invariants \n')
    FILE.writelines('%4i   nfamilies\n'%(nf))
    FILE.writelines('%4i  multiplicit(nfamilies)\n'%(nm))

    for i in range(len(ss)):
        spn = ss[i][0]
        for j in range(len(ss[i])-1):
            FILE.writelines('  %3i%3i%3i   '%(spn[0],spn[1],spn[2]))
            for k in range(3):
                FILE.writelines('%3i'%(ss[i][j+1][k]))
            FILE.writelines('\n')


def fcc(filename = '/mnt/hgfs/10_19_simple_shear/VPSC/sx/hijhiihb.sx',
        hii = 1.0, hij = 1.4, hb = -0.4,
        tau0 = 1.0, tau1 = 0.2, thet0 = 1.0, thet1 = 0.05,
        #tau0 = 1.045e2, tau1= 70., thet0 = 2.6e2, thet1 = 0.95e2,
        hpfac = 0., gndfac = 0., header = '** material info',
        c11 = 205., c12 = 138., c44 = 126.,
        #nmodesx = 2, modes =[1],
        iopsysx = 0,
        dependency = 'indv',
        #dependency = 'iso',
        nrsx = 20,
        twshx = 0., isectw = 0, thres1 = 0., thres2 = 0.,
        #interaction ='Bauschinger'
        interaction ='Bauschinger Latent'
        ):    
    """
    ** INFO
    filename = 'sx/fcc.sx'
    header = '**materialinfo',
    
    ** HARDENING MATRIX
    hii = 1.0, hij = 1.0, hb = 1.0
    hpfac = 0., gndfac = 0.,

    ** VOCE HARDENING PARAMETERS, RATE-DEPENDENT EXPONENT
    tau0 = 1.045e2, tau1 = 70., thet0 = 2.6e2, thet1 = 0.95e2
    nrsx = 20,

    ** SINGLE CUBIC CRYSTAL ELASTIC COEFFICIENTS
    c11 = 205, c12 = 138., c44 = 126.

    ** INTERACTION AND LATENT HARDENING PARAMETERS
    iopsysx = 1, 0
    dependency = 'iso', 'indv'
    interaction = 'Bauschinger', 'Bauschinger Latent', 'Latent'

    ** TWINNING RELATED PARAMTERS
    twshx = 0., isectw = 0, thres1 = 0., thre2 = 0.
    """

    
    f = open(filename, 'w')
    f.writelines(header + '\n')
    f.writelines('cubic \n')
    f.writelines('%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f'%(1.0, 1.0, 1.0, 90., 90., 90.,))
    f.writelines('  unit cell axes and angles \n')

    # SINGLE CRYSTAL ELASTIC STIFFNESS MATRIX :
    f.writelines('** Elastic stiffness (single crystal [GPa]; scaled = 0.85 x interpolated)\n')
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %(c11, c12, c12, 0., 0., 0.))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %(c12, c11, c12, 0., 0., 0.))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %(c12, c12, c11, 0., 0., 0.))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %(0., 0., 0., c44, 0., 0.,))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %(0., 0., 0., 0., c44, 0. ))
    f.writelines(' %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f \n' %(0., 0., 0., 0., 0., c44))

    # THERMAL EXPANSION COEFFICIENTS (single crystal in crystal axis)
    f.writelines('** Thermal expansion coefficients (single crystal axis): \n')
    f.writelines(' %9.1e%9.1e%9.1e%9.1e%9.1e%9.1e \n' %( 10.0e-6,10.0e-6,10.0e-6, 0., 0., 0.))


    ## BASIC SLIP SYSTEM DECLARATION OF FCC STRUCTURE {111}<110>
    sn = []
    sb = []

    sn.append([])
    sb.append([])
    sb[0].append([])
    sb[0].append([])
    sb[0].append([])
    sn[0] = [1,1,1]    
    sb[0][0] = [0,1,-1]
    sb[0][1] = [1, 0, -1]
    sb[0][2] = [1,-1, 0]

    sn.append([])
    sb.append([])
    sb[1].append([])
    sb[1].append([])
    sb[1].append([])
    sn[1] = [-1, 1, 1]
    sb[1][0] = [0, 1, -1]
    sb[1][1] = [1, 0, 1]
    sb[1][2] = [1, 1, 0]

    sn.append([])
    sb.append([])
    sb[2].append([])
    sb[2].append([])
    sb[2].append([])
    sn[2] = [-1, -1, 1]
    sb[2][0] = [0, 1, 1]
    sb[2][1] = [1, 0, 1]
    sb[2][2] = [1, -1, 0]

    sn.append([])
    sb.append([])
    sb[3].append([])
    sb[3].append([])
    sb[3].append([])
    sn[3] = [1, -1, 1]
    sb[3][0] = [0, 1, 1]
    sb[3][1] = [1, 0, -1]
    sb[3][2] = [1, 1, 0]

    #dependency !
    # dependency may be 'iso', 'indv'

    ##

    # INFORMATION ABOUT SLIP $ TWINNING MODES IN THIS FILE:

    #nmodesx = 2, modes =[1],
    modes = []
    if dependency =='iso':
        nmodesx = 1
        modes.append(1)
    elif dependency =='indv':
        if iopsysx == 0:
            nmodesx = len(sn)*len(sb[0][0])  * 2
            pass
        elif iopsysx == 1:
            nmodesx = len(sn)*len(sb[0][0])
            pass
        else: return -1
        for i in range(nmodesx):
            modes.append(i+1)

        
    f.writelines('** Slip and Twinning modes(system) information \n')
    f.writelines( str(len(modes)).ljust(4) + '    nmodesx \n')
    f.writelines( str(len(modes)).ljust(4) + '    nmodes \n')    
    for i in range(len(modes)):
        f.writelines(str(modes[i]).ljust(3) )
    f.writelines('  modes \n')

    if dependency == 'iso' :
        
        """ iso case"""
        nsmx = len(sn)* len(sb[0][0])
        if iopsysx == 0: nsmx = nsmx * 2
       
        f.writelines(' {111}<110> slip \n')
        f.writelines('%5i %5i %5i %5i'  %(1, nsmx, nrsx, iopsysx))
        f.writelines('              modex, nsmx, nrsx, oiopsysx \n')
        f.writelines(' %5.3f   %1i     %5.3f   %5.3f ' %(twshx,isectw,thres1,thres2))
        f.writelines('        twshx,isectw,thres1,thres2   \n')
        f.writelines(' %5.3e %5.3e %5.3e %5.3e %5.1f %5.1f \n' %(tau0, tau1, thet0, thet1, hpfac, gndfac))
        f.writelines(' 1.0                                 h_latent  \n')
        for i in range(len(sn)):
            for j in range(len(sb[i])):
                if iopsysx == 1:
                    f.writelines(' %2i  %2i  %2i     %2i  %2i  %2i \n' %(sn[i][0], sn[i][1], sn[i][2],
                                                                         sb[i][j][0], sb[i][j][1], sb[i][j][2]))
                elif iopsysx == 0:
                    f.writelines(' %2i  %2i  %2i     %2i  %2i  %2i \n' %(sn[i][0], sn[i][1], sn[i][2],
                                                                         sb[i][j][0], sb[i][j][1], sb[i][j][2]))
                    f.writelines(' %2i  %2i  %2i     %2i  %2i  %2i \n' %(sn[i][0], sn[i][1], sn[i][2],
                                                                         -sb[i][j][0],-sb[i][j][1],-sb[i][j][2]))

    if dependency == 'indv':
        """ individual crss evolution for each slip system """
        nsmx = 1

        # interaction:
        #        'Bauschinger'
        #        'Latent'
        #        'Bauschinger Latent'
        hmat = []
        if iopsysx == 1: m=1
        elif iopsysx == 0: m=2
        if iopsysx ==1:
            for i in range (len(sn)*len(sb[0])*m):
                if interaction =='Latent':
                    hmat.append(hij)
                    pass
                elif interaction =='Bauschinger Latent':
                    hmat.append(hij)
                    pass
                elif interaction =='Bauschinger':
                    hmat.append(hii)
                    pass
        kount = 0
        for j in range(len(sn)):
            for k in range(len(sb[j])):
                if iopsysx == 1:  #'Bauschinger impossible'
                    f.writelines(' (%2i%2i%2i) <%2i%2i%2i> system \n' %(sn[j][0], sn[j][1], sn[j][2],
                                                                        sb[j][k][0], sb[j][k][1], sb[j][k][2]))

                    f.writelines('%5i %5i %5i %5i'  %(kount+1, nsmx, nrsx, iopsysx))
                    f.writelines('              modex, nsmx, nrsx, oiopsysx \n')
                    f.writelines(' %5.3f   %1i     %5.3f   %5.3f ' %(twshx,isectw,thres1,thres2))
                    f.writelines('        twshx,isectw,thres1,thres2   \n')
                    f.writelines(' %5.3e %5.3e %5.3e %5.3e %5.1f %5.1f\n'%(tau0, tau1, thet0, thet1, hpfac, gndfac))
                    #hlat
                    if interaction =='Bauschinger':
                        count = 0
                        for i in range(len(sn)*len(sb[0])*m):
                            f.writelines(' %5.2f ' %(hmat[count]))
                            count = count + 1
                        f.writelines('\n')
                    elif interaction =='Latent':
                        count = 0
                        for i in range(len(sn)*len(sb[0])*m):
                            if count == kount:
                                f.writelines(' %5.2f ' %(hii))
                            else:
                                f.writelines(' %5.2f ' %(hmat[count]))
                            count = count + 1
                        f.writelines('\n')
                    elif interaction =='Bauschinger Latent':
                        print
                        print 'Basuchinger Latent mode for iopsysx == 1 is not possible'
                        print 'Therefore, it has been hardwired to be Latent in this case'
                        print
                        count = 0
                        for i in range(len(sn)*len(sb[0])*m):
                            if count == kount:
                                f.writelines(' %5.2f ' %(hii))
                            else:
                                f.writelines(' %5.2f ' %(hmat[count]))
                            count = count + 1
                        f.writelines('\n')
                        
                    
                    #index
                    f.writelines(' %2i  %2i  %2i     %2i  %2i  %2i \n' %(sn[j][0], sn[j][1], sn[j][2],
                                                     sb[j][k][0], sb[j][k][1], sb[j][k][2]))
                    kount = kount + 1
                if iopsysx == 0:   #'Bauschinger Latent both possible'

                    f.writelines(' (%2i%2i%2i) <%2i%2i%2i> system \n' %(sn[j][0], sn[j][1], sn[j][2],
                                                                        sb[j][k][0], sb[j][k][1], sb[j][k][2]))
                    f.writelines('%5i %5i %5i %5i'  %(kount+1, nsmx, nrsx, iopsysx))
                    f.writelines('              modex, nsmx, nrsx, oiopsysx \n')
                    f.writelines(' %5.3f   %1i     %5.3f   %5.3f ' %(twshx,isectw,thres1,thres2))
                    f.writelines('        twshx,isectw,thres1,thres2   \n')
                    f.writelines(' %5.3e %5.3e %5.3e %5.3e %5.1f %5.1f \n' %(tau0, tau1, thet0, thet1, hpfac, gndfac))
                    #hlat
                    count = 0
                    for i in range(len(sn)*len(sb[0])*m):
                        if kount == count:
                            f.writelines(' %5.2f ' %(hii))
                        elif kount+1 == count:
                            f.writelines(' %5.2f ' %(hb))
                        else:
                            f.writelines(' %5.2f ' %(hij))
                        count = count + 1
                    f.writelines('\n')
                    
                    #index
                    f.writelines(' %2i  %2i  %2i     %2i  %2i  %2i \n' %(sn[j][0], sn[j][1], sn[j][2],
                                 sb[j][k][0], sb[j][k][1], sb[j][k][2]))
                    kount = kount + 1

                    f.writelines(' (%2i%2i%2i) <%2i%2i%2i> system \n' %(sn[j][0], sn[j][1], sn[j][2],
                                                                        -sb[j][k][0],-sb[j][k][1],-sb[j][k][2]))
                    f.writelines('%5i %5i %5i %5i'  %(kount+1, nsmx, nrsx, iopsysx))
                    f.writelines('              modex, nsmx, nrsx, oiopsysx \n')
                    f.writelines(' %5.3f   %1i     %5.3f   %5.3f ' %(twshx,isectw,thres1,thres2))
                    f.writelines('        twshx,isectw,thres1,thres2   \n')
                    f.writelines(' %5.3e %5.3e %5.3e %5.3e %5.1f %5.1f \n' %(tau0, tau1, thet0, thet1, hpfac, gndfac))
                    #hlat
                    count = 0
                    for i in range(len(sn)*len(sb[0])*m):
                        if kount  == count:
                            f.writelines(' %5.2f ' %(hii))
                        elif kount - 1== count:
                            f.writelines(' %5.2f ' %(hb))
                        else:
                            f.writelines(' %5.2f ' %(hij))
                        count = count + 1
                    f.writelines('\n')
                    #index
                    f.writelines(' %2i  %2i  %2i     %2i  %2i  %2i \n' %(sn[j][0], sn[j][1], sn[j][2],
                                 -sb[j][k][0], -sb[j][k][1], -sb[j][k][2]))
                    kount = kount + 1

                    
        if iopsysx == 1:      #No directionality in slip direction
            pass
        elif iopsysx == 0:    #Slip system has the directionality
            # modes = len(sn) * 2
            pass
        else :
            print
            print 'Err: Wrong dependency!'
            print
            return -1
    f.close()

    print
    print '***********************************************************'
    print 'Creation of the requested single crystal file has completed'
    print '***********************************************************'
    print




def Hij(interaction, iopsysx, sn, sb, j, k):
    """
    given the arguments; interaction, iopsysx, sn, sb
    Returns the appropriate hardening matrix


    -- youngung comment:(2010-12-02)
    Disposable? or Not?
    Am I going to further work on this? (probably not)
    """
    if iopsysx == 0: m = 2
    elif iopsysx == 1: m = 1
    if interaction =='Bauschinger':
        pass
    elif interaction =='Latent':
        pass
    elif interaction =='Bauschinger Latent':
        pass
    
    return hmat

