"""
STRESS AT THE A CERTAIN PLASTIC WORK.
INTERPOLATE THE DISCRETE RESULTS FROM THE 'str_str_r.py' PYTHON SCRIPT
"""
import glob

files = glob.glob('*.str')

print
print '** LIST OF FILES'
print


def work_out(wrk, sig, w, filename = 'tmp.wk'):
    """
    work and stress columns are wrk, sig
    and a certain level of work is given as w
    Arguments:
       wrk
       sig
       w
    """
    fout = open(filename.split('.')[0]+'.wk', 'w')
    fout.writelines('stress w.r.t plastic work estimated from the post-processed str_str_r file \n')
    fout.writelines('pl_work     stress \n')

    #print 'len(wrk)=', len(wrk)
    #print 'len(sig)=', len(sig)
    #print 'len(w)  =', len(w)
    stress = []
    for i in range(len(w)):
        j = 0
        while True:
            if w[i]<wrk[j]:
                #print 'j, wrk[j], sig[j] =', j, wrk[j], sig[j]
                wrk1 = wrk[j-1]
                wrk2 = wrk[j]
                sig1 = sig[j-1]
                sig2 = sig[j]
                slope = (sig2 - sig1)/(wrk2 - wrk1)
                stress.append(slope * (w[i] - wrk1) + sig1)
                fout.writelines(str(w[i]).rjust(4) + '        ' + str(stress[i]).ljust(3) + '\n')
                break
            j = j + 1
            
    
    fout.close()
    return stress

###
###  Work level input
###   : default work levels are [1,2,4,6,8,10,12,14]

w = []
while True:
    tmp = raw_input('work level >>')
    if len(tmp)==0: break
    else: w.append(float(tmp))
if len(w)==0: w=[1,2,4,6,8,10,12,14]


###
###  work and stress are designated to be in a certain column id
###   id_wrk : 3
###   id_sig : 2
###

id_wrk = 3
id_sig = 2
f = []
sigma = []
for i in range(len(files)):
    work, stress = [], []
    sigma.append([])
    print files[i]
    f.append(open(files[i], 'r'))
    lines = f[i].readlines()
    for j in range(len(lines)-2):
        work.append(float(lines[j+2].split()[id_wrk]))
        stress.append(float(lines[j+2].split()[id_sig]))

    f[i].close()
    sigma[i].append(work_out(wrk = work, sig = stress, w = w, filename = files[i]))

"""
print 'sigma[0]', sigma[0]
print 'sigma[0][0]', sigma[0][0]
print 'sigma[1]', sigma[1]
print 'len(sigma)', len(sigma)
print 'len(sigma[0])', len(sigma[0])
print 'len(sigma[0][0])', len(sigma[0][0])
print 'len(w)', len(w)
print 'len(files)', len(files)
"""
answer = [ ]
for i in range(len(w)):
    f = open('w_'+str(w[i])+'.eqw', 'w')
    f.writelines('pl_work         stress\n')

    for j in range(len(files)):
        print '************'
        print files[j]
        if i == 0:
            iRD = raw_input('RD(y) or TD(n) >> ')            
            if iRD == 'y':
                f.writelines(files[j].split('.')[0] +  '  0  0  1  0  '  + str(sigma[j][0][i])[0:6] + '  0 \n')
            elif iRD =='n':
                f.writelines(files[j].split('.')[0] +  '  0  0  0  1        0  '  + str(sigma[j][0][i])[0:6]+ '\n')
            else:
                print 'Wrong Answer you must choose either y or n'
                break
            answer.append(iRD)
        else:
            if answer[j] == 'y': f.writelines(files[j].split('.')[0] +  '  0  0  1  0  '  + str(sigma[j][0][i])[0:6] + '  0 \n')
            elif answer[j] == 'n': f.writelines(files[j].split('.')[0] +  '  0  0  0  1     0  '  + str(sigma[j][0][i])[0:6]+ '\n')
            else :
                print 'Something wrong!!'
                break

    f.close()








        
