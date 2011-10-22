"""
A handy program for disk compression test
"""


# variables or parameters of interests
"""
1. diameter
   Actually two diameters are of concern
    one along RD the other along TD
    Under assumption that the oval shape is maintained,
    the cross section area can be calculated
2. Load
    Load to which the disk is subjected will be mostly given.
     Note that there's a certain limit of load by the machine.
3. stress
    Stress is actually a post-processed value from load and diamter.
    But at the same time if stress is given along with one of the load or diameter,
    the remainder is supposed to be estimated by this program.




Unit system:
   load : kN
   stress: MPa
   diamter: mm
   area: mm^2
"""

pi = 3.141592

def __stress__(diameter = None, load = None):
    """
    stress will be calculated based on the given information.
    """
    try:
        load
    except:
        print 'Load is not given? you should enter load'
        return -1
    else:
        load = load *  (10**3)   # kN to N
    try:
        diameter[0]*diameter[1]/4. * 3.141592
    except:
        print 'diamter must be given as a list containing two diamters along RD and TD'
        return -1
    else:
        d0 = diameter[0]*(10**-3)   # mm to m
        d1 = diameter[1]*(10**-3)
        print d0,d1
        area = d0*d1/4. * 3.141592
    print 'Stress [MPa] : ', '%10.3f'%(load/area/10**6)
    return load / area / (10**6)

while True:
    d = []
    for i in range(2):
        tmp = raw_input('Diameters in [mm] >>')
        if tmp =='': break
        d.append(float(tmp))
    if tmp =='': break
    l = float(raw_input('Load in kN >>'))
    __stress__(diameter = d, load = l)
    raw_input()
