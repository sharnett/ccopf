from scipy.sparse import lil_matrix, csr_matrix 
from scipy.sparse.linalg import factorized
from numpy import zeros, append

def load(conf, data_dir):
    lines = open(data_dir + '/' + conf).readlines()
    casefile, windfile, costfile = [l.split()[1] for l in lines[0:3]]
    line_nu, gen_nu = [float(l.split()[1]) for l in lines[4:6]]
    bus, branch, Bbus = loadcase(data_dir + '/' + casefile)
    loadwind(data_dir + '/' + windfile, bus)
    Binv = getBinv(Bbus, bus)
    loadcost(data_dir + '/' + costfile, bus)
    return bus, branch, Bbus, Binv, line_nu, gen_nu

def loadcase(casefile):
    lines = open(casefile).readlines()
    i = lines.index('mpc.bus = [\n')+1
    n, bus, branch = 0, {}, []
    while lines[i] != '];\n':
        words = lines[i].split()
        busid = int(words[0])-1
        bustype = int(words[1])
        load = float(words[2])
        bus[busid] = {'type': bustype, 'd': load}
        n += 1
        i += 1
    i = lines.index('mpc.gen = [\n')+1
    while lines[i] != '];\n':
        words = lines[i].split()
        busid = int(words[0])-1
        pmax, pmin = [float(w) for w in words[8:10]]
        bus[busid]['pmax'] = bus[busid].get('pmax', 0) + pmax
        bus[busid]['pmin'] = bus[busid].get('pmin', 0) + pmin
        i += 1
    for k,v in bus.iteritems():
        if 'pmax' in v and v['pmax'] <= 0:
            v['type'] = 1
            del v['pmax']
            del v['pmin']
    Bbus = lil_matrix((n,n))
    i = lines.index('mpc.branch = [\n')+1
    nl = 0
    while lines[i] != '];\n':
        nl += 1
        words = lines[i].split()
        fbus, tbus = [int(w)-1 for w in words[0:2]]
        y, fmax = 1/float(words[3]), float(words[5])
        #branch[(fbus, tbus)] = {'y': y, 'fmax': fmax, 'id': nl}
        branch += [{'ij': (fbus,tbus), 'y': y, 'fmax': fmax, 'id': nl}]
        Bbus[fbus,tbus] = Bbus[fbus,tbus]-y;
        Bbus[tbus,fbus] = Bbus[tbus,fbus]-y;
        i += 1
    Bbus = csr_matrix(Bbus)
    for i in range(n): 
        Bbus[i,i] = -Bbus[i,:].sum()
    return bus, branch, Bbus

def loadwind(windfile, bus):
    lines = open(windfile).readlines()
    nw = int(lines[0])
    for i in range(1, nw+1):
        words = lines[i].split()
        busid = int(words[0])-1
        avg, std = [float(w) for w in words[1:3]]
        bus[busid]['avg'] = avg
        bus[busid]['std'] = std

def loadcost(costfile, bus):
    if costfile.find('none') != -1:
        print 'no cost file indicated, using all ones'
        gens = {b for b in bus if bus[b]['type'] == 2}
        for b in gens:
            bus[b]['cost'] = 1.0
        return
    lines = open(costfile).readlines()
    i = 0
    while lines[i].find('END') == -1:
        words = lines[i].split()
        busid = int(words[0])-1
        cost = float(words[1])
        bus[busid]['cost'] = cost
        i += 1

def getBinv(Bbus, bus):
    n = Bbus.shape[0]
    rhs = zeros(n-1)
    solve = factorized(Bbus[:-1,:-1])
    Binv = {} # dictionary with a key per wind bus, value is column array
    for i in bus:
        if not 'avg' in bus[i]: continue # only compute column if wind bus
        rhs[i] = 1
        Binv[i] = append(solve(rhs), 0)
        rhs[i] = 0
    return Binv
