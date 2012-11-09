# TODO
# separate honest-to-goodness DATA from PARAMETERS
# sort out the assemble/__init__/formulate mess
# should Bbus and Binv be considered auxiliary and not OG data?
import loadfiles 
from cplex import Cplex, SparsePair as sp
from gurobipy import Model, quicksum as qs, GRB
from numpy import array as arr, zeros
from numpy.linalg import norm
from scipy.sparse import coo_matrix
from sys import stderr

def solve(conf='conf9', M=None, data_dir='data', solver='cplex', verbose=True, maxits=35):
    ''' High-level look at the ccopf algorithm. Repeatedly solve the linear
    subproblems until either (1) a subproblem fails (2) ccopf converges, or 
    (3) the maximum number of iterations is reached. 
 
    Inputs (all optional):
    conf - configuration file
    M - model object. Reuse after first run to save on loading time.
    data_dir - directory for configuration files
    solver - cplex or gurobi
    verbose - set to False to suppress output
    maxits - maximum number of linear subproblems before giving up
 
    Outputs:
    p - generator production levels
    a - droop coefficients 
    obj - objective value
    M - model object. Use as input on subsequent runs to save on loading time '''
    if not M:
        M = CCOPF(conf=conf, data_dir=data_dir, solver=solver)
    M.formulate()
    for i in xrange(maxits):
        M.solve_subproblem()
        if M.status not in {'good', 'okay'}: break
        M.compute_max_viol(i)
        if M.converged: break
        M.add_cut(i)
    if not M.converged: 
        print >> stderr, 'Failed. {} status {}'.format(solver, M.statcode)
    return M.p, M.a, M.obj, M

class CCOPF:
    def __init__(self, conf='conf9', data_dir='data', solver='cplex', verbose=True):
        if verbose: print 'loading problem data'
        self.conf = conf
        self.solver = solver
        self.tol = 1e-5
        self.verbose = verbose
        bus, branch, Bbus, Binv, self.line_nu, self.gen_nu = loadfiles.load(
                self.conf, data_dir)
        self.data = [bus, branch, Bbus, Binv]

    def set_aux(self):
        ''' Helper variables built from the input data. Call this after
        changing the data; don't modify these directly '''
        bus, branch, Bbus, Binv = self.data 
        n = Bbus.shape[0] # number of buses
        nl = len(branch) # number of lines
        d = [bus[b]['d'] for b in bus] # demand
        mu = [bus[b].get('avg', 0) for b in bus] # mean wind power
        winds = {b for b in bus if 'std' in bus[b]} # wind source bus IDs
        sigma = {b: bus[b]['std'] for b in bus if b in winds} # std wind power
        gens = {b for b in bus if bus[b]['type'] == 2} # generator bus IDs
        self.aux = [n, nl, d, mu, winds, sigma, gens]

    def formulate(self):
        self.set_aux()
        self.M = Cplex() if self.solver == 'cplex' else Model()
        self.addVars()
        self.addCons()
        self.setObj()
        self.p, self.a, self.obj, self.converged, self.tot_its = None,None,None,False,0
        if self.verbose: 
            print 'solving with', self.solver
            print '%2s %10s %8s %4s %3s %10s' % ('it', 'worst line', 'viol', 'stat', 
                    'its', 'obj')

    def addVars(self):
        bus,branch,_,_, n,nl,_,_,_,_,gens = self.data + self.aux
        if self.verbose: print 'defining variables'
        INF = 1e100
        if self.solver == 'cplex':
            p = ['p_%d'%i for i in gens]
            a = ['a_%d'%i for i in gens]
            D = ['D_%d'%i for i in bus]
            t = ['t_%d'%i for i in bus]
            m = ['m{}'.format(i['id']) for i in branch] 
            s = ['s{}'.format(i['id']) for i in branch]
            self.M.variables.add(names = p + a)
            self.M.variables.add(names = D + t, lb = [-INF]*2*n)
            #self.M.variables.add(names = m, lb = [-INF]*nl)
            #self.M.variables.add(names = s)
            self.M.variables.add(names = m + s, lb = [-INF]*2*nl)
            D, t = arr(D), arr(t)
            self.var = (p, a, D, t, m, s)
        else:
            p = {i: self.M.addVar(name='pbar_%d'%i) for i in gens}
            a = {i: self.M.addVar(name='alpha_%d'%i) for i in gens}
            D = {i: self.M.addVar(lb=-INF, name='delta_%d'%i) for i in bus}
            t = {i: self.M.addVar(lb=-INF, name='theta_%d'%i) for i in bus}
            m = {i['id']: self.M.addVar(lb=-INF, name='fbar{}'.format(i['id'])) for 
                    i in branch}
            s = {i['id']: self.M.addVar(lb=-INF, name='std{}'.format(i['id'])) for 
                    i in branch}
            self.var = (p, a, D, t, m, s)
            self.M.update()

    def addCons(self):
        bus,branch,Bbus,_, line_nu,gen_nu, n,_,d,mu,winds,sigma,gens = self.data + \
                [self.line_nu, self.gen_nu] + self.aux
        if self.verbose: print 'defining constraints'
        p, a, D, t, m, s = self.var 
        sumvar = sum(sigma[b]**2 for b in winds)
        if self.solver == 'cplex':
            ng = len(gens)
            self.M.linear_constraints.add(lin_expr = [sp(ind = a, val = [1.0]*ng)], 
                    rhs = [1.0], names=['sum_alpha'])
            self.M.linear_constraints.add(lin_expr = [sp(ind = p, val = [1.0]*ng)], 
                    rhs = [sum(d)-sum(mu)], names=['power_balance'])
            for i in xrange(n-1):
                Bi = Bbus[i,:-1].tocoo()
                J,V = Bi.col, Bi.data
                if i in gens:
                    pi, ai, Di, ti = 'p_%d'%i, 'a_%d'%i, 'delta_%d'%i, 'theta_%d'%i
                    self.M.linear_constraints.add(lin_expr = 
                            [sp(ind = list(D[J])+[ai], val = list(V)+[-1]),
                            sp(ind = list(t[J])+[pi], val = list(V)+[-1])], 
                            rhs=[0, mu[i]-d[i]], names=[Di, ti], senses='EE')
                    self.M.linear_constraints.add(lin_expr = 
                            [sp(ind = [pi,ai], val = [1, sumvar**.5*gen_nu]), 
                            sp(ind = [pi,ai], val = [1, -sumvar**.5*gen_nu])], 
                            rhs=[bus[i]['pmax'], bus[i]['pmin']], 
                            names=['gen+_%d'%i, 'gen-_%d'%i], senses='LG')
                else:
                    self.M.linear_constraints.add(lin_expr = 
                            [sp(ind = D[J], val = V), sp(ind = t[J], val = V)], 
                            rhs=[0, mu[i]-d[i]], names=['delta_%d'%i, 'theta_%d'%i], 
                            senses='EE')
            self.M.linear_constraints.add(lin_expr = 
                    [sp(ind = [D[n-1]], val = [1.0]), sp(ind = [t[n-1]], val = [1.0])], 
                    rhs = [0,0], names=['delta_%d'%(n-1), 'theta_%d'%(n-1)], senses='EE')
            for line in branch:
                ID = line['id']
                ij = line['ij']
                i,j = ij
                mij, sij = 'm{}'.format(ID), 's{}'.format(ID)
                self.M.linear_constraints.add(lin_expr = 
                        [sp(ind = [mij, t[i], t[j]], val = [1.0, -line['y'], line['y']])], 
                        rhs = [0], names=['line_avg_{}'.format(ID)], senses='E')
                self.M.linear_constraints.add(lin_expr = 
                        [sp(ind = [mij, sij], val = [-1.0, line_nu]),
                        sp(ind = [mij, sij], val = [1.0, line_nu])], 
                        names=['line+_{}'.format(ID), 'line-_{}'.format(ID)],
                        rhs = [line['fmax'], line['fmax']], senses='LL')
        else: # gurobi is much nicer
            self.M.addConstr(qs(a[i] for i in gens) == 1, 'sum_alpha_1')
            self.M.addConstr(qs(p[i] for i in gens) + sum(mu) == sum(d), 'power_balance')
            for i in xrange(n-1):
                Bi = Bbus[i,:-1].tocoo()
                Bi = zip(Bi.col, Bi.data)
                if i in gens:
                    self.M.addConstr(qs(v*D[j] for j,v in Bi) == a[i], 'delta_%d'%i)
                    self.M.addConstr(qs(v*t[j] for j,v in Bi) == p[i]+mu[i]-d[i], 
                            'theta_%d'%i)
                    self.M.addConstr(sumvar**.5*a[i]*gen_nu <= bus[i]['pmax']-p[i], 
                            'gen+_%d'%i)
                    self.M.addConstr(sumvar**.5*a[i]*gen_nu <= -bus[i]['pmin']+p[i], 
                            'gen-_%d'%i)
                else:
                    self.M.addConstr(qs(v*D[j] for j,v in Bi) == 0, 'delta_%d'%i)
                    self.M.addConstr(qs(v*t[j] for j,v in Bi) == mu[i]-d[i], 'theta_%d'%i)
            self.M.addConstr(D[n-1] == 0, 'delta')
            self.M.addConstr(t[n-1] == 0, 'theta')
            for line in branch:
                ID = line['id']
                ij = line['ij']
                i,j = ij
                self.M.addConstr(m[ID] == line['y']*(t[i]-t[j]), 
                        'line_avg_{}'.format(ID))
                self.M.addConstr(s[ID]*line_nu <= line['fmax'] - m[ID], 
                        'line+_{}'.format(ID))
                self.M.addConstr(s[ID]*line_nu <= line['fmax'] + m[ID], 
                        'line-_{}'.format(ID))

    def setObj(self):
        bus,_,_,_, _,_,_,_,winds,sigma,gens = self.data + self.aux
        p, a, D, t, m, s = self.var 
        cost = {b: bus[b]['cost'] for b in bus if b in gens}
        sumvar = sum(sigma[b]**2 for b in winds)
        if self.solver == 'cplex':
            ng = len(gens)
            nv = self.M.variables.get_num()
            cost = [2*cost[i] for i in gens] + [2*sumvar*cost[i] for i in gens] + \
                    [0.0]*(nv-2*ng)
            self.M.objective.set_quadratic(cost)
            self.M.set_results_stream(None)
        else:
            self.M.setObjective(qs(cost[i]*(p[i]*p[i] + sumvar*a[i]*a[i]) for i in gens))

    def solve_subproblem(self):
        if self.solver == 'cplex':
            self.M.solve()
            sol = self.M.solution
            self.obj = sol.get_objective_value()
            stat = sol.get_status()
            if stat == sol.status.optimal: 
                self.status = 'good'
            elif stat == sol.status.num_best:
                self.status = 'okay'
            else:
                self.status = 'fail'
                self.statcode = stat
            self.its = sol.progress.get_num_barrier_iterations()
            self.sol = sol
        else:
            self.M.params.outputFlag = 0
            self.M.optimize()
            self.obj = None
            if self.M.status == GRB.OPTIMAL:
                self.status = 'good'
                self.obj = self.M.ObjVal
            elif self.M.status == GRB.SUBOPTIMAL:
                self.status = 'okay'
                self.obj = self.M.ObjVal
            else:
                self.status = 'fail'
                self.statcode = self.M.status
            self.its = self.M.BarIterCount
        self.tot_its += self.its

    def compute_max_viol(self, it):
        tol, _,branch,_,Binv, line_nu, n,_,_,_,winds,sigma,gens = [self.tol] + \
                self.data + [self.line_nu] + self.aux
        p, a, D, t, m, s = self.var 
        self.worst = {'line': -1, 'val': 0}
        for index, line in enumerate(branch):
            ID = line['id']
            i,j = line['ij']
            y, fmax = line['y'], line['fmax']
            if self.solver == 'cplex':
                Di, Dj, mij = self.sol.get_values([D[i], D[j], 'm{}'.format(ID)])
            else:
                Di, Dj, mij = D[i].x, D[j].x, m[ID].x
            x = arr([y*(Binv[k][i]-Binv[k][j]-Di+Dj)*sigma[k] for k in winds])
            violation = (abs(mij)+line_nu*norm(x)-fmax)/fmax
            if violation > self.worst['val']:
                self.worst = {'line': index, 'val': violation}
        if self.worst['val'] < tol:
            self.converged = True
            self.p, self.a = zeros(n), zeros(n)
            for i in gens:
                if self.solver == 'cplex':
                    pi, ai = 'p_%d'%i, 'a_%d'%i
                    self.p[i], self.a[i] = self.sol.get_values(pi), self.sol.get_values(ai)
                else:
                    self.p[i], self.a[i] = p[i].x, a[i].x
        if self.verbose: 
            i = self.worst['line']
            print '%2d %4d->%-4d %8.2e %4s %3d %10.4e' % (it+1, branch[i]['ij'][0], 
                    branch[i]['ij'][1], self.worst['val'], self.status, self.its, self.obj)

    def add_cut(self, it):
        _,branch,_,Binv, _,_,_,_,winds,sigma,_ = self.data + self.aux
        p, a, D, t, m, s = self.var 
        index = self.worst['line']
        ID = branch[index]['id']
        i,j = branch[index]['ij']
        y = branch[index]['y']
        if self.solver == 'cplex':
            Di, Dj = self.sol.get_values([D[i], D[j]])
            sij = 's{}'.format(ID)
            Bij = {k: Binv[k][i]-Binv[k][j] for k in winds}
            rhs = y**2*sum([sigma[k]**2*Bij[k]*(Bij[k]-Di+Dj) for k in winds])
            coeff = y**2*sum([sigma[k]**2*(Bij[k]-Di+Dj) for k in winds])
            x = arr([y*(Bij[k]-Di+Dj)*sigma[k] for k in winds])
            #self.M.linear_constraints.add(lin_expr = [sp(ind = [sij], val = [1.0])], 
                    #rhs = [0], names=['cut_help%d'%it], senses='G')
            self.M.linear_constraints.add(lin_expr = 
                    [sp(ind = [sij, D[i], D[j]], val = [norm(x), coeff, -coeff])], 
                    rhs = [rhs], names=['cut_%d'%it], senses='G')
        else:
            x = arr([y*(Binv[k][i]-Binv[k][j]-D[i].x+D[j].x)*sigma[k] for k in winds])
            self.M.addConstr(s[ID] >= 0, 'cut_%d_help'%(it+1))
            self.M.addConstr(qs(x[l]*y*(Binv[k][i]-Binv[k][j]-D[i]+D[j])*sigma[k]
                    for l,k in enumerate(winds)) <= s[ID]*norm(x), 'cut_%d'%it)

    def set_mean(self, mu):
        ''' mu is a dictionary {busID: mean wind} '''
        bus = self.data[0] 
        winds = {b for b in bus if 'avg' in bus[b]}
        if winds != set(mu.keys()):
            print 'error: mu is different'
            return None
        for i, mean in mu.iteritems():
            bus[i]['avg'] = mean
        self.set_aux()

    def set_std(self, sigma):
        ''' sigma is a dictionary {busID: std wind} '''
        bus = self.data[0] 
        winds = {b for b in bus if 'std' in bus[b]}
        if winds != set(sigma.keys()):
            print 'error: sigma is different'
            return None
        for i, std in sigma.iteritems():
            bus[i]['std'] = std
        self.set_aux()

    def scale_load(self, alpha):
        bus = self.data[0] 
        for i in bus:
            bus[i]['d'] *= alpha
        self.set_aux()

    def scale_lim(self, alpha):
        branch = self.data[1] 
        for line in branch:
            line['fmax'] *= alpha
        self.set_aux()
