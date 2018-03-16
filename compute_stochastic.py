from bttt import *
from bttt.trees import *
from bttt.model import import_tree_model
from matplotlib import pyplot as plt



class PersistentDeathTree(EventTree):

    def __init__(self, N, p, k=5):
        '''
        p: [p_0, p_1]
        '''

        self.nodes = [(0,)]
        nn = self.nodes[-1]
        self.probas = dict()
        for t in range(0, N-1):
            nnb = nn
            for j in range(k):
                if (t==N-2) or (j>0):
                    pb = 1
                else:
                    pb = p
                nnb1 = nnb + (1,)
                self.nodes.append(nnb1)
                self.probas[(nnb,nnb1)] = pb
                nnb = nnb1
            self.probas[(nnb,nnb)] = 1.0

            # create next point:
            if t<N-2:
                nn0 = nn + (0,)
                self.nodes.append(nn0)
                self.probas[(nn, nn0)] = 1-p
                nn = nn0
            else:
                self.probas[(nn,nn)] = 1.0

        n = self.nodes[-3]
        self.probas[(n,n)] = (1-p)


T = 40

p = 0.2
etree = PersistentDeathTree(p=p, N=T, k=2)
for s in etree.nodes:
    if sum(s)==0:
        etree.values[s] = 0.1
    else:
        etree.values[s] = 0.0

def get_ts(etree, sol, vname, terminal_state=None, ts_ind=0):
    import numpy
    if terminal_state is None:
        terminal_states = [e for e in etree.nodes if len(e)==len(etree)]
        terminal_state = terminal_states[ts_ind]
    history = etree.history(terminal_state)
    his_inds = [etree.nodes.index(e) for e in history]
    vals = [sol["{}_{}".format(vname, h)] for h in his_inds]
    return numpy.array(vals).astype(dtype=float)





model = import_tree_model("models.yaml", tree=etree, key="stochastic")
model.calibration["Rbar"] = 1.0
model.calibration["beta"]=0.8/(0.8+0.15)
model.calibration["a"]=0.8
model.calibration["c"]=0.15

sol = model.solve(verbose=True, linear=True, solver_options={"presteps": 2, 'tol': 1e-5})

from numpy import *
from matplotlib import pyplot as mplt
from pandas import DataFrame
import numpy

# ongoing series
z = get_ts(etree, sol, 'z')
f = get_ts(etree, sol, 'f')
e = get_ts(etree, sol, 'e')
Gamma = get_ts(etree, sol, 'Gamma')
R = get_ts(etree, sol, 'R')
tt = linspace(0,len(f)-1,len(f))


from collections import OrderedDict
data = OrderedDict([
    ('t', tt),
    ('z',z),
    ('f',f),
    ('e',e),
    ('Gamma',Gamma),
    ('R',R)

])
df = DataFrame(data).set_index('t')
df

jsnodes = [etree.nodes.index(n) for n in etree.nodes if sum(n)==1]
zjs, fjs, ejs, Rjs, gammajs = [array([sol['{}_{}'.format(k,ind)] for ind in jsnodes]) for k in ['z','f','e','R','Gamma'] ]
ttjs = tt[1:]
ttjs = ttjs[:len(zjs)]


data_stopped = OrderedDict([
    ('t', ttjs),
    ('z_s',zjs),
    ('f_s',fjs),
    ('e_s',ejs),
    ('Gamma_s',gammajs),
    ('R_s',Rjs)
])
df_stopped = DataFrame(data_stopped).set_index('t')


df = df.join(df_stopped)


prob = (1-p)**tt # probability that the shock continues at t.
E_gamma = df['Gamma']*prob + df['Gamma_s']*prob/(1-p)*p
df['E_gamma'] = df['Gamma']*prob + df['Gamma_s']*prob/(1-p)*p
df['E_e'] = df['e']*prob + df['e_s']*prob/(1-p)*p
df['E_f'] = df['f']*prob + df['f_s']*prob/(1-p)*p


import pickle
with open("precomputed_sims_stoch","wb") as f:
    pickle.dump(df,f)
