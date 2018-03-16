######## Solve
import pandas
import numpy
from collections import OrderedDict
from calibrations import *
from time_consistent import solve as solve_time_consistent
from time_consistent import simulate
import bttt
from bttt.trees import DeterministicTree, get_ts
from bttt.trees import DeathTree
from bttt.model import import_tree_model



N = 25
T = 25



model_types =['optimal','peg','volume','time-consistent']


results = OrderedDict()

for c,v in list_of_calibrations.items():
    calib = v.copy()
    beta = calib['beta']
    a = calib['a']
    p = calib['p']
    zbar = calib['zbar']
    max_R=5
    if c[1]=='time-consistent':

        if beta<0.82:
            max_R=4
        else:
            max_R=4
        r = calib['Rbar']
        # sol = solve_time_consistent(max_R=2, **v)
        sol = solve_time_consistent(max_R=max_R, order=100, **v)
        sim_tc = simulate(r, sol, T)
        results[c] = sim_tc
    else:
        print(p,c[1])
        # if p == 1:
        tree = DeterministicTree(N)
        for s in tree.nodes:
            tree.values[s] = zbar
        model = import_tree_model('models.yaml', key=c[1], tree=tree)
        model.calibration.update(calib)
        sol = model.solve(verbose=True)
        df = numpy.concatenate( [get_ts(tree, sol, varname)[:,None] for varname in ['e','f', 'Gamma']], axis=1 )
        df = pandas.DataFrame(df, columns=['e','f','Gamma'])
        results[c] = df



def Gamma(t, e, parm):
    tot = 0
    estar = parm['estar']
    beta = parm['beta']
    p = parm['p']
    alpha = parm['a']/(parm['a']+parm['c'])
    for s in range(t+1):
        tot+=(beta)**s*alpha**(t-s)*(e[s]-estar)
    return tot


list_of_calibrations
for c,sim in results.items():
    print(c)
    parm = list_of_calibrations[c]
    # add level of reserves
    Rbar = c[2]
    p = parm['p']

    sim['R'] = Rbar - sim['f'].cumsum().shift()
    sim['R'][0] = Rbar # - sim['f'].cumsum().shift()
    # add Gamma
    gg = [p**t*Gamma(t, sim['e'], parm) for t in range(T)]
    if 'Gamma' in sim.columns:
        diff = abs(sim['Gamma'] - gg).max()
        # print(diff)
        # if diff>1e-6:
            # raise Exception("Incorrect computation of Gamma")
    else:
        sim['Gamma'] = gg

#### save results
from dolo import groot

groot()

import pickle
with open("precomputed_simulations.pickle", 'wb') as f:
    pickle.dump({'results': results, 'calibrations': list_of_calibrations}, f)




#################333
####################


decision_rules = OrderedDict()

for p in [1.0, 0.9, 0.8]:

    v = list_of_calibrations[('baseline','optimal',1.0)].copy()
    v['p'] = p
    sol = solve_time_consistent(max_R=7, order=1000, verbose=True,**v)
    sol = sol[:3]
    decision_rules[p] = sol



import pickle
with open("precomputed_decision_rules.pickle",'wb') as f:
    pickle.dump({'decision_rules': decision_rules, 'calibration': v}, f)



#################
#################

from bttt.model import import_tree_model
from collections import OrderedDict

model = import_tree_model('models.yaml', key='moving_target', tree=tree)

v = list_of_calibrations[('baseline', 'optimal', 1.0)].copy()

od = OrderedDict()
for R in [0.01, 1]:
    lamvec = [0.8, 0.85, 0.9, 0.95, 1.0]
    for lam in lamvec:
        vv = v.copy()
        vv['lam'] = lam
        vv['Rbar'] = R
        model.calibration.update(vv)
        sol = model.solve()
        df = numpy.concatenate( [get_ts(tree, sol, varname)[:,None] for varname in ['e','f', 'Gamma', 'target']], axis=1 )
        df = pandas.DataFrame(df, columns=['e','f','Gamma','target'])
        od[(R,lam)] = df


import pickle
with open("precomputed_moving_target.pickle",'wb') as f:
    pickle.dump({'simulations': od, 'calibration': v}, f)
