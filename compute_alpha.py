######## Solve
import pandas
import numpy
from collections import OrderedDict
# from solution import *
from calibrations import *
from time_consistent import solve as solve_time_consistent
from time_consistent import simulate
from bttt.trees import DeterministicTree, get_ts
from matplotlib import pyplot as plt
from bttt.trees import DeathTree
from bttt.model import import_tree_model
from matplotlib import pyplot as plt
from numpy import *



N = 25
T = 25


calib0 = calib.copy()
beta = calib0['beta']
a = calib0['a']
p = calib0['p']
zbar = calib0['zbar']
max_R=5


tree = DeterministicTree(N)
for s in tree.nodes:
    tree.values[s] = zbar
    model = import_tree_model('models.yaml', key='stochastic', tree=tree)

def solve_it(**cc):
    tree = DeterministicTree(N)
    for s in tree.nodes:
        tree.values[s] = zbar
    model = import_tree_model('models.yaml', key='optimal', tree=tree)
    model.calibration.update(cc)
    sol = model.solve(verbose=True)
    df = numpy.concatenate( [get_ts(tree, sol, varname)[:,None] for varname in ['e','f', 'Gamma']], axis=1 )
    df = pandas.DataFrame(df, columns=['e','f','Gamma'])
    return df

from collections import OrderedDict


Rvec0 = linspace(0.01, 3.0, 20)



pvec = [0.8, 0.85, 0.9, 0.95, 1.0]

# Unconstrained solutions

unconstrained_sols = OrderedDict()
for p in pvec:
    dfs = [solve_it(Rbar=i, min_f=0, p=p) for i in Rvec0]
    unconstrained_sols[p] = dfs

all_gammas = OrderedDict()
for p in pvec:
    dfs = unconstrained_sols[p]
    max_gammas = numpy.array([dfs[i]['Gamma'].max() for i,r in enumerate(Rvec0)])
    all_gammas[p] = max_gammas

for p in pvec:
    plt.plot(Rvec0, all_gammas[p])

alpha = 0.5

Rlimit = OrderedDict()

for p in pvec:
    dfs = unconstrained_sols[p]
    max_gammas = numpy.array([dfs[i]['Gamma'].max() for i,r in enumerate(Rvec0)])
    j = numpy.where(max_gammas<alpha)[0][0]
    a0 = max_gammas[j-1]
    a1 = max_gammas[j]
    r0 = Rvec0[j-1]
    r1 = Rvec0[j]
    rmax = r0 + (r1-r0)*(alpha-a0)/(a1-a0)
    Rlimit[p] = rmax

# Constrained solutions

R0 = 0.8

constrained_solutions = OrderedDict()

for p in pvec:
    Rm = Rlimit[p]
    R1 = min([R0, Rm])
    if Rm>0:
        df = solve_it(p=p, Rbar=R1)
    elif Rm<=0:
        df = solve_it(p=p, Rbar=0.001)
    df['R'] = R0 - df['f'].cumsum().shift()
    df['R'][0]  = R0
    constrained_solutions[p] = df

d = {'dfs':constrained_solutions}

import pickle
with open("precomputed_alpha_p.pickle",'wb') as f:
    pickle.dump(d,f)

###
##
###


Rmax = 0.998
nRvec0 = numpy.array([0.01, 0.5, 0.9, Rmax, 2.0, 3.0])
ndfs = [solve_it(Rbar=i, min_f=0) for i in nRvec0]
for i in [4,5]:
    df = ndfs[3].copy()
    ndfs[i] = df
for i,df in enumerate(ndfs):
    df['R'] = nRvec0[i] - df['f'].cumsum()
#
#
# # In[43]:
#
#
# def plot_dfs(dfs, labels):
#     attributes = ['b', 'g', 'r']
#     if not isinstance(dfs, list):
#         dfs = [dfs]
#     fig = plt.figure(figsize=(15,10))
# #     plt.clear()
#     plt.subplot(131)
#     plt.plot(dfs[0].index,dfs[0].index*0+0.5,color='black', linestyle='--')
#     for i,df in enumerate(dfs):
#         plt.plot(df["Gamma"])
#
#     plt.grid()
#     plt.title("Marginal Value of Intervention")
# #     plt.figure()
#     plt.subplot(132)
#     for i,df in enumerate(dfs):
#         plt.plot(df["e"], label=labels[i])
#     plt.legend(loc='lower right')
#     plt.grid()
#     plt.title("Exchange Rate")
#     plt.subplot(133)
#     for i,df in enumerate(dfs):
#         plt.plot(df["R"], label=labels[i])
#     plt.grid()
#     plt.title("Reserves")
#     return fig
#
#
# # In[44]:


import pickle
with open("precomputed_option_value.pickle","wb") as f:
    pickle.dump({"commitment":ndfs},f)

#
# # In[82]:
#
#
# f = plot_dfs(ndfs, labels=['$R_0={}$'.format(i) for i in Rvec0])
# plt.savefig("optimal_choice_noconstraint.png", bbox_inches="tight")
# plt.savefig("optimal_choice_noconstraint.pdf", bbox_inches="tight")
# f
#
#
# # In[11]:
#
#
# Rvec0 = linspace(0.3, 3.0, 10)
# dfs = [solve_it(Rbar=i) for i in Rvec0]
#
#
# # In[12]:
#
#
# f = plot_dfs(dfs, labels=['$R_0={}$'.format(i) for i in Rvec0])
# plt.savefig("optimal_choice_constraint.png", bbox_inches="tight")
# plt.savefig("optimal_choice_constraint.pdf", bbox_inches="tight")
# f
