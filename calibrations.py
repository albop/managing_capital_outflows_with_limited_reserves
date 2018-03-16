
import copy
from collections import OrderedDict

params = dict(
    beta=0.8/(0.8+0.15),
    a=0.8,
    c=0.15,
    estar=-0.0,
    Rbar=0.5,
    min_f=0,
    kappa=1.0,
    N=40,
    zbar=0.1,
    p=1,
    model='optimal'
)


reserve_levels = [0.01, 0.5, 1.0, 1.5, 2.0, 3.0]
policies = ['volume','peg','optimal','time-consistent']
cases = {
    'baseline': {},
    'accumulation': {'min_f': -10000},
    'low_beta': {'beta': 0.8},
    'high_beta': {'beta': 0.9},
    'super_low_a': {'a': 0.01},
    'low_a': {'a': 0.4},
    'high_a': {'a': 1.6},
    'low_c': {'c': 0.075},
    'high_c': {'c': 0.30},
    'p_8': {'p': 0.8},
    'p_85': {'p': 0.85},
    'p_9': {'p': 0.9},
    'p_95': {'p': 0.95},
}

list_of_calibrations = OrderedDict()
for case in cases.keys():
    for pol in policies:
        for r in reserve_levels:
            calib = copy.copy(params)
            calib.update(cases[case])
            calib['Rbar'] = r
            calib['model'] = pol
            list_of_calibrations[(case,pol,r)] = calib.copy()
