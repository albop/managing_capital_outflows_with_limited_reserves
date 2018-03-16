from matplotlib import pyplot as plt
from calibrations import *

import numpy

calib = list_of_calibrations[('baseline', 'optimal',1.0)]
calib['a']
calib['c']
calib['beta']
calib['a']/(calib['a']+calib['c'])
T = 20
calib['N'] = T



def perfect_foresight(calib):

    from calibrations import list_of_calibrations
    calibration = calib.copy()

    from bttt.trees import DeterministicTree, get_ts
    import numpy
    import pandas
    import copy

    from bttt.model import import_tree_model

    # R0 = calibration['Rbar']
    N = calibration['N']
    shock = calibration['zbar']
    model = calibration['model']

    # cc = copy.copy(calibration)
    # cc['Rbar'] = R0

    # solve deterministic case
    etree = DeterministicTree(N)
    for s in etree.nodes:
        etree.values[s] = shock

    model = import_tree_model('models.yaml', key=model, tree=etree)
    model.calibration.update(**calib)
    sol = model.solve( verbose=True, solver_options=dict(eps2=1e-6))

    columns = ['Gamma', 'e', 'f']
    tab  = [get_ts(etree, sol, c) for c in columns]
    tab = numpy.row_stack(tab).T

    df = pandas.DataFrame(tab, columns=columns)

    return df


def get_welfare(sol, params):
    beta = params['beta']
    estar = params['estar']
    p = params['p']
    res = sol['e']**2/2.0
    N = len(res)
    discount = (p*beta)**numpy.linspace(0,N,N)
    return -sum(discount*res)

import numpy

filenames = ['simple_rules_welfares.xlsx','simple_rules_welfares_9.xlsx','simple_rules_welfares_8.xlsx']
pvalues = [1,0.9,0.8]


for i,p in enumerate(pvalues):

    calib['p'] = p
    filename = filenames[i]

    if i==0:
        kvec = numpy.linspace(0.5, 1.0, 50)
    else:
        kvec = numpy.linspace(0.6, 1.0, 50)

    # from models import VolumeModel, PegModel, ConstrainedModel
    import numpy
    import pandas
    import pickle

    from time_consistent import solve as solve_time_consistent
    from time_consistent import simulate

    calib['Rbar'] = 1.0

    sol_dr_tc = solve_time_consistent(max_R=3, order=100, **calib)
    sol_tc = simulate(calib['Rbar'], sol_dr_tc, T)


    calib['model'] = 'optimal'
    commitment_solution = perfect_foresight(calib)

    calibs =  [calib.copy() for k in kvec]
    for i,cc in enumerate(calibs):
        cc['kappa'] = kvec[i]
        cc['model'] = 'volume'
    volume_solutions = []
    for cc in calibs:
        print(cc)
        volume_solutions.append(perfect_foresight(cc))


    calibs =  [calib.copy() for k in kvec]
    for i,cc in enumerate(calibs):
        cc['kappa'] = kvec[i]
        cc['model'] = 'peg'
    peg_solutions = []
    for cc in calibs:
        print(cc)
        peg_solutions.append(perfect_foresight(cc))

    do_nothing = commitment_solution.copy()
    do_nothing['e'] = 0.1/calib['c']

    tc_solutions = [sol_tc]*len(volume_solutions)

    optimal_welfares = [get_welfare(commitment_solution, calib) for k in kvec]
    donothing_welfares = [get_welfare(do_nothing, calib) for k in kvec]
    volume_welfares = [get_welfare(sol, calib) for sol in volume_solutions]
    peg_welfares = [get_welfare(sol, calib) for sol in peg_solutions]
    tc_welfares = [get_welfare(sol, calib) for sol in tc_solutions]

    df_welfares = pandas.DataFrame(
        numpy.row_stack([
                optimal_welfares,
                tc_welfares,
                volume_welfares,
                peg_welfares,
                donothing_welfares
        ]).T,
        columns=['commitment', 'time-consistent', 'volume', 'peg', 'do_nothing'],
        index=kvec
    )

    df_optimal = pandas.DataFrame(numpy.row_stack([commitment_solution['e'], sol_tc['e'][:len(commitment_solution)], do_nothing['e']]).T, columns=['commitment', 'time-consistent', 'do_nothing'])
    df_peg = pandas.DataFrame(
        numpy.row_stack([p['e'] for p in peg_solutions]).T, columns=['{}'.format(k) for k in kvec])
    df_volume = pandas.DataFrame(
        numpy.row_stack([p['e'] for p in volume_solutions]).T, columns=['{}'.format(k) for k in kvec])
    writer = pandas.ExcelWriter(filename, engine='xlsxwriter')
    df_welfares.to_excel(writer, sheet_name="welfares")
    df_optimal.to_excel(writer, sheet_name="optimal rules")
    df_peg.to_excel(writer, sheet_name="peg")
    df_volume.to_excel(writer, sheet_name="volume")
    writer.close()
