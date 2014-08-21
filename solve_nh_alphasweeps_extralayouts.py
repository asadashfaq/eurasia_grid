import numpy as np
import multiprocessing as mp
from solve_flows import solve_flow
from FlowCalculation import FlowCalculation

modes = ['lin', 'sqr']
layouts = ['EU_RU', 'EU_NA', 'EU_ME', 'EU_RU_NA_MEstar']

alphas = [''.join(['aHO', str(a)]) for a in np.linspace(0, 1, 21)]

fc_list = []
for m in modes:
    for l in layouts:
        for a in alphas:
            fc = FlowCalculation(l, a, 'copper', m,\
                    savemode='FCResult', hourly_flowhist=False)
            fc_list.append(fc)

#print [str(fc) for fc in fc_list]
pool = mp.Pool(mp.cpu_count())

pool.map(solve_flow, fc_list)
