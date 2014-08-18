import multiprocessing as mp
import numpy as np
from solve_flows import solve_flow
from FlowCalculation import FlowCalculation

modes = ['lin', 'sqr']
layouts = ['Europe', 'Europe_RU', 'Europe_NA', 'Europe_ME', 'Europe_RU_NA_ME']

fc_list = []
for m in modes:
    for l in layouts:
        fc_list.append(FlowCalculation(l, 'aHE', 'copper', m, \
                basisnetwork='europeplus', savemode='FCResult'))

print [str(fc) for fc in fc_list]

pool = mp.Pool(mp.cpu_count())
pool.map(solve_flow, fc_list)
