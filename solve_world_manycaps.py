import numpy as np
import multiprocessing as mp
from solve_flows import solve_flow
import FlowCalculation as fc

modes = ['capR', 'capM']
capacities = [''.join([str(a), 'q99']) for a in np.linspace(0, 2, 41)]

fc_list = []

for m in modes:
    for c in capacities:
        fc_list.append(fc.FlowCalculation('w', 'aHE', c, m, basisnetwork='w'))
        print fc_list[-1]

pool = mp.Pool(mp.cpu_count())
pool.map(solve_flow, fc_list)
