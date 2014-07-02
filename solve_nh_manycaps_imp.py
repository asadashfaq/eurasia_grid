import numpy as np
import multiprocessing as mp
from solve_flows import solve_flow
from FlowCalculation import FlowCalculation

modes = ['lin', 'sqr', 'lin_imp', 'sqr_imp']
layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']

capacities = [''.join([str(a), 'q99']) for a in np.linspace(0,1.5,31)]

fc_list = []
for m in modes:
    for l in layouts:
        for c in capacities:
            fc = FlowCalculation(l, 'aHE', c, m, savemode='FCResult',\
                    basisnetwork = 'nh')

            fc_list.append(fc)

#print [str(fc) for fc in fc_list]
pool = mp.Pool(mp.cpu_count())
pool.map(solve_flow, fc_list)
