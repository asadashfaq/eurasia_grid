import multiprocessing as mp
import numpy as np
from solve_flows import solve_flow
from FlowCalculation import FlowCalculation

modes = ['lin', 'sqr', 'lin_imp', 'sqr_imp']
layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']

fc_list = []
for m in modes:
    for l in layouts:
        fc_list.append(FlowCalculation(l, 'aHE', 'copper', m,\
                basisnetwork='nh'))


pool = mp.Pool(mp.cpu_count())
pool.map(solve_flow, fc_list)
