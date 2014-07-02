import numpy as np
import multiprocessing as mp
from solve_flows import solve_flow
from FlowCalculation import FlowCalculation

modes = ['lin', 'sqr', 'lin_imp', 'sqr_imp']
layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']

alphas = [''.join(['aHO', str(a)]) for a in np.linspace(0, 1, 21)]

fc_list = []
for m in modes:
    for l in layouts:
        for a in alphas:
            fc = FlowCalculation(l, a, 'copper', m,\
                    savemode='FCResult', hourly_flowhist=True)
            fc_list.append(fc)


#for f in [str(fc) for fc in fc_list]:
#    print f
#print [fc.hourly_flowhist for fc in fc_list]
pool = mp.Pool(mp.cpu_count())

pool.map(solve_flow, fc_list)
