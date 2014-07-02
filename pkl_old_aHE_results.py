import numpy as np
import multiprocessing as mp
from solve_flows import solve_flow
from FlowCalculation import FlowCalculation
from FCResult import FCResult
from nhgrid import nh_Nodes

def convert_to_FCResult_file(fc):
    pkl_filename = str(fc) + '.pkl'
    npz_filename = str(fc) + '.npz'
    flow_filename = str(fc) + '_flows.npy'

    N = nh_Nodes(load_filename=npz_filename)
    F = np.load('./results/'+flow_filename)

    myresult = FCResult(pkl_filename)
    myresult.add_instance(N, F, fc)
    myresult.save_results(pkl_filename)

    return

modes = ['lin', 'sqr', 'lin_imp', 'sqr_imp']
layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']

fc_list = []
for m in modes:
    for l in layouts:
        fc_list.append(FlowCalculation(l, 'aHE', 'copper', m,\
                savemode='FCResult', hourly_flowhist=True))


pool = mp.Pool(mp.cpu_count())

pool.map(convert_to_FCResult_file, fc_list)
