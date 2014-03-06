import numpy as np
import aurespf.solvers as au
from worldgrid import world_Nodes
import FlowCalculation as fc

def extract_BC_vs_TC_data(solvermode, TCscalefactors, savefilename=None):
    """ This function returns the a dataset with total balancing
        capacity (normalized to the total mean load) of the world
        along with the total transmission capacity of the world.
        The balancing capacity found as the maximum balancing over
        a timeseries, as well as the 99% quantile is returned.

        """


    h0 = au.get_quant_caps(filename =\
               ''.join(['./results/w_aHE_copper_', solvermode, '_flows.npy']))
    total_TC_q99 = np.sum(au.biggestpair(h0))

    BC_max = []
    BC_q99 = []
    TC = []

    for a in TCscalefactors:
        flow_calc = fc.FlowCalculation('w', 'aHE',\
                                       ''.join([str(a), 'q99']), solvermode)
        filename = ''.join([str(flow_calc), '.npz'])
        nodes = world_Nodes(load_filename = filename)

        total_mean_load = np.sum([n.mean for n in nodes])

        BC_max.append(np.sum(np.max(n.balancing) for n in nodes)\
                            /total_mean_load)
        BC_q99.append(np.sum(au.get_q(n.balancing, 0.99) for n in nodes)\
                                                        /total_mean_load)

        TC.append(a*total_TC_q99)

    if not savefilename:
        savefilename = ''.join(['./results/BC_TC_data', solvermode, '.npz'])

    np.savez(savefilename, TC=TC, BC_max=BC_max, BC_q99=BC_q99)

