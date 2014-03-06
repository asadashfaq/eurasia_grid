import numpy as np
import multiprocessing as mp

import aurespf.solvers as au
from worldgrid import world_Nodes

modes = ['capped rolando copper', 'capped martin copper']

def solve_copper(mode):
    admat = './settings/wadmat.txt'
    nodes = world_Nodes(admat=admat)
    solved_nodes, flows = au.solve(nodes, mode=mode)

    if 'rolando' in mode:
        shortmodename = 'capR'
    if 'martin' in mode:
        shortmodename = 'capM'
    filename = ''.join(['w_aHE_copper_', shortmodename])
    solved_nodes.save_nodes(filename)
    np.save('./results/' + filename + '_flows', flows)

pool = mp.Pool(mp.cpu_count())

pool.map(solve_copper, modes)



