import sys
import numpy as np
import multiprocessing as mp

import aurespf.solvers as au
from worldgrid import world_Nodes
from FlowCalculation import FlowCalculation # my own class for passing info about calculations

def solve_flow(flow_calc):
    admat = ''.join(['./settings/', flow_calc.layout, 'admat.txt'])
    filename = str(flow_calc)
    copperflow_filename = ''.join(['./results/', flow_calc.layout, '_',
        flow_calc.alphas, '_copper_', flow_calc.solvermode, '_flows.npy'])

    if flow_calc.alphas=='aHE':
        nodes = world_Nodes(admat=admat)
    elif flow_calc.alphas=='aHO0':
        nodes = world_Nodes(admat=admat, alphas=np.zeros(8))
    elif flow_calc.alphas=='aHO1':
        nodes = world_Nodes(admat=admat, alphas=np.ones(8))
    else:
        sys.stderr.write('The object has an distribution of mixes that\
                          is not accounted for.')

    if flow_calc.solvermode=='lin':
        mode = 'linear'
    elif flow_calc.solvermode=='sqr':
        mode = 'square'
    elif flow_calc.solvermode=='capM':
        mode = 'capped martin'
    elif flow_calc.solvermode=='capR':
        mode = 'capped rolando'
    else:
        sys.stderr.write('The solver mode must be "lin", "sqr", "capM" or "capR"')

    flowfilename = ''.join(['./results/', str(flow_calc), '_flows.npy'])
    print flowfilename
    print flow_calc.capacities
    if flow_calc.capacities=='copper':
        solved_nodes, flows = au.solve(nodes, mode=''.join([mode, ' copper']))
        print "Checkpt. copper"
    elif flow_calc.capacities=='q99':
        h0 = au.get_quant_caps(filename=copperflow_filename)
        print h0
        print h0.mean()
        solved_nodes, flows = au.solve(nodes, h0=h0, mode=mode)
        print "Checkpt. q99"
    elif flow_calc.capacities=='hq99': # corresponds to half the capacities
                                         # of the 99% quantile layout
        h0 = 0.5*au.get_quant_caps(filename=copperflow_filename)
        print h0
        print h0.mean()
        solved_nodes, flows = au.solve(nodes, h0=h0, mode=mode)
        print "Checkpt. hq99"
    elif flow_calc.capacities.endswith('q99'):
        scale = float(flow_calc.capacities[0:-3])
        h0 = scale*au.get_quant_caps(filename=copperflow_filename)
        print h0
        solved_nodes, flows = au.solve(nodes, h0=h0, mode=mode)
        print "Checkpt. a*q99"
    else:
        sys.stderr.write('The capacities must be either "copper", "q99",\
                            "hq99", or on the form "<number>q99"')

    solved_nodes.save_nodes(filename)
    print filename
    try:
        flows
    except NameError:
        print "Flows not defined."
    np.save('./results/' + filename + '_flows', flows)

