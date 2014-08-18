import sys
import numpy as np
import multiprocessing as mp

import aurespf.solvers as au
from worldgrid import world_Nodes
from nhgrid import nh_Nodes
from europe_plusgrid import europe_plus_Nodes
from FlowCalculation import FlowCalculation # my own class for passing info about calculations
from FCResult import FCResult # my class for storing results in a space-efficient way

def solve_flow(flow_calc):
    admat = ''.join(['./settings/', flow_calc.layout, 'admat.txt'])
    filename = str(flow_calc)
    copperflow_filename = ''.join(['./results/', flow_calc.layout, '_',
        flow_calc.alphas, '_copper_', flow_calc.solvermode, '_flows.npy'])

    if flow_calc.alphas=='aHE':
        if flow_calc.basisnetwork == 'w':
            nodes = world_Nodes(admat=admat)
        elif flow_calc.basisnetwork == 'nh':
            nodes = nh_Nodes(admat=admat)
        elif flow_calc.basisnetwork == 'europeplus':
            nodes = europe_plus_Nodes(admat=admat)
        else:
            sys.stderr.write('The object has a basisnetwork that\
                          is not accounted for. Use "w" or "nh" or "europeplus".')
    elif flow_calc.alphas.startswith('aHO'):
        alpha = float(flow_calc.alphas[3:]) # expects alphas on the form aHO0.4
        if flow_calc.basisnetwork == 'nh':
            nodes = nh_Nodes(admat=admat, alphas=alpha)
        elif flow_calc.basisnetwork == 'europeplus':
            nodes = europe_plus_Nodes(admat=admat, alphas=alpha)
        else:
            sys.stderr.write('The object has a basisnetwork that\
                          is not accounted for. Use "nh".')
    else:
        sys.stderr.write('The object has an distribution of mixes that\
                          is not accounted for.')

    mode_str_list = []
    if 'lin' in flow_calc.solvermode:
        mode_str_list.append('linear ')
    elif 'sqr' in flow_calc.solvermode:
        mode_str_list.append('square ')
    elif 'cap' in flow_calc.solvermode:
        if 'M' in flow_calc.solvermode:
            mode_str_list.append('capped martin ')
        elif 'R' in flow_calc.solvermode:
            mode_str_list.append('capped rolando ')
    else:
        sys.stderr.write('The solver mode must be "lin", "sqr", "capM" or "capR"')
    if 'imp' in flow_calc.solvermode:
        mode_str_list.append('impedance ')

    mode = ''.join(mode_str_list)


    flowfilename = ''.join(['./results/', str(flow_calc), '_flows.npy'])

    if flow_calc.capacities=='copper':
        solved_nodes, flows = au.solve(nodes, mode=''.join([mode, ' copper']),\
                                msg=str(flow_calc))
    elif flow_calc.capacities=='q99':
        h0 = au.get_quant_caps(filename=copperflow_filename)
        solved_nodes, flows = au.solve(nodes, h0=h0, mode=mode, \
                                         msg=str(flow_calc))
    elif flow_calc.capacities=='hq99': # corresponds to half the capacities
                                         # of the 99% quantile layout
        h0 = 0.5*au.get_quant_caps(filename=copperflow_filename)
        solved_nodes, flows = au.solve(nodes, h0=h0, mode=mode, \
                                        msg=str(flow_calc))
    elif flow_calc.capacities.endswith('q99'):
        scale = float(flow_calc.capacities[0:-3])
        h0 = scale*au.get_quant_caps(filename=copperflow_filename)
        solved_nodes, flows = au.solve(nodes, h0=h0, mode=mode,\
                                        msg=str(flow_calc))
    else:
        sys.stderr.write('The capacities must be either "copper", "q99",\
                            "hq99", or on the form "<number>q99"')

    if flow_calc.savemode == 'full':
        solved_nodes.save_nodes(filename)
        try:
            flows
            np.save('./results/' + filename + '_flows', flows)
        except NameError:
            print "Flows not defined."

    elif flow_calc.savemode == 'FCResult':
        result = FCResult(filename+'.pkl')
        result.add_instance(solved_nodes, flows, flow_calc)
        result.save_results(filename+'.pkl')
    else:
        print "Results not saved, invalid savemod provided"

