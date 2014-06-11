import matplotlib.pyplot as plt
import numpy as np
import aurespf.solvers as au
from nhgrid import nh_Nodes

""" This file looks for incidents in the small network,
    EU + first neighbors where for a region the balancing
    energy is bigger for the unconstrained the the constrained case.
    Look in the ./results/figures/CheckStrangeBE/ folder for results.

    """

regions = ['EU', 'RU', 'NA', 'ME']

Ncopper = nh_Nodes(load_filename="EU_RU_NA_ME_aHE_copper_sqr.npz")
Nconstrained = nh_Nodes(load_filename="EU_RU_NA_ME_aHE_0.6q99_sqr.npz")
copper_flows = np.load("./results/EU_RU_NA_ME_aHE_copper_sqr_flows.npy")
constr_flows = np.load("./results/EU_RU_NA_ME_aHE_0.6q99_sqr_flows.npy")

linklist = [au.AtoKh(Ncopper)[-1][i][0] \
            for i in range(len(au.AtoKh(Ncopper)[-1]))]
print(linklist)

def analyze_incident(hour):
    mismatch = {}
    copper_injection = {}
    constr_injection = {}
    copper_flow_dict = {}
    constr_flow_dict = {}
    for i in xrange(4):
        mismatch[regions[i]] = Ncopper[i].mismatch[hour]
        copper_injection[regions[i]] = mismatch[regions[i]]\
                                        - Ncopper[i].curtailment[hour]\
                                        + Ncopper[i].balancing[hour]
        constr_injection[regions[i]] = mismatch[regions[i]]\
                                        - Nconstrained[i].curtailment[hour]\
                                        + Nconstrained[i].balancing[hour]

    constr_link_capacities = 0.6*au.get_quant_caps(\
                     filename="./results/EU_RU_NA_ME_aHE_copper_sqr_flows.npy")
    constr_link_cap_dict = {}
    for j in xrange(len(linklist)):
        copper_flow_dict[linklist[j]] = copper_flows[j][hour]
        constr_flow_dict[linklist[j]] = constr_flows[j][hour]
        constr_link_cap_dict[linklist[j]] = [constr_link_capacities[2*j], \
                                            constr_link_capacities[2*j+1]]


    print(np.sum(copper_injection.values()))
    print(np.sum(constr_injection.values()))

    return mismatch, copper_injection, constr_injection, copper_flow_dict,\
            constr_flow_dict, constr_link_cap_dict

def get_bal_at_hour(copper_nodes, constr_nodes, hour, region='ME'):
    mean_load = copper_nodes[regions.index(region)].mean
    copper_bal = copper_nodes[regions.index(region)].balancing[hour]/mean_load
    constr_bal = constr_nodes[regions.index(region)].balancing[hour]/mean_load

    return [copper_bal, constr_bal]


def find_incidents(copper_nodes=Ncopper, constr_nodes=Nconstrained):
    successes = 0
    N = 280512
    hours_of_interest = []
    for h in xrange(N):
        copper_bal = get_bal_at_hour(Ncopper, Nconstrained, hour=h, region='ME')[0]
        constr_bal = get_bal_at_hour(Ncopper, Nconstrained, hour=h, region='ME')[1]

        if ((copper_bal - constr_bal)/copper_bal > 1e-4) and \
                copper_bal > 1e-6:
            successes = successes + 1
            print(copper_bal, constr_bal, h)
            hours_of_interest.append(h)


    print('Successrate: ', float(successes)/N)

    plt.ion()
    plt.plot(Ncopper[3].balancing/Ncopper[3].mean, label='copper')
    plt.plot(Nconstrained[3].balancing/Nconstrained[3].mean, label='constrained')
    plt.title('ME, normalized balancing energy')
    plt.legend()
    plt.xlabel('Hour')
    plt.ylabel('Balancing [normalized]')

    return hours_of_interest
