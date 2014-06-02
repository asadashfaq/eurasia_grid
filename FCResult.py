import numpy as np
import cPickle as pickle
import aurespf.solvers as au

HOURS_IN_A_YEAR = 8766

class FCResult:

    def __init__(self, filename, path="./results/"):
        self.cache = []

        try:
            self.__load_results__(filename, path)
        except:
            self.save_results(filename, path)

    def __load_results__(self, filename, path="./results/"):
        with open(path+filename) as currentfile:
            self.cache = pickle.load(currentfile)
            assert(self.cache.__class__ == list),\
                    "Attempt to load wrong data format "\
                    + filename + "may be corrupt."

    def __len__(self):
        return len(self.cache)

    def __str__(self):
        return filename

    def save_results(self, filename, path="./results/"):
        with open(path+filename, 'w') as currentfile:
            pickle.dump(self.cache, currentfile)

    def add_instance(self, N, F, FC=None):

        inst={}
        inst['alphas'] = [n.alpha for n in N]
        inst['gammas'] = [n.gamma for n in N]
        length_of_timeseries = len(N[0].balancing)
        inst['BE'] = [np.sum(n.balancing)/\
                          (length_of_timeseries*n.mean) for n in N]
        inst['BC'] = [au.get_q(n.balancing, 0.99)/n.mean for n in N]
        if FC:
            inst['FlowCalculation'] = FC
            if FC.capacities[-3:len(FC.capacities)] == "q99":
                inst['Total_TC'] = TC_from_FC(FC)

        self.cache.append(inst)

def TC_from_FC(FC):
    copperFC = FC.copy()
    copperFC.capacities = "copper"
    copperflowfile = ''.join(["results/", str(copperFC),\
                         '_flows.npy'])
    h0 = au.get_quant_caps(filename=copperflowfile)
    # expecting capacities of format: 0.35q99
    scalefactor = float(FC.capacities[0:-3])

    total_TC = scalefactor*np.sum(au.biggestpair(h0))
    return total_TC
