import numpy as np
import aurespf.solvers as au
from nhgrid import nh_Nodes
from FCResult import FCResult
from FlowCalculation import FlowCalculation

HOURS_PER_YEAR = 8760
def get_data(filename, field, path='./results/'):
    """ Returns the data in a certain field,
        from an FCResult file.
        See the FCResult class for available fields.

        Example
        -------
        >>> get_data("eurasia_aHE_0.95q99_lin.pkl", 'Total_TC')

        """

    result = FCResult(filename, path=path)
    returnvalue = result.cache[0][field]

    return returnvalue

def filename_from_flowcalc(flowcalc):
    return str(flowcalc) + '.pkl'


def total_annual_BE(flowcalc, datapath='./results/AlphaSweepsCopper/', onlyEU=False):
    filename = filename_from_flowcalc(flowcalc)
    if onlyEU:
        return get_data(filename, 'BE', path=datapath)[0]\
            *HOURS_PER_YEAR
    return np.sum(get_data(filename, 'BE', path=datapath))\
            *HOURS_PER_YEAR


def get_total_BC(flowcalc, datapath='./results/AlphaSweepsCopper/', onlyEU=False):
    filename = filename_from_flowcalc(flowcalc)
    if onlyEU:
        return get_data(filename, 'BC', path=datapath)[0]
    return np.sum(get_data(filename, 'BC', path=datapath))


def get_total_wind_capacity(flowcalc, capacityfactor,
        datapath='./results/AlphaSweepsCopper/', onlyEU=False):
    filename = filename_from_flowcalc(flowcalc)
    alphas = get_data(filename, 'alphas', path=datapath)
    gammas = get_data(filename, 'gammas', path=datapath)
    meanloads = np.load('./results/' + flowcalc.layout + '_meanloads.npy')

    if onlyEU:
        return gammas[0]*alphas[0]*meanloads[0]/capacityfactor
    return np.sum([gammas[i]*alphas[i]*meanloads[i]\
                   /capacityfactor for i in range(len(meanloads))])


def get_total_solar_capacity(flowcalc, capacityfactor,
        datapath='./results/AlphaSweepsCopper/', onlyEU=False):
    filename = filename_from_flowcalc(flowcalc)
    alphas = get_data(filename, 'alphas', path=datapath)
    gammas = get_data(filename, 'gammas', path=datapath)
    meanloads = np.load('./results/' + flowcalc.layout + '_meanloads.npy')

    if onlyEU:
        return gammas[0]*(1-alphas[0])*meanloads[0]/capacityfactor
    return np.sum([gammas[i]*(1-alphas[i])*meanloads[i]\
                   /capacityfactor for i in range(len(meanloads))])


def get_TCs(flowcalc, datapath='./results/AlphaSweepsCopper/'):
    filename = filename_from_flowcalc(flowcalc)
    h0 = get_data(filename, 'TC', path=datapath)
    return au.biggestpair(h0)

def total_annual_energy_consumption(flowcalc, r=4.0, onlyEU=False):
    admat = './settings/' + flowcalc.layout + 'admat.txt'
    N = nh_Nodes(admat=admat)
    if onlyEU:
        return N[0].mean*HOURS_PER_YEAR*au.ann(r)
    return np.sum([n.mean for n in N])*HOURS_PER_YEAR*au.ann(r)






