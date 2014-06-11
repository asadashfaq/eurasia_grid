import matplotlib.pyplot as plt
import numpy as np

import aurespf.solvers as au
from nhgrid import nh_Nodes
from FCResult import FCResult
from FlowCalculation import FlowCalculation

#### Colorschemes ################################

dth = (3.425)
dcolwidth = (2*3.425+0.236)

blue = '#134b7c'
yellow = '#f8ca00'
orange = '#e97f02'
brown = '#876310'
green = '#4a8e05'
red = '#ae1215'
purple = '#4f0a3d'
darkred= '#4f1215'
pink = '#bd157d'
lightpink = '#d89bc2'
aqua = '#37a688'
darkblue = '#09233b'
lightblue = '#8dc1e0'
grayblue = '#4a7fa2'

blue_cycle = [darkblue, blue, grayblue, lightblue]

color_cycle = [blue, red, orange, purple, green, \
               pink, lightblue, darkred, yellow]

def make_layout_flow_hists(flowcalc, interactive=False):
    """ Takes a FlowCalculation object,
        assuming there is a corresponding saved nodes
        file (.npz) and flows (.npy).

        """

    filename = str(flowcalc) + '.npz'
    N = nh_Nodes(load_filename=filename)
    EU_meanload = N[0].mean
    print(EU_meanload)
    print(N.pathadmat)
    flowfilename = str(flowcalc) + '_flows.npy'
    flows = np.load('./results/'+flowfilename)

    listflows = [au.AtoKh(N)[-1][i][0] for i in range(len(au.AtoKh(N)[-1]))]
    print([min(flows[i]) for i in range(len(flows))])
    print([max(flows[i]) for i in range(len(flows))])
    print(listflows)

    xmin = -1.9e5/EU_meanload
    xmax = 1.5e5/EU_meanload
    bins = np.linspace(xmin, xmax, 200)
    if interactive:
        plt.ion()

    for f in listflows:
        index = listflows.index(f)
        plt.figure()
        plt.hist(flows[index]/EU_meanload, bins=bins, normed=True)
        plt.title(f + ': ' + flowcalc.capacities)
        plt.xlabel("Directed power flow [normalized to EU mean load]")
        plt.ylabel("$P(F_l)$")
        figfilename = f.replace(' ','_') + '_' + flowcalc.capacities + '_flowhist.pdf'
        savepath = './results/figures/CheckStrangeBE/'
        if not interactive:
            plt.savefig(savepath + figfilename)
            plt.close()


def make_layout_mismatchhists(flowcalc, interactive=False):
    """ Takes a FlowCalculation object,
        assuming there is a corresponding saved nodes
        file (.npz).

        """

    filename = str(flowcalc) + '.npz'
    N = nh_Nodes(load_filename=filename)
    xmin = -1
    xmax = 2
    bins = np.linspace(xmin,xmax,200)
    if interactive:
        plt.ion()

    for n in N:
        if n.id in get_indices_from_layout('EU_RU_NA_ME'):
            mismatch = n.curtailment - n.balancing
            nonzero_mismatch = []
            for w in mismatch:
                if w>=1 or w<-1:
                    nonzero_mismatch.append(float(w)/n.mean)

            plt.figure()
            plt.hist(nonzero_mismatch, bins=bins, normed=True)
            plt.title(''.join([str(n.label), '(', flowcalc.capacities, ')']))
            plt.xlabel('Mismatch [normalized]')
            plt.ylabel('Mismatch distribution')

            figfilename = ''.join([str(n.label), '_mismatchhist_',\
                    flowcalc.capacities, '.pdf'])
            savepath = "./results/figures/CheckStrangeBE/"
            if not interactive:
                plt.savefig(savepath+figfilename)
                plt.close()


def ion_BE_vs_TC_testplot(layout):
    """ This function plots the normalized balancing energy
        as a function of transmission capacity for all the
        regions in a provided layout. (sqr flow)

        """

    plt.ion()
    plt.figure()
    plt.rcParams['axes.color_cycle'] = color_cycle
    filenames = []
    for a in np.linspace(0,1.5,31):
        capacity = ''.join([str(a), 'q99'])
        filenames.append(''.join([layout,'_aHE_',capacity,'_sqr.pkl']))

    indexlist = get_indices_from_layout(layout)
    N = nh_Nodes()
    for i in indexlist:
        TC = []
        BE = []
        for filename in filenames:
            TC.append(get_data(filename, 'Total_TC',\
                      './results/BalvsTransNoImpedance/')/1e6) # now in TW
            BE.append(get_data(filename, 'BE',\
                     './results/BalvsTransNoImpedance/')[i])

        plt.plot(TC, BE, '-o', label=str(N[i].label) )
    plt.xlabel('Total transmission capacity [TW]')
    plt.ylabel('Balancing energy [normalized]')
    plt.legend()
    plt.title(layout)


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

def make_all_noimp_BEBC_TC():
    """ Running this generates the graphs
        for BE vs TC and BC vs TC for all
        layouts both for EU and for the whole layout together.

        """

    modes = ['lin', 'sqr']
    layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']
    ydatalabels = ['BE', 'BC']

    regions = ['EU', 'all']

    for r in regions:
        for l in layouts:
            fc_list = []
            for m in modes:
                fc = FlowCalculation(l, 'aHE', 'copper', m)
                fc_list.append(fc)
            for ydatalabel in ydatalabels:
                savepath = ''.join(['./results/figures/', ydatalabel, 'vsTC/'])
                make_bal_vs_trans_graph(fc_list, ydatalabel, region=r,\
                                        savepath=savepath)


def make_bal_vs_trans_graph(flowcalcs, ydatalabel, region='EU', \
                    trans_scalerange=np.linspace(0,1.5, 31),\
                    figfilename=None, savepath = './results/figures/', datapath='./results/BalvsTransNoImpedance/', interactive=False):

    """ Example
        -------
        make_bal_vs_trans_graph([flowcalclin, flowcalcsqr], 'BE')

        """

    plt.close()
    plt.rc('lines', lw=2)

    if not flowcalcs.__class__==list:
        flowcalcs = [flowcalcs]

    xmaxlist = []

    if interactive:
        plt.ion()

    for fc in flowcalcs:
        filenames = []
        for a in trans_scalerange:
            capacity = ''.join([str(a), 'q99'])
            filenames.append(''.join([str(FlowCalculation(fc.layout, fc.alphas, \
                        capacity, fc.solvermode)), '.pkl']))

        ydata = []
        TC = []
        if region=='EU':
            for filename in filenames:
                ydata.append(get_data(filename, ydatalabel, path=datapath)[0])
                TC.append(get_data(filename, 'Total_TC', path=datapath)\
                          /1e6) # now in TW

        if region=='all':
            N = nh_Nodes()
            indexlist = get_indices_from_layout(fc.layout)
            mean_load_list = [N[i].mean for i in xrange(9)]
            ### the total mean load for the regions in the current layout
            total_mean_load = np.sum([mean_load_list[i] for i in indexlist])

            for filename in filenames:
                unnormalized_data = [mean_load_list[i]*\
                                     get_data(filename, ydatalabel, \
                                     path=datapath)[i] for i in indexlist]
                ydata.append(np.sum(unnormalized_data)/total_mean_load)
                TC.append(get_data(filename, 'Total_TC', path=datapath)\
                          /1e6) # now in TW


        xmaxlist.append(np.max(TC))
        plt.plot(TC,ydata, label=''.join([fc.layout,': ', fc.solvermode]))

    plt.xlabel('Total transmission capacity [TW]')
    if ydatalabel=='BE':
        plt.ylabel('Balancing energy [normalized]')
    if ydatalabel=='BC':
        plt.ylabel('Balancing capacity [normalized]')

    plt.legend()
    plt.title(region + ': ' + flowcalcs[0].layout)
    plt.xlim((0,np.min(xmaxlist)))
    if not figfilename:
        figfilename = ''.join([flowcalcs[0].layout, '_', ydatalabel, '_', \
                region, '.pdf'])

    if not interactive:
        plt.savefig(savepath+figfilename)
        plt.close()

def get_indices_from_layout(layout):
    """ When given a layout, in the Nothern
        hemisphere network, it returns a list
        of indices of only the regions in said layuout.

        """

    indexlist = []
    if 'EU_RU_NA_ME' in layout:
        indexlist.extend(range(4))

    if 'eurasia' in layout:
        indexlist.extend(range(8))

    if 'US' in layout:
        indexlist.append(8)

    return indexlist


def BC_vs_TC_old(filenames = ['BC_TC_datacapM.npz', 'BC_TC_datacapR.npz'],\
        path='./results/', labels = None, savepath='./results/figures/',\
        figfilename = 'BC_vs_TC_RolMar.pdf', BC_type='both'):
    plt.close()
    plt.rc('lines', lw=2)

    labelcount = 0
    for f in filenames:
        data = np.load(''.join([path, f]))
        BC_max = data['BC_max']
        BC_q99 = data['BC_q99']
        TC = data['TC']
        if not labels:
            label = f
        else:
            label = labels[labelcount]
        if BC_type in ['max', 'both']:
            plt.plot(1e-6*TC, BC_max, 'o', label='BCmax: '+label)
        if BC_type in ['q99', 'both']:
            plt.plot(1e-6*TC, BC_q99, 'o', label='BCq99: '+label)
        else:
            print 'Try "max", "q99" or "both" for BC_type'
            return

        labelcount = labelcount + 1

    plt.xlabel('Transmission capacity [TW]')
    plt.ylabel('Balancing capacity [normalized]')
    plt.legend()

    plt.savefig(savepath+figfilename)
    plt.close()

