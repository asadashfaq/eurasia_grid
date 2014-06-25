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
               pink, lightblue, darkred, yellow,
               darkblue, grayblue, brown]


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

    admat = "./settings/" + layout + "admat.txt"
    N = nh_Nodes(admat=admat)
    for n in N:
        TC = []
        BE = []
        for filename in filenames:
            TC.append(get_data(filename, 'Total_TC',\
                      './results/BalvsTrans/')/1e6) # now in TW
            BE.append(get_data(filename, 'BE',\
                     './results/BalvsTrans/')[n.id])

        plt.plot(TC, BE, '-o', label=str(n.label) )
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

def make_all_BEBC_TC():
    """ Running this generates the graphs
        for BE vs TC and BC vs TC for all
        layouts both for EU and for the whole layout together.

        """
    modes = ['lin', 'lin_imp', 'sqr', 'sqr_imp']

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
                    figfilename=None, savepath = './results/figures/', \
                    datapath='./results/BalvsTrans/', interactive=False):

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

        admat = "./settings/" + fc.layout + "admat.txt"
        N = nh_Nodes(admat=admat)

        filenames = []
        for a in trans_scalerange:
            capacity = ''.join([str(a), 'q99'])
            filenames.append(''.join([str(FlowCalculation(fc.layout, fc.alphas,\
                        capacity, fc.solvermode)), '.pkl']))

        ydata = []
        TC = []
        if region=='EU':
            for filename in filenames:
                ydata.append(get_data(filename, ydatalabel,\
                        path=datapath)[0]/N[0].mean)
                TC.append(get_data(filename, 'Total_TC', path=datapath)\
                          /1e6) # now in TW

        if region=='all':
            ### the total mean load for the regions in the current layout
            total_mean_load = np.sum([n.mean for n in N])

            for filename in filenames:
                unnormalized_data = [get_data(filename, ydatalabel, \
                                     path=datapath)[n.id] for n in N]
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

def get_total_linkflows(filename, path='./results/BalvsTrans/'):
    """ Function that return a list of the total flows over the
        full time series (not directed) for each link in a given
        layout. Takes a filename and path to an FCResults file.
        The total flows are (well) estimated based on the flow
        histogramst that are save in the FCResults class.

        """

    flow_hists = get_data(filename, 'flowhists', path)

    total_flows = []
    for i in range(len(flow_hists)):
        total_flows.append(np.sum(\
                [np.abs(flow_hists[i][0][j])*flow_hists[i][1][j] for j in \
                range(len(flow_hists[i][0]))]))

    return total_flows

def normalize_list(mylist):
    total = np.sum(mylist)
    normalized_list = [float(x)/total for x in mylist]
    return normalized_list

def make_relflow_vs_TC_graph(flow_calcs, \
                             trans_scalerange=np.linspace(0.05,1.5, 30),\
                    figfilename=None, savepath = './results/figures/', \
                    datapath='./results/BalvsTrans/', interactive=False):
    plt.close()
    plt.rc('lines', lw=2)

    plt.rcParams['axes.color_cycle'] = color_cycle

    if interactive:
        plt.ion()

    for fc in flow_calcs:
        plt.subplot(1,len(flow_calcs), flow_calcs.index(fc))
        admat = './settings/' + fc.layout + 'admat.txt'
        N= nh_Nodes(admat=admat)

        filenames = []
        for a in trans_scalerange:
            capacity = ''.join([str(a), 'q99'])
            filenames.append(''.join([str(FlowCalculation(fc.layout, fc.alphas,\
                        capacity, fc.solvermode)), '.pkl']))


        rel_flows = []
        TC = []
        for filename in filenames:
            rel_flows.append(normalize_list(get_total_linkflows(filename,\
                                            path=datapath)))
            TC.append(get_data(filename, 'Total_TC', path=datapath)\
                    /1e6) # Now in TW

        linklist = au.AtoKh(N)[-1]
        for i in range(len(linklist)):
            linkflow = [rel_flows[j][i] for j in range(len(rel_flows))]
            plt.plot(TC, linkflow, label=''.join([linklist[i][0], " ", \
                    fc.solvermode]))

            plt.legend(prop={'size':8})
        plt.xlabel('Total transmission capacity [TW]')
        plt.ylabel('Relative link usage')
        plt.ylim(0, 0.4)
        plt.xlim(min(TC), max(TC))
        plt.title(fc.layout)

    if not figfilename:
        figfilename = ''.join([flow_calcs[0].layout, '_relflows', '.pdf'])

    if not interactive:
        plt.savefig(savepath+figfilename)
        plt.close()

def make_all_relflowgraphs():
    modes = [['lin_imp', 'lin'], ['sqr_imp', 'sqr']]

    layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']

    savepath = './results/figures/RelflowsImpNoImp/'
    for l in layouts:
        for modesublist in modes:
            fclist = [FlowCalculation(l, 'aHE', 'copper', modesublist[0]),
                        FlowCalculation(l, 'aHE', 'copper', modesublist[1])]
            figfilename = l + '_' + modesublist[0] + '.pdf'
            make_relflow_vs_TC_graph(fclist, figfilename=figfilename,\
                                         savepath=savepath)



