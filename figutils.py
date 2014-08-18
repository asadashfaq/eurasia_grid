import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

import aurespf.solvers as au
from nhgrid import nh_Nodes
from europe_plusgrid import europe_plus_Nodes
from FCResult import FCResult
from FlowCalculation import FlowCalculation
import costtools as ct

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
                    datapath='./results/BalvsTrans/', interactive=False,
                    title=True, legend=True, ylim=None):

    """ Example
        -------
        make_bal_vs_trans_graph([flowcalclin, flowcalcsqr], 'BE')

        """

    plt.close()
    plt.rc('lines', lw=2)
    plt.rcParams['axes.color_cycle'] = color_cycle

    if type(flowcalcs)!=list:
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
        plt.plot(TC,ydata, label=fc.pretty_solvermode())

    plt.xlabel('Total transmission capacity [TW]')
    if ydatalabel=='BE':
        plt.ylabel('Balancing energy [normalized]')
    if ydatalabel=='BC':
        plt.ylabel('Balancing capacity [normalized]')

    if legend:
        plt.legend()
    if title:
        plt.title(region + ': ' + flowcalcs[0].layout)

    plt.xlim((0,np.min(xmaxlist)))
    if ylim:
        plt.ylim(ylim)
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


def make_all_y_vs_alpha_graph():
    modes = ['lin', 'lin_imp', 'sqr', 'sqr_imp']

    layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']
    ydatalabels = ['BE', 'BC', 'Total_TC']
    for l in layouts:
        fc_list = []
        for m in modes:
            fc = FlowCalculation(l, 'aHO1', 'copper', m)
            fc_list.append(fc)
        for ydatalabel in ydatalabels:
            savepath = ''.join(['./results/figures/', ydatalabel, 'vsAlpha/'])

            make_y_vs_alpha_graph(fc_list, ydatalabel, savepath=savepath, \
                                   zerotrans=True, showminima=True)

def make_y_vs_alpha_graph(flowcalcs, ydatalabel, alphas=np.linspace(0,1,21), \
                          figfilename=None, savepath='./results/figures/', \
                          datapath='./results/AlphaSweepsCopper/', \
                          interactive=False, zerotrans=False,\
                          showminima=False, hetpoints=True, labels=None,\
                          title=True, small_legend=True):

    """ ydatalabel should be 'BE', 'BC' or 'Total_TC'
        The alphas's field of the FlowCalculation objects are disregared,
        an a graph is generated based on the other fields in the object,
        in a range of alphas.
        Example
        -------
        make_y_vs_alpha_graph([flowcalclin, flowcalcsqr], 'BE')

        """

    plt.close()
    plt.rc('lines', lw=2)
    plt.rcParams['axes.color_cycle'] = color_cycle

    if type(flowcalcs)!=list:
        flowcalcs = [flowcalcs]

    if interactive:
        plt.ion()

    layoutlist = [fc.layout for fc in flowcalcs]
    samelayout = (layoutlist[1:]==layoutlist[:-1]) ## True if all the layouts
                                                    # are the same, False
                                                    # otherwise
    minimum_text = 'Minima:\n'
    label_counter = 0

    for fc in flowcalcs:
        admat = "./settings/" + fc.layout + "admat.txt"
        N = nh_Nodes(admat=admat)

####### Plot the zero transmission case ######################################
        if zerotrans and ydatalabel!='Total_TC':
            if (not samelayout) or flowcalcs.index(fc)==0:
                zerotrans_ydata = []
                for a in alphas:
                    zerotrans_ydata.append(get_zerotrans_data(alpha=a, \
                                         ydatalabel=ydatalabel, fc=fc))

                if not labels:
                    zerotrans_label = ''.join([fc.layout, ': zerotrans'])
                else:
                    zerotrans_label = labels[label_counter]
                    label_counter += 1

                plt.plot(alphas, zerotrans_ydata, \
                        label=zerotrans_label)

####### extract an plot data from homogenous alpha layout ####################
        filenames = []
        for a in alphas:
            alphacode = ''.join(['aHO', str(a)])
            filenames.append(''.join([str(FlowCalculation(fc.layout,\
                         alphacode, fc.capacities, fc.solvermode)), '.pkl']))

        total_mean_load = np.sum([n.mean for n in N])
        ydata = []
        for filename in filenames:
            if ydatalabel in ['BE', 'BC']:
                unnormalized_data = get_data(filename, ydatalabel, \
                                 path=datapath)
                assert(len(unnormalized_data)==len(N))
                ydata.append(np.sum(unnormalized_data)/total_mean_load)
            elif ydatalabel == 'Total_TC':
                ydata.append(get_data(filename, ydatalabel, path=datapath)\
                             /1e6) # now in TW
        if not labels:
            label = ''.join([fc.layout, ': ', fc.solvermode])
        else:
            label = labels[label_counter]
            label_counter += 1

        plt.plot(alphas, ydata, \
                label=label)

####### find and show minima if this option is set True ######################
        if showminima:
            if zerotrans and ydatalabel!='Total_TC':
                if (not samelayout) or flowcalcs.index(fc)==0:
                    alphamin, ymin = find_interp_minimum(alphas,\
                                                 zerotrans_ydata)
                    plt.plot(alphamin, ymin, 'ok')
                    minimum_text = minimum_text + '%s %s: (%f, %f)\n' \
                            %(fc.layout, 'zerotrans', alphamin, ymin)
            alphamin, ymin = find_interp_minimum(alphas, ydata)
            plt.plot(alphamin, ymin, 'ok')
            minimum_text = minimum_text + '%s %s: (%f, %f)\n' \
                    %(fc.layout, fc.solvermode, alphamin, ymin)


####### generate point from heterogeneous alpha layout (optimal mixes) #####
####### if the hetpoints option is set as true #############################
        if hetpoints:
            het_filename = ''.join([str(FlowCalculation(fc.layout, 'aHE',\
                                    'copper', fc.solvermode)), '.pkl'])
            # the folowing i an average of the mixes in the heterogeneous
            # (optimal wrt. balancing energy) layout of alphas. This works
            # because N is loaded with these mixes as default
            avg_het_alpha = np.sum([n.mean*n.alpha for n in N])\
                                /total_mean_load
            if ydatalabel in ['BE', 'BC']:
                het_y_value = np.sum(get_data(het_filename, ydatalabel, \
                                    path=datapath))/total_mean_load
            elif ydatalabel == 'Total_TC':
                het_y_value = get_data(het_filename, ydatalabel,\
                        path=datapath)/1e6 # now in TW

            if not labels:
                het_label = fc.layout + ' ' + r'$\alpha^W_{\mathrm{opt}}$' + ': ' \
                            + fc.solvermode
            else:
                het_label = labels[label_counter]
                label_counter += 1
            plt.plot(avg_het_alpha, het_y_value, 'x', markersize=8, \
                    label=het_label)
#### finish up the plot ####################################################
    if showminima:
        plt.text(0.1, 0, minimum_text)
    plt.xlabel(r'$\alpha^W$')

    if ydatalabel=='BE':
        plt.ylabel('Balancing energy [normalized]')
    elif ydatalabel=='BC':
        plt.ylabel('Balancing capacity [normalized]')
    elif ydatalabel=='Total_TC':
        plt.ylabel('Total transmission capacity [TW]')

    if small_legend:
        plt.legend(prop={'size':7})
    else:
        plt.legend()

    if samelayout:
        if title:
            plt.title(flowcalcs[0].layout)

    plt.xlim((0, 1))
    plt.ylim(ymin=0)
    if not figfilename:
        figfilename = ''.join([flowcalcs[0].layout, '_', ydatalabel, \
                              '_vs_alpha.pdf'])

    if not interactive:
        plt.savefig(savepath+figfilename)
        plt.close()

def get_zerotrans_data(alpha, ydatalabel, fc):
    """ If ydatalabel is 'BE' it returns the normalized
        total backup energy in the layout specified in the
        Flowcalculation object fc, in a homogenous mixing
        layout with mixing alpha. If 'BC' the backup capacity
        is retured.

        """

    admat = "./settings/" + fc.layout + "admat.txt"
    N = nh_Nodes(admat=admat, alphas=alpha)
    total_mean_load = np.sum([n.mean for n in N])
    length_of_timeseries = len(N[0].mismatch)
    for n in N:
        n.balancing = -au.get_negative(n.mismatch)

    if ydatalabel=='BE':
        result = np.sum([np.sum(n.balancing) for n in N])\
                 /(length_of_timeseries*total_mean_load)
    elif ydatalabel=='BC':
        result = np.sum([au.get_q(n.balancing, 0.99) for n in N])\
                /total_mean_load
    else:
        print "ydatalabel must be 'BE' or 'BC'"
        return

    return result


def find_interp_minimum(x, y, rel_tol = 1e-3):
    """ Given an independent variable x and a variable
        dependent on x, y, this function creates a cubic
        interpolation, to estimate the minimum value of
        y as a function of x, with the relative precicion
        rel_tol.

        """

    N = int(1.0/rel_tol)
    xfine = np.linspace(np.min(x), np.max(x), N)
    f = interp1d(x, y, kind='cubic')
    ymin = f(xfine).min()
    xmin = xfine[f(xfine).argmin()]

    return xmin, ymin


def make_all_hourly_flowhists():

    modes = ['lin', 'lin_imp', 'sqr', 'sqr_imp']

    layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
               'US_eurasia_closed', 'US_EU_RU_NA_ME']

    # make plots with one mode in each plot,
    # at 4 different times a day
    for layout in layouts:
        for m in modes:
            fc = FlowCalculation(layout, 'aHO1', 'copper', m)
            savepath = './results/figures/HourlyFlowhists/Every6hours/' \
                         + layout + '/'
            make_hourly_flowhists(fc, hours = [0, 6, 12, 18],
                                  savepath = savepath)

    modes = [['lin', 'lin_imp'], ['sqr', 'sqr_imp']]
    for layout in layouts:
        for m in modes:
            flowcalcs = [FlowCalculation(layout, 'aHO1', 'copper', m[0]),\
                         FlowCalculation(layout, 'aHO1', 'copper', m[1])]
            savepath = './results/figures/HourlyFlowhists/ImpVsNoImp/'\
                         + layout + '/'
            make_hourly_flowhists(flowcalcs, hours = [0, 12],
                                  figfileending=m[0] + '_hflowhist',
                                  savepath = savepath)

    return


def make_hourly_flowhists(flowcalcs, links='all', alphas=[0.0, 0.5, 0.9, 1.0],\
                          hours = [0, 12], figfileending=None, \
                          savepath='./results/figures/', \
                          datapath='./results/AlphaSweepsCopper/', \
                          interactive=False, subplotshape=(2,2)):

    """ The produced plots will be named <link>_<figfileending>.pdf
        The arguments links must be either 'all' or on the form:
        ['EU to RU', ...] (a list!). Make sure that the alphas match
        existing .pkl-files, from the FCResults-class, number of decimals
        matter.

        """

    plt.close()
    plt.rc('lines', lw=2)
    plt.rcParams['axes.color_cycle'] = color_cycle

    if type(flowcalcs)!=list:
        flowcalcs = [flowcalcs]

    if interactive:
        plt.ion()


    layoutlist = [fc.layout for fc in flowcalcs]
    samelayout = (layoutlist[1:]==layoutlist[:-1]) ## True if all the layouts
                                                    # are the same, False
                                                    # otherwise
    only_openclosed = all([('US_eurasia' in l) for l in layoutlist])
    if links=='all' and (samelayout or only_openclosed) :
        admat = "./settings/" + layoutlist[0] + "admat.txt"
        N = nh_Nodes(admat=admat)
        linklist = [link[0] for link in au.AtoKh(N)[-1]]
    elif links=='all':
        print "'all' option for links, is only possible when all flowcalcs\
                have the same layout, or all the layouts are either\
                US_eurasia_open or -closed. No figure produced. "
        return
    else:
        linklist = links

    print linklist

    for link in linklist:
        fig = plt.figure(figsize=(12,12))
        minflow = 0
        maxflow = 0
        if samelayout:
            linkindex = get_link_number_in_layout(link, flowcalcs[0].layout)
        lines = []
        labels = []
        for fc in flowcalcs:
            if not samelayout:
                linkindex = get_link_number_in_layout(link, fc.layout)
            for a in alphas:
                filename = ''.join([fc.layout, '_aHO', str(a),\
                                     '_copper_', fc.solvermode, '.pkl'])
                flowhistdata = get_data(filename, 'hourly_flowhists', \
                                        path=datapath)

                plt.subplot(subplotshape[0], subplotshape[1], alphas.index(a)+1)

                for h in hours:

                    flow = flowhistdata[linkindex][h][0]
                    count = flowhistdata[linkindex][h][1]
                    totalcount = np.sum(count)
                    Pflow = [float(c)/totalcount for c in count]
                    # plotting flow in GW
                    if samelayout:
                        label = ''.join([fc.pretty_solvermode(),': ', r'$h=$', str(h)])
                    else:
                        label = ''.join([fc.pretty_layout(), ' ', fc.pretty_solvermode(),': ',\
                                r'$h=$', str(h)])
                    line = plt.plot(flow/1e3, Pflow, label=label)

                    if alphas.index(a)==0:
                        lines.append(line)
                        labels.append(label)

                    if maxflow < np.max(flow)/1e3:
                        maxflow = np.max(flow)/1e3
                    if minflow > np.min(flow)/1e3:
                        minflow = np.min(flow)/1e3

                plt.legend(prop={'size':10})
                plt.title(r'$\alpha^W$ = ' + str(a))
                plt.xlabel(r'$F_l$' + ' [GW]')
                plt.ylabel(r'$p(F_l | h)$')

        # loop for adjusting all the axis, to be equal
        for a in alphas:
            plt.subplot(subplotshape[0], subplotshape[1], alphas.index(a)+1)
            plt.xlim(minflow, maxflow)
            if 'lin' in [fc.solvermode for fc in flowcalcs]:
                plt.ylim(0, 0.04)
            else:
                plt.ylim(0, 0.024)
        lines = [l[0] for l in lines]
        #fig.legend(lines, labels, 'best')
        if samelayout:
            fig.suptitle(''.join([link, ' (in ', fc.pretty_layout(), ' layout)']),\
                         fontsize='20')
        else:
            fig.suptitle(link, fontsize='20')

        if not figfileending:
            figfilename = ''.join([link.replace(' ', '_'), '_', \
                                   flowcalcs[0].solvermode, '_',\
                                   'hflowhist.pdf'])
        else:
            figfilename = ''.join([link.replace(' ', '_'), '_', \
                                    figfileending, '.pdf'])

        if not interactive:
            plt.savefig(savepath+figfilename, bbox_inches='tight')
            plt.close()

    return


def get_link_number_in_layout(link, layout):
    """ Takes a link and a layout and returns the index if the
        link in the layout.

        Example:
        -------
        get_link_number_in_layout('IN to SE', 'eurasia')

        """

    admat = './settings/' + layout + 'admat.txt'
    N = nh_Nodes(admat=admat)
    linkinfo = au.AtoKh(N)[-1]
    reversedlink = ''.join([link[-2:], link[2:-2], link[0:2]])
    for i in xrange(len(linkinfo)):
        if linkinfo[i][0] in [link, reversedlink]:
            return linkinfo[i][1]

    print "Link not found"
    return


def make_bal_vs_layout_barplot(flowcalcs, ydatalabel,\
                    figfilename=None, savepath = './results/figures/', \
                    datapath='./results/BalvsTrans/', interactive=False,
                    barwidth=0.5, title=None):

    plt.close()
    plt.rcParams['axes.color_cycle'] = color_cycle

    if type(flowcalcs)!=list:
        flowcalcs = [flowcalcs]

    if interactive:
        plt.ion()

    ydata = []
    layoutlist = []
    for fc in flowcalcs:
        admat = "./settings/" + fc.layout + "admat.txt"
        N = nh_Nodes(admat=admat)
        total_mean_load = np.sum([n.mean for n in N])

        filename = str(fc) + '.pkl'
        unnormalized_data = [get_data(filename, ydatalabel, \
                             path=datapath)[n.id] for n in N]
        ydata.append(np.sum(unnormalized_data)/total_mean_load)
        layoutlist.append(fc.pretty_layout())

    index = np.arange(len(ydata))
    left = index+0.5*barwidth
    plt.ion()
    ax = plt.subplot(1,1,1)
    plt.bar(left, ydata, width=barwidth, color=blue)
    plt.xticks(left + 0.5*barwidth, layoutlist)
    if ydatalabel=='BE':
        plt.ylabel('Balancing energy [normalized]')
    elif ydatalabel=='BC':
        plt.ylabel('Balancing capacity [normalized]')

    if title:
        plt.title(title)

    if not figfilename:
        figfilename = ydatalabel + '_vs_layout.pdf'

    if not interactive:
        plt.savefig(savepath+figfilename)

    return


def make_all_LCOE_vs_alpha_graphs():
    modes = ['lin', 'lin_imp', 'sqr', 'sqr_imp']

    layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']

    for m in modes:
        for l in layouts:
            fc = FlowCalculation(l, 'aHO1.0',  'copper', m)
            make_LCOE_vs_alpha_graph(fc, savepath='./results/figures/LCOEvsAlpha/')

    return

def make_LCOE_vs_alpha_graph(masterflowcalc, alphas=np.linspace(0,1,21), \
                          figfilename=None, savepath='./results/figures/', \
                          datapath='./results/AlphaSweepsCopper/', \
                          interactive=False, CFw=0.35, CFs=0.15, title=True):
    plt.close()
    if interactive:
        plt.ion()

    total_energy = ct.total_annual_energy_consumption(masterflowcalc)
    admat = './settings/' + masterflowcalc.layout + 'admat.txt'
    BE_LCOE = []
    BC_LCOE = []
    wind_LCOE = []
    solar_LCOE = []
    TC_LCOE = []
    zerotrans_total_LCOE = []
    zerotrans_datapath = './results/AlphaSweepsZerotrans/'

    for a in alphas:
        alphacode = ''.join(['aHO', str(a)])
        fc = FlowCalculation(masterflowcalc.layout, alphacode, \
                masterflowcalc.capacities, masterflowcalc.solvermode)
        BE_LCOE.append(au.cbe(ct.total_annual_BE(fc, datapath))/total_energy)
        BC_LCOE.append(au.cbc(ct.get_total_BC(fc, datapath))/total_energy)
        wind_LCOE.append(au.cwc(ct.get_total_wind_capacity(fc, CFw, datapath))\
                         /total_energy)
        solar_LCOE.append(au.csc(\
                                ct.get_total_solar_capacity(fc, CFs, datapath))\
                         /total_energy)
        TC_LCOE.append(au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat)\
                       /total_energy)
        zerotrans_fc = FlowCalculation(masterflowcalc.layout, alphacode,\
                                         'zerotrans', 'raw')
        zerotrans_total_LCOE.append(
                (au.cbe(ct.total_annual_BE(zerotrans_fc, zerotrans_datapath))
         + au.cbc(ct.get_total_BC(zerotrans_fc, zerotrans_datapath))
         + au.cwc(\
           ct.get_total_wind_capacity(zerotrans_fc, CFw, zerotrans_datapath))
         + au.csc(\
           ct.get_total_solar_capacity(zerotrans_fc, CFs, zerotrans_datapath)))\
           /total_energy)


    plt.ion()
    plt.fill_between(alphas,
                     np.array(BE_LCOE) + np.array(BC_LCOE) + \
                     np.array(solar_LCOE) + np.array(wind_LCOE) +\
                     np.array(TC_LCOE), label='Backup energy', color=orange,
                     edgecolor='k')
    plt.fill_between(alphas,
                     np.array(BC_LCOE) +
                     np.array(solar_LCOE) + np.array(wind_LCOE) +
                     np.array(TC_LCOE), label='Backup capacity', color=red,
                     edgecolor='k')
    plt.fill_between(alphas,
                     np.array(solar_LCOE) + np.array(wind_LCOE) +
                     np.array(TC_LCOE), label='Solar capacity', color=yellow,
                     edgecolor='k')
    plt.fill_between(alphas,
                     np.array(wind_LCOE) +
                     np.array(TC_LCOE), label='Wind capacity', color=blue,
                     edgecolor='k')
    plt.fill_between(alphas, np.array(TC_LCOE), label='Transmission capacity',\
                      color=green, edgecolor='k')

    plt.plot(alphas, zerotrans_total_LCOE, color=lightblue, lw=2, ls='--',\
             label="Total LCOE, zero transmission")

    colors = [orange, red, yellow, blue, green]
    rectangles = [plt.Rectangle((0,0), 1, 1, fc=c) for c in colors]
    plt.legend(rectangles, ['Backup energy', 'Backup capacity', \
                            'Solar capacity', 'Wind capacity',\
                            'Transmission capacity'])

    plt.ylim(0,250)
    plt.xlabel(r'$\alpha_W$')
    plt.ylabel('LCOE [' + u'\u20AC' + '/MWh]')
    if title:
        plt.title(masterflowcalc.layout + ' ' + masterflowcalc.solvermode)


    if not figfilename:
        figfilename = ''.join([masterflowcalc.layout, '_', \
                        masterflowcalc.solvermode, '_LCOEvsalpha', '.pdf'])

    if not interactive:
        plt.savefig(savepath+figfilename)
        plt.close()

def make_europeplus_barplot(ydatalabel, interactive=False,\
        barwidth=0.5):
    """ Creates a bar plot showing the property specified by
        ydatalabel in the five different layouts in the Europeplus
        configuration.

        """

    plt.close()
    plt.rcParams['axes.color_cycle'] = color_cycle
    if interactive:
        plt.ion()
    datapath = "./results/Europeplus/"
    savepath = "./results/figures/Europeplusfigs/"
    sqr_and_lin = False
    if ydatalabel != 'BE':
        sqr_and_lin = True

    layoutlist = ['Europe',\
                 'Europe-RU', 'Europe-NA', 'Europe-ME', 'Europe-RU-NA-ME']

    ydata_lin = []
    ydata_sqr = []
    ydata_lin_ext = []
    ydata_sqr_ext = []

    for l in layoutlist:
        layout = l.replace('-', '_')
        fc_string_lin = str(FlowCalculation(layout, 'aHE', 'copper', 'lin'))
        fc_string_sqr = str(FlowCalculation(layout, 'aHE', 'copper', 'sqr'))
        admat = "./settings/" + layout + "admat.txt"
        N = europe_plus_Nodes(admat=admat)
        if layout=='Europe':
            internal_links = [au.AtoKh(N)[-1][i][0] \
                    for i in range(len(au.AtoKh(N)[-1]))]

        total_mean_load = np.sum([n.mean for n in N])
        if ydatalabel == 'BE':
            unnormalized_data = get_data(fc_string_lin + '.pkl', \
                                         ydatalabel, path=datapath)
            ydata_lin.append(np.sum(unnormalized_data)/total_mean_load)
            ylabel = "Backup energy [normalized]"
        elif ydatalabel == 'BC':
            unnormalized_data_lin = get_data(fc_string_lin + '.pkl',\
                                             ydatalabel, path=datapath)
            unnormalized_data_sqr = get_data(fc_string_sqr + '.pkl',\
                                             ydatalabel, path=datapath)
            ydata_lin.append(np.sum(unnormalized_data_lin)/total_mean_load)
            ydata_sqr.append(np.sum(unnormalized_data_sqr)/total_mean_load)
            ylabel = "Backup capacity [normalized]"
        elif ydatalabel == 'TC':
            internal_indices, external_indices =\
                    get_internal_external_link_indices(N, internal_links)
            print internal_indices, external_indices
            all_lin_TC = au.biggestpair(\
                    get_data(fc_string_lin + '.pkl', 'TC', path=datapath))
            all_sqr_TC = au.biggestpair(\
                    get_data(fc_string_sqr + '.pkl', 'TC', path=datapath))
            ydata_lin.append(\
                    sum([all_lin_TC[i] for i in internal_indices])/1e3)
            ydata_lin_ext.append(\
                    sum([all_lin_TC[i] for i in external_indices])/1e3)
            ydata_sqr.append(\
                    sum([all_sqr_TC[i] for i in internal_indices])/1e3)
            ydata_sqr_ext.append(\
                    sum([all_sqr_TC[i] for i in external_indices])/1e3)
            ylabel = "Transmission capacity [GW]" ## note the unit, this was
                                                  ## obtained by /1e3

################# implement transmission capacity that only include internal european links######

    index = np.array([0, 1.5, 2.5, 3.5, 5])
    left = index + 0.5*barwidth
    left2 = index + barwidth
    ax = plt.subplot(1,1,1)

    if ydatalabel=='BE':
        plt.bar(left, ydata_lin, width=barwidth, color=blue)
    elif ydatalabel=='BC':
        plt.bar(left, ydata_lin, width=0.5*barwidth,\
                color=blue, label = 'Localized flow')
        plt.bar(left2, ydata_sqr, width=0.5*barwidth,\
                color=red, label = 'Synchronized flow')
        plt.legend()
    elif ydatalabel=='TC':
        plt.bar(left, np.array(ydata_lin)+np.array(ydata_lin_ext), \
                width=0.5*barwidth, color=darkblue,\
                label = 'External capacity: Localized flow')
        plt.bar(left, np.array(ydata_lin), width=0.5*barwidth,\
                color=blue, label = 'Internal capacity: Localized flow')
        plt.bar(left2, np.array(ydata_sqr)+np.array(ydata_sqr_ext), \
                width=0.5*barwidth, color=darkred,\
                label = 'External capacity: Synchronized flow')
        plt.bar(left2, np.array(ydata_sqr), width=0.5*barwidth,\
                color=red, label = 'Internal capacity: Synchronized flow')
        plt.legend(prop={'size':13}, loc=2)
    plt.xticks(left + 0.5*barwidth, layoutlist)
    plt.ylabel(ylabel)
    plt.title('Copper flow')

    figfilename = ydatalabel + "vslayout.pdf"
    if not interactive:
        plt.savefig(savepath + figfilename)

def get_internal_external_link_indices(N, internal_links):
    all_links = [au.AtoKh(N)[-1][i][0] for i in range(len(au.AtoKh(N)[-1]))]

    internal_indices = []
    external_indices = []
    for link in all_links:
        if link in internal_links:
            internal_indices.append(all_links.index(link))
        else:
            external_indices.append(all_links.index(link))

    return internal_indices, external_indices



