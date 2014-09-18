import figutils as fig
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

import aurespf.solvers as au
from nhgrid import nh_Nodes
from europe_plusgrid import europe_plus_Nodes
from FCResult import FCResult
from FCResult import myhist
from FlowCalculation import FlowCalculation
import costtools as ct

savepath = './results/figures/Articlefigs/'

def make_mismatchfig(interactive=True):


    plt.close()
    plt.rcParams['axes.color_cycle'] = fig.color_cycle
    if interactive:
        plt.ion()

    layouts = ['EU_RU_NA_ME', 'US_EU_RU_NA_ME', 'eurasia', 'US_eurasia_closed']
    labels = ['EU-RU-NA-ME', 'US-EU-RU-NA-ME', 'Eurasia', 'US-Eurasia']

    bins = np.linspace(-1.2, 2.3,500)
    N = nh_Nodes(admat='./settings/EU_RU_NA_MEadmat.txt')
    EUmismatch = N[0].mismatch/N[0].mean
    EUvalue, EUcount = myhist(EUmismatch, bins=bins, normed=True)
    ax1 = plt.subplot(2,1,1)
    ax1.plot(EUvalue, EUcount, label='EU', lw=2)

    for layout in layouts:
        admat = './settings/' + layout + 'admat.txt'
        N = nh_Nodes(admat=admat)
        total_mean_load = sum([n.mean for n in N])
        mismatch = sum([n.mismatch for n in N])/total_mean_load
        value, count = myhist(mismatch, bins=bins, normed=True)
        ax1.plot(value, count, label=labels[layouts.index(layout)], lw=2)

    ax1.legend()
    ax1.set_xlim((-1.2,2.3))
    ax1.set_xlabel('Aggregated mismatch [normalized]')
    ax1.set_ylabel('Probability density')
    ax1.text(-1.12,2.1,'(a)', fontdict={'size':20})

    ax2 = plt.subplot(2,1,2)
    ax2.plot(EUvalue, EUcount, label='EU', lw=2)
    for layout in layouts:
        load_filename = layout+'_aHE_copper_lin.npz'
        N = nh_Nodes(load_filename=load_filename)
        normed_res_mismatches = []
        for n in [N[0]]:
            normed_res_mismatches.extend((n.curtailment - n.balancing)/n.mean)
        nonzero_res_mismatch = [d for d in normed_res_mismatches if (d>0.01 or d<-0.01)]
        value, count = myhist(normed_res_mismatches, bins=bins, normed=True)
        nonzero_values = [v for v in value if (v<-0.005 or v>0.01)]
        nonzero_indices = [np.where(value==v)[0][0] for v in value if (v<-0.005 or v>0.01)]
        nonzero_count = count[nonzero_indices]
        ax2.plot(nonzero_values, nonzero_count,label=labels[layouts.index(layout)], lw=2)
    legend2 = ax2.legend(title='EU embedded in:')
    legend2 = ax2.legend(title='EU embedded in:')
    plt.setp(legend2.get_title(),fontsize=14)
    ax2.set_xlim((-1.2,2.3))
    ax2.set_ylim((0,1.2))
    ax2.set_xlabel('Non-zero residual mismatch [normalized]')
    ax2.set_ylabel('Probability density')
    ax2.text(-1.12,1,'(b)', fontdict={'fontsize':20})

    plt.tight_layout()


    figfilename = 'mismatch_hists.pdf'
    if not interactive:
        plt.savefig(savepath + figfilename)



def make_vs_layout_barplot(interactive=True):

    barwidth=0.5
    datapath='./results/BalvsTrans/'
    plt.close()
    plt.rcParams['axes.color_cycle'] = fig.color_cycle

    if interactive:
        plt.ion()

    BE = []
    BC = []
    TC = []

    layouts = ['EU_RU_NA_ME', 'US_EU_RU_NA_ME', 'eurasia', 'US_eurasia_closed']
    flowcalcs = [FlowCalculation(layout, 'aHE', '1.5q99', 'lin')\
                 for layout in layouts]
    print [str(fc) for fc in flowcalcs]
    N = nh_Nodes()
    EU = N[0]
    EU_bal = -au.get_negative(EU.mismatch)
    EU_BE = np.sum(EU_bal)/(EU.mean*len(EU.mismatch))
    BE.append(EU_BE)
    EU_BC = au.get_q(EU_bal, 0.99)/EU.mean
    BC.append(EU_BC)
    TC.append(0.01)

    for fc in flowcalcs:
        print fc.layout
        admat = "./settings/" + fc.layout + "admat.txt"
        N = nh_Nodes(admat=admat)
        total_mean_load = np.sum([n.mean for n in N])
        print total_mean_load

        filename = str(fc) + '.pkl'
        flowfilename = './results/' + fc.layout + '_aHE_copper_lin_flows.npy'
        unnormalized_BE = [fig.get_data(filename, 'BE', \
                             path=datapath)[n.id] for n in N]
        unnormalized_BC = [fig.get_data(filename, 'BC', \
                             path=datapath)[n.id] for n in N]
        BE.append(np.sum(unnormalized_BE)/total_mean_load)
        BC.append(np.sum(unnormalized_BC)/total_mean_load)
        LI = au.linfo(admat)
        energywiseTC = au.biggestpair(au.get_quant_caps(filename=flowfilename))
        print np.sum(energywiseTC)
        TC.append(np.sum([energywiseTC[i]*float(LI[i][2]) \
                            for i in range(len(LI))])/total_mean_load)

        print TC
    index1 = np.arange(len(BE))
    left = index1 + 0.5*barwidth
    plt.ion()
    ax1 = plt.subplot(3,1,1)
    plt.gcf().set_size_inches([7.5,12])
    plt.gcf().set_dpi(400)
    ax1.bar(left, BE, width=barwidth, color=fig.orange)
    layoutlist1 = ['EU', 'EU-RU-NA-ME', 'US-EU-RU-NA-ME', 'Eurasia', 'US-Eurasia']
    ax1.set_xticks(left + 0.5*barwidth)
    ax1.set_xticklabels(layoutlist1)
    ax1.set_ylabel('Backup energy [normalized]')#r'$E_\mathrm{total}^B$' + ' [normalized]')
    ax1.set_ylim((0,0.2))
    ax1.text(0.1*barwidth, 0.17, '(a)',fontdict={'fontsize':20})

    ax2 = plt.subplot(3,1,2)
    ax2.bar(left, BC, width=barwidth, color=fig.red)
    ax2.set_xticks(left + 0.5*barwidth)
    ax2.set_xticklabels(layoutlist1)
    ax2.set_ylabel('Backup capacity [normalized]')#r'$C_\mathrm{total}^B$'+ ' [normalized]')
    ax2.set_ylim((0.5, 0.8))
    ax2.text(0.1*barwidth, 0.748, '(b)',fontdict={'fontsize':20})

    layoutlist2 = ['', 'EU-RU-NA-ME', 'US-EU-RU-NA-ME', 'Eurasia', 'US-Eurasia']
    ax3 = plt.subplot(3,1,3)
    ax3.bar(left, TC, width=barwidth, color=fig.green)
    ax3.set_xticks(left + 0.5*barwidth)
    ax3.set_xticklabels(layoutlist2)
    ax3.set_ylabel('Transmission capacity [km]')#r'$C_\mathrm{total}^T$' + ' [km]')
    ax3.text(0.1*barwidth, 3370, '(c)',fontdict={'fontsize':20})

    plt.tight_layout()

    figfilename = 'BEBCTC_vs_layout.pdf'

    if not interactive:
        plt.savefig(savepath+figfilename)


def make_LCOE_barplot(interactive=True):
    plt.close()
    if interactive:
        plt.ion()

    CFw = 0.35
    CFs = 0.15
    datapath = './results/AlphaSweepsCopper/'
    layouts = ['EU_RU_NA_ME', 'US_EU_RU_NA_ME', 'eurasia', 'US_eurasia_closed']
    total_LCOE = []
    LCOE2 = []
    LCOE3 = []
    LCOE4 = []
    LCOE5 = []

    # find LCOE for EU isolated
    zerotrans_fc = FlowCalculation('eurasia', 'aHE', 'zerotrans', 'raw')
    zerotrans_datapath = './results/AlphaSweepsZerotrans/'
    total_EU_energy = ct.total_annual_energy_consumption(zerotrans_fc, onlyEU=True)
    EU_BE_LCOE = au.cbe(ct.total_annual_BE(zerotrans_fc, datapath=zerotrans_datapath, onlyEU=True))/total_EU_energy

    EU_BC_LCOE = au.cbc(ct.get_total_BC(zerotrans_fc, datapath=zerotrans_datapath, onlyEU=True))/total_EU_energy
    EU_wind_LCOE = au.cwc(ct.get_total_wind_capacity(zerotrans_fc, CFw, zerotrans_datapath, onlyEU=True))\
                         /total_EU_energy
    EU_solar_LCOE = au.csc(ct.get_total_solar_capacity(zerotrans_fc, CFs, datapath, onlyEU=True))\
                         /total_EU_energy

    total_LCOE.append(sum([EU_BE_LCOE, EU_BC_LCOE, EU_wind_LCOE, EU_solar_LCOE]))
    LCOE2.append(sum([EU_BC_LCOE, EU_wind_LCOE, EU_solar_LCOE]))
    LCOE3.append(sum([EU_wind_LCOE, EU_solar_LCOE]))
    LCOE4.append(sum([EU_wind_LCOE]))
    LCOE5.append(0)


    fclist = [FlowCalculation(l, 'aHE', 'copper', 'lin') for l in layouts]
    for fc in fclist:

        admat = './settings/' + fc.layout + 'admat.txt'
        total_energy = ct.total_annual_energy_consumption(fc)
        print fc.layout, total_energy
        BE_LCOE = au.cbe(ct.total_annual_BE(fc, datapath))/total_energy
        BC_LCOE = au.cbc(ct.get_total_BC(fc, datapath))/total_energy
        wind_LCOE = au.cwc(ct.get_total_wind_capacity(fc, CFw, datapath))\
                         /total_energy
        solar_LCOE = au.csc(ct.get_total_solar_capacity(fc, CFs, datapath))\
                         /total_energy
        TC_LCOE = au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat)\
                       /total_energy
        total_LCOE.append(sum([BE_LCOE, BC_LCOE, wind_LCOE, solar_LCOE, TC_LCOE]))
        LCOE2.append(sum([BC_LCOE, wind_LCOE, solar_LCOE, TC_LCOE]))
        LCOE3.append(sum([wind_LCOE, solar_LCOE, TC_LCOE]))
        LCOE4.append(sum([wind_LCOE, TC_LCOE]))
        LCOE5.append(sum([TC_LCOE]))


        print total_LCOE

    barwidth = 0.5
    ax = plt.subplot(1,1,1)
    index = np.arange(len(total_LCOE))
    left = index + 0.5*barwidth
    ax.bar(left, total_LCOE, width=barwidth, color=fig.orange, label='Backup energy')
    ax.bar(left, LCOE2, width=barwidth, color=fig.red, label='Backup capacity')
    ax.bar(left, LCOE3, width=barwidth, color=fig.yellow, label='Solar capacity')
    ax.bar(left, LCOE4, width=barwidth, color=fig.blue, label='Wind capacity')
    ax.bar(left, LCOE5, width=barwidth, color=fig.green, label='Transmission capacity')
    ax.set_xticks(left + 0.5*barwidth)
    layoutlist = ['EU', 'EU-RU-NA-ME', 'US-EU-RU-NA-ME', 'Eurasia', 'US-Eurasia']
    ax.set_xticklabels(layoutlist)
    ax.set_ylabel('LCOE [' + u'\u20AC' + '/MWh]')
    ax.legend(loc=2, prop={'size':12})
    plt.tight_layout()

    if not interactive:
        figfilename = 'LCOEvslayout.pdf'
        plt.savefig(savepath+figfilename)


def LCOE_vs_alpha(interactive=True):
    masterflowcalc = FlowCalculation('US_eurasia_closed', 'aHO1', 'copper', 'lin')
    alphas = np.linspace(0,1,21)
    datapath = './results/AlphaSweepsCopper/'
    CFw = 0.35
    CFs = 0.15

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
    EU_BE_LCOE = []
    EU_BC_LCOE = []
    EU_wind_LCOE = []
    EU_solar_LCOE = []
    total_EU_energy = ct.total_annual_energy_consumption(\
                                masterflowcalc, onlyEU=True)


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

        EU_BE_LCOE.append(au.cbe(ct.total_annual_BE(\
                    zerotrans_fc, zerotrans_datapath, onlyEU=True))\
                    /total_EU_energy)
        EU_BC_LCOE.append(au.cbc(ct.get_total_BC(\
            zerotrans_fc, zerotrans_datapath, onlyEU=True))/total_EU_energy)
        EU_wind_LCOE.append(au.cwc(ct.get_total_wind_capacity(\
                zerotrans_fc, CFw, zerotrans_datapath, onlyEU=True))\
                /total_EU_energy)
        EU_solar_LCOE.append(au.csc(ct.get_total_solar_capacity(\
                zerotrans_fc, CFs, zerotrans_datapath, onlyEU=True))\
                /total_EU_energy)

    print EU_BE_LCOE, EU_BC_LCOE, EU_wind_LCOE, EU_solar_LCOE

    plt.ion()
    ax1 = plt.subplot(2,1,1)
    plt.gcf().set_size_inches([7.5,12])
    plt.gcf().set_dpi(400)
    ax2 = plt.subplot(2,1,2)
    ax1.fill_between(alphas,
                     np.array(EU_BE_LCOE) + np.array(EU_BC_LCOE) + \
                     np.array(EU_solar_LCOE) + np.array(EU_wind_LCOE),
                     label='Backup energy', color=fig.orange, edgecolor='k')
    ax1.fill_between(alphas,
                     np.array(EU_BC_LCOE) + \
                     np.array(EU_solar_LCOE) + np.array(EU_wind_LCOE),
                     label='Backup capacity', color=fig.red,
                     edgecolor='k')
    ax1.fill_between(alphas,
                     np.array(EU_solar_LCOE) + np.array(EU_wind_LCOE),
                     label='Solar capacity', color=fig.yellow,
                     edgecolor='k')
    ax1.fill_between(alphas,
                     np.array(EU_wind_LCOE),
                     label='Wind capacity', color=fig.blue,
                     edgecolor='k')


    ax2.fill_between(alphas,
                     np.array(BE_LCOE) + np.array(BC_LCOE) + \
                     np.array(solar_LCOE) + np.array(wind_LCOE) +\
                     np.array(TC_LCOE), label='Backup energy', color=fig.orange,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(BC_LCOE) +
                     np.array(solar_LCOE) + np.array(wind_LCOE) +
                     np.array(TC_LCOE), label='Backup capacity', color=fig.red,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(solar_LCOE) + np.array(wind_LCOE) +
                     np.array(TC_LCOE), label='Solar capacity', color=fig.yellow,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(wind_LCOE) +
                     np.array(TC_LCOE), label='Wind capacity', color=fig.blue,
                     edgecolor='k')
    ax2.fill_between(alphas, np.array(TC_LCOE), label='Transmission capacity',\
                      color=fig.green, edgecolor='k')

    ax2.plot(alphas, zerotrans_total_LCOE, color='w', lw=2, ls='--',\
             label="Total LCOE, zero transmission")

    colors = [fig.orange, fig.red, fig.yellow, fig.blue, fig.green]
    rectangles = [plt.Rectangle((0,0), 1, 1, fc=c) for c in colors]
    ax1.legend(rectangles, ['Backup energy', 'Backup capacity', \
                            'Solar capacity', 'Wind capacity',\
                            'Transmission capacity'])#, prop={'size':12})

    ax1.set_ylim(0,200)
    ax1.set_xlabel('Wind/solar mix')
    ax1.set_ylabel('LCOE [' + u'\u20AC' + '/MWh]')
    ax2.set_ylim(0,200)
    ax2.set_xlabel('Wind/solar mix')
    ax2.set_ylabel('LCOE [' + u'\u20AC' + '/MWh]')
    ax1.text(0.04, 185, '(a)', fontdict={'fontsize':20})
    ax2.text(0.04, 185, '(b)', fontdict={'fontsize':20})
    plt.tight_layout()

    if not interactive:
        figfilename = 'LCOEvsAlpha.pdf'
        plt.savefig(savepath+figfilename)
