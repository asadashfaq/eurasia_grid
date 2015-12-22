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

green1 = '#7bed08'
green2 = '#5eb406'
green3 = '#407b04'
green4 = '#234202'

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
    BC_sqr = []
    TC = []
    TC_sqr = []

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
    TC.append(1e-8)
    TC_sqr.append(1e-8)

    for fc in flowcalcs:
        fc_sqr = FlowCalculation(fc.layout, 'aHE', '1.5q99', 'sqr')
        print fc.layout
        admat = "./settings/" + fc.layout + "admat.txt"
        N = nh_Nodes(admat=admat)
        total_mean_load = np.sum([n.mean for n in N])
        print total_mean_load

        filename = str(fc) + '.pkl'
        filename_sqr = str(fc_sqr) + '.pkl'
        flowfilename = './results/' + fc.layout + '_aHE_copper_lin_flows.npy'
        flowfilename_sqr = './results/' + fc.layout + '_aHE_copper_sqr_flows.npy'
        unnormalized_BE = [fig.get_data(filename, 'BE', \
                             path=datapath)[n.id] for n in N]
        unnormalized_BC = [fig.get_data(filename, 'BC', \
                             path=datapath)[n.id] for n in N]
        unnormalized_BC_sqr = [fig.get_data(filename_sqr, 'BC', \
                             path=datapath)[n.id] for n in N]
        BE.append(np.sum(unnormalized_BE)/total_mean_load)
        BC.append(np.sum(unnormalized_BC)/total_mean_load)
        BC_sqr.append(np.sum(unnormalized_BC_sqr)/total_mean_load)
        LI = au.linfo(admat)
        energywiseTC = au.biggestpair(au.get_quant_caps(filename=flowfilename))
        energywiseTC_sqr = au.biggestpair(au.get_quant_caps(filename=flowfilename_sqr))
        print np.sum(energywiseTC)
        TC.append(np.sum([energywiseTC[i]*float(LI[i][2]) \
                            for i in range(len(LI))])/(1e3*total_mean_load))
        TC_sqr.append(np.sum([energywiseTC_sqr[i]*float(LI[i][2]) \
                            for i in range(len(LI))])/(1e3*total_mean_load))
        print TC
    index1 = np.arange(len(BE))
    left = index1 + 0.5*barwidth
    left2 = left + 0.5*barwidth
    plt.ion()
    ax1 = plt.subplot(3,1,1)
    if not interactive:
        plt.gcf().set_size_inches([7.5,12])
        plt.gcf().set_dpi(400)
    ax1.bar(left, BE, width=barwidth, color=fig.orange,\
                label='Localized/synchronized flow')
    layoutlist1 = ['EU', 'EU-RU-NA-ME', 'US-EU-RU-NA-ME', 'Eurasia', 'US-Eurasia']
    ax1.set_xticks(left + 0.5*barwidth)
    ax1.set_xticklabels(layoutlist1)
    ax1.set_ylabel('Backup energy [normalized]')#r'$E_\mathrm{total}^B$' + ' [normalized]')
    ax1.set_ylim((0,0.2))
    ax1.text(0.1*barwidth, 0.182, '(a)',fontdict={'fontsize':20})
    ax1.legend()

    ax2 = plt.subplot(3,1,2)
    ax2.bar(left[0], BC[0], width=barwidth, color=fig.red)
    ax2.bar(left[1:], BC[1:], width=0.5*barwidth, color=fig.red,\
            label='Localized flow')
    ax2.bar(left2[1:], BC_sqr, width=0.5*barwidth, color=fig.darkred,\
            label='Synchronized flow')
    ax2.set_xticks(left + 0.5*barwidth)
    ax2.set_xticklabels(layoutlist1)
    ax2.set_ylabel('Backup capacity [normalized]')#r'$C_\mathrm{total}^B$'+ ' [normalized]')
    ax2.set_ylim((0.0, 0.82))
    ax2.text(0.1*barwidth, 0.748, '(b)',fontdict={'fontsize':20})
    ax2.legend(ncol=2)

    layoutlist2 = ['', 'EU-RU-NA-ME', 'US-EU-RU-NA-ME', 'Eurasia', 'US-Eurasia']
    ax3 = plt.subplot(3,1,3)
    ax3.bar(left, TC, width=0.5*barwidth, color=fig.green, \
             label = 'Localized flow')
    ax3.bar(left2, TC_sqr, width=0.5*barwidth, color=green4,
             label = 'Synchronized flow')
    ax3.set_xticks(left + 0.5*barwidth)
    ax3.set_xticklabels(layoutlist2)
    ax3.set_ylabel('Transmission capacity\n [normalized'\
            + r'$\times$' + '1000 km]')#r'$C_\mathrm{total}^T$' + ' [km]')
    ax3.text(0.1*barwidth, 5.46, '(c)',fontdict={'fontsize':20})
    ax3.legend(ncol=2)

    plt.tight_layout()

    figfilename = 'BEBCTC_vs_layout.pdf'

    if not interactive:
        plt.savefig(savepath+figfilename)


def make_LCOE_barplot(interactive=True, solvermode='lin', ax=None):
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
    LCOE050 = []
    LCOE025 = []
    LCOE015 = []

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

    total_LCOE.append(0)
    LCOE2.append(sum([EU_BE_LCOE, EU_BC_LCOE, EU_wind_LCOE, EU_solar_LCOE]))
    LCOE3.append(sum([EU_BC_LCOE, EU_wind_LCOE, EU_solar_LCOE]))
    LCOE4.append(sum([EU_wind_LCOE, EU_solar_LCOE]))
    LCOE5.append(sum([EU_wind_LCOE]))
    LCOE050.append(0)
    LCOE025.append(0)
    LCOE015.append(0)


    fclist = [FlowCalculation(l, 'aHE', 'copper', solvermode) for l in layouts]
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
        TC_LCOE = au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat,\
                scale_factor=1)/total_energy
        TC_LCOE50 = au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat,\
                scale_factor=0.5)/total_energy
        TC_LCOE25 = au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat,\
                scale_factor=0.25)/total_energy
        TC_LCOE15 = au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat,\
                scale_factor=0.15)/total_energy

##################### Make TC_LCOE75 etc variablse
        total_LCOE.append(sum([TC_LCOE, BE_LCOE, BC_LCOE, wind_LCOE, solar_LCOE]))
        LCOE050.append(sum([TC_LCOE50, BE_LCOE, BC_LCOE, wind_LCOE, solar_LCOE]))
        LCOE025.append(sum([TC_LCOE25, BE_LCOE, BC_LCOE, wind_LCOE, solar_LCOE]))
        LCOE015.append(sum([TC_LCOE15, BE_LCOE, BC_LCOE, wind_LCOE, solar_LCOE]))
        LCOE2.append(sum([BE_LCOE, BC_LCOE, wind_LCOE, solar_LCOE]))
        LCOE3.append(sum([BC_LCOE, wind_LCOE, solar_LCOE]))
        LCOE4.append(sum([wind_LCOE, solar_LCOE]))
        LCOE5.append(sum([wind_LCOE]))
        #LCOE2.append(sum([BC_LCOE, wind_LCOE, solar_LCOE, TC_LCOE]))
        #LCOE3.append(sum([wind_LCOE, solar_LCOE, TC_LCOE]))
        #LCOE4.append(sum([wind_LCOE, TC_LCOE]))
        #LCOE5.append(sum([TC_LCOE]))


        print total_LCOE

    barwidth = 0.5
    if ax==None:
        ax = plt.subplot(1,1,1)
    index = np.arange(len(total_LCOE))
    left = index + 0.5*barwidth
    ax.bar(left, total_LCOE, width=barwidth, color=green1, label='Transmission capacity')
    ax.bar(left, LCOE050, width=barwidth, color=green2)
    ax.bar(left, LCOE025, width=barwidth, color=green3)
    ax.bar(left, LCOE015, width=barwidth, color=green4)
    ax.bar(left, LCOE2, width=barwidth, color=fig.orange, label='Backup energy')
    ax.bar(left, LCOE3, width=barwidth, color=fig.red, label='Backup capacity')
    ax.bar(left, LCOE4, width=barwidth, color=fig.yellow, label='Solar capacity')
    ax.bar(left, LCOE5, width=barwidth, color=fig.blue, label='Wind capacity')

    ax.set_xticks(left + 0.5*barwidth)
    layoutlist = ['EU', 'EU-RU-NA-ME', 'US-EU-RU-NA-ME', 'Eurasia', 'US-Eurasia']
    ax.set_xticklabels(layoutlist)
    ax.set_ylabel('LCOE [' + u'\u20AC' + '/MWh]')
    ax.legend(loc=2, prop={'size':12})
    ax.set_ylim(0,100)

    if not interactive:
        plt.tight_layout()
        figfilename = 'LCOEvslayout.pdf'
        plt.savefig(savepath+figfilename)

def LCOE_vs_layout_double():
    plt.figure(figsize=(8,12))
    ax1 = plt.subplot(2,1,1)
    make_LCOE_barplot(solvermode='lin', ax=ax1)
    ax1.text(1.85,92, 'a) Localized flow', size=14)
    ax2 = plt.subplot(2,1,2)
    make_LCOE_barplot(solvermode='sqr', ax=ax2)
    ax2.text(1.85,92, 'b) Synchronized flow', size=14)
    plt.tight_layout()
    plt.savefig(savepath + 'LCOE_vs_layout_double.pdf')



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
    TC050_LCOE = []
    TC025_LCOE = []
    TC015_LCOE = []
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
        TC050_LCOE.append(au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat,\
                scale_factor=0.5)/total_energy)
        TC025_LCOE.append(au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat,\
                scale_factor=0.25)/total_energy)
        TC015_LCOE.append(au.ctc(ct.get_TCs(fc, datapath), pathadmat=admat,\
                scale_factor=0.15)/total_energy)
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
    if not interactive:
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
                     np.array(TC_LCOE) + np.array(BE_LCOE) +\
                     np.array(BC_LCOE) + \
                     np.array(solar_LCOE) + np.array(wind_LCOE),\
                     label='Transmission capacity', color=green1,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(TC050_LCOE) + np.array(BE_LCOE) +\
                     np.array(BC_LCOE) + \
                     np.array(solar_LCOE) + np.array(wind_LCOE),\
                     color=green2,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(TC025_LCOE) + np.array(BE_LCOE) +\
                     np.array(BC_LCOE) + \
                     np.array(solar_LCOE) + np.array(wind_LCOE),\
                     color=green3,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(TC015_LCOE) + np.array(BE_LCOE) +\
                     np.array(BC_LCOE) + \
                     np.array(solar_LCOE) + np.array(wind_LCOE),\
                     color=green4,
                     edgecolor='k')

    ax2.fill_between(alphas,
                     np.array(BE_LCOE) + np.array(BC_LCOE) + \
                     np.array(solar_LCOE) + np.array(wind_LCOE),\
                     label='Backup energy', color=fig.orange,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(BC_LCOE) +
                     np.array(solar_LCOE) + np.array(wind_LCOE),\
                     label='Backup capacity', color=fig.red,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(solar_LCOE) + np.array(wind_LCOE),
                     label='Solar capacity', color=fig.yellow,
                     edgecolor='k')
    ax2.fill_between(alphas,
                     np.array(wind_LCOE),
                     label='Wind capacity', color=fig.blue,
                     edgecolor='k')

    ax2.plot(alphas, zerotrans_total_LCOE, color='w', lw=2, ls='--',\
             label="Total LCOE, zero transmission")

    colors = [green1, fig.orange, fig.red, fig.yellow, fig.blue]
    rectangles = [plt.Rectangle((0,0), 1, 1, fc=c) for c in colors]
    ax1.legend(rectangles, ['Transmission capacity', 'Backup energy',\
                             'Backup capacity', \
                            'Solar capacity', 'Wind capacity'])


    ax1.set_ylim(0,200)
    ax1.set_xlabel('Wind/solar mix')
    ax1.set_ylabel('LCOE [' + u'\u20AC' + '/MWh]')
    ax2.set_ylim(0,200)
    ax2.set_xlabel('Wind/solar mix')
    ax2.set_ylabel('LCOE [' + u'\u20AC' + '/MWh]')
    ax1.text(0.04, 185, '(a) EU isolated', fontdict={'fontsize':20})
    ax2.text(0.04, 185, '(b) US-Eurasia', fontdict={'fontsize':20})
    plt.tight_layout()

    if not interactive:
        figfilename = 'LCOEvsAlpha.pdf'
        plt.savefig(savepath+figfilename)
