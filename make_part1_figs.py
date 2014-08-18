import figutils as fig
from FlowCalculation import FlowCalculation

savepath = './results/figures/FigsOfInterestJul10/'
########### This first section is for hour-dependent flowhistograms #########

## illustrate alphadepencence
fig.make_hourly_flowhists(FlowCalculation('eurasia', 'aHO1', 'copper', 'sqr'),
                      links=['EU to ME', 'SE to CN'], hours=[0, 6, 12, 18],
                      figfileending='alphadependence', savepath=savepath)

## illustrate open vs closed US_eurasiagrid
fig.make_hourly_flowhists(
            [FlowCalculation('US_eurasia_open', 'aHO1', 'copper', 'sqr'),
             FlowCalculation('US_eurasia_closed', 'aHO1', 'copper', 'sqr')],
                      links=['EU to US', 'SE to CN'], hours=[6],
                      figfileending='openvsclosed', savepath=savepath)

## illustrate lin vs sqr flow
fig.make_hourly_flowhists(
             [FlowCalculation('eurasia', 'aHO1', 'copper', 'lin'),
              FlowCalculation('eurasia', 'aHO1', 'copper', 'sqr')],
                      links=['EU to ME', 'SE to CN'], hours=[0, 12],
                      figfileending='linvssqr', savepath=savepath)
## illustrate effect of impedance
fig.make_hourly_flowhists(
             [FlowCalculation('eurasia', 'aHO1', 'copper', 'sqr'),
              FlowCalculation('eurasia', 'aHO1', 'copper', 'sqr_imp')],
                      links=['EU to ME', 'SE to CN', 'EU to RU'],
                      hours=[0, 12],
                      figfileending='impnoimp', savepath=savepath)

"""

########### This second section is for BC, BE, TC vs alpha plots ############

## illustrate lin vs sqr for balancing capacity
fig.make_y_vs_alpha_graph(
             [FlowCalculation('eurasia', 'aHO1', 'copper', 'lin'),
              FlowCalculation('eurasia', 'aHO1', 'copper', 'sqr')],
             ydatalabel='BC', savepath=savepath,
             figfilename='eurasia_BC_vs_alpha_linsqr.pdf', zerotrans=True,
             title=False, small_legend=True,
             labels=['No transmission', 'Localized flow', \
                     'Localized flow: ' + r'$\alpha^W_{\mathrm{opt}}$',
                     'Synchronized flow',\
                     'Synchronized flow: ' + r'$\alpha^W_{\mathrm{opt}}$'])

## Balancing capacity for all layouts
layouts = ['EU_RU_NA_ME', 'US_EU_RU_NA_ME', 'eurasia', 'US_eurasia_open']
alllayoutflowcalcs = [FlowCalculation(l, 'aHO1', 'copper', 'sqr') \
                      for l in layouts]
fig.make_y_vs_alpha_graph(alllayoutflowcalcs,
             ydatalabel='BC', savepath=savepath,
             figfilename='BC_vs_alpha_all_layouts.pdf', zerotrans=False,
             hetpoints=False, small_legend=False,
             labels= ['EU-RU-NA-ME', 'US-EU-RU-NA-ME',\
                      'Eurasia', 'US-Eurasia-open'])

## Balancing energy for eurasia
fig.make_y_vs_alpha_graph(
             FlowCalculation('eurasia', 'aHO1', 'copper', 'sqr'),
             ydatalabel='BE', savepath=savepath,
             figfilename='eurasia_BE_vs_alpha.pdf', zerotrans=True,
             labels=['No transmission', 'Copper flow', 'Copper flow: '\
                     + r'$\alpha^W_{\mathrm{opt}}$'],
             title=False, small_legend=False)

## Balancing energy for all layouts
fig.make_y_vs_alpha_graph(alllayoutflowcalcs,
             ydatalabel='BE', savepath=savepath,
             figfilename='BE_vs_alpha_all_layouts.pdf', zerotrans=False,
             hetpoints=False, small_legend=False,
             labels= ['EU-RU-NA-ME', 'US-EU-RU-NA-ME',\
                      'Eurasia', 'US-Eurasia-open'])

## Transmission capacity lin vs sqr
fig.make_y_vs_alpha_graph(
             [FlowCalculation('eurasia', 'aHO1', 'copper', 'lin'),
              FlowCalculation('eurasia', 'aHO1', 'copper', 'sqr')],
             ydatalabel='Total_TC', savepath=savepath,
             figfilename='eurasia_TC_vs_alpha_linsqr.pdf',
             title=False, small_legend=False,
             labels=['Localized flow', \
                     'Localized flow: ' + r'$\alpha^W_{\mathrm{opt}}$',
                     'Synchronized flow',\
                     'Synchronized flow: ' + r'$\alpha^W_{\mathrm{opt}}$'])

## Transmission capacity open vs closed
fig.make_y_vs_alpha_graph(
            [FlowCalculation('US_eurasia_open', 'aHO1', 'copper', 'sqr'),
             FlowCalculation('US_eurasia_closed', 'aHO1', 'copper', 'sqr')],
             ydatalabel='Total_TC', savepath=savepath,
             figfilename='TC_vs_alpha_openclosed.pdf',
             small_legend=False,
             labels=['US-Eurasia-open',
                    'US-Eurasia-open: ' + r'$\alpha^W_{\mathrm{opt}}$',
                    'US-Eurasia-closed',
                    'US-Eurasia-closed: ' + r'$\alpha^W_{\mathrm{opt}}$'])


######### This section is for BC and BE vs TC graphs ######################
## Balancing capacity vs transmission capacity, lin vs sqr
fig.make_bal_vs_trans_graph(
            [FlowCalculation('eurasia', 'aHE', 'copper', 'lin'),
             FlowCalculation('eurasia', 'aHE', 'copper', 'sqr')],
            ydatalabel='BC', region='all', figfilename='BCvsTClinsqr.pdf',
            savepath=savepath, title=False, ylim=(0,0.9))

## example of balancing energy vs transmission capacity
fig.make_bal_vs_trans_graph(
             FlowCalculation('eurasia', 'aHE', 'copper', 'sqr'),
            ydatalabel='BE', region='all', figfilename='BEvsTC.pdf',
            savepath=savepath, legend=False, title=False, ylim=(0,0.18))

######### This section is for BC and BE vs layout barplots ###############
fig.make_bal_vs_layout_barplot(
            [FlowCalculation('EU_RU_NA_ME', 'aHE', '1.5q99', 'sqr'),
             FlowCalculation('US_EU_RU_NA_ME', 'aHE', '1.5q99', 'sqr'),
             FlowCalculation('eurasia', 'aHE', '1.5q99', 'sqr'),
             FlowCalculation('US_eurasia_open', 'aHE', '1.5q99', 'sqr')],
             ydatalabel='BE', savepath=savepath,
             title='Copper flow')

fig.make_bal_vs_layout_barplot(
            [FlowCalculation('EU_RU_NA_ME', 'aHE', '1.5q99', 'sqr'),
             FlowCalculation('US_EU_RU_NA_ME', 'aHE', '1.5q99', 'sqr'),
             FlowCalculation('eurasia', 'aHE', '1.5q99', 'sqr'),
             FlowCalculation('US_eurasia_open', 'aHE', '1.5q99', 'sqr')],
             ydatalabel='BC', savepath=savepath,
             title='Synchronized copper flow')
"""

fig.make_LCOE_vs_alpha_graph(
        FlowCalculation('eurasia', 'aHO1', 'copper', 'lin'),
        savepath=savepath, title=False)
