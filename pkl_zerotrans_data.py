import numpy as np
from FlowCalculation import FlowCalculation
from FCResult import FCResult
from nhgrid import nh_Nodes

layouts = ['EU_RU_NA_ME', 'eurasia', 'US_eurasia_open', \
           'US_eurasia_closed', 'US_EU_RU_NA_ME']
alphas = np.linspace(0, 1, 21)
savepath = './results/AlphaSweepsZerotrans/'

for l in layouts:
    admat = './settings/' + l + 'admat.txt'
    for a in alphas:
        N = nh_Nodes(admat=admat, alphas=a)
        fc = FlowCalculation(l, 'aHO' + str(a), 'zerotrans', 'raw')
        filename = str(fc) + '.pkl'

        myresult = FCResult(filename, savepath)
        myresult.add_instance(N, [], fc)
        myresult.save_results(filename, savepath)

