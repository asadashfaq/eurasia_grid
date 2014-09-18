import aurespf.solvers as au
import figutils as fig
from FlowCalculation import FlowCalculation
from europe_plusgrid import europe_plus_Nodes

layouts = ['Europe_RU', 'Europe_NA', 'Europe_ME', 'Europe_RU_NA_ME']
modes = ['lin', 'sqr']

datapath = './results/Europeplus/'

Europe = europe_plus_Nodes(admat='./settings/Europeadmat.txt')
internal_links = [au.AtoKh(Europe)[-1][i][0]\
        for i in range(len(au.AtoKh(Europe)[-1]))]

for layout in layouts:
    admat = './settings/' + layout + 'admat.txt'
    print admat
    N = europe_plus_Nodes(admat=admat)
    for mode in modes:
        fc = FlowCalculation(layout, 'aHE', 'copper', mode)
        filename = str(fc) + '.pkl'
        internal_indices, external_indices = fig.get_internal_external_link_indices(N, internal_links)
        all_TC = au.biggestpair(fig.get_data(filename, 'TC', path=datapath))
        print fc.pretty_layout(), fc.pretty_solvermode()
        for link_index in external_indices:
            print au.AtoKh(N)[-1][link_index]
            print all_TC[link_index]


