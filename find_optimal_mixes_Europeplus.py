import numpy as np
import aurespf.solvers as au
regions = ['AUT', 'FIN', 'NLD', 'BIH', 'FRA', 'NOR', 'BEL','GBR', \
                 'POL', 'BGR', 'GRC', 'PRT', 'CHE', 'HRV', 'ROU', 'CZE',\
                 'HUN', 'SRB', 'DEU', 'IRL', 'SWE', 'DNK', 'ITA', 'SVN',\
                 'ESP', 'LUX', 'SVK', 'EST', 'LVA', 'LTU', 'RU', 'NA', 'ME']


optimalmixes = []
for r in regions:
    filename = ''.join(['data/VE_', r, '.npz'])
    data = np.load(filename)
    load = data['L']
    Gw = data['Gw']
    Gs = data['Gs']

    optimalmixes.append(au.optimal_mix_balancing(load, Gw, Gs)[0])

np.save("results/europeplusOptimalmixes.npy", optimalmixes)

