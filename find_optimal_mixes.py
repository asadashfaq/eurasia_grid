import numpy as np
import aurespf.solvers as au

regions = ['EU', 'RU', 'NA', 'ME', 'IN', 'SE', 'CN', 'JK', 'US']

nhOptimalmixes = []
for r in regions:
    filename = ''.join(['data/VE_', r, '.npz'])
    data = np.load(filename)
    load = data['L']
    Gw = data['Gw']
    Gs = data['Gs']

    nhOptimalmixes.append(au.optimal_mix_balancing(load, Gw, Gs)[0])

np.save("results/nhOptimalmixes.npy", nhOptimalmixes)

