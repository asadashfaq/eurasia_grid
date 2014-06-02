import numpy as np
from EUgrid import EU_Nodes

N = EU_Nodes()

wind = np.zeros_like(N[0].get_wind())
solar = np.zeros_like(wind)
load = np.zeros_like(wind)

for n in N:
    wind = wind + n.get_wind()
    solar = solar + n.get_solar()
    load = load + n.load

t = np.linspace(0, len(load)-1, len(load))
wind = wind/np.mean(wind)
solar = solar/np.mean(solar)
np.savez("./data/VE_EU", Gw=wind, Gs=solar, L=load, t=t, datalabel='EU')
