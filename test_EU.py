import numpy as np
import aurespf.solvers as au
from EUgrid import EU_Nodes

EU = EU_Nodes()

solvedEU, solvedflows = au.solve(EU, mode="copper capped rolando verbose")

solvedEU.save_nodes("eu_aHE_copper_capR")
np.save("./results/eu_aHE_copper_capR_flows.npy", solvedflows)
