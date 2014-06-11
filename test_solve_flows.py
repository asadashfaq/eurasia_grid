import FlowCalculation as fc
from solve_flows import solve_flow

flowCalc = fc.FlowCalculation('eurasia', 'aHE', 'copper', 'lin', 'FCResult')

print(str(flowCalc))
print(flowCalc.basisnetwork)

solve_flow(flowCalc)
