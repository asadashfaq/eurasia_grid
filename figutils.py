import matplotlib.pyplot as plt
import numpy as np

def BC_vs_TC(filenames = ['BC_TC_datacapM.npz', 'BC_TC_datacapR.npz'],\
        path='./results/', labels = None, savepath='./results/figures/', figfilename = 'BC_vs_TC_RolMar.pdf', BC_type='both'):
    plt.close()
    plt.rc('lines', lw=2)

    labelcount = 0
    for f in filenames:
        data = np.load(''.join([path, f]))
        BC_max = data['BC_max']
        BC_q99 = data['BC_q99']
        TC = data['TC']
        if not labels:
            label = f
        else:
            label = labels[labelcount]
        if BC_type in ['max', 'both']:
            plt.plot(1e-6*TC, BC_max, 'o', label='BCmax: '+label)
        if BC_type in ['q99', 'both']:
            plt.plot(1e-6*TC, BC_q99, 'o', label='BCq99: '+label)
        else:
            print 'Try "max", "q99" or "both" for BC_type'
            return

        labelcount = labelcount + 1

    plt.xlabel('Transmission capacity [TW]')
    plt.ylabel('Balancing capacity [normalized]')
    plt.legend()

    plt.savefig(savepath+figfilename)
    plt.close()

