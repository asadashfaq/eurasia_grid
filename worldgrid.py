import numpy as np
import aurespf.solvers as au

def world_Nodes(admat=None, load_filename=None, full_load=False, alphas=None):

    if alphas==None:
        alphas = np.load("./data/wOptimalMixes.npy")

    if not admat:
        admat = './settings/wadmat.txt' # the whole world of connected
                                        # superregions, as in Martins
                                        # Frankfurt talk

    regions = ['EU', 'RU', 'NA', 'ME', 'IN', 'SE', 'CN', 'JK']
    files = [''.join([r,'.npz']) for r in regions]

    prefix='VE_'
    admat = admat

    return au.Nodes(admat=admat, path='./data/', prefix=prefix,
                    files=['EU.npz', 'RU.npz', 'NA.npz', 'ME.npz', 'IN.npz',
                        'SE.npz', 'CN.npz', 'JK.npz'], load_filename=load_filename,
                    full_load=full_load, alphas=alphas,
                    gammas=np.ones(8))

