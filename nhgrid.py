import numpy as np
import aurespf.solvers as au

def nh_Nodes(admat=None, load_filename=None, full_load=False, alphas=None):

    if alphas==None:
        alphas = np.load("./results/nhOptimalMixes.npy")

    if not admat and not load_filename:
        admat = './settings/US_eurasia_openadmat.txt' # the whole world of connected
                                        # superregions, as in Martins
                                        # Frankfurt talk, plus the US connected
                                        # to Europe with a link

    regions = ['EU', 'RU', 'NA', 'ME', 'IN', 'SE', 'CN', 'JK', 'US']
    files = [''.join([r,'.npz']) for r in regions]

    prefix='VE_'
    admat = admat

    return au.Nodes(admat=admat, path='./data/', prefix=prefix,
                    files=files, load_filename=load_filename,
                    full_load=full_load, alphas=alphas,
                    gammas=np.ones(9))

