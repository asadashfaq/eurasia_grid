import numpy as np
import aurespf.solvers as au

def nh_Nodes(admat=None, load_filename=None, full_load=False, alphas=None):


    if not admat and not load_filename:
        admat = './settings/US_eurasia_openadmat.txt' # the whole world of connected
                                        # superregions, as in Martins
                                        # Frankfurt talk, plus the US connected
                                        # to Europe with a link

    all_regions = ['EU', 'RU', 'NA', 'ME', 'IN', 'SE', 'CN', 'JK', 'US']
    regions = []
    if admat:
        if 'eurasia' in admat:
            regions.extend(['EU', 'RU', 'NA', 'ME', 'IN', 'SE', 'CN', 'JK'])
        if 'EU_RU_NA_ME' in admat:
            regions.extend(['EU', 'RU', 'NA', 'ME'])
        if 'US' in admat:
            regions.append('US')
    elif load_filename:
        if 'eurasia' in load_filename:
            regions.extend(['EU', 'RU', 'NA', 'ME', 'IN', 'SE', 'CN', 'JK'])
        if 'EU_RU_NA_ME' in load_filename:
            regions.extend(['EU', 'RU', 'NA', 'ME'])
        if 'US' in load_filename:
            regions.append('US')
    else:
        print("You must provide either an admat or a load_filename")

    files = [''.join([r,'.npz']) for r in regions]

    prefix='VE_'

    if alphas==None:
        all_alphas = np.load("./results/nhOptimalMixes.npy")
        alphas = []
        for r in regions:
            alphas.append(all_alphas[all_regions.index(r)])


    return au.Nodes(admat=admat, path='./data/', prefix=prefix,
                    files=files, load_filename=load_filename,
                    full_load=full_load, alphas=alphas,
                    gammas=np.ones(len(regions)))

