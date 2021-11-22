#!/usr/bin/env python3
import logging, os
from linearFEM import LinearFEM
import math
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np

logger = logging.getLogger(__name__)

def get_order(errk, errkP1, Ek, EkP1):
    return 2 * (math.log(errk / errkP1)) / (math.log(EkP1 / Ek))

def evaluate():
    meshes = [f"gd{i}" for i in range(0, 5)]

    trials = []
    for mesh in meshes:
        fem = LinearFEM()
        fem.prepareMesh(mesh)
        faceCount = fem.mesh.num_elements

        f = lambda x, y: (2+math.pi**2 * (1-y)*y) * math.sin(math.pi*x)
        u = lambda x, y: math.sin(math.pi * x) * (1-y) * y

        fem.solve(f)

        _, err_Linf_bary, err_L2_bary = fem.errorAnalysisBarycentric(u, 16)

        (X, Y, Z), err_Linf_uni, err_L2_uni = (1,2,3), 1, 1
        #(X, Y, Z), err_Linf_uni, err_L2_uni = fem.errorAnalysisUniform(u, 100)
        
        trials.append({
            'name': mesh,
            'faceCount': faceCount,
            'err_Linf_bary': err_Linf_bary,
            'err_L2_bary': err_L2_bary,
            'err_Linf_uni': err_Linf_uni,
            'err_L2_uni': err_L2_uni,
        })

        # Visualize
        # fig, ax = plt.subplots()
        # im = ax.imshow(np.abs(Z), interpolation='bilinear', origin='lower')

        # CBI = fig.colorbar(im)

        # plt.show()

    for idx, trial in enumerate(trials):
        if idx == 0:
            EkP1 = trial['faceCount']
            ord_Linf_bary = 'N/A'
            ord_L2_bary = 'N/A'
            ord_Linf_uni = 'N/A'
            ord_L2_uni = 'N/A'
        else:
            Ek = trials[idx - 1]['faceCount']
            EkP1 = trial['faceCount']
            ord_Linf_bary = get_order(
                trials[idx-1]['err_Linf_bary'],
                trials[idx]['err_Linf_bary'],
                Ek,
                EkP1
            )
            ord_L2_bary = get_order(
                trials[idx-1]['err_L2_bary'],
                trials[idx]['err_L2_bary'],
                Ek,
                EkP1
            )
            ord_Linf_uni = get_order(
                trials[idx-1]['err_Linf_uni'],
                trials[idx]['err_Linf_uni'],
                Ek,
                EkP1
            )
            ord_L2_uni = get_order(
                trials[idx-1]['err_L2_bary'],
                trials[idx]['err_L2_bary'],
                Ek,
                EkP1
            )

        print(f"|{EkP1}|{trials[idx]['err_Linf_bary']}|{trials[idx]['err_L2_bary']}|{ord_Linf_bary}|{ord_L2_bary}|")



if __name__ == "__main__":
    evaluate()

if 0:
    logging.basicConfig(level=logging.INFO)

    fem = LinearFEM()
    fem.prepareMesh("gd1")
    #fem.visualize()

    #f = lambda x, y: x * (x-1) * math.sin(math.pi * y)
    f = lambda x, y: (2+math.pi**2 * (1-y)*y) * math.sin(math.pi*x)

    logger.info("Begin FEM solve()")
    fem.solve(f)
    logger.info("FEM solve() finished")
    #fem.visualize()

    u = lambda x, y: math.sin(math.pi * x) * (1-y) * y
    ret = fem.errorAnalysisBarycentric(u, 10000)
    #print(ret)

    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111, projection='3d')

    # ax.scatter([elem[0] for elem in ret[0]], [elem[1] for elem in ret[0]], [elem[2] for elem in ret[0]])
    # plt.show()

    X, Y, Z = fem.errorAnalysisUniform(u, 10)
    surf = ax.plot_surface(X, Y, Z)

    plt.show()
