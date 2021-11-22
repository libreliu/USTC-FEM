#!/usr/bin/env python3
import logging, os
from linearFEM import LinearFEM
import math
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np

logger = logging.getLogger(__name__)

def evaluate():
    meshes = [f"gd{i}" for i in range(0, 5)]

    trials = {}
    for mesh in meshes:
        fem = LinearFEM()
        fem.prepareMesh(mesh)
        faceCount = fem.mesh.num_elements

        f = lambda x, y: (2+math.pi**2 * (1-y)*y) * math.sin(math.pi*x)
        u = lambda x, y: math.sin(math.pi * x) * (1-y) * y

        fem.solve(f)

        _, err_Linf_bary, err_L2_bary = fem.errorAnalysisBarycentric(u, 16)

        (X, Y, Z), err_Linf_uni, err_L2_uni = fem.errorAnalysisUniform(u, 100)
        
        # Visualize
        fig, ax = plt.subplots()
        im = ax.imshow(np.abs(Z), interpolation='bilinear', origin='lower')
        
        CBI = fig.colorbar(im)

        plt.show()



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
