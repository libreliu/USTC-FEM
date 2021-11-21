#!/usr/bin/env python3
import logging, os
from linearFEM import LinearFEM
import math
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    fem = LinearFEM()
    fem.prepareMesh("gd0")
    #fem.visualize()

    #f = lambda x, y: x * (x-1) * math.sin(math.pi * y)
    f = lambda x, y: (2+math.pi**2 * (1-y)*y) * math.sin(math.pi*x)

    logger.info("Begin FEM solve()")
    fem.solve(f)
    logger.info("FEM solve() finished")
    fem.visualize()

    u = lambda x, y: math.sin(math.pi * x) * (1-y) * y
    ret = fem.errorAnalysisBarycentric(u, 8)
    #print(ret)

    plt.scatter(
        [elem[0] for elem in ret[0]],
        [elem[1] for elem in ret[0]],
        color = "green"
    )
    plt.show()