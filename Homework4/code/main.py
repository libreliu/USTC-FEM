#!/usr/bin/env python3
import logging, os
from linearFEM import LinearFEM
import math

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    fem = LinearFEM()
    fem.prepareMesh("gd0")
    fem.visualize()

    #f = lambda x, y: x * (x-1) * math.sin(math.pi * y)
    f = lambda x, y: (2+math.pi**2 * (1-y)*y) * math.sin(math.pi*x)

    logger.info("Begin FEM solve()")
    fem.solve(f)
    logger.info("FEM solve() finished")
    fem.visualize()