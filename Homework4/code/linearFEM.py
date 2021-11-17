import logging, os

import numpy as np
import meshUtil

logger = logging.getLogger(__name__)

class LinearFEM:
    def __init__(self):
        self.mesh = None
        self.free_nodes = None

    def prepare_mesh(self, meshPrefix):
        cwd = os.getcwd()
        self.mesh = meshUtil.Mesh2D.from_easymesh(
            os.path.join(cwd, f"./meshes/{meshPrefix}.n"),
            os.path.join(cwd, f"./meshes/{meshPrefix}.e"),
            os.path.join(cwd, f"./meshes/{meshPrefix}.s")
        )

    def solve(self):
        numNodes = len(self.mesh.nodes)
        K = np.zeros((numNodes, numNodes), dtype=np.float64)
        F = np.zeros((numNodes,), dtype=np.float64)

        for nodeIdx in self.mesh.nodes:
            neighFaces = self.mesh.get_adjacent_elements(nodeIdx)
            for faceIdx in neighFaces:
                # TODO: compute F_i

                