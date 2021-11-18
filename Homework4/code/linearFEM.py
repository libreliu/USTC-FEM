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

    def phiF_intg(self, f, faceIdx, nodeIdx):
        # determine nodeLocalIdx: 0~2
        faceVertices = self.mesh.fv_indices[faceIdx]
        nodeLocalIdx = None
        for i in range(0, 3):
            if faceVertices[i] == nodeIdx:
                nodeLocalIdx = i
        assert(nodeLocalIdx != None)
        fVCoords = []

        for localIdx in range(0, 3):
            fVCoords.append(self.mesh.points[faceVertices[localIdx]])
        
        # (x1-x3)*(y2-y3) - (x2-x3)(y1-y3)
        X = lambda idx: fVCoords[idx-1][0]
        Y = lambda idx: fVCoords[idx-1][1]
        jacobian = (X(1)-X(3))*(Y(2)-Y(3)) - (X(2)-X(3))*(Y(1)-Y(3))

        detA = X(1)*Y(2) - X(1)*Y(3) - X(2)*Y(1) + X(2)*Y(3) + X(3)*Y(1) - X(3)*Y(2)
        area = 0.5 * np.abs(detA)
        # scheme: evaluate on three midpoints (order 1 method)
        # (0.5, 0), (0, 0.5), (0.5, 0.5) on barycentric coordinate
        res = 0
        for xi1, xi2 in [(0.5, 0), (0, 0.5), (0.5, 0.5)]:
            


    def solve(self, f: 'function'):
        numNodes = len(self.mesh.nodes)
        K = np.zeros((numNodes, numNodes), dtype=np.float64)
        F = np.zeros((numNodes,), dtype=np.float64)

        for nodeIdx in self.mesh.nodes:
            neighFaces = self.mesh.get_adjacent_elements(nodeIdx)
            for faceIdx in neighFaces:
                # compute F_i
                F[nodeIdx] += phiF_intg(f, faceIdx, nodeIdx)
                