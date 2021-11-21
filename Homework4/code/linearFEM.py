from ctypes import create_string_buffer
import logging, os
from platform import node

import numpy as np
import math
import meshUtil

logger = logging.getLogger(__name__)

class LinearFEM:
    def __init__(self):
        self.mesh = None
        self.free_map = None
        self.U = None

    def prepareMesh(self, meshPrefix):
        fileDir = os.path.dirname(os.path.abspath(__file__))
        self.mesh = meshUtil.Mesh2D.from_easymesh(
            os.path.join(fileDir, f"./meshes/{meshPrefix}.n"),
            os.path.join(fileDir, f"./meshes/{meshPrefix}.e"),
            os.path.join(fileDir, f"./meshes/{meshPrefix}.s")
        )
        # nodeIdx -> freeNodeIdx
        self.free_map = {}
        # freeNodeIdx -> nodeIdx
        self.free_map_inv = {}
        freeIdx = 0
        for idx in range(0, len(self.mesh.boundary_marks)):
            if self.mesh.boundary_marks[idx] == False:
                self.free_map[idx] = freeIdx
                self.free_map_inv[freeIdx] = idx
                freeIdx += 1

    def visualize(self):
        if self.U is not None:
            pointsZ = np.ndarray((self.mesh.num_points, 1), dtype=np.float64)
            for idx in range(0, self.mesh.num_points):
                if idx in self.free_map:
                    matIdx = self.free_map[idx]
                    pointsZ[idx] = self.U[matIdx]
                else:
                    # Boundary points
                    pointsZ[idx] = 0
            self.mesh.visualize(pointsZ)
        else:
            self.mesh.visualize()

    def phiFIntg(self, f, faceIdx, nodeIdx):
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
        print(f"jacob: {jacobian}")

        detA = X(1)*Y(2) - X(1)*Y(3) - X(2)*Y(1) + X(2)*Y(3) + X(3)*Y(1) - X(3)*Y(2)
        area = 0.5 * np.abs(detA)
        # scheme: evaluate on three midpoints (order 1 method)
        # (0.5, 0), (0, 0.5), (0.5, 0.5) on barycentric coordinate
        res = 0
        for xi1, xi2 in [(0.5, 0), (0, 0.5), (0.5, 0.5)]:
            xi3 = 1-xi1-xi2
            x = X(1)*xi1 + X(2)*xi2 + X(3)*xi3
            y = Y(1)*xi1 + Y(2)*xi2 + Y(3)*xi3
            pointVal = f(x, y) * (xi1, xi2, xi3)[nodeLocalIdx]
            print(f"{faceVertices}: ({xi1}, {xi2}) ({x}, {y}), nodeLocalIdx={nodeLocalIdx} contrib={pointVal}")
            res += pointVal
        
        res *= jacobian * area / 3
        print(f"{res}, area={area}")
        return res

    def phiFIntgDirect(self, f, faceIdx, nodeIdx):
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
        X = lambda idx: fVCoords[idx-1][0]
        Y = lambda idx: fVCoords[idx-1][1]

        detA = X(1)*Y(2) - X(1)*Y(3) - X(2)*Y(1) + X(2)*Y(3) + X(3)*Y(1) - X(3)*Y(2)
        area = 0.5 * np.abs(detA)

        res = 0
        for adjV in range(0, 3):
            if adjV == nodeLocalIdx:
                continue
            res += f((X(adjV+1) + X(nodeLocalIdx+1))/ 2, (Y(adjV+1)+Y(nodeLocalIdx+1)) /2)
        # res = res * area / 3 * 0.5
        res = res * area / 3 * 0.5
        print(f"{faceVertices}: nodeLocalIdx={nodeLocalIdx}, {res}, area={area}")
        return res

    def phiDerivPhiDerivIntg(self, jGIdx, iGIdx, faceIdx):
        jIdx = self.free_map[jGIdx]
        iIdx = self.free_map[iGIdx]
        assert(iIdx >= jIdx)

        faceVertices = self.mesh.fv_indices[faceIdx]
        def getLocalIndex(faceVertices, nodeIdx):
            nodeLocalIdx = None
            for i in range(0, 3):
                if faceVertices[i] == nodeIdx:
                    nodeLocalIdx = i
            assert(nodeLocalIdx != None)
            return nodeLocalIdx
        
        fVCoords = []
        for localIdx in range(0, 3):
            fVCoords.append(self.mesh.points[faceVertices[localIdx]])
        X = lambda idx: fVCoords[idx-1][0]
        Y = lambda idx: fVCoords[idx-1][1]

        jLIdx = getLocalIndex(faceVertices, jGIdx)
        iLIdx = getLocalIndex(faceVertices, iGIdx)

        detA = X(1)*Y(2) - X(1)*Y(3) - X(2)*Y(1) + X(2)*Y(3) + X(3)*Y(1) - X(3)*Y(2)
        area = 0.5 * np.abs(detA)

        res = area / (detA**2) * {
            (1, 1): (Y(2)-Y(3))**2 + (X(3)-X(2))**2,
            (1, 2): (Y(2)-Y(3))*(Y(3)-Y(1)) + (X(3)-X(2))*(X(1)-X(3)),
            (1, 3): (Y(2)-Y(3))*(Y(1)-Y(2)) + (X(3)-X(2))*(X(2)-X(1)),
            (2, 2): (Y(3)-Y(1))**2 + (X(1)-X(3))**2,
            (2, 3): (Y(3)-Y(1))*(Y(1)-Y(2)) + (X(1)-X(3))*(X(2)-X(1)),
            (3, 3): (Y(1)-Y(2))**2 + (X(2)-X(1))**2,
        }[tuple(sorted((jLIdx+1, iLIdx+1)))]
        
        return res


    def solve(self, f: 'function'):
        numNodes = len(self.mesh.nodes)
        K = np.zeros((numNodes, numNodes), dtype=np.float64)
        F = np.zeros((numNodes,), dtype=np.float64)
        #print(self.mesh.nodes)
        for nodeIdx in self.mesh.nodes:
            neighFaces = self.mesh.get_adjacent_elements(nodeIdx)
            for faceIdx in neighFaces:
                # compute F_i
                F[self.free_map[nodeIdx]] += self.phiFIntgDirect(f, faceIdx, nodeIdx)

                adjNodeIdxList = []
                faceVertices = self.mesh.fv_indices[faceIdx]
                for vIdx in faceVertices:
                    if vIdx != nodeIdx:
                        adjNodeIdxList.append(vIdx)
                assert(len(adjNodeIdxList) == 2)

                for adjNodeIdx in adjNodeIdxList:
                    if self.mesh.is_boundary_node(adjNodeIdx):
                        continue

                    if self.free_map[nodeIdx] < self.free_map[adjNodeIdx]:
                        continue

                    K[self.free_map[nodeIdx], self.free_map[adjNodeIdx]] += self.phiDerivPhiDerivIntg(adjNodeIdx, nodeIdx, faceIdx)

                    K[self.free_map[adjNodeIdx], self.free_map[nodeIdx]] = K[self.free_map[nodeIdx], self.free_map[adjNodeIdx]]
                
                K[self.free_map[nodeIdx], self.free_map[nodeIdx]] += self.phiDerivPhiDerivIntg(nodeIdx, nodeIdx, faceIdx)
        
        np.set_printoptions(precision=3, linewidth=150)
        print(K)
        print(F)
        # F = np.asarray([0.4312, 0.3011, 0.3011, 0.4312])

        self.U = np.linalg.solve(K, F)
        print(self.U)

        chk = K @ self.U
        assert(np.allclose(chk, F))

        return self.U

    def getSolvedValueBarycentric(self, faceIdx, xi1, xi2, xi3):
        faceVertices = self.mesh.fv_indices[faceIdx]
        localBasisVal = []
        for localIdx in range(0, 3):
            if faceVertices[localIdx] in self.free_map:
                localBasisVal.append(self.U[self.free_map[faceVertices[localIdx]]])
            else:
                localBasisVal.append(0.0)
        
        return xi1 * localBasisVal[0] + xi2 * localBasisVal[1] + xi3 * localBasisVal[2]

    def errorAnalysisBarycentric(self, u_ref, samplesPerCell=100):
        """Returns error w.r.t L1 norm and Linf norm evaluated on each cell"""
        samplePerEdge = math.ceil(math.sqrt(samplesPerCell))

        err_Linf = 0
        err_L2 = 0
        errs = []

        for faceIdx in range(0, self.mesh.num_elements):
            err_L2_local = 0
            faceVertices = self.mesh.fv_indices[faceIdx]
            fVCoords = []

            for localIdx in range(0, 3):
                fVCoords.append(self.mesh.points[faceVertices[localIdx]])
            
            X = lambda idx: fVCoords[idx-1][0]
            Y = lambda idx: fVCoords[idx-1][1]

            # sample uniformly on barycentric coordinate
            samplePoints = []
            bStep = 1.0 / (samplePerEdge - 1)
            for sampleIdxX in range(0, samplePerEdge):
                for sampleIdxY in range(0, samplePerEdge - sampleIdxX):
                    xi1 = bStep * sampleIdxX
                    xi2 = bStep * sampleIdxY
                    xi3 = 1 - xi1 - xi2
                    assert(xi3 <= 1)
                    samplePoints.append((xi1, xi2, xi3))

            for (xi1, xi2, xi3) in samplePoints:
                x = X(1)*xi1 + X(2)*xi2 + X(3)*xi3
                y = Y(1)*xi1 + Y(2)*xi2 + Y(3)*xi3
                ref_val = u_ref(x, y)
                my_val = self.getSolvedValueBarycentric(faceIdx, xi1, xi2, xi3)                    
        
                errs.append((x, y, ref_val - my_val))
                err_Linf = max(abs(ref_val - my_val), err_Linf)
                err_L2_local += abs(ref_val - my_val)
            
            detA = X(1)*Y(2) - X(1)*Y(3) - X(2)*Y(1) + X(2)*Y(3) + X(3)*Y(1) - X(3)*Y(2)
            area = 0.5 * np.abs(detA)

            err_L2 += area * err_L2_local

        return (errs, err_Linf, err_L2)


