from types import DynamicClassAttribute
import numpy as np

import logging
logger = logging.getLogger(__name__)

class Mesh2D:
    def __init__(self, points: np.ndarray, boundary_marks: np.ndarray, fv_indices: np.ndarray):
        self.points = points
        self.num_points = len(self.points)
        # Boundar markers of points
        self.boundary_marks = boundary_marks

        self.num_elements = len(fv_indices)
        self.fv_indices = fv_indices

        # Calculate free nodes
        self.nodes = []
        for idx in range(0, self.num_points):
            if self.boundary_marks[idx] == False:
                self.nodes.append(idx)

    def get_adjacent_elements(self, pointIdx: int):
        adjElems = []
        for elemIdx in range(0, self.num_elements):
            if pointIdx == self.fv_indices[elemIdx][0] \
                or pointIdx == self.fv_indices[elemIdx][1] \
                or pointIdx == self.fv_indices[elemIdx][2]:
                adjElems.append(elemIdx)
        
        return adjElems

    @staticmethod
    def from_easymesh(nodefilePath, elementFilePath, sideFilePath):
        """
        https://web.mit.edu/easymesh_v1.4/www/output.html
        """
        with open(nodefilePath, "r") as nf:
            numNodes = int(nf.readline())
            points = np.ndarray((numNodes, 2), dtype=np.float64)
            boundary_marks = np.ndarray((numNodes,), dtype=np.bool8)
            for idx in range(0, numNodes):
                lineData = nf.readline().split()
                pointIdx = int(lineData[0])
                pointX = float(lineData[1])
                pointY = float(lineData[2])
                isBoundary = bool(int(lineData[3]))

                assert(pointIdx == idx)
                points[idx][0] = pointX
                points[idx][1] = pointY
                boundary_marks[idx] = isBoundary
                logger.debug(f"- ({pointX}, {pointY}) boundary={isBoundary}")
        
        logger.info(f"numNodes: {numNodes}")

        with open(elementFilePath, "r") as ef:
            numElements = int(ef.readline())
            fv_indices = np.ndarray((numElements, 3), dtype=np.int32)
            for idx in range(0, numElements):
                lineData = ef.readline().split()
                elemIdx = int(lineData[0])
                elemVi, elemVj, elemVk = int(lineData[1]), int(lineData[2]), int(lineData[3])

                assert(elemIdx == idx)
                logger.debug(f"- elemIdx={elemIdx}, {elemVi} {elemVj} {elemVk}")
            logger.info(f"numElements: {numElements}")

        return Mesh2D(points, boundary_marks, fv_indices)