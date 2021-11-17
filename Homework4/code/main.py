#!/usr/bin/env python3
import logging, os
import meshUtil

def prepare_mesh(meshPrefix):
    cwd = os.getcwd()
    meshUtil.Mesh2D.from_easymesh(
        os.path.join(cwd, f"./meshes/{meshPrefix}.n"),
        os.path.join(cwd, f"./meshes/{meshPrefix}.e"),
        os.path.join(cwd, f"./meshes/{meshPrefix}.s")
    )

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)

    prepare_mesh("gd0")