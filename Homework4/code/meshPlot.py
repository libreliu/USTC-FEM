import pyvista as pv
import numpy as np

def gen_vis_polydata(points, fv_indices):
    """ Plot (triangulated) mesh via pyvista """
    # print(points)

    orig_indices = fv_indices
    # print(orig_indices)

    trans_indices = np.zeros((orig_indices.shape[0], orig_indices.shape[1] + 1), dtype=np.int32)
    trans_indices[:,1:4] = orig_indices
    trans_indices[:,0] = 3

    trans_indices = np.hstack(trans_indices)
    # print(trans_indices)

    surf = pv.PolyData(points, trans_indices)
    return surf

def offscreen_combine_plot(filename, *args):
    assert(len(args) >= 1)
    p = pv.Plotter(off_screen=True, window_size=[1920, 1080])
    p.set_position(np.array([4.0, 6.4, 5.6]))
    for idx, arg in enumerate(args):
        p.add_mesh(arg[0], **arg[1])
    
    p.screenshot(filename)

def combine_plot(*args):
    assert(len(args) >= 1)
    p = pv.Plotter()
    for idx, arg in enumerate(args):
        p.add_mesh(arg[0], **arg[1])
    p.show()
