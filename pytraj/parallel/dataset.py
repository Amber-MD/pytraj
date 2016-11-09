from __future__ import absolute_import
import numpy as np
from pytraj.core.c_options import info as compiled_info
from pytraj import matrix
from pytraj import mean_structure
from pytraj import volmap
from pytraj import Frame
from pytraj import ired_vector_and_matrix
from pytraj import rotation_matrix
from pytraj import NH_order_parameters
from pytraj.utils.tools import concat_dict
from pytraj.externals.six import string_types
from pytraj.utils.get_common_objects import get_reference

class PmapDataset(object):
    '''Dataset handlder for both pmap and pmap_mpi

    Parameters
    ----------
    data_collection : List[Tuple[rank, OrdereDict[key, array], n_frames]] or some numpy
    array based on return type
    '''

    def __init__(self, data_collection):
        pass

    def process(self):
        if func in [matrix.dist, matrix.idea, volmap]:
            mat = np.sum((val[1] * val[2] for val in data)) / traj.n_frames
            return mat
        elif func in [ired_vector_and_matrix, ]:
            # data is a list of (rank, (vectors, matrix), n_frames)
            mat = np.sum((val[1][1] * val[2] for val in data)) / traj.n_frames
            vecs = np.column_stack(val[1][0] for val in data)
            return (vecs, mat)
        elif func in [rotation_matrix, ]:
            if 'with_rmsd' in kwd.keys() and kwd['with_rmsd']:
                # data is a list of (rank, (mat, rmsd), n_frames)
                mat = np.row_stack(val[1][0] for val in data)
                rmsd_ = np.hstack(val[1][1] for val in data)
                return OrderedDict(out=(mat, rmsd_))
            else:
                mat = np.row_stack(val[1] for val in data)
                return OrderedDict(mat=mat)
        elif func == mean_structure:
            xyz = np.sum((x[2] * x[1].xyz for x in data)) / traj.n_frames
            frame = Frame(xyz.shape[0])
            frame.xyz[:] = xyz
            return frame
        elif 'hbond' in func.__name__:
            return concat_hbond(data)
        else:
            return concat_dict((x[1] for x in data))

