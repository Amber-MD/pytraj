from __future__ import absolute_import
from .DataSetList import DataSetList
from .cast_dataset import cast_dataset
from .DataSet import DataSet
from .DataSet_1D import DataSet_1D
from .DataSet_2D import DataSet_2D
from .DatasetDouble import DatasetDouble
from .DatasetFloat import DatasetFloat
from .DatasetInteger import DatasetInteger
from .DatasetString  import DatasetString
from .DataSet_Mesh import DataSet_Mesh
from .DatasetMatrixDouble import DatasetMatrixDouble
from .DataSet_MatrixFlt import DataSet_MatrixFlt
from .DatasetGridFloat import DatasetGridFloat
from .DatasetVector import DatasetVector
from .DataSet_Coords import DataSet_Coords
from .DataSet_Coords_REF import DataSet_Coords_REF
from .DataSet_Coords_CRD import DataSet_Coords_CRD
from .DataSet_Coords_TRJ import DataSet_Coords_TRJ


def from_pickle(filename):
    from .. DataSetList import DataSetList
    dslist = DataSetList()
    dslist.from_pickle(filename)
    return dslist

def from_json(filename):
    from .. DataSetList import DataSetList
    dslist = DataSetList()
    dslist.from_json(filename)
    return dslist
