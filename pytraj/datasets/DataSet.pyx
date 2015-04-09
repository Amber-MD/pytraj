# distutils: language = c++
from cpython.array cimport array as pyarray
from pytraj.cpptraj_dict import DataTypeDict, scalarDict, scalarModeDict, get_key
from pytraj.decorators import makesureABC

cdef class DataSet:
    """
    Original doc from cpptraj:
    =========================
    Class: DataSet
        Base class that all DataSet types will inherit.
        DataSets are given certain attributes to make DataSet selection easier; 
        these are name, index, and aspect. Name is typically associated with the
        action that creates the DataSet, e.g. RMSD or distance. Index is used
        when and action outputs subsets of data, e.g. with RMSD it is possible to 
        output per-residue RMSD, where the DataSet index corresponds to the residue
        number. Aspect is used to further subdivide output data type; e.g. with 
        nucleic acid analysis each base pair (divided by index) has shear,
        stagger etc calculated.
    """

    def __cinit__(self):
        pass
        # don't create instance for this abstract class
        #self.baseptr0 = new _DataSet()

    def __dealloc__(self):
        pass
        # let sub-class do this job
        #if self.baseptr0 != NULL:
        #    del self.baseptr0

    def __iter__(self):
        raise NotImplementedError("Must over-write DataSet data attr")

    def __getitem__(self, idx):
        raise NotImplementedError("Must over-write DataSet data attr")

    def __setitem__(self, idx, value):
        raise NotImplementedError("Must over-write DataSet data attr")

    @property
    def size(self):
        return self.baseptr0.Size()

    def set_width(self, int width):
        self.baseptr0.SetWidth(width)

    def set_precision(self, int width , int precision):
        self.baseptr0.SetPrecision(width, precision)

    def set_legend(self, lengend):
        self.baseptr0.SetLegend(lengend.encode())

    def set_dim(self, DimIdxType i, Dimension d):
        self.baseptr0.SetDim(i, d.thisptr[0])

    def set_scalar(self,scalar_mode, scalar_type=None):
        scalar_mode = scalar_mode.upper()
        if scalar_type is None:
            self.baseptr0.SetScalar(scalarModeDict[scalar_mode])
        else:
            scalar_type = scalar_type.upper()
            self.baseptr0.SetScalar(scalarModeDict[scalar_mode], scalarDict[scalar_type])

    def set_format(self, bint leftAlignIn):
        return self.baseptr0.SetDataSetFormat(leftAlignIn)

    def scalar_descr(self):
        from pytraj._utils import set_worl_silent
        set_worl_silent(False)
        self.baseptr0.ScalarDescription()
        set_worl_silent(True)

    def is_empty(self):
        return self.baseptr0.Empty()

    @property
    def legend(self):
        legend = self.baseptr0.Legend()
        return legend.decode()
    
    @property
    def name(self):
        cdef string t
        t = self.baseptr0.Name()
        return t.decode()

    @property
    def idx(self):
        return self.baseptr0.Idx()
    
    @property
    def aspect(self):
        aspect = self.baseptr0.Aspect()
        return aspect.decode()

    @property
    def column_width(self):
        return self.baseptr0.ColumnWidth()
    
    @property
    def dtype(self):
        """Using `dtype` keyword since this is commond term"""
        return get_key(self.baseptr0.Type(), DataTypeDict).lower()

    @property
    def scalar_mode(self):
        return get_key(self.baseptr0.ScalarMode(), scalarModeDict).lower()

    @property
    def scalar_type(self):
        return get_key(self.baseptr0.ScalarType(), scalarDict).lower()

    @property 
    def ndim(self):
        return self.baseptr0.Ndim()

    def __richcmp__(DataSet self, DataSet other, opt):
        if opt == 0:
            # operator "<"
            return self.baseptr0[0] < other.baseptr0[0]
        else:
            raise NotImplemented()

    @property
    def data_format(self):
        return self.baseptr0.DataFormat()

    @property
    def data(self):
        """return 1D python array of `self`
        ABC method, must override
        """
        raise NotImplementedError("Must over-write DataSet data attr")

    def tolist(self):
        return list(self.data)

    def to_ndarray(self):
        """return ndarray view of self.data"""
        from pytraj.utils import _import_numpy
        has_np, np = _import_numpy()
        if has_np:
            return np.asarray(self.data)
        else:
            raise ImportError("require numpy")
