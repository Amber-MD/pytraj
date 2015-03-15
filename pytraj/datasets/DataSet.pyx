# distutils: language = c++
from cpython.array cimport array as pyarray
from pytraj.cpptraj_dict import DataTypeDict, get_key

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

    @property
    def size(self):
        return self.baseptr0.Size()

    def set_width(self,int widthIn):
        self.baseptr0.SetWidth(widthIn)

    def set_precision(self, int widthIn , int precisionIno):
        self.baseptr0.SetPrecision(widthIn, precisionIno)

    #def setup_set(self, string nameIn, int idxIn, string aspectIn):
    #    return self.baseptr0.SetupSet(nameIn, idxIn, aspectIn)

    def set_legend(self, string lIn):
        self.baseptr0.SetLegend(lIn)

    def set_scalar(self, scalarMode mIn):
        self.baseptr0.SetScalar(mIn)

    def set_dim(self,DimIdxType i, Dimension d):
        self.baseptr0.SetDim(i, d.thisptr[0])

    def set_scalar(self,scalarMode mIn, scalarType mT):
        self.baseptr0.SetScalar(mIn, mT)

    def set_data_set_format(self, bint leftAlignIn):
        return self.baseptr0.SetDataSetFormat(leftAlignIn)

    def scalar_descr(self):
        self.baseptr0.ScalarDescription()

    def is_empty(self):
        return self.baseptr0.Empty()

    @property
    def legend(self):
        return self.baseptr0.Legend()
    
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
        return self.baseptr0.Aspect()

    @property
    def column_width(self):
        return self.baseptr0.ColumnWidth()
    
    @property
    def dtype(self):
        """Using `dtype` keyword since this is commond term"""
        return get_key(self.baseptr0.Type(), DataTypeDict).lower()

    @property
    def scalar_mode(self):
        return self.baseptr0.ScalarMode()

    @property
    def scalar_type(self):
        return self.baseptr0.ScalarType()

    @property 
    def ndim(self):
        return self.baseptr0.Ndim()

    def dim(self,int i):
        # TODO: what does this do?
        cdef Dimension dim = Dimension()
        dim.thisptr[0] = self.baseptr0.Dim(i)

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
        """return 1D python array of `self`"""
        cdef pyarray arr = pyarray('d', [])
        cdef int i

        for i in range(self.baseptr0.Size()):
            arr.append(self[i])
        return arr
