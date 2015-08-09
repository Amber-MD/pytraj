# distutils: language = c++
from __future__ import division
from cpython.array cimport array as pyarray
from ..cpptraj_dict import DataTypeDict, scalarDict, scalarModeDict, get_key
from ..decorators import makesureABC, require_having
from ..core.DataFileList import DataFileList
from ..core.DataFile import DataFile
from ..utils import _import_numpy

_, np = _import_numpy()

cdef class DataSet:
    __cpptraj_dataset__ = None
    """
    Original doc from cpptraj
    -------------------------
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

    pytraj doc
    ----------
        Add some other methods: `tolist`, `to_ndarray`, `plot` and attribute `data`

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

    def __str__(self):
        cname = self.class_name
        dname = self.name
        dformat = self.format
        size = self.size
        legend = self.legend
        aspect = self.aspect
        dtype = self.dtype

        msg0 = """<pytraj.datasets.{0}: size={1}, key={2}> """.format(cname, size, legend)
        return msg0

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        raise NotImplementedError("Must over-write DataSet data attr")

    def __getitem__(self, idx):
        raise NotImplementedError("Must over-write DataSet data attr")

    def __setitem__(self, idx, value):
        raise NotImplementedError("Must over-write DataSet data attr")


    def __array__(self):
        """
        Aim: directly use numpy to perform analysis without casting to ndararay again

        Examples
        -------
            d = DataSet_integer()
            d.resize(200)
            d.data[:] = np.arange(200)
            np.mean(d)
        """
        from pytraj.utils import _import_numpy
        _, np = _import_numpy()
        try:
            return np.asarray(self.data)
        except:
            raise NotImplementedError("don't know how to cast to ndarray")

    def __add__(self, value):
        dnew = self.copy()
        dnew.__iadd__(value)
        return dnew

    def __iadd__(self, value):
        if hasattr(value, '_npdata'):
            self._npdata += value._npdata
        else:
            self._npdata += value
        return self

    def __sub__(self, value):
        dnew = self.copy()
        dnew.__isub__(value)
        return dnew

    def __isub__(self, value):
        if hasattr(value, '_npdata'):
            self._npdata -= value._npdata
        else:
            self._npdata -= value
        return self

    def __mul__(self, value):
        dnew = self.copy()
        dnew.__imul__(value)
        return dnew

    def __imul__(self, value):
        if hasattr(value, '_npdata'):
            self._npdata *= value._npdata
        else:
            self._npdata *= value
        return self

    def __div__(self, value):
        dnew = self.copy()
        dnew.__imul__(value)
        return dnew

    def __idiv__(self, value):
        if hasattr(value, '_npdata'):
            self._npdata /= value._npdata
        else:
            self._npdata /= value
        return self

    def __itruediv__(self, value):
        if hasattr(value, '_npdata'):
            self._npdata /= value._npdata
        else:
            self._npdata /= value
        return self

    def copy(self):
        cdef int i
        cdef int size = self.size

        new_ds = self.__class__()
        new_ds.legend = self.legend
        new_ds.set_scalar(self.scalar_mode, self.scalar_type)
        new_ds.set_name_aspect_index_ensemble_num(self.name,
                                                  self.aspect,
                                                  self.idx,
                                                  0)
        try:
            try:
                new_ds.resize(self.size)
                new_ds.data[:] = self.data
            except:
                # try to make copy (Vector, ...)
                # slower
                for i in range(size):
                    new_ds.append(self[i])
            return new_ds
        except:
            raise TypeError("don't know how to copy %s" % self.class_name)

    @property
    def class_name(self):
        return self.__class__.__name__

    @property
    def size(self):
        return self.baseptr0.Size()

    def set_width(self, int width):
        self.baseptr0.SetWidth(width)

    def set_precision(self, int width , int precision):
        self.baseptr0.SetPrecision(width, precision)

    def set_legend(self, legend):
        self.baseptr0.SetLegend(legend.encode())

    def set_scalar(self,scalar_mode, scalar_type=None):
        scalar_mode = scalar_mode.upper()
        if scalar_type is None:
            self.baseptr0.SetScalar(scalarModeDict[scalar_mode])
        else:
            scalar_type = scalar_type.upper()
            self.baseptr0.SetScalar(scalarModeDict[scalar_mode], scalarDict[scalar_type])

    def set_format(self, bint leftAlignIn):
        return self.baseptr0.SetDataSetFormat(leftAlignIn)

    def set_name_aspect_index_ensemble_num(self, name, aspect, idx, num):
        self.baseptr0.SetupSet(name.encode(), idx, aspect.encode(), num)

    def scalar_descr(self):
        from pytraj import set_worl_silent
        set_worl_silent(False)
        self.baseptr0.ScalarDescription()
        set_worl_silent(True)

    def is_empty(self):
        return self.baseptr0.Empty()

    property legend:
        def __get__(self):
            legend = self.baseptr0.Legend()
            return legend.decode()
        def __set__(self, legend):
            self.baseptr0.SetLegend(legend.encode())

    property key:
        # retire self.legend?
        def __get__(self):
            legend = self.baseptr0.Legend()
            return legend.decode()
        def __set__(self, legend):
            self.baseptr0.SetLegend(legend.encode())
    
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
    def format(self):
        my_format = self.baseptr0.DataFormat().decode()
        return my_format.strip()

    @property
    def data(self):
        """return 1D python array of `self`
        ABC method, must override
        """
        raise NotImplementedError("Must over-write DataSet data attr")

    property _npdata:
        def __get__(self):
            """return memoryview as numpy array for self.data"""
            # NOTE: overwrite by using `raise NotImplementedError`
            # for some DataSet subclasses not returning a `view`
            _, np = _import_numpy()
            return np.asarray(self.data)
        def __set__(self, value):
            _, np = _import_numpy()
            arr = np.asarray(self.data)
            arr[:] = value

    def tolist(self):
        return list(self.data)

    def to_pyarray(self):
        type_dict = {'float' : 'f',
                     'double' : 'd',
                     'integer' : 'i',
                     'string' : 's',
                    }
        try:
            return pyarray(type_dict[self.dtype], self.data)
        except:
            msg = "not implemented for %s" % self.__class__.__name__
            raise NotImplementedError(msg)

    @property
    def values(self):
        return self.to_ndarray()

    def to_ndarray(self, copy=False):
        """return ndarray view of self.data"""
        from pytraj.utils import _import_numpy
        has_np, np = _import_numpy()
        if has_np:
            if not copy:
                return np.asarray(self.data)
            else:
                # make copy
                return np.array(self.data)
        else:
            raise ImportError("require numpy. Use `tolist` or `to_pyarray`"
                              " or `to_dict` to access data")

    def to_dict(self, use_numpy=False):
        if np and use_numpy:
            return {self.legend : self.values}
        if not np and use_numpy:
            raise ImportError("require numpy. Set `use_numpy=False`")
        return {self.legend : self.tolist()}

    def hist(self, plot=True, show=True, *args, **kwd):
        """
        Parameters
        ----------
        plot : bool, default True
            if True, use `matplotlib` to plot. 
            if False, return `2D numpy array`
        """
        if not plot:
            import numpy as np
            return np.histogram(self.values)
        else:
            try:
                from matplotlib import pyplot as plt
                ax = plt.hist(self.values, *args, **kwd)
                if show:
                    plt.show()
                return ax
            except ImportError:
                raise ImportError("require matplotlib")

    def split(self, n_chunks_or_array):
        """split `self.data` to n_chunks

        Notes : require numpy (same as `array_split`)
        """
        return np.array_split(self.to_ndarray(), n_chunks_or_array)

    def write_to_cpptraj_format(self, filename):
        dflist = DataFileList()
        d = dflist.add_datafile(filename)
        d.add_dataset(self)
        d.write_data()

    def save(self, filename, format='cpptraj'):
        '''TODO: pickle, json'''
        if format == 'cpptraj':
            self.write_to_cpptraj_format(filename)
        else:
            raise NotImplementedError("not yet, stay tuned")

    def plot(self, show=False, *args, **kwd):
        """return matplotlib object
        Notes
        ----
        Need to over-write this method for subclass if needed.
        """
        from pytraj.utils import _import
        _, plt = _import("matplotlib.pyplot")
        try:
            ax = plt.pyplot.plot(self.data, *args, **kwd)
            if show:
                plt.pyplot.show()
            return ax
        except ImportError:
            raise ImportError("require matplotlib")

    def chunk_average(self, n_chunk):
        import numpy as np
        return np.array(list(map(np.mean, self.split(n_chunk))))

    def std(self, *args, **kwd):
        import numpy as np
        return np.std(self.values, *args, **kwd)

    def sum(self, *args, **kwd):
        return np.sum(self.values, *args, **kwd)

    def topk(self, k):
        """pick top k max-values
        Returns
        -------
        a list with len = k

        # TODO : array?
        """
        return sorted(self.values, reverse=True)[:k]

    def head(self, k=20, restype='ndarray'):
        if restype == 'ndarray':
            return self.values[:k]
        elif restype == 'list':
            return self.tolist()[:k]

    def tail(self, k=20):
        return self.values[-k:]

    def is_(self, DataSet other):
        return self.baseptr0 == other.baseptr0

    def filter(self, func):
        """return a numpy array with all elements that satisfy `func`

        Example
        -------
        >>> d0 = traj.calc_radgyr(dtype='dataset')[0]
        >>> d0.filter(lambda x : 105. < x < 200.)
        """
        import numpy as np
        return np.array(list(filter(func, self.values)))
