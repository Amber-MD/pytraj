from __future__ import absolute_import
from .externals.six import string_types
from .Topology import Topology
from .ArgList import ArgList

def _get_top(traj, top):
    if isinstance(top, string_types):
        _top = Topology(top)
    elif top is None: 
        if hasattr(traj, 'top'):
           _top = traj.top 
        else:
            # list, tuple of traj objects 
            try:
                for tmp in traj:
                    if hasattr(tmp, 'top'):
                        _top = tmp.top 
                        break
            except:
                #print("Topology is None")
                _top = None
    else:
        _top = top
    return _top

def _get_arglist(arg):
    if isinstance(arg, ArgList):
        return arg
    else:
        return ArgList(arg)

def _get_data_from_dtype(d0, dtype='dataset'):
   dtype = dtype.lower()
   if dtype == 'dataset':
       return d0
   elif dtype == 'list':
       return d0.tolist()
   elif dtype == 'ndarray':
       return d0.to_ndarray()
   elif dtype == 'pyarray':
       return d0.to_pyarray()
   elif dtype == 'dict':
       return d0.to_dict()
   else:
       raise NotImplenmentedError()
