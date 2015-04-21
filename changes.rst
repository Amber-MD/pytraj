=======
pytraj change log
======

Lastest change
=============

Features and bugfixs added (from April 2015 - )
-----------------------------------
* support distance-based mask selection (start v0.1.2.dev3 version)
* add "groupby" method for DataSetList (end v0.1.2.dev2 version)
* directly access "calc_..." methods from FrameArray and TrajReadOnly (traj.calc_COM, ...)
* "tolist" in Traj should return 3D list (n_frame, n_atoms, 3)
* "tolist" of Frame should return 2D list (n_atoms, 3)
* add "tolist" for DataSet_Vector
* add "keys()" to DataSetList class. Now we can use this as a list and a dict
* support Grid calculation (calc_volmap, ...)
* enhance speed of AtomMask selection (about 5 times faster)
* add `set_error_silent` to  turn-off (and on) cpptraj' warning
* add attribute `indices` for AtomMask object (shorter than `selected_indices` method's name
* support loading `MDAnalysis` objects (Universe)
* add method `load_pseudo_parm` to load ParmEd and mdtraj Topology objects
* fully compat with cpptraj, user can get CpptrajState from trajin.in file
  (`pytraj.io.load_cpptraj_file`)
* add `to_ndarray` for DataSet, Frame, Vec3,... objects
* introduce memoryview for DataSet
* fix memory for `calculate` method
* add simple wrapper to pysander to be used with pandas' DataFrame
* add calc_multidihedral for quickly scanning dihedral
* improve code for AtomMask: `add_selected_indices` could accept memview
* add `_frame_iter_master` to iterate anything that results Frame objects
* add frame_iter with atom indices (before: only support atom mask): traj(mask=index_array)

Features added (from March 2015 - April 2015)
----------------------------------
* can combine with `pysander` to get energy
* update DataSet class inheritance (following `cpptraj` changes)
* support Analysis classes in cpptraj (including **clustering**)
* introduce `dataframe` to interact with `pandas` (`to_dataframe`)
* support OPENMP
* have option to make faster building (less than 2 minutes) by using multiprocess
* `load` method of FrameArray is much more robust (it can load raw numpy, list, tuple of numbers, any iterable objects that resulting Frame ... and much more)
* support loading `mdtraj` objects (Trajectory and Topology). Support loading array with ndim=1, 2, 3
* support more than 100 kinds of calculations and analyses in `cpptraj`
* many bugfixs to avoid segmentation fault 
* support Python 2.6
* support simple plots (plot_matrix, ...)
* add get_* methods for DataSetList (dslist.get_legends())
* update more pytraj notebooks
* fix DataSet_MatrixDbl segmentation fault, now we can get full matrix (get_full_matrix)
* add traj iter with mask (traj(0, 10, 2, '@CA')
* add method 'fit_to' for traj to get alignment
* clean up, don't keep methods pytraj does not need  (from cpptraj)

Features added (from Jan 2015 - Feb 2015)
----------------------------------------
* partially supported `pip pytraj install` (or `pip install pytraj -d .`)
* support Python version: 2.7, 3.3, 3.4 (pytraj.v0.1.beta)
* support most cpptraj action classes (>70 Actions)
* shorten codes 
* quickly get FrameArray with given mask: traj['@CA :frame']
* faster iterator for reading trajectory
* support reading hd5f file from mdtraj
* better setup.py

Bugs fixed
----------
* (add here)
* fix several segmenation faults

Other stuff
----------
* removing unnecessary cpptraj methods
