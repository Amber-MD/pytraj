=======
pytraj change log
======

Lastest change
=============

Features and bugfixs added (from May 2015 - )
-----------------------------------
* speed up compiling time (~10 times) without using multiprocess
* move all cpptraj's Analysis classes to a single file (CpptrajAnalyses)
* move all cpptraj's Action classes to a single file (CpptrajActions)
* add show_version()
* change name for several DataSet classes (DataSet_integer to DatasetInteger)
* change islice to iterator_slice
* rmsfit_to --> rmsfit
* add load_cpptraj_datafile to read cpptraj's datafile
* use TrajectoryCpptraj as backend for TrajectoryIterator (able to load multiple files)
* add TrajectoryCpptraj
* add seperate command for each dihedral calculation (phi, psi,...)
* top.select("@CA") returns an AtomMask object
* add memoryview for DataSet_GridFlt
* start working on pytraj-dev 0.1.2.dev5
* release pytraj-dev 0.1.2.dev4 (May 26, 2015) (binstar, github)
* add pretty print for DataSetList when having `pandas` 
* enhance speed for frame_iter and chunk_iter
* add 'merge_traj'
* add `apply` method to change traj's internal data
* add basic math for Trajectory 
* add XYZ class as a simple wrapper for numpy ndarray
  traj.xyz is no longer ndarray, just a read-only wrapper
* add TrajectoryMDAnalysisIterator
* add read_data for DataFile and DataSetList
* add frame `swap` for Trajectory
* add atom coords assignment traj['@CA'] = xyz
* add `stack` for DataSetList to join 2 datasets
* rename `set_frame_m` to `set_frame_mass`
* update `to_amber_mask` to convert integer array-like to string mask
* add autoimage and rmsfit_to for chunk_iter and frame_iter
* tune speed for chunk_iter and frame_iter for TrajectoryIterator
* add `plot` method for DataSet
* add `write_to_cpptraj_format`
* add __array__ for Frame 
* introduce array interface to numpy for cpptraj's DataSet 
* add median, std, ... for DataSetList
* update run_all_and_find_fails.py for capturing segmetation fault
* add Trajectory joining (traj += traj)
* optimize slicing for `TrajectoryIterator` 
* optimize strip_atoms (1000 times speed up)
* add `copy` method for DataSet
* add fancy indexing for `DataSetList`
* add `_fast_slice` for Trajectory (1000 times faster)
* add regular expression for `groupby`
* add `to_pickle`, `to_json`, `read_pickle`, `read_json`
* add `_guess_filetype`
* change `fit_to` to `rmsfit_to`
* enhance smart _frame_iter_master
* create alias `rmsd` for `calc_rmsd` in `common_actions` module
* update loading hd5f file without going through `mdtraj`
* correctly handle `box` in `api.Trajectory`
* add `TrajectoryREMDIterator` class to handal REMD
* change `io.load_remd_iterator` to `io.iterload_remd`
* add `io._iterload_from_filelist` and `io._load_from_filelist`
* add `count` for DataSet_integer
* API CHANGE (05-06-2015): 
  * change `FrameArray` to `Trajectory`
  * change `TrajReadOnly` to `TrajectoryIterator` (issue #262)
    (https://github.com/pytraj/pytraj/issues/262)
  * `io.load` will return `Trajectory` instead of `TrajectoryIterator`
  * `io.iterload` will return `TrajectoryIterator`
* add "legend" as property of DataSet
* remove "get_box" from Frame. use "box"
* support `autoimage` when loading mdtraj object
* bugfix for `io.loadpdb_rcsb`
* add 'calc_rdf'
* trick to add reference to Action
  * Example: calc_abcxyz([[ref,], traj], command), where `ref` is a Frame object, `abcxyz` is a specific action.
* add option to load coords from ParmEd object
* add `align_principal_axis`
* fix iteration in api.Trajectory: 2 orders of magniture faster
* add "as_ndmatrix" for Matrix_3x3
* start working on pytraj-dev 0.1.2.dev4

Features and bugfixs added (from April 2015 - )
-----------------------------------
* release pytraj-dev 0.1.2.dev3 (May 1st, 2015) (binstar, github)
* add set_nobox for FrameArray
* add `dtype` in `calculate` method to return `list`, `pyarray`, `dataset` or `ndarray`
* support ParmEd with bonds, dih, angles loading
* add support loading more infor (bonds, box, angles, dihedrals) from MDAnalysis 
  and mdtraj
* add `bonds`, `angles` and `dihedrals` iterators for Topology
* support MPI (with mpi4py)
* add `calc_rmsd` 
* rename "do_rotation" to "rotate", "do_translation" to "translate"
  "do_autoimage" to "autoimage"
* add "box" to Topology
* reorganize folders (04-25-2015)
    * moved Atom, Molecule, Residue, Box to `pytraj.core`
    * moved Trajin_Single to `pytraj.trajs`
    * moved Matrix_3x3, Vec3, Grid ... to `pytraj.math` 
    * moved all unused files away from main folder
    * dont support Energy routine in cpptraj, use pysander instead
    * move all optional package importing to pytraj.externals. 
    * enhance `groupby` method: dslist.groupby("omega", mode='aspect')
    * added pytraj.api.Trajectory as a new Trajectory object with `xyz` is a numpy array
* add universal _get_top method
    _get_top(traj, top)
* add DataSet_Mesh
* better setup.py: automatically install `libcpptraj`
* change API of FrameArray and TrajReadOnly
   traj['CA'] return a stripped-atoms traj (not a xyz coords)
* add calc_temperature
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

Log from very first days of pytraj (from earliest to newest)
----------------------------------
* Oct - Nov 2014: Start writting pycpptraj (then was renamed to `pytraj` following Jason's suggestion)
* Nov-11-2014: Start writing gencode.py to automatically convert cpptraj header file to python 
 (\*.pxd and \*.pyx) files
* Nov-14-2014: Finish importing all (most?) cpptraj functions (most header files) to Cython. 
  Passed 229 tests (./src/tests/test_*.pyx)
* Dec-04-2014: Write setup.py and test installing
* Dec-11-2014: Know how to create wrapper for abstract and sub-class
* Jan-24-2015: Taking 3D numpy array just by traj[:, :, :] or traj[0:20, :, :] (read only)
* Jan-29-2015: working on pip: python setup.py sdist
    * good tutorial: http://peterdowns.com/posts/first-time-with-pypi.html
    * to test:
        * python setup.py register -r pypitest
        * python setup.py sdist upload -r pypitest
    * to upload:
        * register: python setup.py register
        * upload: python setup.py sdist upload -r pypi
* Jan-31-2015: working on porting code to Python3. Stuck
    * Error: "TypeError: expected bytes, str found". tried but don't know how yet.
        Try: filename = filename.encode("UTF-8") (work well for the communication between C++ string and Python3 string)
    * us 2to3 to make compat: 2to3 --output-dir=python3-version/mycode -W -n python2-version/mycode
