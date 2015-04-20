TODO
----

* better way to manage memory for DataSetList
    * calling traj.search_hbonds().to_dict() a thousand times,  the results are appended to 
      single DataSetList object
* add method to check class hiercharchy (in case Dan updates his codes)
* build conda `libcpptraj` with netcdf: DONE
* FIXME: `rmsd` method in Frame class will update frame object coordinate (I do not want this): DONE
* rewrite cpptraj classes for pytraj (more flexibility): ActionList, CpptrajState, TrajinList, AnalysisList: LOW
* `frame_iter` is very slow, enhance this : DONE
* organize codes to subfolder 
* compat with official cpptraj in Ambertools14: DONE
* compat python 2.7 and 3.4: DONE
* write log: what is new in each version?: DOING
* add more examples: DOING
* load pdb file without topology: DONE
* write manual: DOING
* database for pycpptraj: REALLY NEED?
* Sometimes get every large ts[-1:-9:-1][0].n_atoms (example: 30401312 atoms vs Tc5b = 304 atoms): DONE
* Add exception
* Rename "./examples" folder to "tests" and make REAL example script: DONE
* make action, analysis dictionary: DONE
* write script to mine the enum in cpptraj code and convert to dict : DONE
* Write automated script to convert *.h (cpptraj) to *.pxd (pycpptraj) files : DONE
* rename modules to lower case : should we?
* adding more doc
* to know why getting "ImportError: No module named trajs.TrajectoryIO" in `travis-ci`: DONE
    * wrong setup.py (this one did not inlcude pycpptraj.trajs as a package)
* Remove unused *pxd and *pyx files : DONE
* protein view using `ivew` (like mdtraj package)
