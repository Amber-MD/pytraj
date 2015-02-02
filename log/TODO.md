TODO
----

* organize codes to subfolder
* compat with official cpptraj in Ambertools14
* compat python 2.7 and 3.4
* write log: what is new in each version?
* add more examples
* load pdb file without topology: DONE
* write manual: DOING
* database for pycpptraj
* Sometimes get every large ts[-1:-9:-1][0].n_atoms (example: 30401312 atoms vs Tc5b = 304 atoms): DONE
    (ts is FrameArray or Trajin_Single instance)
    * create tmpfarray for FrameArray and Trajin_Single classes
* Add exception
* Rename "./examples" folder to "tests" and make REAL example script: DONE
* make action, analysis dictionary
* write script to mine the enum in cpptraj code and convert to dict : DONE
* Write automated script to convert *.h (cpptraj) to *.pxd (pycpptraj) files : DONE
* rename modules to lower case : should we?
* adding more doc
* to know why getting "ImportError: No module named trajs.TrajectoryIO" in `travis-ci`: DONE
    * wrong setup.py (this one did not inlcude pycpptraj.trajs as a package)
* Remove unused *pxd and *pyx files : DONE
