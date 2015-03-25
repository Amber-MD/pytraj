=======
pytraj change log
======

Lastest change
=============

Features added (from March 2015 - )
----------------------------------
* support Python 2.6
* ported to Ambertools 15
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
