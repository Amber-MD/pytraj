Welcom to pycpptraj!

Goal: Python API for cpptraj library, a data analysis package for biomolecular simulation
-----------------------------------------
- Why using pycpptraj:
    * It's fast
        * it's a wrapper of cpptraj (was written in C++)
        * it has interface with numpy. Data calculation are performed without copying to numpy array
    * It has clean syntax
        * Python/Cython style with fancy indexing 
    * It has been extensively testes
    * It's flexible: 
        * you can write extension modules in either high (Python) or low (C/C++ or Cython) level
