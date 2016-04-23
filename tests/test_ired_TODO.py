#!/usr/bin/env python
from __future__ import print_function
import os
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import eq, aa_eq, cpptraj_test_dir
from pytraj.compat import zip
from pytraj.nmr import _ired

parm_dir = cpptraj_test_dir + '/Test_IRED/1IEE_A_prot.prmtop'
traj_dir = cpptraj_test_dir + '/Test_IRED/1IEE_A_test.mdcrd'

txt = '''
 parm {0}
 trajin {1}
 vector v0 @25 ired @26
 vector v1 @41 ired @42
 vector v2 @61 ired @62
 vector v3 @68 ired @69
 vector v4 @92 ired @93
 vector v5 @102 ired @103
 vector v6 @117 ired @118
 vector v7 @136 ired @137
 vector v8 @146 ired @147
 vector v9 @156 ired @157
 vector v10 @166 ired @167
 vector v11 @183 ired @184
 vector v12 @205 ired @206
 vector v13 @229 ired @230
 vector v14 @246 ired @247
 vector v15 @253 ired @254
 vector v16 @272 ired @273
 vector v17 @284 ired @285
 vector v18 @298 ired @299
 vector v19 @319 ired @320
 vector v20 @343 ired @344
 vector v21 @350 ired @351
 vector v22 @371 ired @372
 vector v23 @382 ired @383
 vector v24 @401 ired @402
 vector v25 @408 ired @409
 vector v26 @422 ired @423
 vector v27 @446 ired @447
 vector v28 @462 ired @463
 vector v29 @472 ired @473
 vector v30 @482 ired @483
 vector v31 @492 ired @493
 vector v32 @514 ired @515
 vector v33 @534 ired @535
 vector v34 @549 ired @550
 vector v35 @560 ired @561
 vector v36 @574 ired @575
 vector v37 @594 ired @595
 vector v38 @608 ired @609
 vector v39 @622 ired @623
 vector v40 @639 ired @640
 vector v41 @649 ired @650
 vector v42 @663 ired @664
 vector v43 @677 ired @678
 vector v44 @701 ired @702
 vector v45 @715 ired @716
 vector v46 @729 ired @730
 vector v47 @741 ired @742
 vector v48 @748 ired @749
 vector v49 @759 ired @760
 vector v50 @773 ired @774
 vector v51 @785 ired @786
 vector v52 @806 ired @807
 vector v53 @813 ired @814
 vector v54 @832 ired @833
 vector v55 @851 ired @852
 vector v56 @868 ired @869
 vector v57 @887 ired @888
 vector v58 @901 ired @902
 vector v59 @912 ired @913
 vector v60 @936 ired @937
 vector v61 @960 ired @961
 vector v62 @984 ired @985
 vector v63 @994 ired @995
 vector v64 @1008 ired @1009
 vector v65 @1020 ired @1021
 vector v66 @1027 ired @1028
 vector v67 @1051 ired @1052
 vector v68 @1079 ired @1080
 vector v69 @1086 ired @1087
 vector v70 @1097 ired @1098
 vector v71 @1121 ired @1122
 vector v72 @1135 ired @1136
 vector v73 @1154 ired @1155
 vector v74 @1164 ired @1165
 vector v75 @1178 ired @1179
 vector v76 @1211 ired @1212
 vector v77 @1221 ired @1222
 vector v78 @1232 ired @1233
 vector v79 @1242 ired @1243
 vector v80 @1261 ired @1262
 vector v81 @1280 ired @1281
 vector v82 @1291 ired @1292
 vector v83 @1302 ired @1303
 vector v84 @1314 ired @1315
 vector v85 @1333 ired @1334
 vector v86 @1347 ired @1348
 vector v87 @1357 ired @1358
 vector v88 @1368 ired @1369
 vector v89 @1384 ired @1385
 vector v90 @1398 ired @1399
 vector v91 @1408 ired @1409
 vector v92 @1418 ired @1419
 vector v93 @1440 ired @1441
 vector v94 @1462 ired @1463
 vector v95 @1481 ired @1482
 vector v96 @1497 ired @1498
 vector v97 @1508 ired @1509
 vector v98 @1520 ired @1521
 vector v99 @1527 ired @1528
 vector v100 @1541 ired @1542
 vector v101 @1548 ired @1549
 vector v102 @1565 ired @1566
 vector v103 @1579 ired @1580
 vector v104 @1589 ired @1590
 vector v105 @1613 ired @1614
 vector v106 @1629 ired @1630
 vector v107 @1639 ired @1640
 vector v108 @1663 ired @1664
 vector v109 @1687 ired @1688
 vector v110 @1701 ired @1702
 vector v111 @1725 ired @1726
 vector v112 @1735 ired @1736
 vector v113 @1757 ired @1758
 vector v114 @1764 ired @1765
 vector v115 @1778 ired @1779
 vector v116 @1790 ired @1791
 vector v117 @1806 ired @1807
 vector v118 @1823 ired @1824
 vector v119 @1833 ired @1834
 vector v120 @1857 ired @1858
 vector v121 @1876 ired @1877
 vector v122 @1900 ired @1901
 vector v123 @1907 ired @1908
 vector v124 @1917 ired @1918
 vector v125 @1941 ired @1942
 matrix ired name matired order 2
 analyze matrix matired vecs 126
 createcrd CRD1
'''.format(parm_dir, traj_dir)


class TestIred(unittest.TestCase):

    def test_simple_for_coverage(self):
        '''
        '''
        traj = pt.iterload(traj_dir, parm_dir)
        h_indices = pt.select_atoms('@H', traj.top)
        n_indices = pt.select_atoms('@H', traj.top) - 1
        nh_indices = list(zip(n_indices, h_indices))
        vecs_and_mat = pt.ired_vector_and_matrix(traj, nh_indices, dtype='tuple')
        vecs_and_mat = pt.ired_vector_and_matrix(traj, nh_indices, dtype='tuple')
        state_vecs = vecs_and_mat[0]
        mat_ired = vecs_and_mat[1]

        # get eigenvalues and eigenvectors
        modes = pt.matrix.diagonalize(mat_ired, n_vecs=len(state_vecs))
        evals, evecs = modes

        data_0 = _ired(state_vecs,
                       modes=(evals, evecs),
                       NHbond=True,
                       tcorr=10000,
                       tstep=1.)

        data_1 = _ired(state_vecs,
                       modes=modes,
                       NHbond=True,
                       tcorr=10000,
                       tstep=1)

        for d0, d1 in zip(data_0, data_1):
            if d0.dtype not in ['modes', ]:
                aa_eq(d0.values, d1.values)
            else:
                # modes
                # values: tuple
                aa_eq(d0.values[0], d1.values[0])
                aa_eq(d0.values[1], d1.values[1])

        # try different dtype
        out_try_new_dtype = pt.ired_vector_and_matrix(traj, nh_indices, dtype='cpptraj_dataset')

    # TODO: how can I get order paramters?
    def test_ired_need_lapack_cpptraj(self):
        state = pt.load_cpptraj_state(txt)
        state.run()
        xyz = state.data['CRD1'].xyz
        top = state.data['CRD1'].top
        traj = pt.Trajectory(xyz=xyz, top=top)
        state_vecs = state.data[1:-3].values

        h_indices = pt.select_atoms('@H', traj.top)
        n_indices = pt.select_atoms('@H', traj.top) - 1
        nh_indices = list(zip(n_indices, h_indices))
        mat_ired = pt.ired_vector_and_matrix(traj,
                                             mask=nh_indices,
                                             order=2)[-1]
        mat_ired /= mat_ired[0, 0]

        # matired: make sure to reproduce cpptraj output
        aa_eq(mat_ired, state.data['matired'].values)

        # get modes
        modes = state.data[-2]
        cpp_eigenvalues = modes.eigenvalues
        cpp_eigenvectors = modes.eigenvectors
        evals, evecs = np.linalg.eigh(mat_ired)

        # need to sort a bit
        evals = evals[::-1]
        # cpptraj's eigvenvalues
        aa_eq(evals, cpp_eigenvalues)

        # cpptraj's eigvenvectors
        # use absolute values to avoid flipped sign
        # from Dan Roe
        # In practice, the "sign" of an eigenvector depends on the math library used to calculate it.
        # This is in fact why the modes command displacement test is disabled for cpptraj.
        # I bet if you use a different math library (e.g. use your system BLAS/LAPACK instead of the one bundled with Amber
        # or vice versa) you will get different signs.
        # Bottom line is that eigenvector sign doesn't matter.

        aa_eq(np.abs(evecs[:, ::-1].T), np.abs(cpp_eigenvectors), decimal=4)
        data = _ired(state_vecs, modes=(cpp_eigenvalues, cpp_eigenvectors))
        order_s2 = data['IRED_00127[S2]']

        # load cpptraj's output and compare to pytraj' values for S2 order paramters
        cpp_order_s2 = np.loadtxt(os.path.join(cpptraj_test_dir, 'Test_IRED',
                                               'orderparam.save')).T[-1]
        aa_eq(order_s2, cpp_order_s2, decimal=5)

    @unittest.skip('do not test now, get nan in some runs')
    def test_ired_lapack_in_numpy(self):
        parmfile = parm_dir
        trajfile = traj_dir

        # load to TrajectoryIterator
        traj = pt.iterload(trajfile, parmfile)

        # create N-H vectors
        h_indices = pt.select_atoms('@H', traj.top)
        n_indices = pt.select_atoms('@H', traj.top) - 1
        nh_indices = list(zip(n_indices, h_indices))

        # compute N-H vectors and ired matrix
        vecs_and_mat = pt.ired_vector_and_matrix(traj,
                                                 mask=nh_indices,
                                                 order=2)
        state_vecs = vecs_and_mat[:-1].values
        mat_ired = vecs_and_mat[-1]
        mat_ired /= mat_ired[0, 0]

        # cpptraj
        data_cpp = pt.matrix.diagonalize(mat_ired, n_vecs=len(state_vecs))[0]
        print(data_cpp.eigenvectors)

        # numpy
        data_np = pt.matrix._diag_np(mat_ired, n_vecs=len(state_vecs))

        def order_(modes):
            data = _ired(state_vecs, modes=modes)
            order_s2_v0 = data['IRED_00127[S2]']
            # make sure the S2 values is 1st array

            # load cpptraj's output and compare to pytraj' values for S2 order paramters
            cpp_order_s2 = np.loadtxt(os.path.join(
                cpptraj_test_dir, 'Test_IRED', 'orderparam.save')).T[-1]
            aa_eq(order_s2_v0.values, cpp_order_s2, decimal=4)

        order_(data_cpp.values)

        def plot_(x, y):
            import seaborn as sb
            sb.heatmap(x - y)
            pt.show()

        print((data_cpp.values[1] - data_np[1]).shape)


if __name__ == "__main__":
    unittest.main()
