# distutils: language = c++

import pytraj.common_actions as pyca
from pytraj._get_common_objects import _get_data_from_dtype
from pytraj.utils import is_int


class ActionTrajectory(object):

    def calc_distance(self, mask="", *args, **kwd):
        return pyca.calc_distance(self, mask, *args, **kwd)

    def calc_radgyr(self, mask="", *args, **kwd):
        return pyca.calc_radgyr(self, mask, *args, **kwd)

    def calc_angle(self, mask="", *args, **kwd):
        return pyca.calc_angle(self, mask)

    def calc_matrix(self, mask="", *args, **kwd):
        return pyca.calc_matrix(self, mask, *args, **kwd)

    def distance_matrix(self, mask="", *args, **kwd):
        from . matrix_analysis import distance_matrix
        return distance_matrix(self, mask, *args, **kwd)

    def calc_dssp(self, mask="", *args, **kwd):
        return pyca.calc_dssp(self, mask, *args, **kwd)

    def calc_dihedral(self, mask="", *args, **kwd):
        return pyca.calc_dihedral(self, mask, *args, **kwd)

    def calc_multidihedral(self, mask="", *args, **kwd):
        return pyca.calc_multidihedral(self, mask, *args, **kwd)

    def calc_molsurf(self, mask="", *args, **kwd):
        return pyca.calc_molsurf(self, mask, *args, **kwd)

    def calc_center_of_mass(self, mask="", *args, **kwd):
        return pyca.calc_center_of_mass(self, mask, *args, **kwd)

    def calc_COM(self, mask="", *args, **kwd):
        return pyca.calc_center_of_mass(self, mask, *args, **kwd)

    def calc_center_of_geometry(self, mask="", *args, **kwd):
        return pyca.calc_center_of_geometry(self, mask, *args, **kwd)

    def calc_COG(self, mask="", *args, **kwd):
        return pyca.calc_center_of_geometry(self, mask, *args, **kwd)

    def calc_vector(self, mask="", dtype='dataset', *args, **kwd):
        """
        dtype = {'dataset', 'list', 'ndarray'}
        """
        from pytraj.actions.CpptrajActions import Action_Vector
        from pytraj.datasets.DataSetList import DataSetList as CpptrajDatasetList
        act = Action_Vector()
        dslist = CpptrajDatasetList()

        act(mask, self, dslist=dslist, *args, **kwd)
        dtype = dtype.lower()
        return _get_data_from_dtype(dslist, dtype)

    def calc_pairwise_rmsd(self, mask="", *args, **kwd):
        return pyca.calc_pairwise_rmsd(self, mask, *args, **kwd)

    def calc_distrmsd(self, mask="", *args, **kwd):
        return pyca.calc_distrmsd(self, mask, *args, **kwd)

    def rmsd(self, ref=None, mask="", mass=False, fit=True, *args, **kwd):
        """"""
        if is_int(ref):
            # index
            ref = self[ref]
        return pyca.calc_rmsd(mask=mask, traj=self,
                              ref=ref, mass=mass, fit=fit, *args, **kwd)

    def calc_rmsd(self, *args, **kwd):
        return self.rmsd(*args, **kwd)

    def calc_bfactors(self, *args, **kwd):
        return pyca.calc_bfactors(*args, **kwd)

    def search_hbonds(self, mask="*", *args, **kwd):
        """return CpptrajDatasetList
        """
        return pyca.search_hbonds(self, mask, *args, **kwd)

    def average(self, mask="", *args, **kwd):
        """return Frame
        """
        return pyca.get_average_frame(self, mask, *args, **kwd)

    @property
    def temperatures(self):
        return pyca.calc_temperatures(self, "frame")
