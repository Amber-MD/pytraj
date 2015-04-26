# distutils: language = c++

from pytraj._common_actions import calculate
import pytraj.common_actions as pyca
from pytraj.utils import is_int

class ActionInTraj(object):
    def calc_distance(self, mask=""):
        return pyca.calc_distance(mask, self)

    def calc_distrmsd(self, mask=""):
        return pyca.calc_distrmsd(mask, self)

    def calc_radgyr(self, mask=""):
        return pyca.calc_radgyr(mask, self)

    def calc_angle(self, mask=""):
        return pyca.calc_angle(mask, self)

    def calc_matrix(self, mask=""):
        return pyca.calc_matrix(mask, self)

    def calc_dssp(self, mask="", *args, **kwd):
        return pyca.calc_dssp(mask, self, *args, **kwd)

    def calc_dihedral(self, mask=""):
        return pyca.calc_dihedral(mask, self)
    
    def calc_multidihedral(self, mask=""):
        return pyca.calc_multidihedral(mask, self)

    def calc_molsurf(self, mask=""):
        return pyca.calc_molsurf(mask, self)

    def calc_center_of_mass(self, mask=""):
        return pyca.calc_center_of_mass(mask, self)

    def calc_COM(self, mask=""):
        return pyca.calc_center_of_mass(mask, self)

    def calc_center_of_geometry(self, mask=""):
        return pyca.calc_center_of_geometry(mask, self)

    def calc_COG(self, mask=""):
        return pyca.calc_center_of_geometry(mask, self)

    def calc_vector(self, mask=""):
        from pytraj.actions.Action_Vector import Action_Vector
        from pytraj.DataSetList import DataSetList
        act = Action_Vector()
        dslist = DataSetList()

        if 'name' not in mask:
            # for some reasons, I got segmentation fault without 'name' keyword
            # need to check cpptraj code
            mask = "myvector " + mask
        act(mask, self, dslist=dslist)
        dslist.set_py_free_mem(False)
        return dslist[0]

    def calc_pairwise_rmsd(self, mask=""):
        return pyca.calc_pairwise_rmsd(mask, self)

    def calc_rmsd(self, mask="", ref=None, mass=False, fit=True):
        """"""
        assert ref is not None
        if is_int(ref):
            # index
            ref = self[ref]
        return pyca.calc_rmsd(mask, traj=self, ref=ref, mass=mass, fit=fit)

    def search_hbonds(self, mask="*"):
        return pyca.search_hbonds(self, mask)

    def get_average_frame(self, mask=""):
        return pyca.get_average_frame(mask, self)

    def calc_temperatures(self, mask=""):
        return pyca.calc_temperatures(mask, self)

    def calc_watershell(self, mask=""):
        return pyca.calc_watershell(mask, self)

    @property
    def temperatures(self):
        return pyca.calc_temperatures(self, "frame")
