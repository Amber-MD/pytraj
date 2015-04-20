# distutils: language = c++

from pytraj._common_actions import calculate
import pytraj.common_actions as pyca

class ActionInTraj(object):
    def calc_distance(self, mask=""):
        return pyca.calc_distance(mask, self).tolist()

    def calc_distrmsd(self, mask=""):
        return pyca.calc_distrmsd(mask, self).tolist()

    def calc_radgyr(self, mask=""):
        return pyca.calc_radgyr(mask, self).tolist()

    def calc_angle(self, mask=""):
        return pyca.calc_angle(mask, self).tolist()

    def calc_matrix(self, mask=""):
        return pyca.calc_matrix(mask, self)

    def calc_dssp(self, mask="", *args, **kwd):
        return pyca.calc_dssp(mask, self, *args, **kwd)

    def calc_dihedral(self, mask=""):
        return pyca.calc_dihedral(mask, self).tolist()
    
    def calc_multidihedral(self, mask=""):
        return pyca.calc_multidihedral(mask, self)

    def calc_molsurf(self, mask=""):
        return pyca.calc_molsurf(mask, self).tolist()

    def calc_center_of_mass(self, mask=""):
        return pyca.calc_center_of_mass(mask, self).tolist()

    def calc_COM(self, mask=""):
        return pyca.calc_center_of_mass(mask, self).tolist()

    def calc_center_of_geometry(self, mask=""):
        return pyca.calc_center_of_geometry(mask, self).tolist()

    def calc_COG(self, mask=""):
        return pyca.calc_center_of_geometry(mask, self).tolist()

    def calc_vector(self, mask=""):
        from pytraj.actions.Action_Vector import Action_Vector
        from pytraj.DataSetList import DataSetList
        act = Action_Vector()
        dslist = DataSetList()

        if 'name' not in mask:
            # for some reasons, I got segmentation fault without 'name' keyword
            # need to check cpptraj code
            mask = "name myvector " + mask
        act(mask, self, dslist=dslist)
        return [d0.tolist() for d0 in dslist]

    def search_hbonds(self, mask="*"):
        return pyca.search_hbonds(self, mask)

    def get_average_frame(self, mask=""):
        return pyca.get_average_frame(mask, self)

    def calc_watershell(self, mask=""):
        return pyca.calc_watershell(mask, self).tolist()
