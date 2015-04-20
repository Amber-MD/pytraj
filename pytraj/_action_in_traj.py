# distutils: language = c++

from pytraj._common_actions import calculate
from pytraj.hbonds import search_hbonds

class ActionInTraj(object):
    def calc_matrix(self, mask=""):
        return calculate("matrix", mask, self)

    def search_hbonds(self, mask="*"):
        return search_hbonds(self, mask)
