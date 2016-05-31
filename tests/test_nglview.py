#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from ipykernel.comm import Comm
import ipywidgets as widgets

from traitlets import TraitError
from ipywidgets import Widget

#-----------------------------------------------------------------------------
# Utility stuff from ipywidgets tests
#-----------------------------------------------------------------------------

class DummyComm(Comm):
    comm_id = 'a-b-c-d'

    def open(self, *args, **kwargs):
        pass

    def send(self, *args, **kwargs):
        pass

    def close(self, *args, **kwargs):
        pass

_widget_attrs = {}
displayed = []
undefined = object()

def setup():
    _widget_attrs['_comm_default'] = getattr(Widget, '_comm_default', undefined)
    Widget._comm_default = lambda self: DummyComm()
    _widget_attrs['_ipython_display_'] = Widget._ipython_display_
    def raise_not_implemented(*args, **kwargs):
        raise NotImplementedError()
    Widget._ipython_display_ = raise_not_implemented


def teardown():
    for attr, value in _widget_attrs.items():
        if value is undefined:
            delattr(Widget, attr)
        else:
            setattr(Widget, attr, value)


def test_nglview():
    def copy_coordinates_dict(view):
        return dict((k, v.copy())
                for (k, v) in view.coordinates_dict.items())

    traj0 = pt.iterload("data/tz2.nc", "data/tz2.parm7")
    traj1 = pt.datafiles.load_trpcage()
    traj2 = pt.datafiles.load_rna()

    view = traj0.view()
    view.add(traj1)
    view.add_trajectory(traj2)

    frame_index = 2
    view.frame = frame_index
    coordinates_dict = copy_coordinates_dict(view)

    aa_eq(coordinates_dict[0], traj0[frame_index].xyz)
    aa_eq(coordinates_dict[1], traj1[frame_index].xyz)
    aa_eq(coordinates_dict[2], traj2[frame_index].xyz)
