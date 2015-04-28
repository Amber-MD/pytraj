from __future__ import absolute_import

from .font_config import * # call rc('font',**font)
from ..utils import _import
from ..utils.check_and_assert import require
_, plt = _import("matplotlib.pyplot")
_, np = _import("numpy")
