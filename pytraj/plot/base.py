from __future__ import absolute_import

from .font_config import *  # call rc('font',**font)
from ..utils import _import

try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = None
