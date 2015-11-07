from __future__ import absolute_import

from .font_config import font, rc

try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = None
