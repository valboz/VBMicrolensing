from .VBMicrolensing import VBMicrolensing

import inspect
import os

_fil_in=inspect.getfile(VBMicrolensing)
VBMicrolensing.SetESPLtablefile(os.path.join(os.path.join(os.path.dirname(_fil_in), 'data'), 'ESPL.tbl'))
VBMicrolensing.SetSuntablefile(os.path.join(os.path.join(os.path.dirname(_fil_in), 'data'), 'SunEphemeris.txt'))

from importlib.metadata import version

__version__ = version("VBMicrolensing")
