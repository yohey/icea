# -*- coding: utf-8 -*-

import os
import sys

_libpath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

if sys.platform == 'linux':
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = _libpath + ';' + os.environ['LD_LIBRARY_PATH']
    else:
        os.environ['LD_LIBRARY_PATH'] = _libpath

elif sys.platform == 'darwin':
    if 'DYLD_LIBRARY_PATH' in os.environ:
        os.environ['DYLD_LIBRARY_PATH'] = _libpath + ';' + os.environ['DYLD_LIBRARY_PATH']
    else:
        os.environ['DYLD_LIBRARY_PATH'] = _libpath

elif sys.platform == 'win32':
    if 'PATH' in os.environ:
        os.environ['PATH'] = _libpath + ';' + os.environ['PATH']
    else:
        os.environ['PATH'] = _libpath


import ctypes
from ctypes.util import find_library

_libceapath = find_library('cea')


if _libceapath is None:
    raise RuntimeError('Library is not found.')


if sys.platform == 'win32':
    _libcea = ctypes.WinDLL(_libceapath)
else:
    _libcea = ctypes.CDLL(_libceapath)
