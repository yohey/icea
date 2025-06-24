# -*- coding: utf-8 -*-

import os
import sys

_libpath = os.path.dirname(os.path.dirname(__file__))

if sys.platform == 'linux':
    if 'LD_LIBRARY_PATH' in os.environ:
        if _libpath not in os.environ['LD_LIBRARY_PATH']:
            os.environ['LD_LIBRARY_PATH'] = _libpath + ':' + os.environ['LD_LIBRARY_PATH']
    else:
        os.environ['LD_LIBRARY_PATH'] = _libpath

elif sys.platform == 'darwin':
    if 'DYLD_LIBRARY_PATH' in os.environ:
        if _libpath not in os.environ['DYLD_LIBRARY_PATH']:
            os.environ['DYLD_LIBRARY_PATH'] = _libpath + ':' + os.environ['DYLD_LIBRARY_PATH']
    else:
        os.environ['DYLD_LIBRARY_PATH'] = _libpath

elif sys.platform == 'win32':
    if 'PATH' in os.environ:
        if _libpath not in os.environ['PATH']:
            os.environ['PATH'] = _libpath + ';' + os.environ['PATH']
    else:
        os.environ['PATH'] = _libpath


import ctypes
from ctypes.util import find_library

_libceapath = find_library('cea')

if _libceapath is None:
    _libceapath = find_library('libcea')

if _libceapath is None:
    raise RuntimeError('Library is not found.')

if not os.path.isabs(_libceapath):
    if sys.platform == 'linux':
        for d in os.environ['LD_LIBRARY_PATH'].split(':'):
            if os.path.exists(os.path.join(d, _libceapath)):
                _libceapath = os.path.join(d, _libceapath)
                break

    elif sys.platform == 'darwin':
        for d in os.environ['DYLD_LIBRARY_PATH'].split(':'):
            if os.path.exists(os.path.join(d, _libceapath)):
                _libceapath = os.path.join(d, _libceapath)
                break

    elif sys.platform == 'win32':
        for d in os.environ['PATH'].split(';'):
            if os.path.exists(os.path.join(d, _libceapath)):
                _libceapath = os.path.join(d, _libceapath)
                break

if sys.platform == 'win32':
    os.add_dll_directory(os.path.dirname(_libceapath))

    _mingw64path = find_library('libgfortran-5')
    if not os.path.isabs(_mingw64path):
        for d in os.environ['PATH'].split(';'):
            if os.path.exists(os.path.join(d, _mingw64path)):
                _mingw64path = os.path.join(d, _mingw64path)

    os.add_dll_directory(os.path.dirname(_mingw64path))

_libcea = ctypes.CDLL(_libceapath)
