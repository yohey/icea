# -*- coding: utf-8 -*-

from ctypes import POINTER, Structure, byref, cast, c_bool, c_char_p, c_double, c_size_t, c_void_p
from typing import Iterable

from ._findlib import _libcea


class FFI_C_Ptr_Array(Structure):
    _fields_ = [('addr', c_void_p),
                ('size', c_size_t)]


class CEA_Problem:

    def __init__(self):
        _libcea.ffi_cea_new_problem.restype = c_void_p
        self._ffi = c_void_p(_libcea.ffi_cea_new_problem())

    def __del__(self):
        _libcea.ffi_cea_del_problem.argtypes = [c_void_p]
        _libcea.ffi_cea_del_problem(self._ffi)

    def set_problem(self, mode, name = None, mole_ratios = None, equilibrium = None, ions = None,
                    frozen = None, frozen_at_throat = None, thermo_lib = None, trans_lib = None):
        _libcea.ffi_cea_set_problem.argtypes = [c_void_p, c_char_p, c_char_p, POINTER(c_bool), POINTER(c_bool), POINTER(c_bool),
                                                POINTER(c_bool), POINTER(c_bool), c_char_p, c_char_p]
        _libcea.ffi_cea_set_problem(self._ffi, mode.encode(), _c_char_or_None(name), _c_bool_or_None(mole_ratios), _c_bool_or_None(equilibrium), _c_bool_or_None(ions),
                                    _c_bool_or_None(frozen), _c_bool_or_None(frozen_at_throat), _c_char_or_None(thermo_lib), _c_char_or_None(trans_lib))

    def set_output_options(self, SI = None, debug_points = None, mass_fractions = None, short = None, trace_tol = None, transport = None):
        _libcea.ffi_cea_set_output_options.argtypes = [c_void_p, POINTER(c_bool), POINTER(FFI_C_Ptr_Array), POINTER(c_bool), POINTER(c_bool), POINTER(c_double), POINTER(c_bool)]
        _libcea.ffi_cea_set_output_options(self._ffi, _c_bool_or_None(SI), _c_size_t_array_or_None(debug_points), _c_bool_or_None(mass_fractions),
                                           _c_bool_or_None(short), _c_double_or_None(trace_tol), _c_bool_or_None(transport))

    def set_chamber_pressures(self, pressures, unit = None):
        _libcea.ffi_cea_set_chamber_pressures.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array), c_char_p]
        _libcea.ffi_cea_set_chamber_pressures(self._ffi, _c_double_array(pressures), _c_char_or_None(unit))

    def set_mixture_ratios(self, ratios, type_ = None):
        _libcea.ffi_cea_set_mixture_ratios.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array), c_char_p]
        _libcea.ffi_cea_set_mixture_ratios(self._ffi, _c_double_array(ratios), _c_char_or_None(type_))

    def set_pressure_ratios(self, ratios):
        _libcea.ffi_cea_set_pressure_ratios.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_pressure_ratios(self._ffi, _c_double_array(ratios))

    def set_subsonic_area_ratios(self, ratios):
        _libcea.ffi_cea_set_subsonic_area_ratios.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_subsonic_area_ratios(self._ffi, _c_double_array(ratios))

    def set_supersonic_area_ratios(self, ratios):
        _libcea.ffi_cea_set_supersonic_area_ratios.argtypes = [c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_supersonic_area_ratios(self._ffi, _c_double_array(ratios))

    def set_finite_area_combustor(self, contraction_ratio = None, mass_flow_ratio = None):
        _libcea.ffi_cea_set_finite_area_combustor.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)]
        _libcea.ffi_cea_set_finite_area_combustor(self._ffi, _c_double_or_None(contraction_ratio), _c_double_or_None(mass_flow_ratio))

    def add_reactant(self, type_, name, formula = None, ratio = None, T = None, rho = None, h = None, u = None,
                     T_unit = None, rho_unit = None, h_unit = None, u_unit = None):
        _libcea.ffi_cea_add_reactant.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p, POINTER(c_double),
                                                 POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                                 c_char_p, c_char_p, c_char_p, c_char_p]
        _libcea.ffi_cea_add_reactant(self._ffi, type_.encode(), name.encode(), _c_char_or_None(formula), _c_double_or_None(ratio),
                                     _c_double_or_None(T), _c_double_or_None(rho), _c_double_or_None(h), _c_double_or_None(u),
                                     _c_char_or_None(T_unit), _c_char_or_None(rho_unit), _c_char_or_None(h_unit), _c_char_or_None(u_unit))

    def set_reactant(self, index, ratio = None, T = None, rho = None, h = None, u = None,
                     T_unit = None, rho_unit = None, h_unit = None, u_unit = None):
        _libcea.ffi_cea_set_reactant.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_double),
                                                 POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                                 c_char_p, c_char_p, c_char_p, c_char_p]
        _libcea.ffi_cea_set_reactant(self._ffi, byref(c_size_t(index + 1)), _c_double_or_None(ratio),
                                     _c_double_or_None(T), _c_double_or_None(rho), _c_double_or_None(h), _c_double_or_None(u),
                                     _c_char_or_None(T_unit), _c_char_or_None(rho_unit), _c_char_or_None(h_unit), _c_char_or_None(u_unit))

    def set_insert_species(self, species):
        _libcea.ffi_cea_set_insert_species.argtypes = [c_void_p, c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_insert_species(self._ffi, *_c_char_array(species))

    def set_omit_species(self, species):
        _libcea.ffi_cea_set_omit_species.argtypes = [c_void_p, c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_omit_species(self._ffi, *_c_char_array(species))

    def set_only_species(self, species):
        _libcea.ffi_cea_set_only_species.argtypes = [c_void_p, c_void_p, POINTER(FFI_C_Ptr_Array)]
        _libcea.ffi_cea_set_only_species(self._ffi, *_c_char_array(species))

    def set_legacy_mode(self, legacy_mode = None):
        _libcea.ffi_cea_set_legacy_mode.argtypes = [c_void_p, POINTER(c_bool)]
        _libcea.ffi_cea_set_legacy_mode(self._ffi, _c_bool_or_None(legacy_mode))

    def run(self, out_filename = None):
        _libcea.ffi_cea_run.argtypes = [c_void_p, c_char_p]
        _libcea.ffi_cea_run(self._ffi, _c_char_or_None(out_filename))

    def get_pressure(self, iOF, ipt):
        _libcea.ffi_cea_get_pressure.restype = c_double
        _libcea.ffi_cea_get_pressure.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_pressure(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_temperature(self, iOF, ipt):
        _libcea.ffi_cea_get_temperature.restype = c_double
        _libcea.ffi_cea_get_temperature.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_temperature(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_chamber_temperature(self):
        _libcea.ffi_cea_get_chamber_temperature.restype = c_double
        _libcea.ffi_cea_get_chamber_temperature.argtypes = [c_void_p]
        return _libcea.ffi_cea_get_chamber_temperature(self._ffi)

    def get_molecular_weight(self, iOF, ipt):
        _libcea.ffi_cea_get_molecular_weight.restype = c_double
        _libcea.ffi_cea_get_molecular_weight.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_molecular_weight(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_specific_heat(self, iOF, ipt):
        _libcea.ffi_cea_get_specific_heat.restype = c_double
        _libcea.ffi_cea_get_specific_heat.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_specific_heat(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_specific_heat_ratio(self, iOF, ipt):
        _libcea.ffi_cea_get_specific_heat_ratio.restype = c_double
        _libcea.ffi_cea_get_specific_heat_ratio.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t)]
        return _libcea.ffi_cea_get_specific_heat_ratio(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)))

    def get_characteristic_velocity(self):
        _libcea.ffi_cea_get_characteristic_velocity.restype = c_double
        _libcea.ffi_cea_get_characteristic_velocity.argtypes = [c_void_p]
        return _libcea.ffi_cea_get_characteristic_velocity(self._ffi)

    def get_specific_impulse(self, iOF, ipt, vacuum):
        _libcea.ffi_cea_get_specific_impulse.restype = c_double
        _libcea.ffi_cea_get_specific_impulse.argtypes = [c_void_p, POINTER(c_size_t), POINTER(c_size_t), POINTER(c_bool)]
        return _libcea.ffi_cea_get_specific_impulse(self._ffi, byref(c_size_t(iOF + 1)), byref(c_size_t(ipt + 1)), _c_bool_or_None(vacuum))

    def calc_frozen_exhaust(self, T):
        _libcea.ffi_cea_calc_frozen_exhaust.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        _cp, _gamma, _mu, _k, _Pr = c_double(), c_double(), c_double(), c_double(), c_double()
        _libcea.ffi_cea_calc_frozen_exhaust(self._ffi, _c_double_or_None(T), byref(_cp), byref(_gamma), byref(_mu), byref(_k), byref(_Pr))
        return _cp.value, _gamma.value, _mu.value, _k.value, _Pr.value

    def get_thermo_reference_properties(self, name: str):
        _libcea.ffi_cea_get_thermo_reference_properties.argtypes = [c_void_p, c_char_p, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        _M, _T_ref, _h0_ref = c_double(), c_double(), c_double()
        _libcea.ffi_cea_get_thermo_reference_properties(self._ffi, name.encode(), byref(_M), byref(_T_ref), byref(_h0_ref))
        return _M.value, _T_ref.value, _h0_ref.value


    def write_debug_output(self, filename: str):
        _libcea.ffi_cea_write_debug_output.argtypes = [c_void_p, c_char_p]
        _libcea.ffi_cea_write_debug_output(self._ffi, filename.encode())


def _c_bool_or_None(b: bool | None):
    if b is None:
        return None
    else:
        return byref(c_bool(b))

def _c_double_or_None(x: float | None):
    if x is None:
        return None
    else:
        return byref(c_double(x))

def _c_char_or_None(s: str | None):
    if s is None:
        return None
    else:
        return s.encode()

def _c_double_array(a: Iterable[float]):
    _a = (c_double * len(a))(*a)
    ptr = FFI_C_Ptr_Array(addr = cast(_a, c_void_p), size = len(a))
    return byref(ptr)

def _c_size_t_array_or_None(a: Iterable[int] | None):
    if a is None:
        return None
    else:
        _a = (c_size_t * len(a))(*a)
        ptr = FFI_C_Ptr_Array(addr = cast(_a, c_void_p), size = len(a))
        return byref(ptr)

def _c_char_array(a: Iterable[str]):
    char_array_ptr = (c_char_p * len(a))(*[s.encode() for s in a])
    char_length_array = (c_size_t * len(a))(*[len(s) for s in a])
    char_length_array_ptr = FFI_C_Ptr_Array(addr = cast(char_length_array, c_void_p), size = len(a))
    return byref(char_array_ptr), byref(char_length_array_ptr)
