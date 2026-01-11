#pragma once

#include <cstddef>


using _FFI_CEA_Problem_Ptr = void*;

struct _FFI_CEA_Problem_Array {
  _FFI_CEA_Problem_Ptr* addr;
  size_t size;
};

struct _FFI_SizeT_Array {
  const size_t* addr;
  size_t size;
};

struct _FFI_Double_Array {
  const double* addr;
  size_t size;
};

#ifdef __cplusplus
extern "C" {
#endif

  _FFI_CEA_Problem_Ptr _ffi_cea_new_problem();
  void _ffi_cea_del_problem(const _FFI_CEA_Problem_Ptr);
  _FFI_CEA_Problem_Ptr _ffi_cea_replica(const _FFI_CEA_Problem_Ptr);

  void _ffi_cea_set_problem(const _FFI_CEA_Problem_Ptr, const char*, const char* = nullptr, const bool& = false, const bool& = false, const bool& = false,
                            const bool& = false, const bool& = false, const bool& = false, const bool& = false, const char* = nullptr, const char* = nullptr);

  void _ffi_cea_set_output_options(const _FFI_CEA_Problem_Ptr, const bool& = true, const _FFI_SizeT_Array& = {},
                                   const bool& = false, const bool& = false, const double& = 0, const bool& = false,
                                   const char* const* = nullptr, const _FFI_SizeT_Array& = {});

  void _ffi_cea_set_chamber_pressures(const _FFI_CEA_Problem_Ptr, const _FFI_Double_Array&, const char* = nullptr);
  void _ffi_cea_set_chamber_temperatures(const _FFI_CEA_Problem_Ptr, const _FFI_Double_Array&, const char* = nullptr);
  void _ffi_cea_set_chamber_densities(const _FFI_CEA_Problem_Ptr, const _FFI_Double_Array&, const char* = nullptr);
  void _ffi_cea_set_chamber_internal_energy(const _FFI_CEA_Problem_Ptr, const double&);

  void _ffi_cea_set_mixture_ratios(const _FFI_CEA_Problem_Ptr, const _FFI_Double_Array&, const char* = nullptr);
  void _ffi_cea_set_pressure_ratios(const _FFI_CEA_Problem_Ptr, const _FFI_Double_Array&);
  void _ffi_cea_set_subsonic_area_ratios(const _FFI_CEA_Problem_Ptr, const _FFI_Double_Array&);
  void _ffi_cea_set_supersonic_area_ratios(const _FFI_CEA_Problem_Ptr, const _FFI_Double_Array&);

  void _ffi_cea_set_finite_area_combustor(const _FFI_CEA_Problem_Ptr, const double* = nullptr, const double* = nullptr);

  void _ffi_cea_set_initial_velocities(const _FFI_CEA_Problem_Ptr, const _FFI_Double_Array&, const bool& = false);

  void _ffi_cea_add_reactant(const _FFI_CEA_Problem_Ptr, const char*, const char*, const char* = nullptr, const double* = nullptr,
                             const double* = nullptr, const double* = nullptr, const double* = nullptr, const double* = nullptr,
                             const char* = nullptr, const char* = nullptr, const char* = nullptr, const char* = nullptr);

  void _ffi_cea_set_reactant(const _FFI_CEA_Problem_Ptr, const size_t&, const double* = nullptr,
                             const double* = nullptr, const double* = nullptr, const double* = nullptr, const double* = nullptr,
                             const char* = nullptr, const char* = nullptr, const char* = nullptr, const char* = nullptr);

  void _ffi_cea_set_insert_species(const _FFI_CEA_Problem_Ptr, const char* const*, const _FFI_SizeT_Array&);
  void _ffi_cea_set_omit_species(const _FFI_CEA_Problem_Ptr, const char* const*, const _FFI_SizeT_Array&);
  void _ffi_cea_set_only_species(const _FFI_CEA_Problem_Ptr, const char* const*, const _FFI_SizeT_Array&);

  void _ffi_cea_set_legacy_mode(const _FFI_CEA_Problem_Ptr, const bool& = true);
  void _ffi_cea_run(const _FFI_CEA_Problem_Ptr, const char* = nullptr, const char* = nullptr);

  double _ffi_cea_get_pressure(const _FFI_CEA_Problem_Ptr, const size_t&, const size_t&);
  double _ffi_cea_get_temperature(const _FFI_CEA_Problem_Ptr, const size_t&, const size_t&);
  double _ffi_cea_get_chamber_temperature(const _FFI_CEA_Problem_Ptr);
  double _ffi_cea_get_molecular_weight(const _FFI_CEA_Problem_Ptr, const size_t&, const size_t&);
  double _ffi_cea_get_specific_heat(const _FFI_CEA_Problem_Ptr, const size_t&, const size_t&);
  double _ffi_cea_get_specific_heat_ratio(const _FFI_CEA_Problem_Ptr, const size_t&, const size_t&);
  double _ffi_cea_get_characteristic_velocity(const _FFI_CEA_Problem_Ptr);
  double _ffi_cea_get_specific_impulse(const _FFI_CEA_Problem_Ptr, const size_t&, const size_t&, const bool& = false);

  void _ffi_cea_calc_frozen_exhaust(const _FFI_CEA_Problem_Ptr, const double&, const double* = nullptr, const double* = nullptr, const double* = nullptr, const double* = nullptr, const double* = nullptr);
  void _ffi_cea_get_thermo_reference_properties(const _FFI_CEA_Problem_Ptr, const char*, double* const, double* const, double* const);

  void _ffi_cea_write_debug_output(const _FFI_CEA_Problem_Ptr, const char*);

  _FFI_CEA_Problem_Array _ffi_cea_read_legacy_input(const char*);
  void _ffi_cea_deallocate_array_ptr(_FFI_CEA_Problem_Array&);
  void _ffi_cea_run_all_cases(_FFI_CEA_Problem_Array&, const char* = nullptr, const char* = nullptr);

  size_t _ffi_cea_sizeof(const _FFI_CEA_Problem_Ptr);

#ifdef __cplusplus
}
#endif
