#pragma once

typedef void* _FFI_CEA_Problem_Ptr;

#ifdef __cplusplus
extern "C" {
#endif

  _FFI_CEA_Problem_Ptr _ffi_cea_new_problem_C_impl();
  void _ffi_cea_del_problem_C_impl(const _FFI_CEA_Problem_Ptr);

  void _ffi_cea_set_problem_C_impl(const _FFI_CEA_Problem_Ptr, const char*, const char*, int, int, int, int, int, int, int);

  void _ffi_cea_set_output_options_C_impl(const _FFI_CEA_Problem_Ptr, int, const int*, size_t, int, int, double, int, const char**, size_t);

  void _ffi_cea_set_chamber_pressures_C_impl(const _FFI_CEA_Problem_Ptr, const double*, size_t, const char*);

  void _ffi_cea_set_mixture_ratios_C_impl(const _FFI_CEA_Problem_Ptr, const double*, size_t, const char*);
  void _ffi_cea_set_pressure_ratios_C_impl(const _FFI_CEA_Problem_Ptr, const double*, size_t);
  void _ffi_cea_set_subsonic_area_ratios_C_impl(const _FFI_CEA_Problem_Ptr, const double*, size_t);
  void _ffi_cea_set_supersonic_area_ratios_C_impl(const _FFI_CEA_Problem_Ptr, const double*, size_t);

  void _ffi_cea_add_reactant_C_impl(const _FFI_CEA_Problem_Ptr, const char*, const char*, const char*, double,
                                    double, double, double, double, const char*, const char*, const char*, const char*);

  void _ffi_cea_set_reactant_C_impl(const _FFI_CEA_Problem_Ptr, int, double, double, double, double, double,
                                    const char*, const char*, const char*, const char*);

  void _ffi_cea_set_legacy_mode_C_impl(const _FFI_CEA_Problem_Ptr, int);
  void _ffi_cea_run_C_impl(const _FFI_CEA_Problem_Ptr, const char*, const char*);

  double _ffi_cea_get_chamber_temperature_C_impl(const _FFI_CEA_Problem_Ptr);

  void _ffi_cea_get_thermo_reference_properties_C_impl(const _FFI_CEA_Problem_Ptr, const char*, double*, double*, double*);

  void _ffi_cea_write_debug_output_C_impl(const _FFI_CEA_Problem_Ptr, const char*);

  int _ffi_cea_sizeof_C_impl(const _FFI_CEA_Problem_Ptr);

  double _ffi_cea_quiet_NaN_C_impl();

#ifdef __cplusplus
}
#endif
