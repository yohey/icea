
#include <limits>
#include <vector>

#include "mod_ffi.h"
#include "ffi_cea_modelica.h"


extern "C" {

  _FFI_CEA_Problem_Ptr _ffi_cea_new_problem_C_impl()
  {
    return _ffi_cea_new_problem();
  }

  void _ffi_cea_del_problem_C_impl(const _FFI_CEA_Problem_Ptr p)
  {
    _ffi_cea_del_problem(p);
  }

  void _ffi_cea_set_problem_C_impl(const _FFI_CEA_Problem_Ptr p, const char* mode, const char* name, int mole_ratios, int equilibrium, int ions,
                                   int frozen, int frozen_at_throat, int incident, int reflected)
  {
    const bool _mole_ratios      = (mole_ratios != 0);
    const bool _equilibrium      = (equilibrium != 0);
    const bool _ions             = (ions != 0);
    const bool _frozen           = (frozen != 0);
    const bool _frozen_at_throat = (frozen_at_throat != 0);
    const bool _incident         = (incident != 0);
    const bool _reflected        = (reflected != 0);
    _ffi_cea_set_problem(p, mode, name, _mole_ratios, _equilibrium, _ions, _frozen, _frozen_at_throat, _incident, _reflected);
  }

  void _ffi_cea_set_output_options_C_impl(const _FFI_CEA_Problem_Ptr p, int SI, const int* debug_points, size_t num_debug_points,
                                          int mass_fractions, int _short, double trace_tol, int transport,
                                          const char** plot, size_t num_plot)
  {
    const bool _SI = (SI != 0);
    const bool _mass_fractions = (mass_fractions != 0);
    const bool __short = (_short != 0);
    const bool _transport = (transport != 0);

    std::vector<size_t> _debug_points_size_t;
    _debug_points_size_t.reserve(num_debug_points);
    for (size_t i = 0; i < num_debug_points; ++i) {
      _debug_points_size_t.push_back(static_cast<size_t>(debug_points[i]));
    }
    _FFI_SizeT_Array _debug_points{_debug_points_size_t.data(), _debug_points_size_t.size()};

    std::vector<size_t> _char_length_array(num_plot);
    std::vector<const char*> _char_ptr_array(num_plot);
    for (size_t i = 0; i < num_plot; ++i) {
      _char_length_array[i] = strlen(plot[i]);
      _char_ptr_array[i] = plot[i];
    }
    _FFI_SizeT_Array _char_length_array_ptr{_char_length_array.data(), _char_length_array.size()};

    _ffi_cea_set_output_options(p, _SI, _debug_points, _mass_fractions, __short, trace_tol, _transport, _char_ptr_array.data(), _char_length_array_ptr);
  }

  void _ffi_cea_set_chamber_pressures_C_impl(const _FFI_CEA_Problem_Ptr p, const double* pressures, size_t num_pressures, const char* unit)
  {
    _FFI_Double_Array _pressures{pressures, num_pressures};
    const char* _unit = strlen(unit) > 0 ? unit : nullptr;
    _ffi_cea_set_chamber_pressures(p, _pressures, _unit);
  }


  void _ffi_cea_set_mixture_ratios_C_impl(const _FFI_CEA_Problem_Ptr p, const double* ratios, size_t num_ratios, const char* type)
  {
    _FFI_Double_Array _ratios{ratios, num_ratios};
    const char* _type = strlen(type) > 0 ? type : nullptr;
    _ffi_cea_set_mixture_ratios(p, _ratios, _type);
  }

  void _ffi_cea_set_pressure_ratios_C_impl(const _FFI_CEA_Problem_Ptr p, const double* ratios, size_t num_ratios)
  {
    _FFI_Double_Array _ratios{ratios, num_ratios};
    _ffi_cea_set_pressure_ratios(p, _ratios);
  }

  void _ffi_cea_set_subsonic_area_ratios_C_impl(const _FFI_CEA_Problem_Ptr p, const double* ratios, size_t num_ratios)
  {
    _FFI_Double_Array _ratios{ratios, num_ratios};
    _ffi_cea_set_subsonic_area_ratios(p, _ratios);
  }

  void _ffi_cea_set_supersonic_area_ratios_C_impl(const _FFI_CEA_Problem_Ptr p, const double* ratios, size_t num_ratios)
  {
    _FFI_Double_Array _ratios{ratios, num_ratios};
    _ffi_cea_set_supersonic_area_ratios(p, _ratios);
  }

  void _ffi_cea_add_reactant_C_impl(const _FFI_CEA_Problem_Ptr p, const char* type, const char* name, const char* formula, double ratio,
                                    double T, double rho, double h, double u, const char* T_unit, const char* rho_unit, const char* h_unit, const char* u_unit)
  {
    double* _ratio = std::isnan(ratio) ? nullptr : &ratio;
    double* _T = std::isnan(T) ? nullptr : &T;
    double* _rho = std::isnan(rho) ? nullptr : &rho;
    double* _h = std::isnan(h) ? nullptr : &h;
    double* _u = std::isnan(u) ? nullptr : &u;
    const char* _T_unit = strlen(T_unit) > 0 ? T_unit : nullptr;
    const char* _rho_unit = strlen(rho_unit) > 0 ? rho_unit : nullptr;
    const char* _h_unit = strlen(h_unit) > 0 ? h_unit : nullptr;
    const char* _u_unit = strlen(u_unit) > 0 ? u_unit : nullptr;
    _ffi_cea_add_reactant(p, type, name, formula, _ratio, _T, _rho, _h, _u, _T_unit, _rho_unit, _h_unit, _u_unit);
  }

  void _ffi_cea_set_reactant_C_impl(const _FFI_CEA_Problem_Ptr p, int index, double ratio,
                                    double T, double rho, double h, double u, const char* T_unit, const char* rho_unit, const char* h_unit, const char* u_unit)
  {
    size_t _index = static_cast<size_t>(index);
    double* _ratio = std::isnan(ratio) ? nullptr : &ratio;
    double* _T = std::isnan(T) ? nullptr : &T;
    double* _rho = std::isnan(rho) ? nullptr : &rho;
    double* _h = std::isnan(h) ? nullptr : &h;
    double* _u = std::isnan(u) ? nullptr : &u;
    const char* _T_unit = strlen(T_unit) > 0 ? T_unit : nullptr;
    const char* _rho_unit = strlen(rho_unit) > 0 ? rho_unit : nullptr;
    const char* _h_unit = strlen(h_unit) > 0 ? h_unit : nullptr;
    const char* _u_unit = strlen(u_unit) > 0 ? u_unit : nullptr;
    _ffi_cea_set_reactant(p, _index, _ratio, _T, _rho, _h, _u, _T_unit, _rho_unit, _h_unit, _u_unit);
  }

  void _ffi_cea_set_legacy_mode_C_impl(const _FFI_CEA_Problem_Ptr p, int legacy_mode)
  {
    const bool _legacy_mode = (legacy_mode != 0);
    _ffi_cea_set_legacy_mode(p, _legacy_mode);
  }

  void _ffi_cea_run_C_impl(const _FFI_CEA_Problem_Ptr p, const char* out_filename, const char* plt_filename)
  {
    const char* _out_filename = strlen(out_filename) > 0 ? out_filename : nullptr;
    const char* _plt_filename = strlen(plt_filename) > 0 ? plt_filename : nullptr;
    _ffi_cea_run(p, _out_filename, _plt_filename);
  }

  double _ffi_cea_get_chamber_temperature_C_impl(const _FFI_CEA_Problem_Ptr p)
  {
    return _ffi_cea_get_chamber_temperature(p);
  }

  void _ffi_cea_write_debug_output_C_impl(const _FFI_CEA_Problem_Ptr p, const char* filename)
  {
    _ffi_cea_write_debug_output(p, filename);
  }

  void _ffi_cea_get_thermo_reference_properties_C_impl(const _FFI_CEA_Problem_Ptr p, const char* name, double* M, double* T_ref, double* h0_ref)
  {
    _ffi_cea_get_thermo_reference_properties(p, name, M, T_ref, h0_ref);
    *M *= 1e-3;
  }

  int _ffi_cea_sizeof_C_impl(const _FFI_CEA_Problem_Ptr p)
  {
    size_t n = _ffi_cea_sizeof(p);

    if (n > static_cast<size_t>(std::numeric_limits<int>::max())) {
      return -1;
    }

    return static_cast<int>(n);
  }

  double _ffi_cea_quiet_NaN_C_impl()
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

}
