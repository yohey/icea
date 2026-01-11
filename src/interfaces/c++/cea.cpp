
#include "cea.h"

#include <cstring>


namespace CEA {

  Problem::Problem(): _ffi(_ffi_cea_new_problem(), [](void* p){ _ffi_cea_del_problem(p); }) {
  }

  Problem::Problem(Problem&& other): _ffi(std::move(other._ffi)) {
  }

  Problem::Problem(const Problem& other): _ffi(_ffi_cea_replica(other._ffi.get()), [](void* p){ _ffi_cea_del_problem(p); }) {
  }

  Problem::Problem(_FFI_CEA_Problem_Ptr ptr): _ffi(ptr, [](void* p){ _ffi_cea_del_problem(p); }) {
  }


  void Problem::set_problem(const char* mode, const std::optional<const char*>& name, const bool& mole_ratios, const bool& equilibrium, const bool& ions,
                            const bool& frozen, const bool& frozen_at_throat, const bool& incident, const bool& reflected,
                            const std::optional<const char*>& thermo_lib, const std::optional<const char*>& trans_lib) {
    const char* name_ffi = name.value_or(nullptr);
    const char* thermo_lib_ffi = thermo_lib.value_or(nullptr);
    const char* trans_lib_ffi = trans_lib.value_or(nullptr);

    _ffi_cea_set_problem(this->_ffi.get(), mode, name_ffi, mole_ratios, equilibrium, ions, frozen, frozen_at_throat, incident, reflected, thermo_lib_ffi, trans_lib_ffi);
  }

  void Problem::set_output_options(const bool& SI, const std::vector<size_t>& debug_points, const bool& mass_fractions, const bool& _short,
                                   const double& trace_tol, const bool& transport, const std::vector<std::string>& plot) {
    std::vector<size_t> debug_points_fidx;
    debug_points_fidx.reserve(debug_points.size());
    for (size_t i: debug_points) { debug_points_fidx.push_back(i + 1); }
    _FFI_SizeT_Array array{debug_points_fidx.data(), debug_points_fidx.size()};

    std::vector<size_t> char_length_array(plot.size());
    std::vector<const char*> char_ptr_array(plot.size());

    for (size_t i = 0; i < plot.size(); ++i) {
      char_length_array[i] = plot[i].length();
      char_ptr_array[i] = plot[i].c_str();
    }

    _FFI_SizeT_Array char_length_array_ptr{char_length_array.data(), char_length_array.size()};

    _ffi_cea_set_output_options(this->_ffi.get(), SI, array, mass_fractions, _short, trace_tol, transport, char_ptr_array.data(), char_length_array_ptr);
  }

  void Problem::set_chamber_pressures(const std::vector<double>& pressure_list, const std::optional<const char*>& unit) {
    const _FFI_Double_Array array{pressure_list.data(), pressure_list.size()};
    _ffi_cea_set_chamber_pressures(this->_ffi.get(), array, unit.value_or(nullptr));
  }

  void Problem::set_chamber_temperatures(const std::vector<double>& temperature_list, const std::optional<const char*>& unit) {
    const _FFI_Double_Array array{temperature_list.data(), temperature_list.size()};
    _ffi_cea_set_chamber_temperatures(this->_ffi.get(), array, unit.value_or(nullptr));
  }

  void Problem::set_chamber_densities(const std::vector<double>& density_list, const std::optional<const char*>& unit) {
    const _FFI_Double_Array array{density_list.data(), density_list.size()};
    _ffi_cea_set_chamber_densities(this->_ffi.get(), array, unit.value_or(nullptr));
  }

  void Problem::set_chamber_internal_energy(const double& uByR) {
    _ffi_cea_set_chamber_internal_energy(this->_ffi.get(), uByR);
  }

  void Problem::set_mixture_ratios(const std::vector<double>& ratio_list, const std::optional<const char*>& type) {
    _FFI_Double_Array array{ratio_list.data(), ratio_list.size()};
    _ffi_cea_set_mixture_ratios(this->_ffi.get(), array, type.value_or(nullptr));
  }

  void Problem::set_pressure_ratios(const std::vector<double>& ratio_list) {
    _FFI_Double_Array array{ratio_list.data(), ratio_list.size()};
    _ffi_cea_set_pressure_ratios(this->_ffi.get(), array);
  }

  void Problem::set_subsonic_area_ratios(const std::vector<double>& ratio_list) {
    _FFI_Double_Array array{ratio_list.data(), ratio_list.size()};
    _ffi_cea_set_subsonic_area_ratios(this->_ffi.get(), array);
  }

  void Problem::set_supersonic_area_ratios(const std::vector<double>& ratio_list) {
    _FFI_Double_Array array{ratio_list.data(), ratio_list.size()};
    _ffi_cea_set_supersonic_area_ratios(this->_ffi.get(), array);
  }

  void Problem::set_finite_area_combustor(const std::optional<double>& contraction_ratio, const std::optional<double>& mass_flow_ratio) {
    const double* contraction_ratio_ptr = (contraction_ratio) ? &contraction_ratio.value() : nullptr;
    const double* mass_flow_ratio_ptr = (mass_flow_ratio) ? &mass_flow_ratio.value() : nullptr;
    _ffi_cea_set_finite_area_combustor(this->_ffi.get(), contraction_ratio_ptr, mass_flow_ratio_ptr);
  }

  void Problem::set_initial_velocities(const std::vector<double>& velocity_list, const bool& is_mach) {
    _FFI_Double_Array array{velocity_list.data(), velocity_list.size()};
    _ffi_cea_set_initial_velocities(this->_ffi.get(), array, is_mach);
  }

  void Problem::add_reactant(const char* type, const char* name, const std::optional<const char*>& formula, const std::optional<double>& ratio,
                             const std::optional<double>& T, const std::optional<double>& rho,
                             const std::optional<double>& h, const std::optional<double>& u,
                             const std::optional<const char*>& T_unit, const std::optional<const char*>& rho_unit,
                             const std::optional<const char*>& h_unit, const std::optional<const char*>& u_unit) {
    const char* formula_ffi = formula.value_or(nullptr);
    const double* ratio_ptr = ratio.has_value() ? &ratio.value() : nullptr;
    const double* T_ptr = T.has_value() ? &T.value() : nullptr;
    const double* rho_ptr = rho.has_value() ? &rho.value() : nullptr;
    const double* h_ptr = h.has_value() ? &h.value() : nullptr;
    const double* u_ptr = u.has_value() ? &u.value() : nullptr;
    const char* T_unit_ffi = T_unit.value_or(nullptr);
    const char* rho_unit_ffi = rho_unit.value_or(nullptr);
    const char* h_unit_ffi = h_unit.value_or(nullptr);
    const char* u_unit_ffi = u_unit.value_or(nullptr);

    _ffi_cea_add_reactant(this->_ffi.get(), type, name, formula_ffi, ratio_ptr, T_ptr, rho_ptr, h_ptr, u_ptr, T_unit_ffi, rho_unit_ffi, h_unit_ffi, u_unit_ffi);
  }

  void Problem::set_reactant(const size_t& index, const std::optional<double>& ratio,
                             const std::optional<double>& T, const std::optional<double>& rho,
                             const std::optional<double>& h, const std::optional<double>& u,
                             const std::optional<const char*>& T_unit, const std::optional<const char*>& rho_unit,
                             const std::optional<const char*>& h_unit, const std::optional<const char*>& u_unit) {
    const double* ratio_ptr = ratio.has_value() ? &ratio.value() : nullptr;
    const double* T_ptr = T.has_value() ? &T.value() : nullptr;
    const double* rho_ptr = rho.has_value() ? &rho.value() : nullptr;
    const double* h_ptr = h.has_value() ? &h.value() : nullptr;
    const double* u_ptr = u.has_value() ? &u.value() : nullptr;
    const char* T_unit_ffi = T_unit.value_or(nullptr);
    const char* rho_unit_ffi = rho_unit.value_or(nullptr);
    const char* h_unit_ffi = h_unit.value_or(nullptr);
    const char* u_unit_ffi = u_unit.value_or(nullptr);

    _ffi_cea_set_reactant(this->_ffi.get(), index, ratio_ptr, T_ptr, rho_ptr, h_ptr, u_ptr, T_unit_ffi, rho_unit_ffi, h_unit_ffi, u_unit_ffi);
  }


  void Problem::set_insert_species(const std::vector<std::string>& species) {
    std::vector<size_t> char_length_array(species.size());
    std::vector<const char*> char_ptr_array(species.size());

    for (size_t i = 0; i < species.size(); ++i) {
      char_length_array[i] = species[i].length();
      char_ptr_array[i] = species[i].c_str();
    }

    _FFI_SizeT_Array char_length_array_ptr{char_length_array.data(), char_length_array.size()};

    _ffi_cea_set_insert_species(this->_ffi.get(), char_ptr_array.data(), char_length_array_ptr);
  }


  void Problem::set_omit_species(const std::vector<std::string>& species) {
    std::vector<size_t> char_length_array(species.size());
    std::vector<const char*> char_ptr_array(species.size());

    for (size_t i = 0; i < species.size(); ++i) {
      char_length_array[i] = species[i].length();
      char_ptr_array[i] = species[i].c_str();
    }

    _FFI_SizeT_Array char_length_array_ptr{char_length_array.data(), char_length_array.size()};

    _ffi_cea_set_omit_species(this->_ffi.get(), char_ptr_array.data(), char_length_array_ptr);
  }


  void Problem::set_only_species(const std::vector<std::string>& species) {
    std::vector<size_t> char_length_array(species.size());
    std::vector<const char*> char_ptr_array(species.size());

    for (size_t i = 0; i < species.size(); ++i) {
      char_length_array[i] = species[i].length();
      char_ptr_array[i] = species[i].c_str();
    }

    _FFI_SizeT_Array char_length_array_ptr{char_length_array.data(), char_length_array.size()};

    _ffi_cea_set_only_species(this->_ffi.get(), char_ptr_array.data(), char_length_array_ptr);
  }


  void Problem::set_legacy_mode(const bool& legacy_mode) {
    _ffi_cea_set_legacy_mode(this->_ffi.get(), legacy_mode);
  }


  void Problem::run(const std::optional<const char*>& out_filename, const std::optional<const char*>& plt_filename) {
    _ffi_cea_run(this->_ffi.get(), out_filename.value_or(nullptr), plt_filename.value_or(nullptr));
  }


  double Problem::get_pressure(const size_t& iOF, const size_t& ipt) {
    return _ffi_cea_get_pressure(this->_ffi.get(), iOF + 1, ipt + 1);
  }


  double Problem::get_temperature(const size_t& iOF, const size_t& ipt) {
    return _ffi_cea_get_temperature(this->_ffi.get(), iOF + 1, ipt + 1);
  }

  double Problem::get_chamber_temperature() {
    return _ffi_cea_get_chamber_temperature(this->_ffi.get());
  }


  double Problem::get_molecular_weight(const size_t& iOF, const size_t& ipt) {
    return _ffi_cea_get_molecular_weight(this->_ffi.get(), iOF + 1, ipt + 1);
  }


  double Problem::get_specific_heat(const size_t& iOF, const size_t& ipt) {
    return _ffi_cea_get_specific_heat(this->_ffi.get(), iOF + 1, ipt + 1);
  }


  double Problem::get_specific_heat_ratio(const size_t& iOF, const size_t& ipt) {
    return _ffi_cea_get_specific_heat_ratio(this->_ffi.get(), iOF + 1, ipt + 1);
  }


  double Problem::get_characteristic_velocity() {
    return _ffi_cea_get_characteristic_velocity(this->_ffi.get());
  }


  double Problem::get_specific_impulse(const size_t& iOF, const size_t& ipt, const bool& vacuum) {
    return _ffi_cea_get_specific_impulse(this->_ffi.get(), iOF + 1, ipt + 1, vacuum);
  }


  void Problem::calc_frozen_exhaust(const double& T, const std::optional<const double*>& cp_ptr, const std::optional<const double*>& gamma_ptr,
                                    const std::optional<const double*>& mu_ptr, const std::optional<const double*>& k_ptr, const std::optional<const double*>& Pr_ptr) {
    _ffi_cea_calc_frozen_exhaust(this->_ffi.get(), T, cp_ptr.value_or(nullptr), gamma_ptr.value_or(nullptr),
                                 mu_ptr.value_or(nullptr), k_ptr.value_or(nullptr), Pr_ptr.value_or(nullptr));
  }


  void Problem::get_thermo_reference_properties(const char* name, double* const M_ptr, double* const T_ref_ptr, double* const h0_ref_ptr) {
    _ffi_cea_get_thermo_reference_properties(this->_ffi.get(), name, M_ptr, T_ref_ptr, h0_ref_ptr);
  }


  void Problem::write_debug_output(const char* filename) {
    _ffi_cea_write_debug_output(this->_ffi.get(), filename);
  }


  std::vector<Problem> read_legacy_input(const char* filename) {
    _FFI_CEA_Problem_Array array = _ffi_cea_read_legacy_input(filename);

    std::vector<Problem> problems;
    problems.reserve(array.size);

    for (size_t i = 0; i < array.size; ++i) {
      problems.emplace_back(array.addr[i]);
    }

    _ffi_cea_deallocate_array_ptr(array);

    return problems;
  }


  void run_all_cases(std::vector<Problem>& problems, const std::optional<const char*>& out_filename, const std::optional<const char*>& plt_filename) {
    _FFI_CEA_Problem_Array array;
    array.size = problems.size();

    bool _need_aligned_buffer = ffi_need_aligned(problems);
    size_t obj_size = _ffi_cea_sizeof(problems[0]._ffi.get());

    if (_need_aligned_buffer) {
      array.addr = static_cast<_FFI_CEA_Problem_Ptr*>(malloc(problems.size() * obj_size));

      for (size_t i = 0; i < problems.size(); ++i) {
        _FFI_CEA_Problem_Ptr ptr = reinterpret_cast<_FFI_CEA_Problem_Ptr>((reinterpret_cast<uintptr_t>(array.addr) + i * obj_size));
        std::memcpy(ptr, problems[i]._ffi.get(), obj_size);
      }

    } else {
      array.addr = static_cast<_FFI_CEA_Problem_Ptr*>(problems[0]._ffi.get());
    }

    const char* out_filename_ffi = out_filename.value_or(nullptr);
    const char* plt_filename_ffi = plt_filename.value_or(nullptr);

    _ffi_cea_run_all_cases(array, out_filename_ffi, plt_filename_ffi);

    if (_need_aligned_buffer) {
      for (size_t i = 0; i < problems.size(); ++i) {
        _FFI_CEA_Problem_Ptr ptr = reinterpret_cast<_FFI_CEA_Problem_Ptr>((reinterpret_cast<uintptr_t>(array.addr) + i * obj_size));
        std::memcpy(problems[i]._ffi.get(), ptr, obj_size);
      }

      free(array.addr);
    }
  }


  bool ffi_need_aligned(const std::vector<Problem>& problems) {
    if (problems.size() < 2) return false;

    size_t obj_size = _ffi_cea_sizeof(problems[0]._ffi.get());
    uintptr_t head_addr = reinterpret_cast<uintptr_t>(problems[0]._ffi.get());
    uintptr_t prev_tail_addr = head_addr + obj_size;

    for (size_t i = 1; i < problems.size(); ++i) {
      head_addr = reinterpret_cast<uintptr_t>(problems[i]._ffi.get());
      if (head_addr != prev_tail_addr) {
        return true;
      } else {
        prev_tail_addr = head_addr + obj_size;
      }
    }

    return false;
  }
}
