
#include "cea.h"

#include <iostream>
#include <cstring>


namespace CEA {

  Problem::Problem() {
    this->_ffi = ffi_cea_new_problem();
  }

  Problem::Problem(FFI_CEA_Problem_Ptr ptr) {
    this->_ffi = ptr;
  }

  Problem::~Problem() {
    ffi_cea_del_problem(this->_ffi);
  }


  void Problem::set_problem(const std::string& mode, const std::optional<std::string>& name, const bool& mole_ratios, const bool& equilibrium, const bool& ions,
                            const bool& frozen, const bool& frozen_at_throat, const std::optional<std::string>& thermo_lib, const std::optional<std::string>& trans_lib) {
    const FFI_String mode_ffi{mode.c_str(), mode.length()};
    const FFI_String* name_ffi = (name) ? new FFI_String{name.value().c_str(), name.value().length()} : nullptr;
    const FFI_String* thermo_lib_ffi = (thermo_lib) ? new FFI_String{thermo_lib.value().c_str(), thermo_lib.value().length()} : nullptr;
    const FFI_String* trans_lib_ffi = (trans_lib) ? new FFI_String{trans_lib.value().c_str(), trans_lib.value().length()} : nullptr;

    ffi_cea_set_problem(this->_ffi, mode_ffi, name_ffi, mole_ratios, equilibrium, ions, frozen, frozen_at_throat, thermo_lib_ffi, trans_lib_ffi);

    if (name_ffi) delete name_ffi;
    if (thermo_lib_ffi) delete thermo_lib_ffi;
    if (trans_lib_ffi) delete trans_lib_ffi;
  }

  void Problem::set_output_options(const bool& SI, const std::vector<size_t>& debug_points, const bool& mass_fractions, const bool& _short,
                                   const double& trace_tol, const bool& transport) {
    FFI_SizeT_Array array{debug_points.data(), debug_points.size()};
    ffi_cea_set_output_options(this->_ffi, SI, array, mass_fractions, _short, trace_tol, transport);
  }

  void Problem::set_chamber_pressures(const std::vector<double>& pressure_list, const std::optional<std::string>& unit) {
    const FFI_Double_Array array{pressure_list.data(), pressure_list.size()};
    const FFI_String* unit_ffi = (unit) ? new FFI_String{unit.value().c_str(), unit.value().length()} : nullptr;

    ffi_cea_set_chamber_pressures(this->_ffi, array, unit_ffi);

    if (unit_ffi) delete unit_ffi;
  }

  void Problem::set_mixture_ratios(const std::vector<double>& ratio_list, const std::optional<std::string>& type) {
    FFI_Double_Array array{ratio_list.data(), ratio_list.size()};
    const FFI_String* type_ffi = (type) ? new FFI_String{type.value().c_str(), type.value().length()} : nullptr;

    ffi_cea_set_mixture_ratios(this->_ffi, array, type_ffi);

    if (type_ffi) delete type_ffi;
  }

  void Problem::set_pressure_ratios(const std::vector<double>& ratio_list) {
    FFI_Double_Array array{ratio_list.data(), ratio_list.size()};
    ffi_cea_set_pressure_ratios(this->_ffi, array);
  }

  void Problem::set_subsonic_area_ratios(const std::vector<double>& ratio_list) {
    FFI_Double_Array array{ratio_list.data(), ratio_list.size()};
    ffi_cea_set_subsonic_area_ratios(this->_ffi, array);
  }

  void Problem::set_supersonic_area_ratios(const std::vector<double>& ratio_list) {
    FFI_Double_Array array{ratio_list.data(), ratio_list.size()};
    ffi_cea_set_supersonic_area_ratios(this->_ffi, array);
  }

  void Problem::set_finite_area_combustor(const std::optional<double>& contraction_ratio, const std::optional<double>& mass_flow_ratio) {
    const double* contraction_ratio_ptr = (contraction_ratio) ? &(*contraction_ratio) : nullptr;
    const double* mass_flow_ratio_ptr = (mass_flow_ratio) ? &(*mass_flow_ratio) : nullptr;
    ffi_cea_set_finite_area_combustor(this->_ffi, contraction_ratio_ptr, mass_flow_ratio_ptr);
  }

  void Problem::add_reactant(const std::string& type, const std::string& name, const std::optional<std::string>& formula, const std::optional<double>& ratio,
                             const std::optional<double>& T, const std::optional<double>& rho,
                             const std::optional<double>& h, const std::optional<double>& u,
                             const std::optional<std::string>& T_unit, const std::optional<std::string>& rho_unit,
                             const std::optional<std::string>& h_unit, const std::optional<std::string>& u_unit) {
    const FFI_String type_ffi{type.c_str(), type.length()};
    const FFI_String name_ffi{name.c_str(), name.length()};
    const FFI_String* formula_ffi = (formula) ? new FFI_String{formula.value().c_str(), formula.value().length()} : nullptr;
    const double* ratio_ptr = (ratio) ? &(*ratio) : nullptr;
    const double* T_ptr = (T) ? &(*T) : nullptr;
    const double* rho_ptr = (rho) ? &(*rho) : nullptr;
    const double* h_ptr = (h) ? &(*h) : nullptr;
    const double* u_ptr = (u) ? &(*u) : nullptr;
    const FFI_String* T_unit_ffi = (T_unit) ? new FFI_String{T_unit.value().c_str(), T_unit.value().length()} : nullptr;
    const FFI_String* rho_unit_ffi = (rho_unit) ? new FFI_String{rho_unit.value().c_str(), rho_unit.value().length()} : nullptr;
    const FFI_String* h_unit_ffi = (h_unit) ? new FFI_String{h_unit.value().c_str(), h_unit.value().length()} : nullptr;
    const FFI_String* u_unit_ffi = (u_unit) ? new FFI_String{u_unit.value().c_str(), u_unit.value().length()} : nullptr;

    ffi_cea_add_reactant(this->_ffi, type_ffi, name_ffi, formula_ffi, ratio_ptr, T_ptr, rho_ptr, h_ptr, u_ptr, T_unit_ffi, rho_unit_ffi, h_unit_ffi, u_unit_ffi);

    if (formula_ffi) delete formula_ffi;
    if (T_unit_ffi) delete T_unit_ffi;
    if (rho_unit_ffi) delete rho_unit_ffi;
    if (h_unit_ffi) delete h_unit_ffi;
    if (u_unit_ffi) delete u_unit_ffi;
  }


  void Problem::insert_species(const std::vector<std::string>& species) {
    std::vector<size_t> char_length_array(species.size());
    std::vector<const char*> char_ptr_array(species.size());

    for (size_t i = 0; i < species.size(); ++i) {
      char_length_array[i] = species[i].length();
      char_ptr_array[i] = species[i].c_str();
    }

    ffi_cea_insert_species(this->_ffi, char_ptr_array.data(), char_length_array.data(), species.size());
  }


  void Problem::set_omit_species(const std::vector<std::string>& species) {
    std::vector<size_t> char_length_array(species.size());
    std::vector<const char*> char_ptr_array(species.size());

    for (size_t i = 0; i < species.size(); ++i) {
      char_length_array[i] = species[i].length();
      char_ptr_array[i] = species[i].c_str();
    }

    ffi_cea_set_omit_species(this->_ffi, char_ptr_array.data(), char_length_array.data(), species.size());
  }


  void Problem::set_only_species(const std::vector<std::string>& species) {
    std::vector<size_t> char_length_array(species.size());
    std::vector<const char*> char_ptr_array(species.size());

    for (size_t i = 0; i < species.size(); ++i) {
      char_length_array[i] = species[i].length();
      char_ptr_array[i] = species[i].c_str();
    }

    ffi_cea_set_only_species(this->_ffi, char_ptr_array.data(), char_length_array.data(), species.size());
  }


  void Problem::set_legacy_mode(const bool& legacy_mode) {
    ffi_cea_set_legacy_mode(this->_ffi, legacy_mode);
  }


  void Problem::run(const std::optional<std::string>& out_filename) {
    const FFI_String* filename_ffi = (out_filename) ? new FFI_String{out_filename.value().c_str(), out_filename.value().length()} : nullptr;

    ffi_cea_run(this->_ffi, filename_ffi);

    if (filename_ffi) delete filename_ffi;
  }


  void Problem::write_debug_output(const std::string& filename) {
    const FFI_String filename_ffi{filename.c_str(), filename.length()};
    ffi_cea_write_debug_output(this->_ffi, filename_ffi);
  }


  std::vector<Problem> read_legacy_input(const std::string& filename) {
    const FFI_String filename_ffi{filename.c_str(), filename.length()};
    FFI_CEA_Problem_Array array = ffi_cea_read_legacy_input(filename_ffi);

    std::vector<Problem> problems;
    problems.reserve(array.size);

    for (size_t i = 0; i < array.size; ++i) {
      problems.emplace_back(array.addr[i]);
    }

    ffi_cea_deallocate_array_ptr(array);

    return problems;
  }


  void run_all_cases(std::vector<Problem>& problems, const std::optional<std::string>& out_filename, const std::optional<std::string>& plt_filename) {
    FFI_CEA_Problem_Array array;
    array.size = problems.size();

    bool _need_aligned_buffer = ffi_need_aligned(problems);
    size_t obj_size = ffi_cea_sizeof(problems[0]._ffi);

    if (_need_aligned_buffer) {
      array.addr = static_cast<FFI_CEA_Problem_Ptr*>(malloc(problems.size() * obj_size));

      for (size_t i = 0; i < problems.size(); ++i) {
        FFI_CEA_Problem_Ptr ptr = reinterpret_cast<FFI_CEA_Problem_Ptr>((reinterpret_cast<uintptr_t>(array.addr) + i * obj_size));
        std::memcpy(ptr, problems[i]._ffi, obj_size);
      }

    } else {
      array.addr = static_cast<FFI_CEA_Problem_Ptr*>(problems[0]._ffi);
    }

    const FFI_String* out_filename_ffi = (out_filename) ? new FFI_String{out_filename.value().c_str(), out_filename.value().length()} : nullptr;
    const FFI_String* plt_filename_ffi = (plt_filename) ? new FFI_String{plt_filename.value().c_str(), plt_filename.value().length()} : nullptr;

    ffi_cea_run_all_cases(array, out_filename_ffi, plt_filename_ffi);

    if (out_filename_ffi) delete out_filename_ffi;
    if (plt_filename_ffi) delete plt_filename_ffi;

    if (_need_aligned_buffer) {
      for (size_t i = 0; i < problems.size(); ++i) {
        FFI_CEA_Problem_Ptr ptr = reinterpret_cast<FFI_CEA_Problem_Ptr>((reinterpret_cast<uintptr_t>(array.addr) + i * obj_size));
        std::memcpy(problems[i]._ffi, ptr, obj_size);
      }

      free(array.addr);
    }
  }


  bool ffi_need_aligned(const std::vector<Problem>& problems) {
    if (problems.size() < 2) return false;

    size_t obj_size = ffi_cea_sizeof(problems[0]._ffi);
    uintptr_t head_addr = reinterpret_cast<uintptr_t>(problems[0]._ffi);
    uintptr_t prev_tail_addr = head_addr + obj_size;

    for (size_t i = 1; i < problems.size(); ++i) {
      head_addr = reinterpret_cast<uintptr_t>(problems[i]._ffi);
      if (head_addr != prev_tail_addr) {
        return true;
      } else {
        prev_tail_addr = head_addr + obj_size;
      }
    }

    return false;
  }
}
