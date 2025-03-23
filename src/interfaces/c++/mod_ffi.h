#pragma once

#include <cstddef>
#include <functional>
#include <memory>


namespace CEA {

  using FFI_CEA_Problem_Ptr = void*;
  using FFI_CEA_Problem_Unique_Ptr = std::unique_ptr<void, std::function<void(const FFI_CEA_Problem_Ptr&)>>;

  struct FFI_CEA_Problem_Array {
    FFI_CEA_Problem_Ptr* addr;
    size_t size;
  };

  struct FFI_SizeT_Array {
    const size_t* addr;
    size_t size;
  };

  struct FFI_Double_Array {
    const double* addr;
    size_t size;
  };

  extern "C" {
    FFI_CEA_Problem_Ptr ffi_cea_new_problem();
    void ffi_cea_del_problem(const FFI_CEA_Problem_Ptr&);

    void ffi_cea_set_problem(const FFI_CEA_Problem_Unique_Ptr&, const char*, const char* = nullptr, const bool& = false, const bool& = false, const bool& = false,
                             const bool& = false, const bool& = false, const char* = nullptr, const char* = nullptr);

    void ffi_cea_set_output_options(const FFI_CEA_Problem_Unique_Ptr&, const bool& = true, const FFI_SizeT_Array& = {},
                                    const bool& = false, const bool& = false, const double& = 0, const bool& = false);

    void ffi_cea_set_chamber_pressures(const FFI_CEA_Problem_Unique_Ptr&, const FFI_Double_Array&, const char* = nullptr);

    void ffi_cea_set_mixture_ratios(const FFI_CEA_Problem_Unique_Ptr&, const FFI_Double_Array&, const char* = nullptr);
    void ffi_cea_set_pressure_ratios(const FFI_CEA_Problem_Unique_Ptr&, const FFI_Double_Array&);
    void ffi_cea_set_subsonic_area_ratios(const FFI_CEA_Problem_Unique_Ptr&, const FFI_Double_Array&);
    void ffi_cea_set_supersonic_area_ratios(const FFI_CEA_Problem_Unique_Ptr&, const FFI_Double_Array&);

    void ffi_cea_set_finite_area_combustor(const FFI_CEA_Problem_Unique_Ptr&, const double* = nullptr, const double* = nullptr);

    void ffi_cea_add_reactant(const FFI_CEA_Problem_Unique_Ptr&, const char*, const char*, const char* = nullptr, const double* = nullptr,
                              const double* = nullptr, const double* = nullptr, const double* = nullptr, const double* = nullptr,
                              const char* = nullptr, const char* = nullptr, const char* = nullptr, const char* = nullptr);

    void ffi_cea_set_insert_species(const FFI_CEA_Problem_Unique_Ptr&, const char* const*, const FFI_SizeT_Array&);
    void ffi_cea_set_omit_species(const FFI_CEA_Problem_Unique_Ptr&, const char* const*, const FFI_SizeT_Array&);
    void ffi_cea_set_only_species(const FFI_CEA_Problem_Unique_Ptr&, const char* const*, const FFI_SizeT_Array&);

    void ffi_cea_set_legacy_mode(const FFI_CEA_Problem_Unique_Ptr&, const bool& = true);
    void ffi_cea_run(const FFI_CEA_Problem_Unique_Ptr&, const char* = nullptr);

    double ffi_cea_get_temperature(const FFI_CEA_Problem_Unique_Ptr&, const size_t& iOF, const size_t& ipt);
    double ffi_cea_get_chamber_temperature(const FFI_CEA_Problem_Unique_Ptr&);
    double ffi_cea_get_molecular_weight(const FFI_CEA_Problem_Unique_Ptr&, const size_t& iOF, const size_t& ipt);
    double ffi_cea_get_specific_heat_ratio(const FFI_CEA_Problem_Unique_Ptr&, const size_t& iOF, const size_t& ipt);

    void ffi_cea_calc_frozen_exhaust(const FFI_CEA_Problem_Unique_Ptr&, const double&, const double* = nullptr, const double* = nullptr, const double* = nullptr);

    void ffi_cea_write_debug_output(const FFI_CEA_Problem_Unique_Ptr&, const char*);

    FFI_CEA_Problem_Array ffi_cea_read_legacy_input(const char*);
    void ffi_cea_deallocate_array_ptr(FFI_CEA_Problem_Array&);
    void ffi_cea_run_all_cases(FFI_CEA_Problem_Array&, const char* = nullptr, const char* = nullptr);

    size_t ffi_cea_sizeof(const FFI_CEA_Problem_Unique_Ptr&);
  }

}
