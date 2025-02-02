#pragma once

#include <cstddef>


namespace CEA {

  using FFI_CEA_Problem_Ptr = void*;

  struct FFI_CEA_Problem_Array {
    FFI_CEA_Problem_Ptr* addr;
    size_t size;
  };

  extern "C" {
    FFI_CEA_Problem_Ptr ffi_cea_new_problem();
    void ffi_cea_del_problem(const FFI_CEA_Problem_Ptr&);

    void ffi_cea_print_case(const FFI_CEA_Problem_Ptr&);

    FFI_CEA_Problem_Array ffi_cea_read_legacy_input(const char*);
    void ffi_cea_run_all_cases(FFI_CEA_Problem_Array&, const char* = nullptr, const char* = nullptr);

    size_t ffi_cea_sizeof(const FFI_CEA_Problem_Ptr&);
  }

}
