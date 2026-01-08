#pragma once

typedef void* _FFI_CEA_Problem_Ptr;

#ifdef __cplusplus
extern "C" {
#endif

  _FFI_CEA_Problem_Ptr _ffi_cea_new_problem();
  void _ffi_cea_del_problem(_FFI_CEA_Problem_Ptr);

  int _ffi_cea_sizeof_int(_FFI_CEA_Problem_Ptr);

#ifdef __cplusplus
}
#endif
