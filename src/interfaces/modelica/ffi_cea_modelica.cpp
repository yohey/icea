
#include <limits>

#include "mod_ffi.h"


#ifdef __cplusplus
extern "C" {
#endif

  int _ffi_cea_sizeof_int(void* p)
  {
    size_t n = _ffi_cea_sizeof(p);

    if (n > static_cast<size_t>(std::numeric_limits<int>::max())) {
      return -1;
    }

    return static_cast<int>(n);
  };

#ifdef __cplusplus
}
#endif
