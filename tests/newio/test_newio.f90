program test_newio
  use mod_cea
  use mod_types
  implicit none

  call test_orig

  stop
end program test_newio


subroutine test_orig
  use mod_cea
  use mod_types
  use mod_io
  use mod_legacy_io
  implicit none

  type(CEA_Problem), allocatable:: cea(:)
  character(*), parameter:: inp_filename = '../../../tests/orig/cea2.inp'
  character(*), parameter:: out_filename = 'orig.out'

  call read_legacy_input(cea, inp_filename)

!!$  cea(8)%legacy_mode = .false.
  call run_all_cases(cea(8:8), out_filename)

  call write_debug_output(cea(8), 'orig-08.debug.out')

  deallocate(cea)

  return
end subroutine test_orig
