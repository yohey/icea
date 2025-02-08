module mod_ffi
  use, intrinsic:: iso_c_binding
  implicit none

  private:: get_string

  type, bind(C):: CEA_Problem_Array
     type(c_ptr):: addr
     integer(c_size_t):: size
  end type CEA_Problem_Array

contains

  function new_problem() result(ptr) bind(C, name = "ffi_cea_new_problem")
    use mod_cea_types, only: CEA_Problem

    type(c_ptr):: ptr
    type(CEA_Problem), pointer:: cea
    allocate(cea)

    ptr = c_loc(cea)

    return
  end function new_problem


  subroutine del_problem(ptr) bind(C, name = "ffi_cea_del_problem")
    use mod_cea_types, only: CEA_Problem

    type(c_ptr), intent(in):: ptr

    type(CEA_Problem), pointer:: cea
    call c_f_pointer(ptr, cea)

    if (associated(cea)) then
!!$       write(0, *) '[DEBUG] CEA_Problem (f90) destructor is called: ', trim(cea%Case)
       deallocate(cea)
    end if

    return
  end subroutine del_problem


  subroutine print_case(ptr) bind(C, name = "ffi_cea_print_case")
    use mod_cea_types, only: CEA_Problem

    type(c_ptr), intent(in):: ptr

    type(CEA_Problem), pointer:: cea
    call c_f_pointer(ptr, cea)

    if (associated(cea)) then
       write(0, '("Case = ", a)') trim(cea%Case)
    end if

    return
  end subroutine print_case


  function ffi_read_legacy_input(filename) result(array) bind(C, name = "ffi_cea_read_legacy_input")
    use mod_cea_types, only: CEA_Problem, init_case, IOINP
    use mod_legacy_io, only: count_cases, read_legacy_case

    type(CEA_Problem_Array):: array
    character(len = 1, kind = c_char), intent(in):: filename(*)
    character(len = :, kind = c_char), allocatable:: filename_

    integer:: icase, num_cases
    type(CEA_Problem), pointer:: cea, cea_prev
    type(c_ptr), pointer:: ptr(:)
    logical:: file_exists

    filename_ = get_string(filename)

    inquire(file = filename_, exist = file_exists)
    if (.not. file_exists) then
       print *, filename_, ' DOES NOT EXIST'
       error stop
    end if

    call count_cases(filename_, num_cases)

    allocate(ptr(num_cases))

    open(newunit = IOINP, file = filename_, status = 'old', form = 'formatted', action = 'read')

    do icase = 1, num_cases
       allocate(cea)

       call init_case(cea)

       !! TEMPORARY WORK AROUND TO REPRODUCE KNOWN BUG !!
       if (icase > 1 .and. associated(cea_prev)) then
          cea%Dens(:) = cea_prev%Dens(:)
       end if
       !!!!!!!!!!!!!!!!! TO BE DELETED !!!!!!!!!!!!!!!!!!

       call read_legacy_case(cea)

       ptr(icase) = c_loc(cea)

       cea_prev => cea
       nullify(cea)
    end do

    nullify(cea_prev)

    close(IOINP)

    array%addr = c_loc(ptr(1))
    array%size = num_cases

    return
  end function ffi_read_legacy_input


  subroutine ffi_run_all_cases(array, out_filename, plt_filename) bind(C, name = "ffi_cea_run_all_cases")
    use mod_cea_types, only: CEA_Problem
    use mod_cea_core, only: run_all_cases

    type(CEA_Problem_Array):: array
    character(len = 1, kind = c_char), intent(in), optional:: out_filename(*)
    character(len = 1, kind = c_char), intent(in), optional:: plt_filename(*)

    type(CEA_Problem), pointer:: cea(:)

    call c_f_pointer(array%addr, cea, [array%size])

    if (present(out_filename) .and. present(plt_filename)) then
       call run_all_cases(cea, get_string(out_filename), get_string(plt_filename))
    else if (present(out_filename)) then
       call run_all_cases(cea, get_string(out_filename))
    else if (present(plt_filename)) then
       call run_all_cases(cea, plt_filename = get_string(plt_filename))
    else
       call run_all_cases(cea)
    end if

    return
  end subroutine ffi_run_all_cases


  function ffi_sizeof(ptr) result(size) bind(C, name = "ffi_cea_sizeof")
    use mod_cea_types, only: CEA_Problem

    integer(c_size_t):: size
    type(c_ptr), intent(in):: ptr

    type(CEA_Problem), pointer:: cea
    call c_f_pointer(ptr, cea)

    if (associated(cea)) then
       size = storage_size(cea) / 8  ! bit -> byte
    end if

    return
  end function ffi_sizeof


  function get_string(cstring) result(fstring)
    character(len = 1, kind = c_char), intent(in):: cstring(*)
    character(len = :, kind = c_char), allocatable:: fstring

    integer:: i, length = 1

    do while (cstring(length) /= c_null_char)
       length = length + 1
    end do

    allocate(character(len = length, kind = c_char):: fstring)

    do concurrent (i = 1:length)
       fstring(i:i) = cstring(i)
    end do

    return
  end function get_string

end module mod_ffi
