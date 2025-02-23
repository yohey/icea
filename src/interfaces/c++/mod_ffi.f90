module mod_ffi
  use, intrinsic:: iso_c_binding
  implicit none

  private:: get_string

  type, bind(C):: FFI_C_Ptr_Array
     type(c_ptr):: addr
     integer(c_size_t):: size
  end type FFI_C_Ptr_Array

  type, bind(C):: FFI_C_String
     type(c_ptr):: addr
     integer(c_size_t):: size
  end type FFI_C_String

contains

  function new_problem() result(ptr) bind(C, name = "ffi_cea_new_problem")
    use cea, only: CEA_Problem

    type(c_ptr):: ptr
    type(CEA_Problem), pointer:: cea
    allocate(cea)

    ptr = c_loc(cea)

    return
  end function new_problem


  subroutine del_problem(ptr) bind(C, name = "ffi_cea_del_problem")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr

    type(CEA_Problem), pointer:: cea
    call c_f_pointer(ptr, cea)

    if (associated(cea)) then
       deallocate(cea)
    end if

    return
  end subroutine del_problem


  subroutine set_problem(ptr, mode, name, mole_ratios, equilibrium, ions, &
       frozen, frozen_at_throat, thermo_lib, trans_lib) bind(C, name = "ffi_cea_set_problem")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_String), intent(in):: mode
    type(FFI_C_String), intent(in), optional:: name
    logical(c_bool), intent(in), optional:: mole_ratios
    logical(c_bool), intent(in), optional:: equilibrium
    logical(c_bool), intent(in), optional:: ions
    logical(c_bool), intent(in), optional:: frozen
    logical(c_bool), intent(in), optional:: frozen_at_throat
    type(FFI_C_String), intent(in), optional:: thermo_lib
    type(FFI_C_String), intent(in), optional:: trans_lib

    character(len = :, kind = c_char), allocatable:: mode_fstr
    character(len = :, kind = c_char), allocatable:: name_fstr
    character(len = :, kind = c_char), allocatable:: thermo_lib_fstr
    character(len = :, kind = c_char), allocatable:: trans_lib_fstr

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    mode_fstr = get_string(mode)
    if (present(name)) name_fstr = get_string(name)
    if (present(thermo_lib)) thermo_lib_fstr = get_string(thermo_lib)
    if (present(trans_lib)) trans_lib_fstr = get_string(trans_lib)

    if (associated(this)) then
       call this%set_problem(mode_fstr, name_fstr, logical(mole_ratios), logical(equilibrium), logical(ions), &
            logical(frozen), logical(frozen_at_throat), thermo_lib_fstr, trans_lib_fstr)
    end if

    if (allocated(mode_fstr)) deallocate(mode_fstr)
    if (allocated(name_fstr)) deallocate(name_fstr)
    if (allocated(thermo_lib_fstr)) deallocate(thermo_lib_fstr)
    if (allocated(trans_lib_fstr)) deallocate(trans_lib_fstr)

    return
  end subroutine set_problem


  subroutine set_output_options(ptr, SI, debug_points, mass_fractions, short, trace_tol, transport) &
       bind(C, name = "ffi_cea_set_output_options")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    logical(c_bool), intent(in), optional:: SI
    type(FFI_C_Ptr_Array), intent(in), optional:: debug_points
    logical(c_bool), intent(in), optional:: mass_fractions
    logical(c_bool), intent(in), optional:: short
    real(c_double), intent(in), optional:: trace_tol
    logical(c_bool), intent(in), optional:: transport

    type(CEA_Problem), pointer:: this
    integer(c_size_t), pointer:: debug_points_array(:)

    call c_f_pointer(ptr, this)
    call c_f_pointer(debug_points%addr, debug_points_array, [debug_points%size])

    if (associated(this)) then
       call this%set_output_options(logical(SI), int(debug_points_array), logical(mass_fractions), logical(short), &
            trace_tol, logical(transport))
    end if

    return
  end subroutine set_output_options


  subroutine set_chamber_pressures(ptr, pressures_ptr, unit) bind(C, name = "ffi_cea_set_chamber_pressures")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: pressures_ptr
    type(FFI_C_String), intent(in), optional:: unit

    type(CEA_Problem), pointer:: this
    real(c_double), pointer:: pressure_list(:)
    character(len = :, kind = c_char), allocatable:: unit_fstr

    call c_f_pointer(ptr, this)
    call c_f_pointer(pressures_ptr%addr, pressure_list, [pressures_ptr%size])

    if (present(unit)) unit_fstr = get_string(unit)

    if (associated(this)) then
       call this%set_chamber_pressures(pressure_list, unit_fstr)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_chamber_pressure)'
#endif
    end if

    return
  end subroutine set_chamber_pressures


  subroutine set_mixture_ratios(ptr, ratios_ptr, type) bind(C, name = "ffi_cea_set_mixture_ratios")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: ratios_ptr
    type(FFI_C_String), intent(in), optional:: type

    type(CEA_Problem), pointer:: this
    real(c_double), pointer:: ratio_list(:)
    character(len = :, kind = c_char), allocatable:: type_fstr

    call c_f_pointer(ptr, this)
    call c_f_pointer(ratios_ptr%addr, ratio_list, [ratios_ptr%size])

    if (present(type)) type_fstr = get_string(type)

    if (associated(this)) then
       call this%set_mixture_ratios(ratio_list, type_fstr)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_mixture_ratio)'
#endif
    end if

    return
  end subroutine set_mixture_ratios


  subroutine set_pressure_ratios(ptr, ratios_ptr) bind(C, name = "ffi_cea_set_pressure_ratios")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: ratios_ptr

    type(CEA_Problem), pointer:: this
    real(c_double), pointer:: ratio_list(:)

    call c_f_pointer(ptr, this)
    call c_f_pointer(ratios_ptr%addr, ratio_list, [ratios_ptr%size])

    if (associated(this)) then
       call this%set_pressure_ratios(ratio_list)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_pressure_ratios)'
#endif
    end if

    return
  end subroutine set_pressure_ratios


  subroutine set_subsonic_area_ratios(ptr, ratios_ptr) bind(C, name = "ffi_cea_set_subsonic_area_ratios")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: ratios_ptr

    type(CEA_Problem), pointer:: this
    real(c_double), pointer:: ratio_list(:)

    call c_f_pointer(ptr, this)
    call c_f_pointer(ratios_ptr%addr, ratio_list, [ratios_ptr%size])

    if (associated(this)) then
       call this%set_subsonic_area_ratios(ratio_list)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_subsonic_area_ratios)'
#endif
    end if

    return
  end subroutine set_subsonic_area_ratios


  subroutine set_supersonic_area_ratios(ptr, ratios_ptr) bind(C, name = "ffi_cea_set_supersonic_area_ratios")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: ratios_ptr

    type(CEA_Problem), pointer:: this
    real(c_double), pointer:: ratio_list(:)

    call c_f_pointer(ptr, this)
    call c_f_pointer(ratios_ptr%addr, ratio_list, [ratios_ptr%size])

    if (associated(this)) then
       call this%set_supersonic_area_ratios(ratio_list)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_supersonic_area_ratios)'
#endif
    end if

    return
  end subroutine set_supersonic_area_ratios


  subroutine set_finite_area_combustor(ptr, contraction_ratio, mass_flow_ratio) bind(C, name = "ffi_cea_set_finite_area_combustor")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    real(c_double), intent(in), optional:: contraction_ratio
    real(c_double), intent(in), optional:: mass_flow_ratio

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    if (associated(this)) then
       call this%set_finite_area_combustor(contraction_ratio, mass_flow_ratio)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_finite_area_combustor)'
#endif
    end if

    return
  end subroutine set_finite_area_combustor


  subroutine add_reactant(ptr, type, name, ratio, T, rho, T_unit, rho_unit) bind(C, name = "ffi_cea_add_reactant")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_String), intent(in):: type
    type(FFI_C_String), intent(in):: name
    real(c_double), intent(in), optional:: ratio
    real(c_double), intent(in), optional:: T
    real(c_double), intent(in), optional:: rho
    type(FFI_C_String), intent(in), optional:: T_unit
    type(FFI_C_String), intent(in), optional:: rho_unit

    character(len = :, kind = c_char), allocatable:: type_fstr
    character(len = :, kind = c_char), allocatable:: name_fstr
    character(len = :, kind = c_char), allocatable:: T_unit_fstr
    character(len = :, kind = c_char), allocatable:: rho_unit_fstr

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    type_fstr = get_string(type)
    name_fstr = get_string(name)
    if (present(T_unit)) T_unit_fstr = get_string(T_unit)
    if (present(rho_unit)) rho_unit_fstr = get_string(rho_unit)

    if (associated(this)) then
       call this%add_reactant(type_fstr, name_fstr, ratio, T, rho, T_unit_fstr, rho_unit_fstr)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (add_reactant)'
#endif
    end if

    return
  end subroutine add_reactant


  subroutine insert_species(ptr, species_ptr, char_length_array, array_size) bind(C, name = "ffi_cea_insert_species")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    integer(c_size_t), intent(in):: array_size
    integer(c_size_t), dimension(array_size), intent(in):: char_length_array
    type(c_ptr), dimension(array_size), intent(in):: species_ptr

    character(len = :, kind = c_char), dimension(:), allocatable:: species

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    call get_string_array(species_ptr, char_length_array, array_size, species)

    call this%insert_species(species)

    deallocate(species)

    return
  end subroutine insert_species


  subroutine set_omit_species(ptr, species_ptr, char_length_array, array_size) bind(C, name = "ffi_cea_set_omit_species")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    integer(c_size_t), intent(in):: array_size
    integer(c_size_t), dimension(array_size), intent(in):: char_length_array
    type(c_ptr), dimension(array_size), intent(in):: species_ptr

    character(len = :, kind = c_char), dimension(:), allocatable:: species

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    call get_string_array(species_ptr, char_length_array, array_size, species)

    call this%set_omit_species(species)

    deallocate(species)

    return
  end subroutine set_omit_species


  subroutine set_only_species(ptr, species_ptr, char_length_array, array_size) bind(C, name = "ffi_cea_set_only_species")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    integer(c_size_t), intent(in):: array_size
    integer(c_size_t), dimension(array_size), intent(in):: char_length_array
    type(c_ptr), dimension(array_size), intent(in):: species_ptr

    character(len = :, kind = c_char), dimension(:), allocatable:: species

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    call get_string_array(species_ptr, char_length_array, array_size, species)

    call this%set_only_species(species)

    deallocate(species)

    return
  end subroutine set_only_species


  subroutine set_legacy_mode(ptr, legacy_mode) bind(C, name = "ffi_cea_set_legacy_mode")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    logical(c_bool), intent(in), optional:: legacy_mode

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    if (associated(this)) then
       if (present(legacy_mode)) then
          call this%set_legacy_mode(logical(legacy_mode))
       else
          call this%set_legacy_mode()
       end if
    end if

    return
  end subroutine set_legacy_mode


  subroutine run(ptr, out_filename) bind(C, name = "ffi_cea_run")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_String), intent(in), optional:: out_filename

    character(len = :, kind = c_char), allocatable:: out_filename_fstr

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    if (present(out_filename)) out_filename_fstr = get_string(out_filename)

    if (associated(this)) then
       call this%run(out_filename_fstr)
    end if

    return
  end subroutine run


!!$  function ffi_read_legacy_input(filename) result(array) bind(C, name = "ffi_cea_read_legacy_input")
!!$    use cea, only: CEA_Problem
!!$    use mod_types, only: init_case, IOINP
!!$    use mod_legacy_io, only: count_cases, read_legacy_case
!!$    use mod_io, only: set_library_paths
!!$
!!$    type(FFI_C_Ptr_Array):: array
!!$    type(FFI_C_String), intent(in):: filename
!!$    character(len = :, kind = c_char), allocatable:: filename_fstr
!!$
!!$    integer:: icase, num_cases
!!$    type(CEA_Problem), pointer:: cea, cea_prev
!!$    type(c_ptr), pointer:: ptr(:)
!!$    logical:: file_exists
!!$
!!$    filename_fstr = get_string(filename)
!!$
!!$    inquire(file = filename_fstr, exist = file_exists)
!!$    if (.not. file_exists) then
!!$       print *, filename_fstr, ' DOES NOT EXIST'
!!$       error stop
!!$    end if
!!$
!!$    call count_cases(filename_fstr, num_cases)
!!$
!!$    allocate(ptr(num_cases))
!!$
!!$    open(newunit = IOINP, file = filename_fstr, status = 'old', form = 'formatted', action = 'read')
!!$
!!$    do icase = 1, num_cases
!!$       allocate(cea)
!!$
!!$       call set_library_paths(cea)
!!$
!!$       call init_case(cea)
!!$
!!$       if (icase > 1 .and. associated(cea_prev)) then
!!$          call read_legacy_case(cea, cea_prev)
!!$       else
!!$          call read_legacy_case(cea)
!!$       end if
!!$
!!$       ptr(icase) = c_loc(cea)
!!$
!!$       cea_prev => cea
!!$       nullify(cea)
!!$    end do
!!$
!!$    nullify(cea_prev)
!!$
!!$    close(IOINP)
!!$
!!$    array%addr = c_loc(ptr(1))
!!$    array%size = num_cases
!!$
!!$    return
!!$  end function ffi_read_legacy_input
!!$
!!$
!!$  subroutine ffi_run_all_cases(array, out_filename, plt_filename) bind(C, name = "ffi_cea_run_all_cases")
!!$    use cea, only: CEA_Problem
!!$    use mod_cea, only: run_all_cases
!!$
!!$    type(FFI_C_Ptr_Array):: array
!!$    type(FFI_C_String), intent(in), optional:: out_filename
!!$    type(FFI_C_String), intent(in), optional:: plt_filename
!!$
!!$    character(len = :, kind = c_char), allocatable:: out_filename_fstr
!!$    character(len = :, kind = c_char), allocatable:: plt_filename_fstr
!!$
!!$    type(CEA_Problem), pointer:: cea(:)
!!$
!!$    call c_f_pointer(array%addr, cea, [array%size])
!!$
!!$    if (present(out_filename)) out_filename_fstr = get_string(out_filename)
!!$    if (present(plt_filename)) plt_filename_fstr = get_string(plt_filename)
!!$
!!$    call run_all_cases(cea, out_filename_fstr, plt_filename_fstr)
!!$
!!$    return
!!$  end subroutine ffi_run_all_cases


  subroutine write_debug_output(ptr, filename) bind(C, name = "ffi_cea_write_debug_output")
    use cea, only: CEA_Problem

    type(c_ptr), intent(in):: ptr
    type(FFI_C_String), intent(in):: filename

    type(CEA_Problem), pointer:: this
    call c_f_pointer(ptr, this)

    if (associated(this)) then
       call this%write_debug_output(get_string(filename))
    end if

    return
  end subroutine write_debug_output


  function ffi_sizeof(ptr) result(size) bind(C, name = "ffi_cea_sizeof")
    use cea, only: CEA_Problem

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
    type(FFI_C_String), intent(in):: cstring
    character(len = :, kind = c_char), allocatable:: fstring
    character(len = cstring%size, kind = c_char), pointer:: fptr

    call c_f_pointer(cstring%addr, fptr)

    fstring = fptr

    nullify(fptr)

    return
  end function get_string


  subroutine get_string_array(cstr_ptr, char_length_array, array_size, fstr_array)
    integer(c_size_t), intent(in):: array_size
    integer(c_size_t), dimension(array_size), intent(in):: char_length_array
    type(c_ptr), dimension(array_size), intent(in):: cstr_ptr
    character(len = :, kind = c_char), dimension(:), allocatable, intent(out):: fstr_array

    integer(c_size_t):: i, j, max_length
    character(len = 1, kind = c_char), dimension(:), pointer:: cstr

    max_length = maxval(char_length_array)

    allocate(character(len = max_length, kind = c_char):: fstr_array(array_size))
    fstr_array(:) = ''

    do i = 1, array_size
       call c_f_pointer(cstr_ptr(i), cstr, [char_length_array(i)])

       do concurrent (j = 1:char_length_array(i))
          fstr_array(i)(j:j) = cstr(j)
       end do
    end do

    return
  end subroutine get_string_array

end module mod_ffi
