module mod_ffi
  use, intrinsic:: iso_c_binding
  implicit none

  private:: get_string, get_string_array

  type, bind(C):: FFI_C_Ptr_Array
     type(c_ptr):: addr
     integer(c_size_t):: size
  end type FFI_C_Ptr_Array

contains

  function new_problem() result(ptr) bind(C, name = "ffi_cea_new_problem")
    use cea, only: CEA_Problem

    type(c_ptr):: ptr
    type(CEA_Problem), pointer:: cea => null()
    allocate(cea)

    ptr = c_loc(cea)

    return
  end function new_problem


  subroutine del_problem(ptr) bind(C, name = "ffi_cea_del_problem")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr

    type(CEA_Problem), pointer:: cea => null()
    call c_f_pointer(ptr, cea)

    if (associated(cea)) then
#ifndef NDEBUG
       write(0, *) '[DEBUG] CEA_Problem (mod_ffi.f90) destructor is called.'
#endif
       deallocate(cea)
    end if

    return
  end subroutine del_problem


  function replica(ptr) result(ptr_replica) bind(C, name = "ffi_cea_replica")
    use cea, only: CEA_Problem

    type(c_ptr):: ptr_replica
    type(c_ptr), value, intent(in):: ptr

    type(CEA_Problem), pointer:: cea => null()
    type(CEA_Problem), pointer:: cea_replica => null()

    call c_f_pointer(ptr, cea)

    if (associated(cea)) then
       allocate(cea_replica)
       call cea%copy_problem(cea_replica)
    end if

    ptr_replica = c_loc(cea_replica)

    return
  end function replica


  subroutine set_problem(ptr, mode, name, mole_ratios, equilibrium, ions, &
       frozen, frozen_at_throat, thermo_lib, trans_lib) bind(C, name = "ffi_cea_set_problem")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    character(len = 1, kind = c_char), intent(in):: mode
    character(len = 1, kind = c_char), intent(in), optional:: name
    logical(c_bool), intent(in), optional:: mole_ratios
    logical(c_bool), intent(in), optional:: equilibrium
    logical(c_bool), intent(in), optional:: ions
    logical(c_bool), intent(in), optional:: frozen
    logical(c_bool), intent(in), optional:: frozen_at_throat
    character(len = 1, kind = c_char), intent(in), optional:: thermo_lib
    character(len = 1, kind = c_char), intent(in), optional:: trans_lib

    character(len = :, kind = c_char), allocatable:: mode_fstr
    character(len = :, kind = c_char), allocatable:: name_fstr
    character(len = :, kind = c_char), allocatable:: thermo_lib_fstr
    character(len = :, kind = c_char), allocatable:: trans_lib_fstr
    logical, pointer:: mole_ratios_ptr
    logical, pointer:: equilibrium_ptr
    logical, pointer:: ions_ptr
    logical, pointer:: frozen_ptr
    logical, pointer:: frozen_at_throat_ptr

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    mode_fstr = get_string(mode)
    if (present(name)) name_fstr = get_string(name)

    if (present(mole_ratios)) then
       allocate(mole_ratios_ptr)
       mole_ratios_ptr = logical(mole_ratios)
    else
       nullify(mole_ratios_ptr)
    end if

    if (present(equilibrium)) then
       allocate(equilibrium_ptr)
       equilibrium_ptr = logical(equilibrium)
    else
       nullify(equilibrium_ptr)
    end if

    if (present(ions)) then
       allocate(ions_ptr)
       ions_ptr = logical(ions)
    else
       nullify(ions_ptr)
    end if

    if (present(frozen)) then
       allocate(frozen_ptr)
       frozen_ptr = logical(frozen)
    else
       nullify(frozen_ptr)
    end if

    if (present(frozen_at_throat)) then
       allocate(frozen_at_throat_ptr)
       frozen_at_throat_ptr = logical(frozen_at_throat)
    else
       nullify(frozen_at_throat_ptr)
    end if

    if (present(thermo_lib)) then
       thermo_lib_fstr = get_string(thermo_lib)
       if (len_trim(thermo_lib_fstr) == 0) deallocate(thermo_lib_fstr)
    end if

    if (present(trans_lib)) then
       trans_lib_fstr = get_string(trans_lib)
       if (len_trim(trans_lib_fstr) == 0) deallocate(trans_lib_fstr)
    end if

    if (associated(this)) then
       call this%set_problem(mode_fstr, name_fstr, mole_ratios_ptr, equilibrium_ptr, ions_ptr, &
            frozen_ptr, frozen_at_throat_ptr, thermo_lib_fstr, trans_lib_fstr)
    end if

    if (allocated(mode_fstr)) deallocate(mode_fstr)
    if (allocated(name_fstr)) deallocate(name_fstr)
    if (allocated(thermo_lib_fstr)) deallocate(thermo_lib_fstr)
    if (allocated(trans_lib_fstr)) deallocate(trans_lib_fstr)
    if (associated(mole_ratios_ptr)) deallocate(mole_ratios_ptr)
    if (associated(equilibrium_ptr)) deallocate(equilibrium_ptr)
    if (associated(ions_ptr)) deallocate(ions_ptr)
    if (associated(frozen_ptr)) deallocate(frozen_ptr)
    if (associated(frozen_at_throat_ptr)) deallocate(frozen_at_throat_ptr)

    return
  end subroutine set_problem


  subroutine set_output_options(ptr, SI, debug_points, mass_fractions, short, trace_tol, transport, &
       plot_ptr, char_length_array_ptr) bind(C, name = "ffi_cea_set_output_options")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    logical(c_bool), intent(in), optional:: SI
    type(FFI_C_Ptr_Array), intent(in), optional:: debug_points
    logical(c_bool), intent(in), optional:: mass_fractions
    logical(c_bool), intent(in), optional:: short
    real(c_double), intent(in), optional:: trace_tol
    logical(c_bool), intent(in), optional:: transport
    type(c_ptr), intent(in), optional:: plot_ptr
    type(FFI_C_Ptr_Array), intent(in), optional:: char_length_array_ptr

    type(CEA_Problem), pointer:: this => null()
    integer(c_size_t), pointer:: debug_points_array(:)
    integer, allocatable:: debug_points_array_fint(:)
    logical, pointer:: SI_ptr
    logical, pointer:: mass_fractions_ptr
    logical, pointer:: short_ptr
    logical, pointer:: transport_ptr
    integer(c_size_t), dimension(:), pointer:: char_length_array
    character(len = :, kind = c_char), dimension(:), allocatable:: plot

    call c_f_pointer(ptr, this)

    if (present(debug_points)) then
       call c_f_pointer(debug_points%addr, debug_points_array, [debug_points%size])
       allocate(debug_points_array_fint(debug_points%size))
       debug_points_array_fint(:) = int(debug_points_array(:))
    end if

    if (present(SI)) then
       allocate(SI_ptr)
       SI_ptr = logical(SI)
    else
       nullify(SI_ptr)
    end if

    if (present(mass_fractions)) then
       allocate(mass_fractions_ptr)
       mass_fractions_ptr = logical(mass_fractions)
    else
       nullify(mass_fractions_ptr)
    end if

    if (present(short)) then
       allocate(short_ptr)
       short_ptr = logical(short)
    else
       nullify(short_ptr)
    end if

    if (present(transport)) then
       allocate(transport_ptr)
       transport_ptr = logical(transport)
    else
       nullify(transport_ptr)
    end if

    if (present(plot_ptr) .and. present(char_length_array_ptr)) then
       call c_f_pointer(char_length_array_ptr%addr, char_length_array, [char_length_array_ptr%size])
       plot = get_string_array(plot_ptr, char_length_array, char_length_array_ptr%size)
    end if

    if (associated(this)) then
       call this%set_output_options(SI_ptr, debug_points_array_fint, mass_fractions_ptr, short_ptr, &
            trace_tol, transport_ptr, plot)
    end if

    if (present(debug_points)) deallocate(debug_points_array_fint)

    if (present(plot_ptr) .and. present(char_length_array_ptr)) deallocate(plot)

    return
  end subroutine set_output_options


  subroutine set_chamber_pressures(ptr, pressures_ptr, unit) bind(C, name = "ffi_cea_set_chamber_pressures")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: pressures_ptr
    character(len = 1, kind = c_char), intent(in), optional:: unit

    type(CEA_Problem), pointer:: this => null()
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

    type(c_ptr), value, intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: ratios_ptr
    character(len = 1, kind = c_char), intent(in), optional:: type

    type(CEA_Problem), pointer:: this => null()
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

    type(c_ptr), value, intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: ratios_ptr

    type(CEA_Problem), pointer:: this => null()
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

    type(c_ptr), value, intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: ratios_ptr

    type(CEA_Problem), pointer:: this => null()
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

    type(c_ptr), value, intent(in):: ptr
    type(FFI_C_Ptr_Array), intent(in):: ratios_ptr

    type(CEA_Problem), pointer:: this => null()
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

    type(c_ptr), value, intent(in):: ptr
    real(c_double), intent(in), optional:: contraction_ratio
    real(c_double), intent(in), optional:: mass_flow_ratio

    type(CEA_Problem), pointer:: this => null()
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


  subroutine add_reactant(ptr, type, name, formula, ratio, T, rho, h, u, T_unit, rho_unit, h_unit, u_unit) &
       bind(C, name = "ffi_cea_add_reactant")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    character(len = 1, kind = c_char), intent(in):: type
    character(len = 1, kind = c_char), intent(in):: name
    character(len = 1, kind = c_char), intent(in), optional:: formula
    real(c_double), intent(in), optional:: ratio
    real(c_double), intent(in), optional:: T
    real(c_double), intent(in), optional:: rho
    real(c_double), intent(in), optional:: h
    real(c_double), intent(in), optional:: u
    character(len = 1, kind = c_char), intent(in), optional:: T_unit
    character(len = 1, kind = c_char), intent(in), optional:: rho_unit
    character(len = 1, kind = c_char), intent(in), optional:: h_unit
    character(len = 1, kind = c_char), intent(in), optional:: u_unit

    character(len = :, kind = c_char), allocatable:: type_fstr
    character(len = :, kind = c_char), allocatable:: name_fstr
    character(len = :, kind = c_char), allocatable:: formula_fstr
    character(len = :, kind = c_char), allocatable:: T_unit_fstr
    character(len = :, kind = c_char), allocatable:: rho_unit_fstr
    character(len = :, kind = c_char), allocatable:: h_unit_fstr
    character(len = :, kind = c_char), allocatable:: u_unit_fstr

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    type_fstr = get_string(type)
    name_fstr = get_string(name)

    if (present(formula)) formula_fstr = get_string(formula)
    if (present(T_unit)) T_unit_fstr = get_string(T_unit)
    if (present(rho_unit)) rho_unit_fstr = get_string(rho_unit)
    if (present(h_unit)) h_unit_fstr = get_string(h_unit)
    if (present(u_unit)) u_unit_fstr = get_string(u_unit)

    if (associated(this)) then
       call this%add_reactant(type_fstr, name_fstr, formula_fstr, ratio, T, rho, h, u, &
                              T_unit_fstr, rho_unit_fstr, h_unit_fstr, u_unit_fstr)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (add_reactant)'
#endif
    end if

    return
  end subroutine add_reactant


  subroutine set_reactant(ptr, index, ratio, T, rho, h, u, T_unit, rho_unit, h_unit, u_unit) &
       bind(C, name = "ffi_cea_set_reactant")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    integer(c_size_t), intent(in):: index
    real(c_double), intent(in), optional:: ratio
    real(c_double), intent(in), optional:: T
    real(c_double), intent(in), optional:: rho
    real(c_double), intent(in), optional:: h
    real(c_double), intent(in), optional:: u
    character(len = 1, kind = c_char), intent(in), optional:: T_unit
    character(len = 1, kind = c_char), intent(in), optional:: rho_unit
    character(len = 1, kind = c_char), intent(in), optional:: h_unit
    character(len = 1, kind = c_char), intent(in), optional:: u_unit

    character(len = :, kind = c_char), allocatable:: T_unit_fstr
    character(len = :, kind = c_char), allocatable:: rho_unit_fstr
    character(len = :, kind = c_char), allocatable:: h_unit_fstr
    character(len = :, kind = c_char), allocatable:: u_unit_fstr

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (present(T_unit)) T_unit_fstr = get_string(T_unit)
    if (present(rho_unit)) rho_unit_fstr = get_string(rho_unit)
    if (present(h_unit)) h_unit_fstr = get_string(h_unit)
    if (present(u_unit)) u_unit_fstr = get_string(u_unit)

    if (associated(this)) then
       call this%set_reactant(int(index), ratio, T, rho, h, u, &
                              T_unit_fstr, rho_unit_fstr, h_unit_fstr, u_unit_fstr)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_reactant)'
#endif
    end if

    return
  end subroutine set_reactant


  subroutine set_insert_species(ptr, species_ptr, char_length_array_ptr) bind(C, name = "ffi_cea_set_insert_species")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    type(c_ptr), intent(in):: species_ptr
    type(FFI_C_Ptr_Array), intent(in):: char_length_array_ptr

    type(CEA_Problem), pointer:: this => null()
    integer(c_size_t), dimension(:), pointer:: char_length_array
    character(len = :, kind = c_char), dimension(:), allocatable:: species

    call c_f_pointer(ptr, this)
    call c_f_pointer(char_length_array_ptr%addr, char_length_array, [char_length_array_ptr%size])

    species = get_string_array(species_ptr, char_length_array, char_length_array_ptr%size)

    call this%set_insert_species(species)

    deallocate(species)

    return
  end subroutine set_insert_species


  subroutine set_omit_species(ptr, species_ptr, char_length_array_ptr) bind(C, name = "ffi_cea_set_omit_species")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    type(c_ptr), intent(in):: species_ptr
    type(FFI_C_Ptr_Array), intent(in):: char_length_array_ptr

    type(CEA_Problem), pointer:: this => null()
    integer(c_size_t), dimension(:), pointer:: char_length_array
    character(len = :, kind = c_char), dimension(:), allocatable:: species

    call c_f_pointer(ptr, this)
    call c_f_pointer(char_length_array_ptr%addr, char_length_array, [char_length_array_ptr%size])

    species = get_string_array(species_ptr, char_length_array, char_length_array_ptr%size)

    call this%set_omit_species(species)

    deallocate(species)

    return
  end subroutine set_omit_species


  subroutine set_only_species(ptr, species_ptr, char_length_array_ptr) bind(C, name = "ffi_cea_set_only_species")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    type(c_ptr), intent(in):: species_ptr
    type(FFI_C_Ptr_Array), intent(in):: char_length_array_ptr

    type(CEA_Problem), pointer:: this => null()
    integer(c_size_t), dimension(:), pointer:: char_length_array
    character(len = :, kind = c_char), dimension(:), allocatable:: species

    call c_f_pointer(ptr, this)
    call c_f_pointer(char_length_array_ptr%addr, char_length_array, [char_length_array_ptr%size])

#ifndef NDEBUG
    if (associated(char_length_array)) then
#endif
       species = get_string_array(species_ptr, char_length_array, char_length_array_ptr%size)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_only_species, char_length_array)'
    end if

    if (associated(this)) then
#endif
       call this%set_only_species(species)
#ifndef NDEBUG
    else
       write(0, *) '[DEBUG] pointer is not associated. (set_only_species, this)'
    end if
#endif

    deallocate(species)

    return
  end subroutine set_only_species


  subroutine set_legacy_mode(ptr, legacy_mode) bind(C, name = "ffi_cea_set_legacy_mode")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    logical(c_bool), intent(in), optional:: legacy_mode

    type(CEA_Problem), pointer:: this => null()
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


  subroutine run(ptr, out_filename, plt_filename) bind(C, name = "ffi_cea_run")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    character(len = 1, kind = c_char), intent(in), optional:: out_filename
    character(len = 1, kind = c_char), intent(in), optional:: plt_filename

    character(len = :, kind = c_char), allocatable:: out_filename_fstr
    character(len = :, kind = c_char), allocatable:: plt_filename_fstr

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (present(out_filename)) out_filename_fstr = get_string(out_filename)
    if (present(plt_filename)) plt_filename_fstr = get_string(plt_filename)

    if (associated(this)) then
       call this%run(out_filename_fstr, plt_filename_fstr)
    end if

    return
  end subroutine run


  function ffi_get_pressure(ptr, iOF, ipt) result(P) bind(C, name = "ffi_cea_get_pressure")
    use cea, only: CEA_Problem

    real(c_double):: P
    type(c_ptr), value, intent(in):: ptr
    integer(c_size_t), intent(in):: iOF
    integer(c_size_t), intent(in):: ipt

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) P = this%get_pressure(int(iOF), int(ipt))

    return
  end function ffi_get_pressure


  function ffi_get_temperature(ptr, iOF, ipt) result(T) bind(C, name = "ffi_cea_get_temperature")
    use cea, only: CEA_Problem

    real(c_double):: T
    type(c_ptr), value, intent(in):: ptr
    integer(c_size_t), intent(in):: iOF
    integer(c_size_t), intent(in):: ipt

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) T = this%get_temperature(int(iOF), int(ipt))

    return
  end function ffi_get_temperature

  function ffi_get_chamber_temperature(ptr) result(T) bind(C, name = "ffi_cea_get_chamber_temperature")
    use cea, only: CEA_Problem

    real(c_double):: T
    type(c_ptr), value, intent(in):: ptr

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) T = this%get_chamber_temperature()

    return
  end function ffi_get_chamber_temperature


  function ffi_get_molecular_weight(ptr, iOF, ipt) result(M) bind(C, name = "ffi_cea_get_molecular_weight")
    use cea, only: CEA_Problem

    real(c_double):: M
    type(c_ptr), value, intent(in):: ptr
    integer(c_size_t), intent(in):: iOF
    integer(c_size_t), intent(in):: ipt

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) M = this%get_molecular_weight(int(iOF), int(ipt))

    return
  end function ffi_get_molecular_weight


  function ffi_get_specific_heat(ptr, iOF, ipt) result(cp) bind(C, name = "ffi_cea_get_specific_heat")
    use cea, only: CEA_Problem

    real(c_double):: cp
    type(c_ptr), value, intent(in):: ptr
    integer(c_size_t), intent(in):: iOF
    integer(c_size_t), intent(in):: ipt

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) cp = this%get_specific_heat(int(iOF), int(ipt))

    return
  end function ffi_get_specific_heat


  function ffi_get_specific_heat_ratio(ptr, iOF, ipt) result(gamma) bind(C, name = "ffi_cea_get_specific_heat_ratio")
    use cea, only: CEA_Problem

    real(c_double):: gamma
    type(c_ptr), value, intent(in):: ptr
    integer(c_size_t), intent(in):: iOF
    integer(c_size_t), intent(in):: ipt

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) gamma = this%get_specific_heat_ratio(int(iOF), int(ipt))

    return
  end function ffi_get_specific_heat_ratio


  function ffi_get_characteristic_velocity(ptr) result(c_star) bind(C, name = "ffi_cea_get_characteristic_velocity")
    use cea, only: CEA_Problem

    real(c_double):: c_star
    type(c_ptr), value, intent(in):: ptr

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) c_star = this%get_characteristic_velocity()

    return
  end function ffi_get_characteristic_velocity


  function ffi_get_specific_impulse(ptr, iOF, ipt, vacuum) result(Isp) bind(C, name = "ffi_cea_get_specific_impulse")
    use cea, only: CEA_Problem

    real(c_double):: Isp
    type(c_ptr), value, intent(in):: ptr
    integer(c_size_t), intent(in):: iOF
    integer(c_size_t), intent(in):: ipt
    logical(c_bool), intent(in), optional:: vacuum

    logical, pointer:: vacuum_ptr

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (present(vacuum)) then
       allocate(vacuum_ptr)
       vacuum_ptr = logical(vacuum)
    else
       nullify(vacuum_ptr)
    end if

    if (associated(this)) Isp = this%get_specific_impulse(int(iOF), int(ipt), vacuum_ptr)

    if (associated(vacuum_ptr)) deallocate(vacuum_ptr)

    return
  end function ffi_get_specific_impulse


  function ffi_read_legacy_input(filename) result(array) bind(C, name = "ffi_cea_read_legacy_input")
    use cea, only: CEA_Problem
    use mod_types, only: init_case, IOINP
    use mod_legacy_io, only: count_cases, read_legacy_case
    use mod_io, only: set_library_paths

    type(FFI_C_Ptr_Array):: array
    character(len = 1, kind = c_char), intent(in):: filename
    character(len = :, kind = c_char), allocatable:: filename_fstr

    integer:: icase, num_cases
    type(CEA_Problem), pointer:: cea => null()
    type(CEA_Problem), pointer:: cea_prev => null()
    type(c_ptr), pointer:: ptr(:)
    logical:: file_exists

    filename_fstr = get_string(filename)

    inquire(file = filename_fstr, exist = file_exists)
    if (.not. file_exists) then
       print *, filename_fstr, ' DOES NOT EXIST'
       error stop
    end if

    call count_cases(filename_fstr, num_cases)

    allocate(ptr(num_cases))

    open(newunit = IOINP, file = filename_fstr, status = 'old', form = 'formatted', action = 'read')

    deallocate(filename_fstr)

    do icase = 1, num_cases
       allocate(cea)

       call set_library_paths(cea)

       call init_case(cea)

       if (icase > 1 .and. associated(cea_prev)) then
          call read_legacy_case(cea, cea_prev)
       else
          call read_legacy_case(cea)
       end if

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


  subroutine ffi_deallocate_array_ptr(array) bind(C, name = "ffi_cea_deallocate_array_ptr")
    use cea, only: CEA_Problem

    type(FFI_C_Ptr_Array), intent(in), target:: array

    type(c_ptr), pointer:: ptr(:)

    call c_f_pointer(array%addr, ptr, [array%size])

    if (associated(ptr)) deallocate(ptr)

    return
  end subroutine ffi_deallocate_array_ptr


  subroutine ffi_run_all_cases(array, out_filename, plt_filename) bind(C, name = "ffi_cea_run_all_cases")
    use cea, only: CEA_Problem
    use mod_cea, only: run_all_cases

    type(FFI_C_Ptr_Array), intent(in):: array
    character(len = 1, kind = c_char), intent(in), optional:: out_filename
    character(len = 1, kind = c_char), intent(in), optional:: plt_filename

    character(len = :, kind = c_char), allocatable:: out_filename_fstr
    character(len = :, kind = c_char), allocatable:: plt_filename_fstr

    type(CEA_Problem), pointer:: cea(:)

    call c_f_pointer(array%addr, cea, [array%size])

    if (present(out_filename)) out_filename_fstr = get_string(out_filename)
    if (present(plt_filename)) plt_filename_fstr = get_string(plt_filename)

    call run_all_cases(cea, out_filename_fstr, plt_filename_fstr)

    if (allocated(out_filename_fstr)) deallocate(out_filename_fstr)
    if (allocated(plt_filename_fstr)) deallocate(plt_filename_fstr)

    return
  end subroutine ffi_run_all_cases


  subroutine ffi_calc_frozen_exhaust(ptr, T, cp, gamma, mu, k, Pr) bind(C, name = "ffi_cea_calc_frozen_exhaust")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    real(c_double), intent(in):: T
    real(c_double), intent(out), optional:: cp
    real(c_double), intent(out), optional:: gamma
    real(c_double), intent(out), optional:: mu
    real(c_double), intent(out), optional:: k
    real(c_double), intent(out), optional:: Pr

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) then
       call this%calc_frozen_exhaust(T, cp, gamma, mu, k, Pr)
    end if

    return
  end subroutine ffi_calc_frozen_exhaust


  subroutine ffi_get_thermo_reference_properties(ptr, name, M, T_ref, h0_ref) &
       bind(C, name = "ffi_cea_get_thermo_reference_properties")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    character(len = 1, kind = c_char), intent(in):: name
    real(c_double), intent(out):: M
    real(c_double), intent(out):: T_ref
    real(c_double), intent(out):: h0_ref

    type(CEA_Problem), pointer:: this => null()
    character(len = :, kind = c_char), allocatable:: name_fstr

    call c_f_pointer(ptr, this)

    name_fstr = get_string(name)

    M = 0
    T_ref = 0
    h0_ref = 0

    if (associated(this)) then
       call this%get_thermo_reference_properties(name_fstr, M, T_ref, h0_ref)
    end if

    return
  end subroutine ffi_get_thermo_reference_properties


  subroutine write_debug_output(ptr, filename) bind(C, name = "ffi_cea_write_debug_output")
    use cea, only: CEA_Problem

    type(c_ptr), value, intent(in):: ptr
    character(len = 1, kind = c_char), intent(in):: filename

    type(CEA_Problem), pointer:: this => null()
    call c_f_pointer(ptr, this)

    if (associated(this)) then
       call this%write_debug_output(get_string(filename))
    end if

    return
  end subroutine write_debug_output


  function ffi_sizeof(ptr) result(size) bind(C, name = "ffi_cea_sizeof")
    use cea, only: CEA_Problem

    integer(c_size_t):: size
    type(c_ptr), value, intent(in):: ptr

    type(CEA_Problem), pointer:: cea => null()
    call c_f_pointer(ptr, cea)

    if (associated(cea)) then
       size = storage_size(cea) / 8  ! bit -> byte
    end if

    return
  end function ffi_sizeof


  function get_string(cstring) result(fstring)
    character(len = 1, kind = c_char), intent(in):: cstring(*)
    character(len = :, kind = c_char), allocatable:: fstring

    integer:: i, length

    length = 1
    do while (cstring(length) /= c_null_char)
       length = length + 1
    end do
    length = length - 1

    allocate(character(len = length, kind = c_char):: fstring)

    do concurrent (i = 1:length)
       fstring(i:i) = cstring(i)
    end do

    return
  end function get_string


  function get_string_array(cstr_ptr, char_length_array, array_size) result(fstr_array)
    character(len = :, kind = c_char), dimension(:), allocatable:: fstr_array

    integer(c_size_t), intent(in):: array_size
    integer(c_size_t), dimension(array_size), intent(in):: char_length_array
    type(c_ptr), dimension(array_size), intent(in):: cstr_ptr

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
  end function get_string_array

end module mod_ffi
