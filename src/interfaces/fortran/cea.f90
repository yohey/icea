module cea
  use mod_types, only: CEA_Core_Problem
  implicit none
  private

  type, extends(CEA_Core_Problem), public:: CEA_Problem
   contains
     procedure, public, pass:: set_problem
     procedure, public, pass:: set_output_options
     procedure, public, pass:: set_chamber_pressures
     procedure, public, pass:: set_mixture_ratios
     procedure, public, pass:: set_pressure_ratios
     procedure, public, pass:: set_subsonic_area_ratios
     procedure, public, pass:: set_supersonic_area_ratios
     procedure, public, pass:: set_finite_area_combustor
     procedure, public, pass:: add_reactant
     procedure, public, pass:: insert_species
     procedure, public, pass:: set_omit_species
     procedure, public, pass:: set_only_species
     procedure, public, pass:: run
     procedure, public, pass:: write_debug_output
     final:: del_problem
  end type CEA_Problem

contains

  subroutine del_problem(this)
    type(CEA_Problem), intent(inout):: this

    write(0, *) '[DEBUG] CEA_Problem (cea.f90) destructor is called: ', trim(this%Case)

    return
  end subroutine del_problem

  subroutine set_problem(this, mode, name, mole_ratios, equilibrium, ions, &
       frozen, frozen_at_throat, thermo_lib, trans_lib)
    use mod_constants, only: stderr
    use mod_types, only: init_case
    use mod_io, only: set_library_paths

    class(CEA_Problem), intent(inout):: this
    character(*), intent(in):: mode
    character(*), intent(in), optional:: name
    logical, intent(in), optional:: mole_ratios
    logical, intent(in), optional:: equilibrium
    logical, intent(in), optional:: ions
    logical, intent(in), optional:: frozen
    logical, intent(in), optional:: frozen_at_throat
    character(*), intent(in), optional:: thermo_lib
    character(*), intent(in), optional:: trans_lib

    call init_case(this)

    select case (mode)
    case ('rocket')
       this%Rkt = .true.
       this%Sp = .true.
       this%Tp = .false.
       this%Hp = .false.
       this%Detn = .false.
       this%Shock = .false.

       this%Nt = 1

       this%Nfz = 1

       if (present(frozen)) then
          this%Froz = frozen

          if (present(frozen_at_throat) .and. frozen_at_throat) then
             if (this%Fac) then
                this%Nfz = 3
             else
                this%Nfz = 2
             end if
          end if
       end if

    case ('thermp')
       write(stderr, *) '[ERROR] Not implemented yet.'
       return

    case ('deton')
       write(stderr, *) '[ERROR] Not implemented yet.'
       return

    case ('shock')
       write(stderr, *) '[ERROR] Not implemented yet.'
       return

    case default
       write(stderr, *) '[ERROR] Invalid mode: ', trim(mode)
       return

    end select

    if (present(name)) then
       this%Case = trim(name)
    end if

    if (present(equilibrium)) then
       this%Eql_in = equilibrium
    end if

    if (present(ions)) then
       this%Ions = ions
    end if

    if (present(mole_ratios)) then
       this%Moles = mole_ratios
    end if

    call set_library_paths(this, thermo_lib, trans_lib)

    return
  end subroutine set_problem

  subroutine set_output_options(this, SI, debug_points, mass_fractions, short, trace_tol, transport)
    use mod_types, only: maxMix

    class(CEA_Problem), intent(inout):: this
    logical, intent(in), optional:: SI
    integer, intent(in), optional:: debug_points(:)
    logical, intent(in), optional:: mass_fractions
    logical, intent(in), optional:: short
    real(8), intent(in), optional:: trace_tol
    logical, intent(in), optional:: transport
    integer:: i, j

    if (present(SI)) then
       this%SIunit = SI
    end if

    if (present(debug_points)) then
       do concurrent (i = 1:maxMix, j = 1:size(debug_points), debug_points(j) <= this%max_points)
          this%points(i, debug_points(j))%Debug = .true.
       end do
    end if

    if (present(mass_fractions)) then
       this%Massf = mass_fractions
    end if

    if (present(short)) then
       this%Short = short
    end if

    if (present(trace_tol)) then
       this%Trace = trace_tol
    end if

    if (present(transport)) then
       this%Trnspt = transport
    end if

    return
  end subroutine set_output_options

  subroutine set_chamber_pressures(this, pressure_list, unit)
    use mod_constants, only: stderr, g0, lb_to_kg, in_to_m

    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: pressure_list(:)
    character(*), intent(in), optional:: unit
    real(8), parameter:: legacy_psi_factor = 1.01325d0 / 14.696006d0
    real(8):: factor = 1

    if (present(unit)) then
       if (trim(unit) == 'psi') then
          factor = 1d-5 * g0 * lb_to_kg / in_to_m**2
       else if (trim(unit) == 'legacy-psi') then
          factor = legacy_psi_factor
       else
          write(stderr, *) '[ERROR] Unsupported unit: ', trim(unit)
          return
       end if
    end if

    this%Np = size(pressure_list)
    this%P(1:this%Np) = pressure_list * factor

    return
  end subroutine set_chamber_pressures

  subroutine set_mixture_ratios(this, ratio_list, type)
    use mod_constants, only: stderr

    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: ratio_list(:)
    character(*), intent(in), optional:: type
    integer:: iof

    this%Nof = size(ratio_list)
    this%Oxf(1:this%Nof) = ratio_list

    if (present(type)) then
       if (trim(type) == '%fuel') then
          do concurrent (iof = 1:this%Nof, this%Oxf(iof) > 0)
             this%Oxf(iof) = (100 - this%Oxf(iof)) / this%Oxf(iof)
          end do
       else
          write(stderr, *) '[ERROR] Unsupported type: ', trim(type)
          return
       end if
    end if

    return
  end subroutine set_mixture_ratios

  subroutine set_pressure_ratios(this, ratio_list)
    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: ratio_list(:)

    this%Npp_in = size(ratio_list)
    this%Pcp(1:this%Npp_in) = ratio_list

    return
  end subroutine set_pressure_ratios

  subroutine set_subsonic_area_ratios(this, ratio_list)
    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: ratio_list(:)

    this%Nsub_in = size(ratio_list)
    this%Subar(1:this%Nsub_in) = ratio_list

    return
  end subroutine set_subsonic_area_ratios

  subroutine set_supersonic_area_ratios(this, ratio_list)
    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: ratio_list(:)

    this%Nsup_in = size(ratio_list)
    this%Supar(1:this%Nsup_in) = ratio_list

    return
  end subroutine set_supersonic_area_ratios

  subroutine set_finite_area_combustor(this, contraction_ratio, mass_flow_ratio)
    use mod_constants, only: stderr

    class(CEA_Problem), intent(inout):: this
    real(8), intent(in), optional:: contraction_ratio
    real(8), intent(in), optional:: mass_flow_ratio

    if (present(contraction_ratio)) then
       this%AcAt_in = contraction_ratio
    else if (present(mass_flow_ratio)) then
       this%mdotByAc_in = mass_flow_ratio
    else
       write(stderr, *) '[ERROR] Either contraction_ratio or mass_flow_ratio must be given.'
       return
    end if

    this%Fac = .true.

    if (this%Froz .and. this%Nfz == 2) then
       this%Nfz = 3
    end if

    return
  end subroutine set_finite_area_combustor


  subroutine add_reactant(this, type, name, ratio, T, rho, T_unit, rho_unit)
    use mod_constants, only: stderr

    class(CEA_Problem), intent(inout):: this
    character(*), intent(in):: type
    character(*), intent(in):: name
    real(8), intent(in), optional:: ratio
    real(8), intent(in), optional:: T
    real(8), intent(in), optional:: rho
    character(*), intent(in), optional:: T_unit
    character(*), intent(in), optional:: rho_unit

    real(8):: T_factor = 1
    real(8):: rho_factor = 1

    if (.not. (type(:2) == 'fu' .or. type(:2) == 'ox' .or. type(:2) == 'na')) then
       write(stderr, *) '[ERROR] Invalid reactant type: ', trim(type)
       return
    end if

    this%Nreac = this%Nreac + 1

    this%Fox(this%Nreac) = trim(type)
    this%Rname(this%Nreac) = trim(name)

    if (present(ratio)) then
       this%Pecwt(this%Nreac) = ratio
    else
       this%Pecwt(this%Nreac) = 1
    end if

    if (present(T_unit)) then
       write(stderr, *) '[ERROR] Not implemented yet: T_unit = ', trim(T_unit)
       return
    end if

    if (present(T)) this%Rtemp(this%Nreac) = T * T_factor

    if (present(rho_unit)) then
       write(stderr, *) '[ERROR] Not implemented yet: rho_unit = ', trim(rho_unit)
       return
    end if

    if (present(rho)) this%Dens(this%Nreac) = rho * rho_factor

    this%Energy(this%Nreac) = 'lib'

    return
  end subroutine add_reactant

  subroutine insert_species(this, species)
    class(CEA_Problem), intent(inout):: this
    character(*), intent(in):: species(:)

    this%Nsert = min(20, size(species))
    this%ensert(1:this%Nsert) = species(1:this%Nsert)

    return
  end subroutine insert_species

  subroutine set_omit_species(this, species)
    use mod_types, only: maxNgc

    class(CEA_Problem), intent(inout):: this
    character(*), intent(in):: species(:)

    this%Nomit = min(maxNgc, size(species))
    this%Omit(1:this%Nomit) = species(1:this%Nomit)

    return
  end subroutine set_omit_species

  subroutine set_only_species(this, species)
    use mod_types, only: maxNgc

    class(CEA_Problem), intent(inout):: this
    character(*), intent(in):: species(:)

    this%Nonly_in = min(maxNgc, size(species))
    this%Prod_in(:) = ''
    this%Prod_in(1:this%Nonly_in) = species(1:this%Nonly_in)
    this%Nonly = this%Nonly_in
    this%Prod(:) = this%Prod_in(:)

    return
  end subroutine set_only_species

  subroutine run(this, out_filename)
    use mod_cea, only: run_case

    class(CEA_Problem), intent(inout):: this
    character(*), intent(in), optional:: out_filename

    call run_case(this, out_filename)

    return
  end subroutine run

  subroutine write_debug_output(this, filename)
    use mod_io, only: mod_io_write_debug_output => write_debug_output

    class(CEA_Problem), intent(in):: this
    character(*), intent(in):: filename

    call mod_io_write_debug_output(this, filename)

    return
  end subroutine write_debug_output

end module cea
