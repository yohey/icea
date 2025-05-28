module cea
  use mod_constants, only: stderr
  use mod_types, only: CEA_Core_Problem
  implicit none
  private

  type, extends(CEA_Core_Problem), public:: CEA_Problem
     private
   contains
     final:: del_problem
     procedure, public, pass:: copy_problem
     procedure, public, pass:: set_problem
     procedure, public, pass:: set_output_options
     procedure, public, pass:: set_chamber_pressures
     procedure, public, pass:: set_chamber_temperatures
     procedure, public, pass:: set_chamber_densities
     procedure, public, pass:: set_chamber_internal_energy
     procedure, public, pass:: set_mixture_ratios
     procedure, public, pass:: set_pressure_ratios
     procedure, public, pass:: set_subsonic_area_ratios
     procedure, public, pass:: set_supersonic_area_ratios
     procedure, public, pass:: set_finite_area_combustor
     procedure, public, pass:: add_reactant
     procedure, public, pass:: set_reactant
     procedure, public, pass:: set_insert_species
     procedure, public, pass:: set_omit_species
     procedure, public, pass:: set_only_species
     procedure, public, pass:: set_legacy_mode
     procedure, public, pass:: run
     procedure, public, pass:: get_pressure
     procedure, public, pass:: get_temperature
     procedure, public, pass:: get_chamber_temperature
     procedure, public, pass:: get_molecular_weight
     procedure, public, pass:: get_specific_heat
     procedure, public, pass:: get_specific_heat_ratio
     procedure, public, pass:: get_characteristic_velocity
     procedure, public, pass:: get_specific_impulse
     procedure, public, pass:: calc_frozen_exhaust
     procedure, public, pass:: get_thermo_reference_properties
     procedure, public, pass:: write_debug_output
  end type CEA_Problem

contains

  subroutine del_problem(this)
    type(CEA_Problem), intent(inout):: this

#ifndef NDEBUG
    write(stderr, *) '[DEBUG] CEA_Problem (cea.f90) destructor is called: ', trim(this%Case)
#endif

    return
  end subroutine del_problem


  subroutine copy_problem(this, cea_copy)
    class(CEA_Problem), intent(in):: this
    type(CEA_Problem), intent(out):: cea_copy

    cea_copy = this

    nullify(cea_copy%thermo_properties)
    nullify(cea_copy%points)

    if (associated(this%thermo_properties)) then
       allocate(cea_copy%thermo_properties(size(this%thermo_properties)))
       cea_copy%thermo_properties = this%thermo_properties
    end if

    if (associated(this%points)) then
       allocate(cea_copy%points(size(this%points, 1), size(this%points, 2)))
       cea_copy%points = this%points
    end if

    return
  end subroutine copy_problem


  subroutine set_problem(this, mode, name, mole_ratios, equilibrium, ions, &
       frozen, frozen_at_throat, thermo_lib, trans_lib)
    use mod_types, only: init_case
    use mod_io, only: set_library_paths, read_thermo_lib

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

    case ('tp', 'pt')
       this%Tp = .true.

    case ('hp', 'ph')
       this%Hp = .true.

    case ('uv', 'vu')
       this%Hp = .true.
       this%Vol = .true.

    case ('tv', 'vt')
       this%Tp = .true.
       this%Vol = .true.

    case ('thermp')
       write(stderr, *) '[ERROR] Not implemented yet.'
       return

    case ('deton')
       this%Detn = .true.

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

    if (.not. associated(this%thermo_properties)) call read_thermo_lib(this)

    return
  end subroutine set_problem


  subroutine set_output_options(this, SI, debug_points, mass_fractions, short, trace_tol, transport, plot)
    use mod_types, only: maxMix

    class(CEA_Problem), intent(inout):: this
    logical, intent(in), optional:: SI
    integer, intent(in), optional:: debug_points(:)
    logical, intent(in), optional:: mass_fractions
    logical, intent(in), optional:: short
    real(8), intent(in), optional:: trace_tol
    logical, intent(in), optional:: transport
    character(*), intent(in), optional:: plot(:)
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

    if (present(plot)) then
       this%Nplt = min(20, size(plot))
       this%Pltvar(:) = ''
       this%Pltvar(1:this%Nplt) = plot(1:this%Nplt)
    end if

    return
  end subroutine set_output_options


  subroutine set_chamber_pressures(this, pressure_list, unit)
    use mod_constants, only: g0, lb_to_kg, in_to_m

    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: pressure_list(:)
    character(*), intent(in), optional:: unit
    real(8), parameter:: legacy_psi_factor = 1.01325d0 / 14.696006d0
    real(8):: factor

    factor = 10 ! default: MPa

    if (present(unit)) then
       select case (trim(unit))
       case ('MPa')
          factor = 10
       case ('kPa')
          factor = 1d-2
       case ('Pa')
          factor = 1d-5
       case ('atm')
          factor = 1.01325d0
       case ('bar')
          factor = 1
       case ('psi')
          factor = 1d-5 * g0 * lb_to_kg / in_to_m**2
       case ('legacy-psi')
          factor = legacy_psi_factor
       case default
          write(stderr, *) '[ERROR] Unsupported unit: ', trim(unit)
          return
       end select
    end if

#ifndef NDEBUG
    if (size(pressure_list) == 0) write(stderr, *) '[ERROR] empty pressure list (set_chamber_pressures).'
#endif

    this%Np = size(pressure_list)
    this%P(1:this%Np) = pressure_list * factor

    return
  end subroutine set_chamber_pressures


  subroutine set_chamber_temperatures(this, temperature_list, unit)
    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: temperature_list(:)
    character(*), intent(in), optional:: unit
    real(8):: factor

    factor = 1 ! default: K

    if (present(unit)) then
       select case (trim(unit))
       case ('K')
          factor = 1
       case default
          write(stderr, *) '[ERROR] Unsupported unit: ', trim(unit)
          return
       end select
    end if

    this%Nt = size(temperature_list)
    this%T(1:this%Nt) = temperature_list * factor

    return
  end subroutine set_chamber_temperatures


  subroutine set_chamber_densities(this, density_list, unit)
    use mod_types, only: maxPv

    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: density_list(:)
    character(*), intent(in), optional:: unit
    real(8):: factor

    factor = 1d5 ! default: kg/m^3

    if (present(unit)) then
       select case (trim(unit))
       case ('kg/m^3')
          factor = 1d5
       case ('g/cc')
          factor = 1d2
       case default
          write(stderr, *) '[ERROR] Unsupported unit: ', trim(unit)
          return
       end select
    end if

    this%Np = min(maxPv, size(density_list))
    this%V(1:this%Np) = factor / density_list

    return
  end subroutine set_chamber_densities


  subroutine set_chamber_internal_energy(this, uByR)
    use mod_types, only: maxPv

    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: uByR

    this%Hsub0 = uByR

    return
  end subroutine set_chamber_internal_energy


  subroutine set_mixture_ratios(this, ratio_list, type)
    class(CEA_Problem), intent(inout):: this
    real(8), intent(in):: ratio_list(:)
    character(*), intent(in), optional:: type
    integer:: iof

    this%Nof = size(ratio_list)
    this%eqrats_in = .false.

    if (allocated(this%Oxf_in)) then
       if (size(this%Oxf_in) /= this%Nof) then
          deallocate(this%Oxf_in)
          allocate(this%Oxf_in(this%Nof))
       end if
    else
       allocate(this%Oxf_in(this%Nof))
    end if

    this%Oxf_in(1:this%Nof) = ratio_list

    if (present(type)) then
       if (trim(type) == 'eq.ratio') then
          this%eqrats_in = .true.
       else if (trim(type) == '%fuel') then
          do concurrent (iof = 1:this%Nof, this%Oxf_in(iof) > 0)
             this%Oxf_in(iof) = (100 - this%Oxf_in(iof)) / this%Oxf_in(iof)
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


  subroutine add_reactant(this, type, name, formula, ratio, T, rho, h, u, T_unit, rho_unit, h_unit, u_unit)
    use mod_constants, only: R0, cal_to_J

    class(CEA_Problem), intent(inout):: this
    character(*), intent(in):: type
    character(*), intent(in):: name
    character(*), intent(in), optional:: formula
    real(8), intent(in), optional:: ratio
    real(8), intent(in), optional:: T
    real(8), intent(in), optional:: rho
    real(8), intent(in), optional:: h
    real(8), intent(in), optional:: u
    character(*), intent(in), optional:: T_unit
    character(*), intent(in), optional:: rho_unit
    character(*), intent(in), optional:: h_unit
    character(*), intent(in), optional:: u_unit

    integer:: istart, iend, iatom
    logical:: atom_found

    real(8):: T_factor
    real(8):: rho_factor
    real(8):: h_factor
    real(8):: u_factor

    T_factor = 1 ! default: K
    rho_factor = 1d-3 ! default: kg/m^3
    h_factor = 1000 / R0 ! default: J/mol
    u_factor = 1000 / R0 ! default: J/mol

    if (.not. (type(:2) == 'fu' .or. type(:2) == 'ox' .or. type(:2) == 'na')) then
       write(stderr, *) '[ERROR] Invalid reactant type: ', trim(type)
       return
    end if

    this%Nreac_in = this%Nreac_in + 1
    this%Nreac = this%Nreac_in

    this%Energy(this%Nreac) = 'lib'

    this%Fox(this%Nreac) = trim(type)
    this%Rname(this%Nreac) = trim(name)

    if (present(formula)) then
       this%Ratom(this%Nreac, :) = ''
       this%Rnum(this%Nreac, :) = 0
       istart = 1
       iatom = 1
       atom_found = .false.

       do iend = 1, len_trim(formula)
          if (formula(iend:iend) == ' ' .or. iend == len_trim(formula)) then

             if (atom_found) then
                read(formula(istart:iend), *) this%Rnum(this%Nreac, iatom)
                iatom = iatom + 1
                atom_found = .false.
             else
                this%Ratom(this%Nreac, iatom) = formula(istart:iend)
                atom_found = .true.
             end if

             istart = iend + 1
          end if
       end do

       this%Nfla(this%Nreac) = iatom - 1
    end if

    if (present(ratio)) then
       this%Pecwt(this%Nreac) = ratio
    else
       this%Pecwt(this%Nreac) = 1
    end if

    if (present(T)) then
       if (present(T_unit)) then
          select case (trim(T_unit))
          case ('K')
             T_factor = 1
          case default
             write(stderr, *) '[ERROR] Unsupported unit: ', trim(T_unit)
             return
          end select
       end if

       this%Rtemp(this%Nreac) = T * T_factor
    end if

    if (present(rho)) then
       if (present(rho_unit)) then
          select case (trim(rho_unit))
          case ('kg/m^3', 'kg/m3')
             rho_factor = 1d-3
          case ('g/cc', 'g/cm^3', 'g/cm3')
             rho_factor = 1
          case default
             write(stderr, *) '[ERROR] Unsupported unit: ', trim(rho_unit)
             return
          end select
       end if

       this%Dens(this%Nreac) = rho * rho_factor
    end if

    if (present(h)) then
       if (.not. present(formula)) then
          write(stderr, *) '[ERROR] Exploded formula is required when h is given.'
          return
       end if

       if (present(h_unit)) then
          select case (trim(h_unit))
          case ('J/mol')
             h_factor = 1d3 / R0
          case ('kJ/mol')
             h_factor = 1d6 / R0
          case ('cal/mol')
             h_factor = 1d3 * cal_to_J / R0
          case default
             write(stderr, *) '[ERROR] Unsupported unit: ', trim(h_unit)
             return
          end select
       end if

       this%Energy(this%Nreac) = 'h,j'
       this%Enth(this%Nreac) = h * h_factor
    end if

    if (present(u)) then
       if (.not. present(formula)) then
          write(stderr, *) '[ERROR] Exploded formula is required when u is given.'
          return
       end if

       if (present(u_unit)) then
          select case (trim(u_unit))
          case ('J/mol')
             u_factor = 1d3 / R0
          case ('kJ/mol')
             h_factor = 1d6 / R0
          case ('cal/mol')
             h_factor = 1d3 * cal_to_J / R0
          case default
             write(stderr, *) '[ERROR] Unsupported unit: ', trim(u_unit)
             return
          end select
       end if

       this%Energy(this%Nreac) = 'u,j'
       this%Enth(this%Nreac) = u * u_factor
    end if

    return
  end subroutine add_reactant


  subroutine set_reactant(this, index, ratio, T, rho, h, u, T_unit, rho_unit, h_unit, u_unit)
    use mod_constants, only: R0
    use mod_legacy_io, only: REACT

    class(CEA_Problem), intent(inout):: this
    integer, intent(in):: index
    real(8), intent(in), optional:: ratio
    real(8), intent(in), optional:: T
    real(8), intent(in), optional:: rho
    real(8), intent(in), optional:: h
    real(8), intent(in), optional:: u
    character(*), intent(in), optional:: T_unit
    character(*), intent(in), optional:: rho_unit
    character(*), intent(in), optional:: h_unit
    character(*), intent(in), optional:: u_unit

    real(8):: T_factor
    real(8):: rho_factor
    real(8):: h_factor
    real(8):: u_factor

    T_factor = 1 ! default: K
    rho_factor = 1d-3 ! default: kg/m^3
    h_factor = 1000 / R0 ! default: J/mol
    u_factor = 1000 / R0 ! default: J/mol

    if (index > this%Nreac) then
       write(stderr, *) '[ERROR] Invalid reactant index: ', index, ' (< ', this%Nreac, ')'
       return
    end if

    this%Energy(index) = 'lib'
    this%Nreac = this%Nreac_in

    if (present(ratio)) then
       this%Pecwt(index) = ratio
    else
       this%Pecwt(index) = 1
    end if

    if (present(T)) then
       if (present(T_unit)) then
          select case (trim(T_unit))
          case ('K')
             T_factor = 1
          case default
             write(stderr, *) '[ERROR] Unsupported unit: ', trim(T_unit)
             return
          end select
       end if

       this%Rtemp(index) = T * T_factor
    end if

    if (present(rho)) then
       if (present(rho_unit)) then
          select case (trim(rho_unit))
          case ('kg/m^3', 'kg/m3')
             rho_factor = 1d-3
          case ('g/cc', 'g/cm^3', 'g/cm3')
             rho_factor = 1
          case default
             write(stderr, *) '[ERROR] Unsupported unit: ', trim(rho_unit)
             return
          end select
       end if

       this%Dens(index) = rho * rho_factor
    end if

    if (present(h)) then
       if (present(h_unit)) then
          select case (trim(h_unit))
          case ('J/mol')
             h_factor = 1000 / R0
          case default
             write(stderr, *) '[ERROR] Unsupported unit: ', trim(h_unit)
             return
          end select
       end if

       this%Energy(index) = 'h,j'
       this%Enth(index) = h * h_factor
    end if

    if (present(u)) then
       if (present(u_unit)) then
          select case (trim(u_unit))
          case ('J/mol')
             u_factor = 1000 / R0
          case default
             write(stderr, *) '[ERROR] Unsupported unit: ', trim(u_unit)
             return
          end select
       end if

       this%Energy(index) = 'u,j'
       this%Enth(index) = u * u_factor
    end if

    this%Nlm = 0
    this%Size = 0
    this%Hsub0 = 1d30

    return
  end subroutine set_reactant


  subroutine set_insert_species(this, species)
    class(CEA_Problem), intent(inout):: this
    character(*), intent(in):: species(:)

    this%Nsert = min(20, size(species))
    this%ensert(1:this%Nsert) = species(1:this%Nsert)

    return
  end subroutine set_insert_species


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


  subroutine set_legacy_mode(this, legacy_mode)
    class(CEA_Problem), intent(inout):: this
    logical, intent(in), optional:: legacy_mode

    if (present(legacy_mode)) then
       this%legacy_mode = legacy_mode
    else
       this%legacy_mode = .true.
    end if

    return
  end subroutine set_legacy_mode


  subroutine run(this, out_filename, plt_filename)
    use mod_cea, only: run_case

    class(CEA_Problem), intent(inout):: this
    character(*), intent(in), optional:: out_filename
    character(*), intent(in), optional:: plt_filename

    call run_case(this, out_filename, plt_filename)

    return
  end subroutine run


  function get_pressure(this, iOF, ipt) result(P)
    real(8):: P
    class(CEA_Problem), intent(in):: this
    integer, intent(in):: iOF
    integer, intent(in):: ipt

    P = this%points(iOF, ipt)%Ppp * 1d-1

    return
  end function get_pressure


  function get_temperature(this, iOF, ipt) result(T)
    real(8):: T
    class(CEA_Problem), intent(in):: this
    integer, intent(in):: iOF
    integer, intent(in):: ipt

    T = this%points(iOF, ipt)%Ttt

    return
  end function get_temperature

  function get_chamber_temperature(this) result(T)
    real(8):: T
    class(CEA_Problem), intent(in):: this

    integer:: ipt

    if (.not. this%Rkt) then
       write(stderr, *) '[ERROR] This function is only for Rocket problem.'
       T = -1
       return
    end if

    if (this%Fac) then
       ipt = 2
    else
       ipt = 1
    end if

    T = this%points(1, ipt)%Ttt

    return
  end function get_chamber_temperature


  function get_molecular_weight(this, iOF, ipt) result(M)
    real(8):: M
    class(CEA_Problem), intent(in):: this
    integer, intent(in):: iOF
    integer, intent(in):: ipt

    M = this%points(iOF, ipt)%Wm

    return
  end function get_molecular_weight


  function get_specific_heat(this, iOF, ipt) result(cp)
    use mod_constants, only: R0
    real(8):: cp
    class(CEA_Problem), intent(in):: this
    integer, intent(in):: iOF
    integer, intent(in):: ipt

    cp = this%points(iOF, ipt)%Cpr * R0 / 1000

    return
  end function get_specific_heat


  function get_specific_heat_ratio(this, iOF, ipt) result(gamma)
    real(8):: gamma
    class(CEA_Problem), intent(in):: this
    integer, intent(in):: iOF
    integer, intent(in):: ipt

    gamma = this%points(iOF, ipt)%Gammas

    return
  end function get_specific_heat_ratio


  function get_characteristic_velocity(this) result(c_star)
    real(8):: c_star
    class(CEA_Problem), intent(in):: this

    c_star = this%Cstr

    return
  end function get_characteristic_velocity


  function get_specific_impulse(this, iOF, ipt, vacuum) result(Isp)
    use mod_constants, only: g0, R0
    use mod_types, only: CEA_Point

    real(8):: Isp
    class(CEA_Problem), intent(in):: this
    integer, intent(in):: iOF
    integer, intent(in):: ipt
    logical, intent(in), optional:: vacuum

    type(CEA_Point), pointer:: p
    p => this%points(iOF, ipt)

    Isp = p%Spim / g0

    if (present(vacuum)) then
       if (vacuum) then
          if (p%Wm * p%Spim == 0) then
             Isp = 0
          else
             Isp = (p%Spim + R0 * p%Ttt / (p%Wm * p%Spim)) / g0
          end if
       end if
    end if

    return
  end function get_specific_impulse


  subroutine calc_frozen_exhaust(this, T, cp, gamma, mu, k, Pr)
    use mod_ext, only: mod_ext_calc_frozen_exhaust => calc_frozen_exhaust

    class(CEA_Problem), intent(in):: this
    real(8), intent(in):: T  !< [K]
    real(8), intent(out), optional:: cp  !< [kJ/(kg·K)]
    real(8), intent(out), optional:: gamma  !< [-]
    real(8), intent(out), optional:: mu  !< [µPa·s]
    real(8), intent(out), optional:: k  !< [W/(m·K)]
    real(8), intent(out), optional:: Pr  !< [-]

    call mod_ext_calc_frozen_exhaust(this, T, cp, gamma, mu, k, Pr)

    return
  end subroutine calc_frozen_exhaust


  subroutine get_thermo_reference_properties(this, name, M, T_ref, h0_ref)
    use mod_ext, only: mod_ext_get_thermo_reference_properties => get_thermo_reference_properties

    class(CEA_Problem), intent(in):: this
    character(*), intent(in):: name
    real(8), intent(out):: M
    real(8), intent(out):: T_ref
    real(8), intent(out):: h0_ref

    call mod_ext_get_thermo_reference_properties(this, name, M, T_ref, h0_ref)

    return
  end subroutine get_thermo_reference_properties


  subroutine write_debug_output(this, filename)
    use mod_io, only: mod_io_write_debug_output => write_debug_output

    class(CEA_Problem), intent(in):: this
    character(*), intent(in):: filename

    call mod_io_write_debug_output(this, filename)

    return
  end subroutine write_debug_output

end module cea
