module mod_io
  implicit none

contains

  subroutine set_problem(cea, mode, name, mole_ratios, equilibrium, ions, &
       frozen, frozen_at_throat)
    use mod_constants, only: stderr
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    character(*), intent(in):: mode
    character(*), intent(in), optional:: name
    logical, intent(in), optional:: mole_ratios
    logical, intent(in), optional:: equilibrium
    logical, intent(in), optional:: ions
    logical, intent(in), optional:: frozen
    logical, intent(in), optional:: frozen_at_throat

    if (mode == 'rocket') then
       cea%Rkt = .true.
       cea%Sp = .true.
       cea%Tp = .false.
       cea%Hp = .false.
       cea%Detn = .false.
       cea%Shock = .false.

       cea%Nt = 1

       cea%Nfz = 1

       if (present(frozen)) then
          cea%Froz = frozen

          if (present(frozen_at_throat) .and. frozen_at_throat) then
             if (cea%Fac) then
                cea%Nfz = 3
             else
                cea%Nfz = 2
             end if
          end if
       end if

    else if (mode == 'thermp') then
       write(stderr, *) '[ERROR] Not implemented yet.'
       return
    else if (mode == 'deton') then
       write(stderr, *) '[ERROR] Not implemented yet.'
       return
    else if (mode == 'shock') then
       write(stderr, *) '[ERROR] Not implemented yet.'
       return
    else
       write(stderr, *) '[ERROR] Invalid mode: ', mode
       return
    end if

    if (present(name)) then
       cea%Case = trim(name)
    end if

    if (present(equilibrium)) then
       cea%Eql_in = equilibrium
    end if

    if (present(ions)) then
       cea%Ions = ions
    end if

    if (present(mole_ratios)) then
       cea%Moles = mole_ratios
    end if

    return
  end subroutine set_problem

  subroutine set_output_options(cea, SI, debug_points, mass_fractions, short, trace_tol, transport)
    use mod_types, only: CEA_Problem, maxMix

    type(CEA_Problem), intent(inout):: cea
    logical, intent(in), optional:: SI
    integer, intent(in), optional:: debug_points(:)
    logical, intent(in), optional:: mass_fractions
    logical, intent(in), optional:: short
    real(8), intent(in), optional:: trace_tol
    logical, intent(in), optional:: transport
    integer:: i, j

    if (present(SI)) then
       cea%SIunit = SI
    end if

    if (present(debug_points)) then
       do concurrent (i = 1:maxMix, j = 1:size(debug_points), debug_points(j) <= cea%max_points)
          cea%points(i, debug_points(j))%Debug = .true.
       end do
    end if

    if (present(mass_fractions)) then
       cea%Massf = mass_fractions
    end if

    if (present(short)) then
       cea%Short = short
    end if

    if (present(trace_tol)) then
       cea%Trace = trace_tol
    end if

    if (present(transport)) then
       cea%Trnspt = transport
    end if

    return
  end subroutine set_output_options

  subroutine set_chamber_pressures(cea, pressure_list, unit)
    use mod_constants, only: stderr, g0, lb_to_kg, in_to_m
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
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

    cea%Np = size(pressure_list)
    cea%P(1:cea%Np) = pressure_list * factor

    return
  end subroutine set_chamber_pressures

  subroutine set_mixture_ratios(cea, ratio_list, type)
    use mod_constants, only: stderr
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    real(8), intent(in):: ratio_list(:)
    character(*), intent(in), optional:: type
    integer:: iof

    cea%Nof = size(ratio_list)
    cea%Oxf(1:cea%Nof) = ratio_list

    if (present(type)) then
       if (trim(type) == '%fuel') then
          do concurrent (iof = 1:cea%Nof, cea%Oxf(iof) > 0)
             cea%Oxf(iof) = (100 - cea%Oxf(iof)) / cea%Oxf(iof)
          end do
       else
          write(stderr, *) '[ERROR] Unsupported type: ', trim(type)
          return
       end if
    end if

    return
  end subroutine set_mixture_ratios

  subroutine set_pressure_ratios(cea, ratio_list)
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    real(8), intent(in):: ratio_list(:)

    cea%Npp_in = size(ratio_list)
    cea%Pcp(1:cea%Npp_in) = ratio_list

    return
  end subroutine set_pressure_ratios

  subroutine set_subsonic_area_ratios(cea, ratio_list)
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    real(8), intent(in):: ratio_list(:)

    cea%Nsub_in = size(ratio_list)
    cea%Subar(1:cea%Nsub_in) = ratio_list

    return
  end subroutine set_subsonic_area_ratios

  subroutine set_supersonic_area_ratios(cea, ratio_list)
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    real(8), intent(in):: ratio_list(:)

    cea%Nsup_in = size(ratio_list)
    cea%Supar(1:cea%Nsup_in) = ratio_list

    return
  end subroutine set_supersonic_area_ratios

  subroutine set_finite_area_combustor(cea, contraction_ratio, mass_flow_ratio)
    use mod_constants, only: stderr
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    real(8), intent(in), optional:: contraction_ratio
    real(8), intent(in), optional:: mass_flow_ratio

    if (present(contraction_ratio)) then
       cea%AcAt_in = contraction_ratio
    else if (present(mass_flow_ratio)) then
       cea%mdotByAc_in = mass_flow_ratio
    else
       write(stderr, *) '[ERROR] Either contraction_ratio or mass_flow_ratio must be given.'
       return
    end if

    cea%Fac = .true.

    if (cea%Froz .and. cea%Nfz == 2) then
       cea%Nfz = 3
    end if

    return
  end subroutine set_finite_area_combustor


  subroutine add_reactant(cea, type, name, ratio, T, rho, T_unit, rho_unit)
    use mod_constants, only: stderr
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
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

    cea%Nreac = cea%Nreac + 1

    cea%Fox(cea%Nreac) = trim(type)
    cea%Rname(cea%Nreac) = trim(name)

    if (present(ratio)) then
       cea%Pecwt(cea%Nreac) = ratio
    else
       cea%Pecwt(cea%Nreac) = 1
    end if

    if (present(T)) cea%Rtemp(cea%Nreac) = T * T_factor

    if (present(rho)) cea%Dens(cea%Nreac) = rho * rho_factor

    cea%Energy(cea%Nreac) = 'lib'

    return
  end subroutine add_reactant

  subroutine insert_species(cea, species)
    use mod_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    character(*), intent(in):: species(:)

    cea%Nsert = min(20, size(species))
    cea%ensert(1:cea%Nsert) = species(1:cea%Nsert)

    return
  end subroutine insert_species

  subroutine set_omit_species(cea, species)
    use mod_types, only: CEA_Problem, maxNgc

    type(CEA_Problem), intent(inout):: cea
    character(*), intent(in):: species(:)

    cea%Nomit = min(maxNgc, size(species))
    cea%Omit(1:cea%Nomit) = species(1:cea%Nomit)

    return
  end subroutine set_omit_species

  subroutine set_only_species(cea, species)
    use mod_types, only: CEA_Problem, maxNgc

    type(CEA_Problem), intent(inout):: cea
    character(*), intent(in):: species(:)

    cea%Nonly_in = min(maxNgc, size(species))
    cea%Prod_in(:) = ''
    cea%Prod_in(1:cea%Nonly_in) = species(1:cea%Nonly_in)
    cea%Nonly = cea%Nonly_in
    cea%Prod(:) = cea%Prod_in(:)

    return
  end subroutine set_only_species


  subroutine write_debug_output(cea, filename)
    use mod_constants
    use mod_types, only: CEA_Problem, CEA_Point

    type(CEA_Problem), intent(in):: cea
    character(*), intent(in):: filename

    integer:: io_out, iof, ipt, j, max_points
    real(8), allocatable:: X(:, :)
    character(12), allocatable:: labels(:)
    type(CEA_Point), pointer:: p


    if (cea%Shock) then
       max_points = cea%Nsk
    else
       max_points = cea%Nt * cea%Np
       if (cea%Rkt) then
          max_points = max_points * (cea%Npp + cea%Nsub + cea%Nsup)
       end if
    end if

    write(0, *) '[DEBUG] max_points = ', max_points

    open(newunit = io_out, file = trim(filename), status = 'unknown', form = 'formatted')

    write(io_out, '("Case = ", a)') trim(cea%Case)

    write(io_out, '(/"OPTIONS:")')
    write(io_out, '("Nof  = ", i5)') cea%Nof
    write(io_out, '("Nt   = ", i5)') cea%Nt
    write(io_out, '("Np   = ", i5)') cea%Np
    write(io_out, '("Npp  = ", i5)') cea%Npp
    write(io_out, '("Nsub = ", i5)') cea%Nsub
    write(io_out, '("Nsup = ", i5)') cea%Nsup

    write(io_out, '("TP = ", l1)') cea%Tp
    write(io_out, '("HP = ", l1)') (cea%Hp .and. .not. cea%Vol)
    write(io_out, '("SP = ", l1)') cea%Sp
    write(io_out, '("TV = ", l1)') (cea%Tp .and. cea%Vol)
    write(io_out, '("UV = ", l1)') (cea%Hp .and. cea%Vol)
    write(io_out, '("SV = ", l1)') (cea%Sp .and. cea%Vol)

    write(io_out, '("DETN   = ", l1)') cea%Detn
    write(io_out, '("SHOCK  = ", l1)') cea%Shock
    write(io_out, '("REFL   = ", l1)') (cea%Shock .and. (cea%Reflfz .or. cea%Refleq))
    write(io_out, '("INCD   = ", l1)') (cea%Shock .and. (cea%Incdfz .or. cea%Incdeq))
    write(io_out, '("Rkt    = ", l1)') cea%Rkt
    write(io_out, '("FROZ   = ", l1)') cea%Froz
    write(io_out, '("EQL    = ", l1)') cea%Eql
    write(io_out, '("IONS   = ", l1)') cea%Ions
    write(io_out, '("SIUNIT = ", l1)') cea%SIunit
    write(io_out, '("DEBUGF = ", l1)') cea%Debugf
    write(io_out, '("SHKDBG = ", l1)') cea%Shkdbg
    write(io_out, '("DETDBG = ", l1)') cea%Detdbg
    write(io_out, '("TRNSPT = ", l1)') cea%Trnspt

    if (max_points == 0) return

    allocate(X(20, cea%Nof * max_points))

    if (cea%Rkt) then
       allocate(labels(cea%Nof * max_points))
       labels(:) = 'EXIT'
       if (cea%Fac) then
          labels(1) = 'INJECTOR'
          labels(2) = 'COMB END'
          labels(3) = 'THROAT'
       else
          labels(1) = 'CHAMBER'
          labels(2) = 'THROAT'
       end if
    end if

    do concurrent (iof = 1:cea%Nof, ipt = 1:max_points)
       p => cea%points(iof, ipt)
       j = (iof - 1) * max_points + ipt
       X(1, j) = cea%Oxf(iof)
       X(2, j) = p%Ppp * 0.1d0
       X(3, j) = p%Ttt
       if (p%Vlm == 0) then
          X(4, j) = 0
       else
          X(4, j) = 1d5 / p%Vlm
       end if
       X(5, j) = p%Hsum * cea%R
       X(6, j) = (p%Hsum - p%Ppp * p%Vlm / R0) * cea%R
       X(7, j) = (p%Hsum - p%Ttt * p%Ssum) * cea%R
       X(8, j) = p%Ssum * cea%R
       X(9, j) = p%Wm
       X(10, j) = p%Dlvpt
       X(11, j) = p%Dlvtp
       X(12, j) = p%Cpr * cea%R
       X(13, j) = p%Gammas
       X(14, j) = p%Sonvel
       X(15, j) = p%Vmoc
    end do

    do concurrent (iof = 1:cea%Nof, ipt = 2:max_points)
       p => cea%points(iof, ipt)
       j = (iof - 1) * max_points + ipt
       X(16, j) = p%AeAt
       X(17, j) = cea%Cstr
       if (cea%Cstr == 0) then
          X(18, j) = 0
       else
          X(18, j) = p%Spim / cea%Cstr
       end if
       if (p%Wm * p%Spim == 0) then
          X(19, j) = 0
       else
          X(19, j) = (p%Spim + R0 * p%Ttt / (p%Wm * p%Spim)) / g0
       end if
       X(20, j) = p%Spim / g0
    end do

    if (.not. cea%SIunit) then
       X(5:8, :) = X(5:8, :) * cal_to_J
       X(12, :) = X(12, :) * cal_to_J
    end if

    write(io_out, *)
    if (cea%Rkt) write(io_out, '(17x, *(a12))') adjustr(labels(:))
    write(io_out, '("O/F            = ", *(f12.5))') X(1, :)
    write(io_out, '("P [MPa]        = ", *(f12.7))') X(2, :)
    write(io_out, '("T [K]          = ", *(f12.2))') X(3, :)
    write(io_out, '("ρ [kg/m^3]     = ", *(f12.6))') X(4, :)
    write(io_out, '("H [kJ/kg]      = ", *(f12.2))') X(5, :)
    write(io_out, '("U [kJ/kg]      = ", *(f12.2))') X(6, :)
    write(io_out, '("G [kJ/kg]      = ", *(f12.1))') X(7, :)
    write(io_out, '("S [kJ/(kg·K)]  = ", *(f12.4))') X(8, :)

    write(io_out, '("M [g/mol]      = ", *(f12.3))') X(9, :)
    write(io_out, '("(dlnV/dlnP)T   = ", *(f12.4))') X(10, :)
    write(io_out, '("(dlnV/dlnT)P   = ", *(f12.4))') X(11, :)
    write(io_out, '("cp [kJ/(kg·K)] = ", *(f12.4))') X(12, :)
    write(io_out, '("γ [-]          = ", *(f12.4))') X(13, :)
    write(io_out, '("a [m/s]        = ", *(f12.2))') X(14, :)
    write(io_out, '("Ma [-]         = ", *(f12.4))') X(15, :)

    write(io_out, '("Ae/At [-]      = ", 12x, *(f12.4))') X(16, 2:)
    write(io_out, '("C* [m/s]       = ", 12x, *(f12.2))') X(17, 2:)
    write(io_out, '("CF [-]         = ", 12x, *(f12.4))') X(18, 2:)
    write(io_out, '("Ivac [s]       = ", 12x, *(f12.3))') X(19, 2:)
    write(io_out, '("Isp [s]        = ", 12x, *(f12.3))') X(20, 2:)

    deallocate(X)
    if (cea%Rkt) deallocate(labels)

    close(io_out)

    return
  end subroutine write_debug_output

end module mod_io
