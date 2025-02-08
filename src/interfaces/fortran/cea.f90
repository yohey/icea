module cea
  implicit none

contains

  subroutine set_problem(cea, mode, name, mole_ratios, equilibrium, ions, &
       frozen, frozen_at_throat)
    use mod_constants, only: stderr
    use mod_cea_types, only: CEA_Problem

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
    use mod_cea_types, only: CEA_Problem, maxMix

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
    use mod_cea_types, only: CEA_Problem

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
    use mod_cea_types, only: CEA_Problem

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
    use mod_cea_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    real(8), intent(in):: ratio_list(:)

    cea%Npp_in = size(ratio_list)
    cea%Pcp(1:cea%Npp_in) = ratio_list

    return
  end subroutine set_pressure_ratios

  subroutine set_subsonic_area_ratios(cea, ratio_list)
    use mod_cea_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    real(8), intent(in):: ratio_list(:)

    cea%Nsub_in = size(ratio_list)
    cea%Subar(1:cea%Nsub_in) = ratio_list

    return
  end subroutine set_subsonic_area_ratios

  subroutine set_supersonic_area_ratios(cea, ratio_list)
    use mod_cea_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    real(8), intent(in):: ratio_list(:)

    cea%Nsup_in = size(ratio_list)
    cea%Supar(1:cea%Nsup_in) = ratio_list

    return
  end subroutine set_supersonic_area_ratios

  subroutine set_finite_area_combustor(cea, contraction_ratio, mass_flow_ratio)
    use mod_constants, only: stderr
    use mod_cea_types, only: CEA_Problem

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
    use mod_cea_types, only: CEA_Problem

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
    use mod_cea_types, only: CEA_Problem

    type(CEA_Problem), intent(inout):: cea
    character(*), intent(in):: species(:)

    cea%Nsert = min(20, size(species))
    cea%ensert(1:cea%Nsert) = species(1:cea%Nsert)

    return
  end subroutine insert_species

  subroutine set_omit_species(cea, species)
    use mod_cea_types, only: CEA_Problem, maxNgc

    type(CEA_Problem), intent(inout):: cea
    character(*), intent(in):: species(:)

    cea%Nomit = min(maxNgc, size(species))
    cea%Omit(1:cea%Nomit) = species(1:cea%Nomit)

    return
  end subroutine set_omit_species

  subroutine set_only_species(cea, species)
    use mod_cea_types, only: CEA_Problem, maxNgc

    type(CEA_Problem), intent(inout):: cea
    character(*), intent(in):: species(:)

    cea%Nonly_in = min(maxNgc, size(species))
    cea%Prod_in(:) = ''
    cea%Prod_in(1:cea%Nonly_in) = species(1:cea%Nonly_in)
    cea%Nonly = cea%Nonly_in
    cea%Prod(:) = cea%Prod_in(:)

    return
  end subroutine set_only_species

end module cea
