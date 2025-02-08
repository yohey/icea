program test_newio
  implicit none

  write(0, *) '[DEBUG] call test_orig'
  call test_orig

  write(0, *) '[DEBUG] call test_rocket'
  call test_rocket

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

  call run_all_cases(cea(8:13), out_filename)

  call write_debug_output(cea(8), 'orig-08.debug.out')
  call write_debug_output(cea(9), 'orig-09.debug.out')
  call write_debug_output(cea(10), 'orig-10.debug.out')
  call write_debug_output(cea(11), 'orig-11.debug.out')
  call write_debug_output(cea(12), 'orig-12.debug.out')
  call write_debug_output(cea(13), 'orig-13.debug.out')

  deallocate(cea)

  return
end subroutine test_orig


subroutine test_rocket
  use mod_cea
  use mod_types
  use mod_io
  use mod_legacy_io, only: REACT

  integer:: icase
  type(CEA_Problem), allocatable:: cea(:)

  allocate(cea(8:13))

  do icase = 8, 13
     call init_case(cea(icase))
  end do

  ! Examples 8 mod
  call set_problem(cea(8), mode = 'rocket', name = 'Example-8-mod', equilibrium = .true., mole_ratios = .false.)
  call set_output_options(cea(8), SI = .true.)
  call set_chamber_pressures(cea(8), [53.3172d0])
  call set_mixture_ratios(cea(8), [5.55157d0])
  call set_pressure_ratios(cea(8), [1d1, 1d2, 1d3])
  call set_subsonic_area_ratios(cea(8), [1.58d0])
  call set_supersonic_area_ratios(cea(8), [25d0, 50d0, 75d0])
  call add_reactant(cea(8), 'fuel', 'H2(L)', 100d0, T = 20.27d0)
  call add_reactant(cea(8), 'oxyd', 'O2(L)', 100d0, T = 90.17d0)

  ! Examples 9 mod
  call set_problem(cea(9), mode = 'rocket', name = 'Example-9-mod', equilibrium = .false., mole_ratios = .false.)
  call set_output_options(cea(9), SI = .true.)
  call set_finite_area_combustor(cea(9), contraction_ratio = 1.58d0)
  call set_chamber_pressures(cea(9), [53.3172d0])
  call set_mixture_ratios(cea(9), [5.55157d0])
  call set_pressure_ratios(cea(9), [1d1, 1d2, 1d3])
  call set_supersonic_area_ratios(cea(9), [25d0, 50d0, 75d0])
  call add_reactant(cea(9), 'fuel', 'H2(L)', 100d0, T = 20.27d0)
  call add_reactant(cea(9), 'oxyd', 'O2(L)', 100d0, T = 90.17d0)

  ! Examples 10 mod
  call set_problem(cea(10), mode = 'rocket', name = 'Example-10-mod', equilibrium = .false., mole_ratios = .false.)
  call set_output_options(cea(10), SI = .true., short = .true.)
  call set_finite_area_combustor(cea(10), mass_flow_ratio = 1333.9d0)
  call set_chamber_pressures(cea(10), [53.3172d0])
  call set_mixture_ratios(cea(10), [5.55157d0])
  call set_pressure_ratios(cea(10), [1d1, 1d2, 1d3])
  call set_supersonic_area_ratios(cea(10), [25d0, 50d0, 75d0])
  call add_reactant(cea(10), 'fuel', 'H2(L)', 100d0, T = 20.27d0)
  call add_reactant(cea(10), 'oxyd', 'O2(L)', 100d0, T = 90.17d0)

  ! Examples 11 mod
  call set_problem(cea(11), mode = 'rocket', name = 'Example-11-mod', equilibrium = .true., mole_ratios = .true., ions = .true.)
  call set_output_options(cea(11), SI = .true., transport = .true.)
  call set_chamber_pressures(cea(11), [1000d0], unit = 'legacy-psi')
  call set_pressure_ratios(cea(11), [68.0457d0])
  call set_subsonic_area_ratios(cea(11), [10d0])
  call set_supersonic_area_ratios(cea(11), [10d0, 20d0, 100d0])
  call add_reactant(cea(11), 'fuel', 'Li(cr)', 1.0000d0, T = 298.15d0)
  call add_reactant(cea(11), 'oxyd', 'F2(L)',  0.5556d0, T =  85.02d0)

  ! Examples 12 mod
  call set_problem(cea(12), mode = 'rocket', name = 'Example-12-mod', &
       equilibrium = .true., mole_ratios = .false., ions = .false., frozen = .true., frozen_at_throat = .true.)
  call set_output_options(cea(12), SI = .true., mass_fractions = .true.)
  call set_chamber_pressures(cea(12), [1000d0], unit = 'legacy-psi')
  call set_mixture_ratios(cea(12), [2.5d0])
  call set_pressure_ratios(cea(12), [68.0457d0])
  call set_supersonic_area_ratios(cea(12), [5d0, 10d0, 25d0, 50d0, 75d0, 100d0, 150d0, 200d0])
  call add_reactant(cea(12), 'fuel', 'CH6N2(L)', rho = 0.874d0)
  call add_reactant(cea(12), 'oxyd', 'N2O4(L)',  rho = 1.431d0)
  call set_only_species(cea(12), &
       ['CO    ', 'CO2   ', 'H     ', 'HNO   ', 'HNO2  ', 'HO2   ', 'H2    ', 'H2O   ', 'H2O2  ', &
        'N     ', 'NO    ', 'NO2   ', 'N2    ', 'N2O   ', 'O     ', 'OH    ', 'O2    ', 'HCO   ', &
        'NH    ', 'CH4   ', 'NH2   ', 'NH3   ', 'H2O(L)', 'C(gr) '])

  ! Examples 13 mod
  call set_problem(cea(13), mode = 'rocket', name = 'Example-13-mod', equilibrium = .true., mole_ratios = .false., ions = .false.)
  call set_output_options(cea(13), SI = .false., trace_tol = 1d-10)
  call set_chamber_pressures(cea(13), [3000d0], unit = 'psi')
  call set_mixture_ratios(cea(13), [67d0], type = '%fuel')
  call set_pressure_ratios(cea(13), [3d0, 10d0, 30d0, 300d0])
  call add_reactant(cea(13), 'fuel', 'N2H4(L)', 0.8d0, T = 298.15d0, rho = 0.874d0)
  call add_reactant(cea(13), 'fuel', 'Be(a)',   0.2d0, T = 298.15d0, rho = 1.431d0)
  call add_reactant(cea(13), 'oxyd', 'H2O2(L)', 1.0d0, T = 298.15d0)
  call insert_species(cea(13), ['BeO(L)'])

  do icase = 8, 13
     cea(icase)%legacy_mode = .true.
  end do

  call write_debug_output(cea(8), 'rocket-08.debug-01.out')
  call write_debug_output(cea(9), 'rocket-09.debug-01.out')
  call write_debug_output(cea(10), 'rocket-10.debug-01.out')
  call write_debug_output(cea(11), 'rocket-11.debug-01.out')
  call write_debug_output(cea(12), 'rocket-12.debug-01.out')
  call write_debug_output(cea(13), 'rocket-13.debug-01.out')

  call run_all_cases(cea(8:13), 'test-rocket-1.out')

  call write_debug_output(cea(8), 'rocket-08.debug-02.out')
  call write_debug_output(cea(9), 'rocket-09.debug-02.out')
  call write_debug_output(cea(10), 'rocket-10.debug-02.out')
  call write_debug_output(cea(11), 'rocket-11.debug-02.out')
  call write_debug_output(cea(12), 'rocket-12.debug-02.out')
  call write_debug_output(cea(13), 'rocket-13.debug-02.out')

  call run_all_cases(cea(8:13), 'test-rocket-2.out')

  call write_debug_output(cea(8), 'rocket-08.debug-03.out')
  call write_debug_output(cea(9), 'rocket-09.debug-03.out')
  call write_debug_output(cea(10), 'rocket-10.debug-03.out')
  call write_debug_output(cea(11), 'rocket-11.debug-03.out')
  call write_debug_output(cea(12), 'rocket-12.debug-03.out')
  call write_debug_output(cea(13), 'rocket-13.debug-03.out')

  return
end subroutine test_rocket
