subroutine run_example_12()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'rocket', name = 'Example 12', equilibrium = .true., frozen = .true., frozen_at_throat = .true.)
  call prob%set_output_options(SI = .true., mass_fractions = .true.)

  call prob%add_reactant('fuel', 'CH6N2(L)', rho = 0.874d0)
  call prob%add_reactant('oxyd', 'N2O4(L)',  rho = 1.431d0)

  call prob%set_chamber_pressures([1000d0], unit = 'legacy-psi')
  call prob%set_mixture_ratios([2.5d0])
  call prob%set_pressure_ratios([68.0457d0])
  call prob%set_supersonic_area_ratios([5d0, 10d0, 25d0, 50d0, 75d0, 100d0, 150d0, 200d0])

  call prob%set_only_species(['CO    ', 'CO2   ', 'H     ', 'HNO   ', 'HNO2  ', 'HO2   ', 'H2    ', 'H2O   ', 'H2O2  ', &
                              'N     ', 'NO    ', 'NO2   ', 'N2    ', 'N2O   ', 'O     ', 'OH    ', 'O2    ', 'HCO   ', &
                              'NH    ', 'CH4   ', 'NH2   ', 'NH3   ', 'H2O(L)', 'C(gr) '])

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-12.out')

  return
end subroutine run_example_12
