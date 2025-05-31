subroutine run_example_07()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'shock', name = 'Example 7', &
                        mole_ratios = .true., equilibrium = .true., frozen = .true., incident = .true.)

  call prob%add_reactant('name', 'H2', ratio = 0.05d0, T = 300d0)
  call prob%add_reactant('name', 'O2', ratio = 0.05d0, T = 300d0)
  call prob%add_reactant('name', 'Ar', ratio = 0.90d0, T = 300d0)

  call prob%set_chamber_pressures([10d0, 20d0], unit = 'mmHg')

  call prob%set_initial_velocities([1000d0, 1100d0, 1200d0, 1250d0, 1300d0, 1350d0, 1400d0])

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-07.out')

  return
end subroutine run_example_07
