module mod_functions
  implicit none

contains

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @brief Calculate specific heat from temperature
  !! NASA RP-1311-1, p. 20, eq. (4.9)
  !---------------------------------------------------------------------------------------------------------------------------------
  pure real(8) function calc_specific_heat(coeffs, T)
    real(8), intent(in):: coeffs(9), T

    calc_specific_heat = coeffs(7) * T**4 + coeffs(6) * T**3 + coeffs(5) * T**2 + coeffs(4) * T &
         + coeffs(3) + coeffs(2) / T + coeffs(1) / T**2

    return
  end function calc_specific_heat

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @brief Calculate enthalpy from temperature
  !! NASA RP-1311-1, p. 20, eq. (4.10)
  !---------------------------------------------------------------------------------------------------------------------------------
  pure real(8) function calc_enthalpy(coeffs, T, lnT, legacy_mode)
    real(8), intent(in):: coeffs(9), T, lnT
    logical, intent(in):: legacy_mode

    if (legacy_mode) then
       calc_enthalpy = coeffs(7) * (T**4 / 5) + coeffs(6) * (T**3 / 4) + (0.6666666666666667d0 / 2) * coeffs(5) * T**2 &
            + coeffs(4) * (T / 2) - coeffs(1) * (1 / T**2) + coeffs(2) * (lnT / T) + coeffs(3) + coeffs(8) / T
    else
       calc_enthalpy = -coeffs(1) / T**2 + coeffs(2) * lnT / T + coeffs(3) &
             + coeffs(4) * T / 2 + coeffs(5) * T**2 / 3 + coeffs(6) * T**3 / 4 + coeffs(7) * T**4 / 5 + coeffs(8) / T
    end if

    return
  end function calc_enthalpy

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @brief Calculate entropy from temperature
  !! NASA RP-1311-1, p. 20, eq. (4.11)
  !---------------------------------------------------------------------------------------------------------------------------------
  pure real(8) function calc_entropy(coeffs, T, lnT, legacy_mode)
    real(8), intent(in):: coeffs(9), T, lnT
    logical, intent(in):: legacy_mode

    if (legacy_mode) then
       calc_entropy = (coeffs(7) * T**4 / 4 + coeffs(6) * T**3 / 3) + coeffs(5) * T**2 / 2 + coeffs(4) * T &
            - coeffs(1) / (2 * T**2) - coeffs(2) / T + coeffs(3) * lnT + coeffs(9)
    else
       calc_entropy = -coeffs(1) / (2 * T**2) - coeffs(2) / T + coeffs(3) * lnT &
             + coeffs(4) * T + coeffs(5) * T**2 / 2 + coeffs(6) * T**3 / 3 + coeffs(7) * T**4 / 4 + coeffs(9)
    end if

    return
  end function calc_entropy

end module mod_functions
