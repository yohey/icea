module mod_ext
  implicit none
  private

  public:: calc_frozen_exhaust

contains

  subroutine calc_frozen_exhaust(cea, T, P, cp, gamma)
    use mod_constants, only: R0
    use mod_types, only: CEA_Core_Problem, CEA_Point

    class(CEA_Core_Problem), intent(in):: cea
    real(8), intent(in):: T  !< [K]
    real(8), intent(out), optional:: P  !< [MPa]
    real(8), intent(out), optional:: cp  !< [kJ/(kg*K)]
    real(8), intent(out), optional:: gamma

    integer:: j
    real(8):: cp_tmp, s, lnwm(cea%Ng), lnPsum
    real(8):: cp_array(cea%Ng), h0_array(cea%Ng), s_array(cea%Ng)
    type(CEA_Point), pointer:: pfz

    if (.not. cea%Rkt) then
       write(0, *) '[ERROR] This function is only for Rocket problem.'
       if (present(P)) P = 0
       if (present(cp)) cp = 0
       if (present(gamma)) gamma = 0
       return
    end if

    pfz => cea%points(cea%iOF, cea%Nfz)

    call calc_thermo_properties(cea, T, cp_array, h0_array, s_array)

    cp_tmp = sum(pfz%En(1:cea%Ng) * cp_array(1:cea%Ng))

    if (cea%Npr /= 0) then
       cp_tmp = cp_tmp + sum(pfz%En(cea%Jcond(1:cea%Npr)) * cp_array(cea%Jcond(1:cea%Npr)))
    end if

    if (present(cp)) cp = cp_tmp * R0 * 1d-3

    if (present(gamma)) gamma = cp_tmp / (cp_tmp - 1 / pfz%Wm)

    if (present(P)) then
       s = pfz%Ssum

       lnwm(:) = 0
       do concurrent (j = 1:cea%Ng, pfz%En(j) /= 0)
          lnwm(j) = -log(pfz%En(j) * pfz%Wm)
       end do

       lnPsum = sum(pfz%En(1:cea%Ng) * (s_array(1:cea%Ng) + lnwm(1:cea%Ng))) - s

       if (cea%Npr /= 0) then
          lnPsum = lnPsum + sum(pfz%En(cea%Jcond(1:cea%Npr)) * s_array(cea%Jcond(1:cea%Npr)))
       end if

       P = 0.1d0 * exp(lnPsum / sum(pfz%En(1:cea%Ng)))
    end if

    return
  end subroutine calc_frozen_exhaust


  subroutine calc_thermo_properties(cea, T, cp, h0, s)
    use mod_types, only: CEA_Core_Problem
    use mod_functions, only: calc_enthalpy, calc_entropy, calc_specific_heat

    class(CEA_Core_Problem), intent(in):: cea
    real(8), intent(in):: T
    real(8), intent(out):: cp(cea%Ng), h0(cea%Ng), s(cea%Ng)

    integer:: ij, j, jj, k
    real(8):: lnT

    lnT = log(T)

    k = 1
    if (T > cea%Tg(2)) k = 2
    if (T > cea%Tg(3)) k = 3

    do concurrent (j = 1:cea%Ng)
       s(j) = calc_entropy(cea%Coef(:, j, k), T, lnT, cea%legacy_mode)
       h0(j) = calc_enthalpy(cea%Coef(:, j, k), T, lnT, cea%legacy_mode)
    end do

    if (.not. cea%Tp .or. cea%Convg) then
       do concurrent (j = 1:cea%Ng)
          cp(j) = calc_specific_heat(cea%Coef(:, j, k), T)
       end do
    end if

    if (cea%Npr /= 0 .and. k /= 3 .and. cea%Ng /= cea%Ngc) then
       do concurrent (ij = 1:cea%Npr)
          j = cea%Jcond(ij)
          jj = cea%Jcond(ij) - cea%Ng

          s(j) = calc_entropy(cea%Cft(:, jj), T, lnT, cea%legacy_mode)
          h0(j) = calc_enthalpy(cea%Cft(:, jj), T, lnT, cea%legacy_mode)
          cp(j) = calc_specific_heat(cea%Cft(:, jj), T)
       end do
    end if

    return
  end subroutine calc_thermo_properties

end module mod_ext
