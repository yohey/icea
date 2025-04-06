module mod_ext
  implicit none
  private

  public:: calc_frozen_exhaust, get_thermo_reference_properties

contains

  subroutine calc_frozen_exhaust(cea, T, P, cp, gamma, mu, k, Pr)
    use mod_constants, only: R0
    use mod_types, only: CEA_Core_Problem, CEA_Point
    use mod_cea, only: calc_thermo_properties

    class(CEA_Core_Problem), intent(in):: cea
    real(8), intent(in):: T  !< [K]
    real(8), intent(out), optional:: P  !< [MPa]
    real(8), intent(out), optional:: cp  !< [kJ/(kg·K)]
    real(8), intent(out), optional:: gamma  !< [-]
    real(8), intent(out), optional:: mu  !< [µPa·s]
    real(8), intent(out), optional:: k  !< [W/(m·K)]
    real(8), intent(out), optional:: Pr  !< [-]

    integer:: j
    real(8):: cp_tmp, s, lnwm(cea%Ng), lnPsum
    real(8):: cp_array(cea%Ng), h0_array(cea%Ng), s_array(cea%Ng)
    real(8):: mu_tmp, k_tmp, Pr_tmp
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

    if (present(mu) .or. present(k) .or. present(Pr)) then
       call calc_frozen_transport_properties(cea, T, mu_tmp, k_tmp, Pr_tmp)
       if (present(mu)) mu = mu_tmp
       if (present(k)) k = k_tmp
       if (present(Pr)) Pr = Pr_tmp
    end if

    return
  end subroutine calc_frozen_exhaust


  subroutine get_thermo_reference_properties(cea, name, M, T_ref, h0_ref)
    use mod_types, only: CEA_Core_Problem, ThermoProperty

    class(CEA_Core_Problem), intent(in):: cea
    character(*), intent(in):: name
    real(8), intent(out):: M
    real(8), intent(out):: T_ref
    real(8), intent(out):: h0_ref

    integer:: i
    type(ThermoProperty), pointer:: th

    M = 0
    T_ref = 0
    h0_ref = 0

    do i = 1, size(cea%thermo_properties)
       th => cea%thermo_properties(i)

       if (th%name == name .or. th%name == '*' // name) then
          if (th%type == 3) then
             M = th%mwt
             T_ref = th%tl(1)
             h0_ref = th%thermo(1, 1)
          else
             write(0, *) '[ERROR] Not implemented yet.'
          end if
       end if
    end do

    return
  end subroutine get_thermo_reference_properties


  subroutine calc_frozen_transport_properties(cea, T, mu, k, Pr)
    use mod_constants, only: stderr, R0
    use mod_types, only: CEA_Core_Problem, maxTr, maxNgc
    use mod_cea, only: TRANIN_core

    class(CEA_Core_Problem), intent(in):: cea
    real(8), intent(in):: T
    real(8), intent(out):: mu, k, Pr

    integer:: ipt, Nm, idummy1(2), idummy2(2, maxTr)
    real(8):: Cprr(maxTr), Con(maxTr), Wmol(maxTr), Xs(maxTr), Eta(maxTr, maxTr)
    real(8):: dummy1(2, maxTr), dummy2(maxTr, maxTr), dummy3(maxNgc)

    integer:: i, j, i1
    real(8):: phi(maxTr, maxTr), psi(maxTr, maxTr), rtpd(maxTr, maxTr), sumc, sumv
    real(8):: cp

    ipt = cea%Nfz + 1

    call TRANIN_core(cea, ipt, T, Nm, idummy1(1), idummy1(2), idummy2(1, :), idummy2(1, :), &
         Cprr, Con, Wmol, Xs, dummy1(1, :), dummy1(2, :), dummy2(:, :), Eta, dummy3(:))

    rtpd(:, :) = 0
    phi(:, :) = 1
    psi(:, :) = 1

    mu = 0
    k = 0

    do i = 1, Nm - 1
       i1 = i + 1
       do j = i1, Nm
          sumc = 2 / (Eta(i, j) * (Wmol(i) + Wmol(j)))
          phi(i, j) = sumc * Wmol(j) * Eta(i, i)
          phi(j, i) = sumc * Wmol(i) * Eta(j, j)
          sumc = (Wmol(i) + Wmol(j))**2
          psi(i, j) = phi(i, j) * (1 + 2.41d0 * (Wmol(i) - Wmol(j)) * (Wmol(i) - 0.142d0 * Wmol(j)) / sumc)
          psi(j, i) = phi(j, i) * (1 + 2.41d0 * (Wmol(j) - Wmol(i)) * (Wmol(j) - 0.142d0 * Wmol(i)) / sumc)
       end do
    end do

    do i = 1, Nm
       sumc = 0
       sumv = 0
       do j = 1, Nm
          sumc = sumc + psi(i, j) * Xs(j)
          sumv = sumv + phi(i, j) * Xs(j)
       end do
       k = k + Con(i) * Xs(i) / sumc
       mu = mu + Eta(i, i) * Xs(i) / sumv
    end do

    cp = R0 * sum(Xs(1:cea%Nm) * Cprr(1:cea%Nm)) / sum(Xs(1:cea%Nm) * Wmol(1:cea%Nm))

    cp = cp * 1d-3  ! [kJ/(kg·K)]
    k = k * 1d-4  ! [W/(m·K)]
    mu = mu / 10  ! [µPa·s]

    Pr = 1d-3 * mu * cp / k

    return
  end subroutine calc_frozen_transport_properties

end module mod_ext
