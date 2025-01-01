!-------------------------------------------------------------------------------
!> @mainpage CEA: Chemical Equlibrium with Applications
!!
!-------------------------------------------------------------------------------
module mod_cea
  implicit none

contains

  subroutine run_all_cases(cea, out_filename, plt_filename)
    use mod_constants, only: stderr
    use mod_types, only: CEA_Problem, MAX_FILENAME, IOOUT
    use mod_legacy_io
    implicit none

    type(CEA_Problem), intent(inout):: cea(:)
    character(MAX_FILENAME), intent(in), optional:: out_filename
    character(MAX_FILENAME), intent(in), optional:: plt_filename
    integer:: icase, num_cases

    num_cases = size(cea)

    do icase = 1, num_cases
       call read_libraries(cea(icase))
    end do

    if (present(out_filename)) then
       call open_legacy_output(IOOUT, out_filename)
    else
       call open_legacy_output(IOOUT)
    end if

    do icase = 1, num_cases
       !! TEMPORARY WORK AROUND TO REPRODUCE KNOWN BUG !!
       if (icase >= 2) then
          cea(icase)%Mu(:) = cea(icase-1)%Mu(:)
       end if
       !!!!!!!!!!!!!!!!! TO BE DELETED !!!!!!!!!!!!!!!!!!

       call write_input_log(cea(icase)%io_log, IOOUT)

       if (cea(icase)%invalid_case) cycle

       if (cea(icase)%Rkt) then
          call ROCKET(cea(icase))
       else if (cea(icase)%Tp .or. cea(icase)%Hp .or. cea(icase)%Sp) then
          call THERMP(cea(icase))
       else if (cea(icase)%Detn) then
          call DETON(cea(icase))
       else if (cea(icase)%Shock) then
          call SHCK(cea(icase))
       end if

#ifndef NDEBUG
       write(stderr, *)
#endif
    end do

    close(IOOUT)

    if (present(plt_filename)) then
       call write_plt_file(cea(1:num_cases), plt_filename)
    end if

    return
  end subroutine run_all_cases


  subroutine CPHS(cea)
    !***********************************************************************
    ! CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
    !***********************************************************************
    use mod_types
    use mod_functions
    implicit none

    type(CEA_Problem), intent(inout):: cea

    integer:: ij, j, jj, k

    k = 1
    if (cea%Tt > cea%Tg(2)) k = 2
    if (cea%Tt > cea%Tg(3)) k = 3

    do concurrent (j = 1:cea%Ng)
       cea%S(j) = calc_entropy(cea%Coef(:, j, k), cea%Tt, cea%Tln, .true.)
       cea%H0(j) = calc_enthalpy(cea%Coef(:, j, k), cea%Tt, cea%Tln, .true.)
    end do

    if (.not. cea%Tp .or. cea%Convg) then
       do concurrent (j = 1:cea%Ng)
          cea%Cp(j) = calc_specific_heat(cea%Coef(:, j, k), cea%Tt)
       end do
    end if

    if (cea%Npr /= 0 .and. k /= 3 .and. cea%Ng /= cea%Ngc) then
       do concurrent (ij = 1:cea%Npr)
          j = cea%Jcond(ij)
          jj = cea%Jcond(ij) - cea%Ng

          cea%S(j) = calc_entropy(cea%Cft(:, jj), cea%Tt, cea%Tln, .true.)
          cea%H0(j) = calc_enthalpy(cea%Cft(:, jj), cea%Tt, cea%Tln, .true.)
          cea%Cp(j) = calc_specific_heat(cea%Cft(:, jj), cea%Tt)
       end do
    end if

    return
  end subroutine CPHS


  subroutine ALLCON(cea)
    !***********************************************************************
    ! CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
    !***********************************************************************
    use mod_types
    use mod_functions
    implicit none

    type(CEA_Problem), intent(inout):: cea

    integer:: j, jj

    do concurrent (jj = 1:cea%Nc)
       j = jj + cea%Ng
       cea%S(j) = calc_entropy(cea%Cft(:, jj), cea%Tt, cea%Tln, .true.)
       cea%H0(j) = calc_enthalpy(cea%Cft(:, jj), cea%Tt, cea%Tln, .true.)
       cea%Cp(j) = calc_specific_heat(cea%Cft(:, jj), cea%Tt)
    end do

    return
  end subroutine ALLCON



  subroutine DETON(cea)
    !***********************************************************************
    ! CHAPMAN-JOUGUET DETONATIONS.
    !***********************************************************************
    use mod_types
    use mod_legacy_io
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(15):: ft1  = 'T1, K'
    character(15):: fh1  = 'H1, CAL/G'
    character(15):: fhs1 = 'H1, KJ/KG'
    character(15):: fm1  = 'M1, (1/n) '
    character(15):: fg1  = 'GAMMA1'
    character(15):: fpp1 = 'P/P1'
    character(15):: ftt1 = 'T/T1'
    character(15):: fmm1 = 'M/M1'
    character(15):: frr1 = 'RHO/RHO1'
    character(15):: fdv  = 'DET VEL,M/SEC'
    character(3):: unit
    integer:: i, ii, iof, itr, j, mdv, mgam, mh, mmach, mp, mson, mt
    integer:: ip, it
    real(8):: a11, a12, a21, a22, alam, alpha, amm, b1, b2, d, gam, &
         p1, pp1, rk, rr1, T1, tem, tt1, ud, x1, x2

    integer:: i_col_start
    real(8), dimension(cea%max_points):: cpl, gm1, h1, pub, rrho, tub

    cea%Eql = .true.

    if (cea%T(1) == 0) then
       cea%T(1) = cea%Rtemp(1)
       cea%Nt = 1
    end if

    iofLoop: do iof = 1, cea%Nof
       cea%Tt = cea%T(1)

       cea%Oxfl = cea%Oxf(iof)

       call NEWOF(cea)

       ! BEGIN T LOOP.
       do it = 1, cea%Nt
          T1 = cea%T(it)

          ! BEGIN P LOOP.
          do ip = 1, cea%Np
             cea%ipt = (it - 1) * cea%Np + ip
             cea%Npt = mod(cea%ipt - 1, Ncol) + 1

             i_col_start = int((cea%ipt - 1) / Ncol) * Ncol + 1

             if (cea%ipt > 1) then
                if (cea%Isv == 1) cea%Isv = -1
                call SETEN(cea)
             end if

             p1 = cea%P(ip)
             cea%Tt = T1
             cea%Pp = p1

             call HCALC(cea)

             if (cea%Tt == 0) return
             if (cea%Detdbg) call OUT1(cea)

             h1(cea%ipt) = cea%Hsub0 * cea%R

             tub(cea%ipt) = T1
             pub(cea%ipt) = p1
             cpl(cea%ipt) = cea%Cpmix * cea%R
             cea%Tt = 3800
             pp1 = 15
             cea%Pp = pp1 * p1

             ! CALCULATE ENTHALPY FOR INITIAL ESTIMATE OF T2(TT AFTER EQLBRM)
             cea%Hsub0 = h1(cea%ipt) / cea%R + 0.75 * T1 * pp1 / cea%Wmix
             cea%Tp = .false.
             cea%Hp = .true.

             call EQLBRM(cea)

             cea%Hsub0 = h1(cea%ipt) / cea%R
             cea%Hp = .false.

             if (cea%Tt /= 0) then
                gam = cea%Gammas(cea%Npt)
                tt1 = cea%Tt / T1
                ii = 0
                tem = tt1 - 0.75 * pp1 / (cea%Cpr(cea%Npt) * cea%Wmix)
                amm = cea%Wm(cea%Npt) / cea%Wmix

                if (cea%Detdbg) write(IOOUT, '(/" T EST.=", F8.2/11X, "P/P1", 17X, "T/T1")') cea%Tt

                ! LOOP FOR IMPROVING T2/T1 AND P2/P1 INITIAL ESTIMATE.
                do ii = 1, 3
                   alpha = amm / tt1
                   pp1 = (1 + gam) * (1 + sqrt(1 - 4 * gam * alpha / (1 + gam)**2)) / (2 * gam * alpha)
                   rk = pp1 * alpha
                   tt1 = tem + 0.5 * pp1 * gam * (rk**2 - 1) / (cea%Wmix * cea%Cpr(cea%Npt) * rk)
                   if (cea%Detdbg) write(IOOUT, '(i5, 2e20.8)') ii, pp1, tt1
                end do

                cea%Tp = .true.
                cea%Tt = T1 * tt1
                rr1 = pp1 * amm / tt1

                ! BEGIN MAIN ITERATION LOOP.
                do itr = 1, 8
                   cea%Pp = p1 * pp1

                   call EQLBRM(cea)

                   if (cea%Npt == 0) exit iofLoop

                   if (cea%Tt == 0) exit

                   gam = cea%Gammas(cea%Npt)
                   amm = cea%Wm(cea%Npt) / cea%Wmix
                   rr1 = pp1 * amm / tt1
                   a11 = 1 / pp1 + gam * rr1 * cea%Dlvpt(cea%Npt)
                   a12 = gam * rr1 * cea%Dlvtp(cea%Npt)
                   a21 = 0.5 * gam * (rr1**2 - 1 - cea%Dlvpt(cea%Npt) * (1 + rr1**2)) + cea%Dlvtp(cea%Npt) - 1
                   a22 = -0.5 * gam * cea%Dlvtp(cea%Npt) * (rr1**2 + 1) - cea%Wm(cea%Npt) * cea%Cpr(cea%Npt)
                   b1 = 1 / pp1 - 1 + gam * (rr1 - 1)
                   b2 = cea%Wm(cea%Npt) * (cea%Hsum(cea%Npt) - h1(cea%ipt) / cea%R) / cea%Tt - 0.5 * gam * (rr1**2 - 1)
                   d = a11 * a22 - a12 * a21
                   x1 = (a22 * b1 - a12 * b2) / d
                   x2 = (a11 * b2 - a21 * b1) / d
                   alam = 1
                   tem = x1

                   if (tem < 0) tem = -tem
                   if (x2 > tem) tem = x2
                   if (-x2 > tem) tem = -x2
                   if (tem > 0.4054652) alam = 0.4054652 / tem

                   pp1 = pp1 * exp(x1 * alam)
                   tt1 = tt1 * exp(x2 * alam)

                   cea%Tt = T1 * tt1
                   ud = rr1 * sqrt(R0 * gam * cea%Tt/cea%Wm(cea%Npt))

                   if (cea%Detdbg) write(IOOUT, '(/" ITER =", i2, 5x, "P/P1 =", e15.8, /7x, "T/T1 =", e15.8, 5x, &
                        & "RHO/RHO1 =", e15.8, /7x, "DEL LN P/P1 =", e15.8, 5x, &
                        & "DEL LN T/T1 =", e15.8)') itr, pp1, tt1, rr1, x1, x2

                   ! CONVERGENCE TEST
                   if (tem <= 0.5E-04) exit
                end do

                if (cea%Tt /= 0) then
                   if (itr < 8) then
                      rrho(cea%ipt) = rr1
                      if (cpl(cea%ipt) == 0) then
                         gm1(cea%ipt) = 0
                         cea%Vmoc(cea%Npt) = 0
                      else
                         gm1(cea%ipt) = cpl(cea%ipt) / (cpl(cea%ipt) - cea%R / cea%Wmix)
                         cea%Vmoc(cea%Npt) = ud / sqrt(R0 * gm1(cea%ipt) * T1 / cea%Wmix)
                      end if
                   else
                      write(IOOUT, '(/" CONSERVATION EQNS NOT SATISFIED IN 8 ITERATIONS (DETON)")')
                      cea%Tt = 0
                   end if

                   if (cea%Trnspt) call TRANP(cea)

                   cea%Isv = 0

                   if (ip /= cea%Np .or. it /= cea%Nt .and. cea%Tt /= 0) then
                      cea%Isv = cea%Npt
                      if (cea%Npt /= Ncol) cycle
                   end if
                end if

                ! OUTPUT
                write(IOOUT, '(//, 21X, "DETONATION PROPERTIES OF AN IDEAL REACTING GAS")')
                call OUT1(cea)

                ! SET MXX ARRAY FOR PLOTTING PARAMETERS
                mp    = 0
                mt    = 0
                mgam  = 0
                mh    = 0
                mdv   = 0
                mson  = 0
                mmach = 0

                do i = 1, cea%Nplt
                   if (index(cea%Pltvar(i)(2:), '1') /= 0) then
                      if (cea%Pltvar(i)(:3) == 'son') then
                         mson = i
                      else if (cea%Pltvar(i)(:3) == 'gam') then
                         mgam = i
                      else if (cea%Pltvar(i)(:1) == 'h') then
                         mh = i
                      else if (cea%Pltvar(i)(:1) == 't') then
                         mt = i
                      else if (cea%Pltvar(i)(:1) == 'p') then
                         mp = i
                      end if
                   else if (index(cea%Pltvar(i), 'vel') /= 0) then
                      mdv = i
                   else if (index(cea%Pltvar(i), 'mach') /= 0) then
                      mmach = i
                   end if
                end do

                write(IOOUT, '(/" UNBURNED GAS"/)')

                cea%fmt(4) = '13'
                cea%fmt(5) = ' '
                cea%fmt(7) = '4,'

                do i = i_col_start, cea%ipt
                   j = i - i_col_start + 1
                   if (cea%SIunit) then
                      cea%V(j) = pub(i)
                      unit = 'BAR'
                   else
                      cea%V(j) = pub(i) / 1.01325d0
                      unit = 'ATM'
                   end if
                   if (mp > 0) cea%Pltout(j+cea%Iplt, mp) = cea%V(j)
                end do

                write(IOOUT, cea%fmt) 'P1, ' // unit // '        ', (cea%V(j), j = 1, cea%Npt)

                cea%fmt(7) = '2,'
                write(IOOUT, cea%fmt) ft1, (tub(i), i = i_col_start, cea%ipt)

                if (.not. cea%SIunit) write(IOOUT, cea%fmt) fh1, (h1(i), i = i_col_start, cea%ipt)
                if (cea%SIunit) write(IOOUT, cea%fmt) fhs1, (h1(i), i = i_col_start, cea%ipt)

                do concurrent (i = i_col_start:cea%ipt)
                   j = i - i_col_start + 1
                   cea%V(j) = cea%Wmix
                   cea%Sonvel(j) = sqrt(R0 * gm1(i) * tub(i) / cea%Wmix)
                end do

                cea%fmt(7) = '3,'
                write(IOOUT, cea%fmt) fm1, (cea%V(j), j = 1, cea%Npt)
                cea%fmt(7) = '4,'
                write(IOOUT, cea%fmt) fg1, (gm1(i), i = i_col_start, cea%ipt)
                cea%fmt(7) = '1,'
                write(IOOUT, cea%fmt) 'SON VEL1,M/SEC ', (cea%Sonvel(j), j = 1, cea%Npt)

                if (cea%Nplt > 0) then
                   do i = i_col_start, cea%ipt
                      j = i - i_col_start + 1
                      if (mt > 0)   cea%Pltout(j+cea%Iplt, mt) = tub(i)
                      if (mgam > 0) cea%Pltout(j+cea%Iplt, mgam) = gm1(i)
                      if (mh > 0)   cea%Pltout(j+cea%Iplt, mh) = h1(i)
                      if (mson > 0) cea%Pltout(j+cea%Iplt, mson) = cea%Sonvel(j)
                   end do
                end if

                write(IOOUT, '(/" BURNED GAS"/)')

                cea%fmt(4) = cea%fmt(6)
                call OUT2(cea, cea%Npt)

                if (cea%Trnspt) call OUT4(cea, cea%Npt)

                write(IOOUT, '(/" DETONATION PARAMETERS"/)')

                cea%fmt(7) = '3,'

                do i = i_col_start, cea%ipt
                   j = i - i_col_start + 1
                   cea%V(j) = cea%Ppp(j) / pub(i)
                   cea%Pcp(j) = cea%Ttt(j) / tub(i)
                   cea%Sonvel(j) = cea%Sonvel(j) * rrho(i)
                   if (mmach > 0) cea%Pltout(j+cea%Iplt, mmach) = cea%Vmoc(j)
                   if (mdv > 0)   cea%Pltout(j+cea%Iplt, mdv) = cea%Sonvel(j)
                end do

                write(IOOUT, cea%fmt) fpp1, (cea%V(j), j = 1, cea%Npt)
                write(IOOUT, cea%fmt) ftt1, (cea%Pcp(j), j = 1, cea%Npt)

                do concurrent (i = 1:cea%Npt)
                   cea%V(i) = cea%Wm(i) / cea%Wmix
                end do

                cea%fmt(7) = '4,'
                write(IOOUT, cea%fmt) fmm1, (cea%V(j), j = 1, cea%Npt)
                write(IOOUT, cea%fmt) frr1, (rrho(i), i = i_col_start, cea%ipt)
                write(IOOUT, cea%fmt) 'DET MACH NUMBER', (cea%Vmoc(j), j = 1, cea%Npt)

                cea%fmt(7) = '1,'
                write(IOOUT, cea%fmt) fdv, (cea%Sonvel(j), j = 1, cea%Npt)

                cea%Eql = .true.

                call OUT3(cea, cea%Npt)

                cea%Iplt = min(cea%Iplt+cea%Npt, 500)

                if (cea%Isv == 0 .and. iof == cea%Nof) exit iofLoop
                if (cea%Np == 1 .and. cea%Nt == 1) cycle iofLoop

                write(IOOUT, '(///)')
             end if
          end do
       end do

       cea%Iplt = min(cea%Iplt + cea%Npt - 1, 500)

    end do iofLoop

    cea%Tp = .false.

    return
  end subroutine DETON


  subroutine EQLBRM(cea)
    !***********************************************************************
    ! CALCULATE EQUILIBRIUM COMPOSITION AND PROPERTIES.
    !***********************************************************************
    use mod_types
    use mod_general
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(12):: ae, cmp(maxEl)
    character(16):: amb
    logical, save:: cpcalc, i2many, newcom, reduce
    integer:: i, il, ilamb, ilamb1, inc, ipr, iq2, iter, ixsing, iz, j, ja, jb, &
         jbx, jc, jcondi, jcons, jdelg, jkg, jneg, jsw, k, kc, kg, kk, &
         kmat, kneg, l, lc, lcs(maxEl), le, lelim, lk, ll, lncvg, ls, lsing, &
         lz, maxitn, ncvg, njc, nn, numb
    real(8), save:: aa, ambda, ambda1, bigen, bigneg, delg, dlnt, dpie, ensol, esize, &
         gap, gasfrc, pie, pisave(maxMat-2), siz9, sizeg, &
         sum0, sum1, szgj, tem, tmelt, tsize, ween, xi, xln, xsize, xx(maxMat)
    real(8):: smalno = 1e-6, smnol = -13.815511
    logical:: mask(cea%Ng)

    type(CEA_Point), pointer:: p !< current point
    p => cea%points(cea%iOF, cea%Npt)

    ixsing = 0
    lsing = 0
    jsw = 0
    jdelg = 0
    maxitn = 50
    ncvg = 0
    lncvg = 3 * cea%Nlm
    reduce = .false.
    siz9 = cea%Size - 9.2103404d0
    tsize = cea%Size
    xsize = cea%Size + 6.90775528d0

    if (cea%Trace /= 0.) then
       maxitn = maxitn + cea%Ngc / 2
       xsize = -log(cea%Trace)
       if (xsize < cea%Size) xsize = cea%Size + 0.1
    end if

    if (xsize > 80) xsize = 80

    esize = min(80d0, xsize + 6.90775528d0)
    jcons = 0
    pie = 0
    i2many = .false.
    cea%Pderiv = .false.
    cea%Convg = .false.
    numb = 0
    cpcalc = .true.

    if (cea%Tp) cpcalc = .false.

    if (cea%Tt /= 0) then
       if (cea%Npr == 0 .or. (cea%Tt /= cea%T(1) .and. .not. cea%Tp)) go to 400
       k = 1
    else
       cea%Tt = 3800
       if (cea%Npr == 0) go to 400
       k = 1
    end if

100 j = cea%Jcond(k)
    jc = j - cea%Ng
    kg = -cea%Ifz(jc)

    do i = 1, 9
       kg = kg + 1
       kc = jc + kg

       if (cea%Tt <= cea%Temp(2, kc)) then
          if (kg /= 0) then
             cea%Jcond(k) = j + kg
             cea%En(j+kg, cea%Npt) = cea%En(j, cea%Npt)
             cea%En(j, cea%Npt) = 0

             if (cea%Prod(j) /= cea%Prod(j+kg) .and. .not. cea%Short) &
                  & write(IOOUT, '(" PHASE CHANGE, REPLACE ", A16, "WITH ", A16)') cea%Prod(j), cea%Prod(j+kg)
          end if
          go to 300
       else if (kc >= cea%Nc .or. cea%Ifz(kc+1) <= cea%Ifz(kc)) then
          exit
       end if
    end do

    if (.not. cea%Tp) then
       cea%Tt = cea%Temp(2, kc) - 10
       k = 1
       go to 100
    end if

    write(IOOUT, '(" REMOVE ", A16)') cea%Prod(j)

    cea%En(j, cea%Npt) = 0
    cea%Enln(j) = 0
    cea%Deln(j) = 0

    do i = k, cea%Npr
       cea%Jcond(i) = cea%Jcond(i+1)
    end do

    cea%Npr = cea%Npr - 1

300 k = k + 1

    if (k <= cea%Npr) go to 100

400 cea%Tln = log(cea%Tt)

    if (cea%Vol) cea%Pp = R0 * cea%Enn * cea%Tt / cea%Vv

    call CPHS(cea)

    cea%Tm = log(cea%Pp / cea%Enn)
    le = cea%Nlm

    if (cea%Lsave /= 0 .and. cea%Nlm /= cea%Lsave) then
       tem = exp(-tsize)
       do i = cea%Lsave + 1, cea%Nlm
          do concurrent (j = 1:cea%Ng, cea%A(i, j) /= 0)
             cea%En(j, cea%Npt) = tem
             cea%Enln(j) = -tsize
          end do
       end do
    end if

    ls = cea%Nlm
    lelim = 0
    lz = ls

    if (cea%Ions) lz = ls - 1

    if (cea%Npt == 1 .and. .not. cea%Shock .and. .not. cea%Short) then
       write(IOOUT, '(/" POINT ITN", 6X, "T", 10X, 4(A4, 8X)/(18X, 5(A4, 8X)))') (cea%Elmt(i), i = 1, cea%Nlm)
    end if

    if (cea%Debug(cea%Npt)) then
       do i = 1, cea%Nlm
          cmp(i) = cea%Elmt(i)
       end do
    end if

    ! BEGIN ITERATION
500 if (cpcalc) then
       cea%Cpsum = sum(cea%En(1:cea%Ng, cea%Npt) * cea%Cp(1:cea%Ng))

       if (cea%Npr /= 0) then
          cea%Cpsum = cea%Cpsum + sum(cea%En(cea%Jcond(1:cea%Npr), cea%Npt) * cea%Cp(cea%Jcond(1:cea%Npr)))
          cpcalc = .false.
       end if
    end if

    numb = numb + 1

    call MATRIX(cea)

    iq2 = cea%Iq1 + 1

    if (cea%Convg) cea%Imat = cea%Imat - 1

    if (cea%Debug(cea%Npt)) then
       if (.not. cea%Convg) then
          write(IOOUT, '(/" ITERATION", I3, 6X, "MATRIX ")') numb
       else
          if (.not. cea%Pderiv) write(IOOUT, '(/" T DERIV MATRIX")')
          if (cea%Pderiv) write(IOOUT, '(/" P DERIV MATRIX")')
       end if

       kmat = cea%Imat + 1

       do i = 1, cea%Imat
          write(IOOUT, '(3X, 5E15.6)') (cea%G(i, k), k = 1, kmat)
       end do
    end if

    call gauss_elimination(cea%G(1:cea%Imat, 1:cea%Imat+1), cea%X(1:cea%Imat), cea%Msing)

    if (cea%Msing == 0) then
       if (cea%Debug(cea%Npt)) then
          write(IOOUT, '(/" SOLUTION VECTOR", /, 6x, 5A15/8X, 5A15)') (cmp(k), k = 1, le)
          write(IOOUT, '(3X, 5E15.6)') (cea%X(i), i = 1, cea%Imat)
       end if

       if (.not. cea%Convg) then
          ! OBTAIN CORRECTIONS TO THE ESTIMATES
          if (cea%Vol) cea%X(iq2) = cea%X(cea%Iq1)
          if (cea%Tp) cea%X(iq2) = 0
          dlnt = cea%X(iq2)
          sum0 = cea%X(cea%Iq1)

          if (cea%Vol) then
             cea%X(cea%Iq1) = 0
             sum0 = -dlnt
          end if

          outerLoop0: do j = 1, cea%Ng
             if (lelim /= 0) then
                cea%Deln(j) = 0
                do i = lelim, ls
                   if (cea%A(i, j) /= 0) cycle outerLoop0
                end do
             end if

             cea%Deln(j) = -cea%Mu(j) + cea%H0(j) * dlnt + sum0


!!$           cea%Deln(j) = cea%Deln(j) + sum(A(1:cea%Nlm, j) * X(1:cea%Nlm)) ???
             do k = 1, cea%Nlm
                cea%Deln(j) = cea%Deln(j) + cea%A(k, j) * cea%X(k)
             end do

             if (pie /= 0) cea%Deln(j) = cea%Deln(j) + cea%A(ls, j) * pie
          end do outerLoop0

          if (cea%Npr /= 0) then
             do concurrent (k = 1:cea%Npr)
                cea%Deln(cea%Jcond(k)) = cea%X(cea%Nlm+k)
             end do
          end if

          ! CALCULATE CONTROL FACTOR, AMBDA
          ambda = 1
          ambda1 = 1
          ilamb = 0
          ilamb1 = 0
          sum0 = 5 * max(abs(cea%X(cea%Iq1)), abs(dlnt))

          do j = 1, cea%Ng
             if (cea%Deln(j) > 0) then
                if ((cea%Enln(j) - cea%Ennl + cea%Size) <= 0) then

                   sum1 = abs(cea%Deln(j)-cea%X(cea%Iq1))

                   if (sum1 >= siz9) then
                      sum1 = abs(-9.2103404d0 - cea%Enln(j) + cea%Ennl) / sum1

                      if (sum1 < ambda1) then
                         ambda1 = sum1
                         ilamb1 = j
                      end if
                   end if

                else if (cea%Deln(j) > sum0) then
                   sum0 = cea%Deln(j)
                   ilamb = j
                end if
             end if
          end do

          if (sum0 > 2) ambda = 2 / sum0

          if (ambda1 <= ambda) then
             ambda = ambda1
             ilamb = ilamb1
          end if

          if (cea%Debug(cea%Npt)) then
             ! INTERMEDIATE OUTPUT
             write(IOOUT, '(/" T=", E15.8, " ENN=", E15.8, " ENNL=", E15.8, " PP=", E15.8, &
                  & /" LN P/N=", E15.8, " AMBDA=", E15.8)') cea%Tt, cea%Enn, cea%Ennl, cea%Pp, cea%Tm, ambda

             if (ambda /= 1) then
                amb = 'ENN'
                if (abs(cea%X(iq2)) > abs(cea%X(cea%Iq1))) amb = 'TEMP'
                if (ilamb /= 0) amb = cea%Prod(ilamb)
                write(IOOUT, '(/" AMBDA SET BY ", A16)') amb
             end if

             if (cea%Vol) write(IOOUT, '(" VOLUME=", E15.8, "CC/G")') cea%Vv * 0.001d0

             write(IOOUT, '(/24X, "Nj", 12X, "LN Nj", 8X, "DEL LN Nj", 6X, "H0j/RT", /, 41X, &
                  & "S0j/R", 10X, " G0j/RT", 8X, " Gj/RT")')

             do j = 1, cea%Ngc
                write(IOOUT, '(1X, A16, 4E15.6, /35x, 3E15.6)') &
                     cea%Prod(j), cea%En(j, cea%Npt), cea%Enln(j), cea%Deln(j), cea%H0(j), cea%S(j), cea%H0(j) - cea%S(j), cea%Mu(j)
             end do
          end if

          ! APPLY CORRECTIONS TO ESTIMATES
          cea%Totn(cea%Npt) = 0

          do concurrent (j = 1:cea%Ng)
             cea%Enln(j) = cea%Enln(j) + ambda * cea%Deln(j)
          end do

          cea%En(1:cea%Ng, cea%Npt) = 0

          if (lelim == 0) then
             do concurrent (j = 1:cea%Ng, (cea%Enln(j) - cea%Ennl + tsize) > 0)
                cea%En(j, cea%Npt) = exp(cea%Enln(j))
             end do
          else
             do concurrent (j = 1:cea%Ng, all(cea%A(lelim:ls, j) == 0))
                cea%En(j, cea%Npt) = exp(cea%Enln(j))
             end do
          end if

          cea%Totn(cea%Npt) = cea%Totn(cea%Npt) + sum(cea%En(1:cea%Ng, cea%Npt))

          if (cea%Ions .and. cea%Elmt(cea%Nlm) == 'E') then
             mask = .false.
             do concurrent (j = 1:cea%Ng, cea%A(ls, j) /= 0 .and. cea%En(j, cea%Npt) == 0 .and. (cea%Enln(j) - cea%Ennl + esize) > 0)
                cea%En(j, cea%Npt) = exp(cea%Enln(j))
                mask(j) = .true.
             end do
             cea%Totn(cea%Npt) = cea%Totn(cea%Npt) + sum(cea%En(1:cea%Ng, cea%Npt), mask=mask)
          end if

          cea%Sumn = cea%Totn(cea%Npt)

          if (cea%Npr /= 0) then
             do concurrent (k = 1:cea%Npr)
                cea%En(cea%Jcond(k), cea%Npt) = cea%En(cea%Jcond(k), cea%Npt) + ambda * cea%Deln(cea%Jcond(k))
             end do
             cea%Totn(cea%Npt) = cea%Totn(cea%Npt) + sum(cea%En(cea%Jcond(1:cea%Npr), cea%Npt))
          end if

          if (.not. cea%Tp) then
             cea%Tln = cea%Tln + ambda * dlnt
             cea%Tt = exp(cea%Tln)
             cpcalc = .true.
             call CPHS(cea)
          end if

          if (cea%Vol) then
             cea%Enn = cea%Sumn
             cea%Ennl = log(cea%Enn)
             if (cea%Vol) cea%Pp = R0 * cea%Tt * cea%Enn / cea%Vv
          else
             cea%Ennl = cea%Ennl + ambda * cea%X(cea%Iq1)
             cea%Enn = exp(cea%Ennl)
          end if

          cea%Tm = log(cea%Pp / cea%Enn)

          if (cea%Elmt(cea%Nlm) == 'E') then
             ! CHECK ON REMOVING IONS
             if (all(cea%A(cea%Nlm, 1:cea%Ngc) == 0 .or. cea%En(1:cea%Ngc, cea%Npt) <= 0)) then
                pie = cea%X(cea%Nlm)
                lelim = cea%Nlm
                cea%Nlm = cea%Nlm - 1
                go to 500
             end if
          end if

          ! TEST FOR CONVERGENCE
          if (numb > maxitn) then
             write(IOOUT, '(/, I4, " ITERATIONS DID NOT SATISFY CONVERGENCE", /, 15x, &
                  & " REQUIREMENTS FOR THE POINT", I5, " (EQLBRM)")') maxitn, cea%Npt

             if (cea%Nc == 0 .or. i2many) go to 1500

             i2many = .true.

             if (.not. cea%Hp .or. cea%Npt /= 1 .or. cea%Tt > 100.) then
                if (cea%Npr /= 1 .or. cea%Enn > 1.E-4) go to 1500
                ! HIGH TEMPERATURE, INCLUDED CONDENSED CONDITION
                write(IOOUT, '(/" TRY REMOVING CONDENSED SPECIES (EQLBRM)")')

                cea%Enn = 0.1
                cea%Ennl = -2.3025851
                cea%Sumn = cea%Enn
                xi = cea%Ng
                xi = cea%Enn/xi
                xln = log(xi)

                do concurrent (j = 1:cea%Ng)
                   cea%En(j, cea%Npt) = xi
                   cea%Enln(j) = xln
                end do

                j = cea%Jcond(1)
                k = 1
                go to 1000

             else
                write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", &
                     & "VE BEEN INSERTED,", &
                     & /" RESTART WITH insert DATASET (EQLBRM)")')
                go to 1500
             end if

          else
             if (abs(cea%X(cea%Iq1) * cea%Enn / cea%Totn(cea%Npt)) > 0.5E-5 .or. &
                  any(abs(cea%Deln(1:cea%Ng)) * cea%En(1:cea%Ng, cea%Npt) / cea%Totn(cea%Npt) > 0.5d-5) .or. &
                  abs(dlnt) > 1.d-04) go to 500

             if (cea%Npr /= 0) then
                do k = 1, cea%Npr
                   j = cea%Jcond(k)
                   if (abs(cea%Deln(j)/cea%Totn(cea%Npt)) > 0.5d-5) go to 500
                   if (cea%En(j, cea%Npt) < 0) go to 700
                end do
             end if

             le = cea%Nlm

             if (any(abs(p%B0(1:cea%Nlm)) >= 1.d-06 .and. &
                  abs(p%B0(1:cea%Nlm) - sum(spread(cea%En(1:cea%Ngc, cea%Npt), 1, cea%Nlm) * cea%A(1:cea%Nlm, 1:cea%Ngc), DIM=2)) &
                  > cea%Bcheck)) go to 500

             if (cea%Trace /= 0.) then
                tsize = xsize
                tem = 1

                if (numb /= 1) then
                   lk = lz
                   if (cea%Nlm < lz) lk = cea%Nlm
                   do i = 1, lk
                      if (i /= lsing) then
                         tem = 0
                         if (cea%X(i) /= 0) then
                            tem = abs((pisave(i) - cea%X(i)) / cea%X(i))
                            if (tem > 0.001) exit
                         end if
                      end if
                   end do
                end if

                do concurrent (i = 1:cea%Nlm)
                   pisave(i) = cea%X(i)
                end do

                if (tem > 0.001) go to 500

                if (cea%Ions) then
                   ! CHECK ON ELECTRON BALANCE
                   iter = 1

                   if (pie /= 0) then
                      cea%X(cea%Nlm+1) = pie
                   end if

566                sum1 = 0
                   sum0 = 0
                   pie = cea%X(le)

                   do j = 1, cea%Ng
                      if (cea%A(ls, j) /= 0) then
                         cea%En(j, cea%Npt) = 0
                         tem = 0

                         if (cea%Enln(j) > -87) tem = exp(cea%Enln(j))

                         if ((cea%Enln(j)-cea%Ennl+tsize) > 0 .and. cea%Elmt(cea%Nlm) == 'E') then
                            pie = 0
                            cea%En(j, cea%Npt) = tem
                         end if

                         aa = cea%A(ls, j) * tem
                         sum0 = sum0 + aa
                         sum1 = sum1 + aa * cea%A(ls, j)
                      end if
                   end do

                   if (sum1 /= 0) then
                      dpie = -sum0 / sum1

                      do concurrent (j = 1:cea%Ng, cea%A(ls, j) /= 0)
                         cea%Enln(j) = cea%Enln(j) + cea%A(ls, j) * dpie
                      end do

                      if (cea%Debug(cea%Npt)) write(IOOUT, '(/" ELECTRON BALANCE ITER NO. =", i4, "  DELTA PI =", e14.7)') iter, dpie

                      if (abs(dpie) > 0.0001) then
                         cea%X(le) = cea%X(le) + dpie
                         iter = iter + 1

                         if (iter <= 80) go to 566
                         write(IOOUT, '(/" DID NOT CONVERGE ON ELECTRON BALANCE (EQLBRM)")')
                         go to 1500

                      else if (cea%Elmt(cea%Nlm) == 'E' .and. pie /= 0) then
                         cea%Nlm = cea%Nlm - 1
                         newcom = .true.
                      end if
                   end if
                end if
             end if
          end if

       else if (.not. cea%Pderiv) then
          ! TEMPERATURE DERIVATIVES--CONVG=T, PDERIV=F
          cea%Dlvtp(cea%Npt) = 1. - cea%X(cea%Iq1)
          cea%Cpr(cea%Npt) = cea%G(iq2, iq2)

          cea%Cpr(cea%Npt) = cea%Cpr(cea%Npt) - sum(cea%G(iq2, 1:cea%Iq1) * cea%X(1:cea%Iq1))

          ! PRESSURE DERIVATIVE--CONVG=T, PDERIV=T
          cea%Pderiv = .true.
          go to 500

       else
          cea%Dlvpt(cea%Npt) = -1 + cea%X(cea%Iq1)
          if (cea%Jliq == 0) then
             cea%Gammas(cea%Npt) = -1 / (cea%Dlvpt(cea%Npt) + (cea%Dlvtp(cea%Npt)**2) * cea%Enn / cea%Cpr(cea%Npt))
          else
             cea%En(cea%Jsol, cea%Npt) = ensol
             cea%Hsum(cea%Npt) = cea%Hsum(cea%Npt) + cea%En(cea%Jliq, cea%Npt) * (cea%H0(cea%Jliq) - cea%H0(cea%Jsol))
             cea%Gammas(cea%Npt) = -1. / cea%Dlvpt(cea%Npt)
             cea%Npr = cea%Npr + 1
             cea%Jcond(cea%Npr) = cea%Jliq
          end if
          go to 1400
       end if

       ! SINGULAR MATRIX
    else
       if (cea%Convg) then
          write(IOOUT, '(/" DERIVATIVE MATRIX SINGULAR (EQLBRM)")')
          cea%Dlvpt(cea%Npt) = -1
          cea%Dlvtp(cea%Npt) = 1
          cea%Cpr(cea%Npt) = cea%Cpsum
          cea%Gammas(cea%Npt) = -1 / (cea%Dlvpt(cea%Npt) + cea%Dlvtp(cea%Npt)**2 * cea%Enn / cea%Cpr(cea%Npt))
          go to 1400

       else
          write(IOOUT, '(/" SINGULAR MATRIX, ITERATION", I3, "  VARIABLE", I3, "(EQLBRM)")') numb, cea%Msing
          lsing = cea%Msing
          ixsing = ixsing + 1
          if (ixsing <= 8) then
             xsize = 80
             tsize = xsize
             if (cea%Msing > cea%Nlm .and. numb < 1 .and. cea%Npr > 1 .and. jdelg > 0) then
                ween = 1000
                j = 0

                do i = 1, cea%Npr
                   jcondi = cea%Jcond(i)
                   if (jcondi /= jdelg) then
                      do ll = 1, cea%Nlm
                         if (cea%A(ll, jdelg) /= 0 .and. cea%A(ll, jcondi) /= 0) then
                            if (cea%En(jcondi, cea%Npt) <= ween) then
                               ween = cea%En(jcondi, cea%Npt)
                               j = jcondi
                               k = i
                            end if
                            exit
                         end if
                      end do
                   end if
                end do

                if (j > 0) then
                   write(IOOUT, '(/" TRY REMOVING CONDENSED SPECIES (EQLBRM)")')
                   go to 1000
                end if

             else if (.not. cea%Hp .or. cea%Npt /= 1 .or. cea%Nc == 0 .or. cea%Tt > 100) then
                if (ixsing >= 3) then
                   if (cea%Msing < cea%Iq1) then
                      if (reduce .and. cea%Msing <= cea%Nlm) then
                         if (cea%Nlm < lelim) go to 1500
                         write(IOOUT, '(/" WARNING!! POINT", I3, &
                              & " USES A REDUCED SET OF COMPONENTS", / &
                              & " SPECIES CONTAINING THE ELIMINATED COMPONENT ARE OMITTED.", &
                              & / &
                              & " IT MAY BE NECESSARY TO RERUN WITH INSERTED CONDENSED SPECIES", &
                              & /" CONTAINING COMPONENT ", A8, "(EQLBRM)")') cea%Npt, cea%Elmt(cea%Nlm)
                         cea%Nlm = cea%Nlm - 1
                         go to 500

                      else if (cea%Msing <= cea%Nlm) then
                         ! FIND NEW COMPONENTS
                         if (.not. cea%Ions) go to 1100
                         if (cea%Elmt(cea%Nlm) /= 'E') go to 1100

                         do concurrent (j = 1:cea%Ng, cea%A(cea%Nlm, j) /= 0)
                            cea%En(j, cea%Npt) = 0
                         end do

                         pie = cea%X(cea%Nlm)
                         cea%Nlm = cea%Nlm - 1
                         if (cea%Msing > cea%Nlm) go to 500
                         go to 1100
                      else
                         ! REMOVE CONDENSED SPECIES TO CORRECT SINGULARITY
                         k = cea%Msing - cea%Nlm
                         j = cea%Jcond(k)

                         if (j /= jcons) then
                            jcons = j
                            go to 1000
                         end if
                      end if
                   end if
                end if

                do concurrent (j = 1:cea%Ng, &
                     .not. (cea%Ions .and. cea%Elmt(cea%Nlm) /= 'E' .and. cea%A(ls, j) /= 0) .and. cea%En(j, cea%Npt) == 0)
                   cea%En(j, cea%Npt) = smalno
                   cea%Enln(j) = smnol
                end do

                go to 500
             else
                write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", &
                     & "VE BEEN INSERTED,", &
                     & /" RESTART WITH insert DATASET (EQLBRM)")')
             end if
          end if
       end if

       go to 1500
    end if

    ! CALCULATE ENTROPY, CHECK ON DELTA S FOR SP PROBLEMS
600 cea%Ssum(cea%Npt) = sum(cea%En(1:cea%Ng, cea%Npt) * (cea%S(1:cea%Ng) - cea%Enln(1:cea%Ng) - cea%Tm))

    if (cea%Npr > 0) then
       cea%Ssum(cea%Npt) = cea%Ssum(cea%Npt) + sum(cea%En(cea%Jcond(1:cea%Npr), cea%Npt) * cea%S(cea%Jcond(1:cea%Npr)))
    end if

    if (.not. cea%Sp) then
       cea%Convg = .true.
    else
       tem = cea%Ssum(cea%Npt) - cea%S0
       if (abs(tem) > 0.0005) go to 500
       if (cea%Debug(cea%Npt)) write(IOOUT, '(/" DELTA S/R =", e15.8)') tem
       cea%Convg = .true.
    end if

    ! CONVERGENCE TESTS ARE SATISFIED, TEST CONDENSED SPECIES.
700 ncvg = ncvg + 1

    if (ncvg > lncvg) then
       ! ERROR, SET TT=0
       write(IOOUT, '(/, I3, " CONVERGENCES FAILED TO ESTABLISH SET OF CONDENSED", " SPECIES (EQLBRM)")') lncvg
       go to 1500
    else
       if (.not. cea%Shock) then
          do concurrent (il = 1:le)
             xx(il) = cea%X(il)
          end do

          if (.not. cea%Short) then
             if (newcom) write(IOOUT, '(/" POINT ITN", 6x, "T", 10x, 4a12/(18x, 5a12))') (cmp(k), k = 1, le)
             write(IOOUT, '(i4, i5, 5f12.3, /(12x, 5f12.3))') cea%Npt, numb, cea%Tt, (xx(il), il = 1, le)
          end if

          if (.not. cea%Tp .and. cea%Npr == 0 .and. cea%Tt <= cea%Tg(1) * 0.2d0) then
             write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", "VE BEEN INSERTED,", &
                  & /" RESTART WITH insert DATASET (EQLBRM)")')
             go to 1500
          end if

          newcom = .false.
       end if

       if (cea%Npr /= 0) then
          bigneg = 0
          jneg = 0

          do k = 1, cea%Npr
             j = cea%Jcond(k)
             if (cea%En(j, cea%Npt) * cea%Cp(j) <= bigneg) then
                bigneg = cea%En(j, cea%Npt) * cea%Cp(j)
                jneg = j
                kneg = k
             end if
          end do

          if (jneg /= 0) then
             j = jneg
             k = kneg
             if (j == cea%Jsol .or. j == cea%Jliq) then
                cea%Jsol = 0
                cea%Jliq = 0
             end if
             go to 1000
          end if
       end if

       if (cea%Ngc /= cea%Ng .or. cea%Tp) then
          cea%Ng = cea%Ngc

          call CPHS(cea)

          cea%Ng = cea%Ngp1 - 1
          cpcalc = .true.

          if (cea%Ngc == cea%Ng) go to 750

          call ALLCON(cea)

          if (cea%Npr /= 0 .and. .not. cea%Tp) then
             gap = 50

             outerLoop1: do ipr = 1, cea%Npr
                j = cea%Jcond(ipr)
                if (j /= cea%Jsol .and. j /= cea%Jliq) then
                   inc = j - cea%Ng
                   kg = -cea%Ifz(inc)

                   do iz = 1, 20
                      kg = kg + 1
                      kc = inc + kg

                      if (cea%Tt <= cea%Temp(2, kc)) then
                         if (kg /= 0) then
                            jkg = j + kg
                            if (abs(kg) > 1 .or. cea%Prod(j) == cea%Prod(jkg)) go to 740
                            if (jkg == jsw) go to 720
                            if (cea%Tt < cea%Temp(1, inc) - gap .or. cea%Tt > cea%Temp(2, inc) + gap) go to 740
                            go to 720
                         end if
                         cycle outerLoop1
                      else if (cea%Ifz(kc+1) <= cea%Ifz(kc)) then
                         cycle outerLoop1
                      end if
                   end do
                   if (cea%Tt > cea%Temp(2, kc) * 1.2d0) go to 1000
                end if
             end do outerLoop1
          end if

          sizeg = 0
          szgj = 0

          do inc = 1, cea%Nc
             j = inc + cea%Ng

             if (cea%Debug(cea%Npt)) then
                write(IOOUT, '(/1x, a15, 2f10.3, 3x, e15.7)') cea%Prod(j), cea%Temp(1, inc), cea%Temp(2, inc), cea%En(j, cea%Npt)
             end if

             if (cea%En(j, cea%Npt) <= 0) then
                if (cea%Tt > cea%Temp(1, inc) .or. cea%Temp(1, inc) == cea%Tg(1)) then
                   if (cea%Tt <= cea%Temp(2, inc)) then
                      delg = (cea%H0(j) - cea%S(j) - sum(cea%A(1:cea%Nlm, j) * cea%X(1:cea%Nlm))) / cea%Mw(j)

                      if (delg < sizeg .and. delg < 0) then
                         if (j /= jcons) then
                            sizeg = delg
                            jdelg = j
                         else
                            szgj = delg
                         end if
                         ipr = ipr - 1
                      end if

                      if (cea%Debug(cea%Npt)) then
                         write(IOOUT, '(" [G0j-SUM(Aij*PIi)]/Mj =", E15.7, 9X, "MAX NEG DELTA G =", E15.7)') delg, sizeg
                      end if
                   end if
                end if
             end if
          end do

          if (sizeg == 0 .and. szgj == 0) go to 750

          if (sizeg /= 0) then
             j = jdelg
             go to 800
          else
             write(IOOUT, '(/" REINSERTION OF ", A16, " LIKELY TO CAUSE SINGULARITY, ", "(EQLBRM)")') cea%Prod(jcons)
             go to 1500
          end if

720       kk = max(0, kg)
          tmelt = cea%Temp(kk+1, inc)
          cea%Tt = tmelt
          cea%Tln = log(cea%Tt)
          cea%Jsol = min(j, jkg)
          cea%Jliq = cea%Jsol + 1
          cea%En(jkg, cea%Npt) = 0.5d0 * cea%En(j, cea%Npt)
          cea%En(j, cea%Npt) = cea%En(jkg, cea%Npt)
          j = jkg
          go to 800

          ! WRONG PHASE INCLUDED FOR T INTERVAL, SWITCH EN
740       cea%En(jkg, cea%Npt) = cea%En(j, cea%Npt)
          cea%Jcond(ipr) = jkg
          cea%En(j, cea%Npt) = 0
          jsw = j

          if (cea%Prod(j) /= cea%Prod(jkg) .and. .not. cea%Short) &
               write(IOOUT, '(" PHASE CHANGE, REPLACE ", A16, "WITH ", A16)') cea%Prod(j), cea%Prod(jkg)

          j = jkg
          go to 900
       end if

       ! CONVERGED WITH NO CONDENSED CHANGES.  IF BOTH SOLID & LIQ PRESENT, 
       ! TEMPORARILY REMOVE LIQ TO PREVENT SINGULAR DERIVATIVE MATRIX.
750    cea%Sumn = cea%Enn
       if (cea%Jsol /= 0) then
          ensol = cea%En(cea%Jsol, cea%Npt)
          cea%En(cea%Jsol, cea%Npt) = cea%En(cea%Jsol, cea%Npt) + cea%En(cea%Jliq, cea%Npt)
          cea%Dlvtp(cea%Npt) = 0
          cea%Cpr(cea%Npt) = 0
          cea%Gammas(cea%Npt) = 0
          cea%Pderiv = .true.

          do k = 1, cea%Npr
             if (cea%Jcond(k) == cea%Jliq) exit
          end do

          cea%Jcond(k:cea%Npr) = cea%Jcond(k+1:cea%Npr+1)

          cea%Npr = cea%Npr - 1
       end if

       go to 500
    end if

    ! ADD CONDENSED SPECIES
800 cea%Npr = cea%Npr + 1

    cea%Jcond(2:cea%Npr) = cea%Jcond(1:cea%Npr-1)
    cea%Jcond(1) = j

    if (.not. cea%Short) write(IOOUT, '(" ADD ", a16)') cea%Prod(j)

900 inc = j - cea%Ng
    cea%Convg = .false.
    if (cea%Tp) cpcalc = .false.
    numb = -1
    go to 500

    ! REMOVE CONDENSED SPECIES
1000 cea%En(j, cea%Npt) = 0
    cea%Deln(j) = 0
    cea%Enln(j) = 0

    cea%Jcond(k:cea%Npr) = cea%Jcond(k+1:cea%Npr+1)

    if (.not. cea%Short) write(IOOUT, '(" REMOVE ", A16)') cea%Prod(j)

    cea%Npr = cea%Npr - 1
    do i = 1, cea%Nlm
       if (cmp(i) == cea%Prod(j)) then
          numb = -1
          cea%Convg = .false.
          if (cea%Tp) cpcalc = .false.
          go to 1100
       end if
    end do

    go to 900

1100 newcom = .false.
    nn = cea%Nlm

    if (cea%Elmt(cea%Nlm) == 'E') nn = cea%Nlm - 1

    ! FIND ORDER OF SPECIES FOR COMPONENTS - BIGGEST TO SMALLEST
    njc = 0
    lcs(1:nn) = 0

1200 bigen = -1d-35

    do j = 1, cea%Ng
       if (cea%En(j, cea%Npt) > bigen) then
          if (.not. cea%Ions .or. cea%A(ls, j) == 0) then
             bigen = cea%En(j, cea%Npt)
             jbx = j
          end if
       end if
    end do

    if (bigen > 0.) then
       do lc = 1, nn
          if (jbx == 0) jbx = cea%Jx(lc)

          if (cea%A(lc, jbx) > smalno) then
             if (njc /= 0) then
                do i = 1, njc
                   l = lcs(i)
                   if (l == lc) cycle
                   if (l == 0) exit
                   j = cea%Jcm(l)
                   if (all(cea%A(1:nn, jbx) == cea%A(1:nn, j))) cycle
                end do
             end if

             do i = 1, nn
                if (i /= lc .and. abs(cea%A(lc, jbx) * cea%A(i, cea%jx(i)) - cea%A(lc, cea%jx(i)) * cea%A(i, jbx)) <= smalno) then
                   cycle
                end if
             end do

             njc = njc + 1
             if (jbx /= cea%Jcm(lc)) newcom = .true.
             cea%Jcm(lc) = jbx
             lcs(njc) = lc
             exit
          end if
       end do

       cea%En(jbx, cea%Npt) = -cea%En(jbx, cea%Npt)
       if (njc < nn) go to 1200
    end if

    do concurrent (j = 1:cea%Ng)
       cea%En(j, cea%Npt) = abs(cea%En(j, cea%Npt))
    end do

    if (newcom) then
       ! SWITCH COMPONENTS
       do lc = 1, nn
          jb = cea%Jcm(lc)

          if (cea%A(lc, jb) == 0) then
             jb = cea%Jx(lc)
             cea%Jcm(lc) = jb
          end if

          tem = cea%A(lc, jb)

          if (tem /= 0) then
             pisave(lc) = cea%H0(jb) - cea%S(jb)

             if (jb <= cea%Ng) pisave(lc) = pisave(lc) + cea%Enln(jb) + cea%Tm
             cmp(lc) = trim(cea%Prod(jb))

             ! CALCULATE NEW COEFFICIENTS
             if (tem /= 1) then
                p%B0(lc) = p%B0(lc) / tem
                p%B0p(lc, 1) = p%B0p(lc, 1) / tem
                p%B0p(lc, 2) = p%B0p(lc, 2) / tem

                do concurrent (j = 1:cea%Nspx)
                   cea%A(lc, j) = cea%A(lc, j) / tem
                end do
             end if

             do i = 1, nn
                if (cea%A(i, jb) /= 0. .and. i /= lc) then
                   tem = cea%A(i, jb)

                   do concurrent (j = 1:cea%Nspx)
                      cea%A(i, j) = cea%A(i, j) - cea%A(lc, j) * tem
                   end do

                   do concurrent (j = 1:cea%Nspx, abs(cea%A(i, j)) < 1.E-5)
                      cea%A(i, j) = 0
                   end do

                   p%B0(i) = p%B0(i) - p%B0(lc) * tem
                   p%B0p(i, 1) = p%B0p(i, 1) - p%B0p(lc, 1) * tem
                   p%B0p(i, 2) = p%B0p(i, 2) - p%B0p(lc, 2) * tem
                end if
             end do
          end if
       end do

       if (cea%Debug(cea%Npt)) then
          write(IOOUT, '(/" NEW COMPONENTS")')
          write(IOOUT, '(/2x, 6A12)') (cmp(k), k = 1, nn)
       end if
    end if

    if (cea%Msing /= 0) then
       ! SWITCH ORDER OF MSING AND NLM COMPONENTS
       reduce = .true.
       lelim = cea%Nlm
       lsing = cea%Nlm

       if (cea%Msing /= cea%Nlm) then

          do j = 1, cea%Nspx
             aa = cea%A(cea%Msing, j)
             cea%A(cea%Msing, j) = cea%A(cea%Nlm, j)
             cea%A(cea%Nlm, j) = aa
          end do

          ja = cea%Jcm(cea%Msing)
          cea%Jcm(cea%Msing) = cea%Jcm(cea%Nlm)
          cea%Jcm(cea%Nlm) = ja
          ae = cmp(cea%Msing)
          cmp(cea%Msing) = cmp(cea%Nlm)
          cmp(cea%Nlm) = ae
          ae = cea%Elmt(cea%Msing)
          cea%Elmt(cea%Msing) = cea%Elmt(cea%Nlm)
          cea%Elmt(cea%Nlm) = trim(ae)
          ja = cea%Jx(cea%Msing)
          cea%Jx(cea%Msing) = cea%Jx(cea%Nlm)
          cea%Jx(cea%Nlm) = ja
          aa = cea%Atwt(cea%Msing)
          cea%Atwt(cea%Msing) = cea%Atwt(cea%Nlm)
          cea%Atwt(cea%Nlm) = aa
          aa = p%B0(cea%Msing)
          p%B0(cea%Msing) = p%B0(cea%Nlm)
          p%B0(cea%Nlm) = aa
          aa = pisave(cea%Msing)
          pisave(cea%Msing) = pisave(cea%Nlm)
          pisave(cea%Nlm) = aa

          do i = 1, 2
             aa = p%B0p(cea%Msing, i)
             p%B0p(cea%Msing, i) = p%B0p(cea%Nlm, i)
             p%B0p(cea%Nlm, i) = aa
          end do
       end if
    else if (.not. newcom .and. cea%Trace == 0.) then
       go to 600
    end if

    cea%Msing = 0
    tsize = xsize
    go to 500

1400 cea%Ttt(cea%Npt) = cea%Tt
    cea%Ppp(cea%Npt) = cea%Pp
    cea%Vlm(cea%Npt) = R0 * cea%Enn * cea%Tt / cea%Pp
    cea%Hsum(cea%Npt) = cea%Hsum(cea%Npt) * cea%Tt
    cea%Wm(cea%Npt) = 1. / cea%Enn
    gasfrc = cea%Enn/cea%Totn(cea%Npt)

    if (gasfrc < 0.0001) write(IOOUT, '(/" WARNING!  RESULTS MAY BE WRONG FOR POINT", i3, " DUE TO", &
         & /" LOW MOLE FRACTION OF GASES (", e15.8, ") (EQLBRM)")') cea%Npt, gasfrc

    if (cea%Trace /= 0) then
       if (lelim == 0) then
          do concurrent (j = 1:cea%Ng, cea%Enln(j) > -87)
             cea%En(j, cea%Npt) = exp(cea%Enln(j))
          end do
       else
          do concurrent (j = 1:cea%Ng, all(cea%A(lelim:ls, j) == 0))
             cea%En(j, cea%Npt) = exp(cea%Enln(j))
          end do
       end if
    end if

    if (cea%Debug(cea%Npt)) write(IOOUT, '(/" POINT=", i3, 3x, "P=", e13.6, 3x, "T=", e13.6, /3x, "H/R=", &
         & e13.6, 3x, "S/R=", e13.6, /3x, "M=", e13.6, 3x, "CP/R=", e13.6, 3x, &
         & "DLVPT=", e13.6, /3x, "DLVTP=", e13.6, 3x, "GAMMA(S)=", e13.6, 3x, "V=", e13.6)') &
         cea%Npt, cea%Pp, cea%Tt, cea%Hsum(cea%Npt), cea%Ssum(cea%Npt), cea%Wm(cea%Npt), cea%Cpr(cea%Npt), &
         cea%Dlvpt(cea%Npt), cea%Dlvtp(cea%Npt), cea%Gammas(cea%Npt), cea%Vlm(cea%Npt)

    if (cea%Tt >= cea%Tg(1) .and. cea%Tt <= cea%Tg(4)) go to 1600

    if (cea%Shock) go to 1600

    write(IOOUT, '(" THE TEMPERATURE=", e12.4, " IS OUT OF RANGE FOR POINT", i5, "(EQLBRM)")') cea%Tt, cea%Npt

    if (cea%Tt >= cea%Tg(1) * 0.8d0 .and. cea%Tt <= cea%Tg(4) * 1.1d0) go to 1600

    cea%Npt = cea%Npt + 1

1500 cea%Tt = 0
    cea%Npt = cea%Npt - 1
    write(IOOUT, '(/" CALCULATIONS STOPPED AFTER POINT", I3, "(EQLBRM)")') cea%Npt

1600 cea%Lsave = cea%Nlm
    cea%Nlm = ls

    if (cea%Npr > 0) cea%Gonly = .false.

    return
  end subroutine EQLBRM



  subroutine FROZEN(cea)
    !***********************************************************************
    ! CALCULATE PROPERTIES WITH FROZEN COMPOSITION AT ASSIGNED ENTROPY
    ! AND PRESSURE.  CALLED FROM ROCKET.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: iter, j
    real(8):: dlnt, dlpm

    cea%Convg = .false.
    cea%Tln = log(cea%Tt)
    dlpm = log(cea%Pp * cea%Wm(cea%Nfz))

    do concurrent (j = 1:cea%Ng, cea%En(j, cea%Nfz) /= 0)
       cea%Deln(j) = -(log(cea%En(j, cea%Nfz)) + dlpm)
    end do

    do iter = 1, 8
       call CPHS(cea)

       cea%Cpsum = sum(cea%En(1:cea%Ng, cea%Nfz) * cea%Cp(1:cea%Ng))
       cea%Ssum(cea%Npt) = sum(cea%En(1:cea%Ng, cea%Nfz) * (cea%S(1:cea%Ng) + cea%Deln(1:cea%Ng)))

       if (cea%Npr /= 0) then
          cea%Cpsum = cea%Cpsum + sum(cea%En(cea%Jcond(1:cea%Npr), cea%Nfz) * cea%Cp(cea%Jcond(1:cea%Npr)))
          cea%Ssum(cea%Npt) = cea%Ssum(cea%Npt) + sum(cea%En(cea%Jcond(1:cea%Npr), cea%Nfz) * cea%S(cea%Jcond(1:cea%Npr)))
       end if

       if (cea%Convg) then
          cea%Hsum(cea%Npt) = sum(cea%En(1:cea%Ngc, cea%Nfz) * cea%H0(1:cea%Ngc)) * cea%Tt

          cea%Ttt(cea%Npt) = cea%Tt
          cea%Gammas(cea%Npt) = cea%Cpsum / (cea%Cpsum - 1 / cea%Wm(cea%Nfz))
          cea%Vlm(cea%Npt) = R0 * cea%Tt / (cea%Wm(cea%Nfz) * cea%Pp)
          cea%Wm(cea%Npt) = cea%Wm(cea%Nfz)
          cea%Dlvpt(cea%Npt) = -1
          cea%Dlvtp(cea%Npt) = 1
          cea%Totn(cea%Npt) = cea%Totn(cea%Nfz)
          cea%Ppp(cea%Npt) = cea%Pp
          cea%Cpr(cea%Npt) = cea%Cpsum

          if (cea%Tt >= cea%Tg(1) * 0.8d0) then
             if (all(cea%En(cea%Ngp1:cea%Ngc, cea%Nfz) == 0)) return
             if (all(cea%Temp(1, cea%Ngp1-cea%Ng:cea%Ngc-cea%Ng) - 50 <= cea%Tt &
                  .and. cea%Tt <= cea%Temp(2, cea%Ngp1-cea%Ng:cea%Ngc-cea%Ng) + 50)) return
          end if

          cea%Tt = 0
          cea%Npt = cea%Npt - 1
          return

       else
          dlnt = (cea%Ssum(cea%Nfz) - cea%Ssum(cea%Npt)) / cea%Cpsum
          cea%Tln = cea%Tln + dlnt
          if (abs(dlnt) < 0.5d-4) cea%Convg = .true.
          cea%Tt = exp(cea%Tln)
       end if
    end do

    write(IOOUT, '(/" FROZEN DID NOT CONVERGE IN 8 ITERATIONS (FROZEN)")')

    cea%Tt = 0
    cea%Npt = cea%Nfz - 1

    return
  end subroutine FROZEN


  subroutine HCALC(cea)
    !***********************************************************************
    ! CALCULATE PROPERTIES FOR TOTAL REACTANT USING THERMO DATA FOR
    ! ONE OR MORE REACTANTS. USED ONLY FOR SHOCK AND DETON PROBLEMS.
    !***********************************************************************
    use mod_types
    use mod_functions
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(6):: date(maxNgc)
    character(2):: el(5)
    character(15):: sub
    integer, save:: i, icf, ifaz, itot, j, k, l, m, n, nall, nint, ntgas, ntot
    real(8):: bb(5), enj, er, sj, t1, t2, tem, thermo(9, 3), Tsave
    integer:: io_thermo

    Tsave = cea%Tt

    cea%Tm = 0
    if (cea%Pp > 0) cea%Tm = log(cea%Pp * cea%Wmix)

    cea%Ssum(cea%Npt) = 0
    cea%Hpp(1) = 0
    cea%Hpp(2) = 0
    cea%Hsub0 = 0
    cea%Cpmix = 0
    tem = (1 + cea%Oxfl)

    ! LOOP ON REACTANTS.
    ! if oxidant, k = 1
    ! if fuel,    k = 2
    cea%Nspr = cea%Nspx
    do n = 1, cea%Nreac
       k = 2
       if (cea%Fox(n)(:1) == 'O' .or. cea%Fox(n)(:1) == 'o') k = 1
       if (cea%Tt == 0) cea%Tt = cea%Rtemp(n)

       j = cea%Jray(n)

       if (j == 0) then
          ! SEARCH FOR REACTANT IN STORED THERMO SPECIES. STORE INDEX IN JRAY(N).
          ifaz = 0
          do j = 1, cea%Ngc
             if (cea%Rname(n) == cea%Prod(j) .or. '*' // cea%Rname(n) == cea%Prod(j)) then
                cea%Jray(n) = j
                if (j > cea%Ng) then
                   write(IOOUT, '(/" REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)")')
                   cea%Tt = 0
                   cea%Cpmix = 0
                   return
                end if
                go to 50
             end if
          end do

          ! SEARCH THERMO.LIB FOR SPECIES.
          open(newunit = io_thermo, file = cea%filename_thermo_lib, status = 'old', form = 'unformatted', action = 'read')

          read(io_thermo) cea%Tg, ntgas, ntot, nall

          cea%Nspr = cea%Nspr + 1
          do itot = 1, nall
             if (itot <= ntot) then
                icf = 3
                if (itot > ntgas) icf = 1
                read(io_thermo) sub, nint, date(cea%Nspr), (el(j), bb(j), j = 1, 5), ifaz, &
                     T1, T2, cea%Mw(cea%Nspr), ((thermo(l, m), l = 1, 9), m = 1, icf)
             else
                read(io_thermo) sub, nint, date(cea%Nspr), (el(j), bb(j), j = 1, 5), ifaz, T1, T2, cea%Mw(cea%Nspr), er
                if (nint /= 0) then
                   read(io_thermo) ((thermo(i, j), i = 1, 9), j = 1, nint)
                   icf = nint
                end if
             end if

             if (sub == cea%Rname(n) .or. sub == '*' // cea%Rname(n)) then
                if (ifaz <= 0 .and. nint > 0) then
                   do j = 1, 5
                      if (bb(j) == 0) exit
                      cea%Nfla(n) = j
                      cea%Ratom(n, j) = el(j)
                      cea%Rnum(n, j) = bb(j)
                   end do

                   cea%Jray(n) = cea%Nspr
                   j = cea%Nspr

                   do concurrent (l = 1:icf, m = 1:9)
                      cea%Coef(m, j, l) = thermo(m, l)
                   end do

                   go to 50

                else
                   if (ifaz > 0) write(IOOUT, '(/" REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)")')
                   if (nint == 0) write(IOOUT, '(/" COEFFICIENTS FOR ", a15, " ARE NOT AVAILABLE (HCALC)")') cea%Rname(n)
                   cea%Tt = 0
                   cea%Cpmix = 0
                   close(io_thermo)
                   return
                end if
             end if
          end do

          close(io_thermo)

          cea%Nspr = cea%Nspr - 1
          write(IOOUT, '(/" ERROR IN DATA FOR ", a15, " CHECK NAME AND TEMPERATURE", &
               & " RANGE IN", /, " thermo.inp (HCALC)")') cea%Rname(n)

          cea%Energy(n) = ' '
          cea%Tt = 0
          cea%Cpmix = 0

          return
       end if

       ! CALCULATE EN FOR REACTANT AND CALCULATE PROPERTIES.
50     if (cea%Moles) then
          enj = cea%Pecwt(n) / cea%Wp(k)
       else
          enj = cea%Pecwt(n) / cea%Rmw(n)
       end if
       enj = enj / tem
       if (k == 1) enj = enj * cea%Oxfl

       cea%Tln = log(cea%Tt)
       cea%En(j, cea%Npt) = enj
       l = 1

       if (ifaz <= 0) then
          if (cea%Tt > cea%Tg(2)) l = 2
          if (cea%Tt > cea%Tg(3) .and. ifaz < 0) l = 3
       end if

       cea%S(j) = calc_entropy(cea%Coef(:, j, l), cea%Tt, cea%Tln, .false.)

       cea%H0(j) = calc_enthalpy(cea%Coef(:, j, l), cea%Tt, cea%Tln, .false.)

       cea%Cp(j) = calc_specific_heat(cea%Coef(:, j, l), cea%Tt)

       if (abs(cea%H0(j)) < 0.01) cea%H0(j) = 0

       ! ADD CONTRIBUTION TO CP, H, AND S OF TOTAL REACTANT.
       cea%Cpmix = cea%Cpmix + cea%Cp(j) * enj

       ! FOR CONDENSED SPECIES:  SJ = S(J)
       sj = cea%S(j) - log(enj) - cea%Tm
       cea%Ssum(cea%Npt) = cea%Ssum(cea%Npt) + enj * sj
       er = cea%H0(j) * enj * cea%Tt
       cea%Hsub0 = cea%Hsub0 + er
       cea%Hpp(k) = cea%Hpp(k) + er
    end do

    if (Tsave /= 0) cea%Tt = Tsave

    return
  end subroutine HCALC


  subroutine MATRIX(cea)
    !***********************************************************************
    ! SET UP ITERATION OR DERIVATIVE MATRIX.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, iq, iq2, iq3, isym, j, k, kk, kmat
    real(8):: energyl, f, h, ss, sss, term, term1

    type(CEA_Point), pointer:: p !< current point
    p => cea%points(cea%iOF, cea%Npt)

    iq = cea%Nlm + cea%Npr
    cea%Iq1 = iq + 1
    iq2 = cea%Iq1 + 1
    iq3 = iq2 + 1
    kmat = iq3
    if (.not. cea%Convg .and. cea%Tp) kmat = iq2
    cea%Imat = kmat - 1
    cea%G = 0
    sss = 0
    cea%Hsum(cea%Npt) = 0

    ! BEGIN SET-UP OF ITERATION OR DERIVATIVE MATRIX
    do j = 1, cea%Ng
       cea%Mu(j) = cea%H0(j) - cea%S(j) + cea%Enln(j) + cea%Tm
       if (cea%En(j, cea%Npt) /= 0) then
          h = cea%H0(j) * cea%En(j, cea%Npt)
          f = cea%Mu(j) * cea%En(j, cea%Npt)
          ss = h - f
          term1 = h
          if (kmat == iq2) term1 = f

          do i = 1, cea%Nlm
             if (cea%A(i, j) /= 0) then
                term = cea%A(i, j) * cea%En(j, cea%Npt)
                do concurrent (k = i:cea%Nlm)
                   cea%G(i, k) = cea%G(i, k) + cea%A(k, j) * term
                end do
                cea%G(i, cea%Iq1) = cea%G(i, cea%Iq1) + term
                cea%G(i, iq2) = cea%G(i, iq2) + cea%A(i, j) * term1
                if (.not. (cea%Convg .or. cea%Tp)) then
                   cea%G(i, iq3) = cea%G(i, iq3) + cea%A(i, j) * f
                   if (cea%Sp) cea%G(iq2, i) = cea%G(iq2, i) + cea%A(i, j) * ss
                end if
             end if
          end do

          if (kmat /= iq2) then
             if (cea%Convg .or. cea%Hp) then
                cea%G(iq2, iq2) = cea%G(iq2, iq2) + cea%H0(j) * h
                if (.not. cea%Convg) then
                   cea%G(iq2, iq3) = cea%G(iq2, iq3) + cea%H0(j) * f
                   cea%G(cea%Iq1, iq3) = cea%G(cea%Iq1, iq3) + f
                end if
             else
                cea%G(iq2, cea%Iq1) = cea%G(iq2, cea%Iq1) + ss
                cea%G(iq2, iq2) = cea%G(iq2, iq2) + cea%H0(j) * ss
                cea%G(iq2, iq3) = cea%G(iq2, iq3) + cea%Mu(j) * ss
                cea%G(cea%Iq1, iq3) = cea%G(cea%Iq1, iq3) + f
             end if
          end if
          cea%G(cea%Iq1, iq2) = cea%G(cea%Iq1, iq2) + term1
       end if
    end do

    ! CONDENSED SPECIES
    if (cea%Npr /= 0) then
       do k = 1, cea%Npr
          j = cea%Jcond(k)
          kk = cea%Nlm + k
          cea%Mu(j) = cea%H0(j) - cea%S(j)

          do concurrent (i = 1:cea%Nlm)
             cea%G(i, kk) = cea%A(i, j)
             cea%G(i, kmat) = cea%G(i, kmat) - cea%A(i, j) * cea%En(j, cea%Npt)
          end do

          cea%G(kk, iq2) = cea%H0(j)
          cea%G(kk, kmat) = cea%Mu(j)
          cea%Hsum(cea%Npt) = cea%Hsum(cea%Npt) + cea%H0(j) * cea%En(j, cea%Npt)

          if (cea%Sp) then
             sss = sss + cea%S(j) * cea%En(j, cea%Npt)
             cea%G(iq2, kk) = cea%S(j)
          end if
       end do
    end if

    sss = sss + cea%G(iq2, cea%Iq1)
    cea%Hsum(cea%Npt) = cea%Hsum(cea%Npt) + cea%G(cea%Iq1, iq2)
    cea%G(cea%Iq1, cea%Iq1) = cea%Sumn - cea%Enn

    ! REFLECT SYMMETRIC PORTIONS OF THE MATRIX
    isym = cea%Iq1
    if (cea%Hp .or. cea%Convg) isym = iq2

    do i = 1, isym
       !DIR$ IVDEP
       do concurrent (j = i:isym)
          cea%G(j, i) = cea%G(i, j)
       end do
    end do

    ! COMPLETE THE RIGHT HAND SIDE
    if (.not. cea%Convg) then
       do concurrent (i = 1:cea%Nlm)
          cea%G(i, kmat) = cea%G(i, kmat) + p%B0(i) - cea%G(i, cea%Iq1)
       end do
       cea%G(cea%Iq1, kmat) = cea%G(cea%Iq1, kmat) + cea%Enn - cea%Sumn

       ! COMPLETE ENERGY ROW AND TEMPERATURE COLUMN
       if (kmat /= iq2) then
          if (cea%Sp) energyl = cea%S0 + cea%Enn - cea%Sumn - sss
          if (cea%Hp) energyl = cea%Hsub0/cea%Tt - cea%Hsum(cea%Npt)
          cea%G(iq2, iq3) = cea%G(iq2, iq3) + energyl
          cea%G(iq2, iq2) = cea%G(iq2, iq2) + cea%Cpsum
       end if

    else
       if (cea%Pderiv) then
          ! PDERIV = .true.-- SET UP MATRIX TO SOLVE FOR DLVPT
          cea%G(cea%Iq1, iq2) = cea%Enn
          do concurrent (i = 1:iq)
             cea%G(i, iq2) = cea%G(i, cea%Iq1)
          end do
       end if
       cea%G(iq2, iq2) = cea%G(iq2, iq2) + cea%Cpsum
    end if

    if (cea%Vol .and. .not. cea%Convg) then
       ! CONSTANT VOLUME MATRIX
       if (kmat == iq2) then
          do concurrent (i = 1:iq)
             cea%G(i, cea%Iq1) = cea%G(i, iq2)
          end do
       else

          !DIR$ IVDEP
          do concurrent (i = 1:iq)
             cea%G(cea%Iq1, i) = cea%G(iq2, i) - cea%G(cea%Iq1, i)
             cea%G(i, cea%Iq1) = cea%G(i, iq2) - cea%G(i, cea%Iq1)
             cea%G(i, iq2) = cea%G(i, iq3)
          end do

          cea%G(cea%Iq1, cea%Iq1) = cea%G(iq2, iq2) - cea%G(cea%Iq1, iq2) - cea%G(iq2, cea%Iq1)
          cea%G(cea%Iq1, iq2) = cea%G(iq2, iq3) - cea%G(cea%Iq1, iq3)
          if (cea%Hp) cea%G(cea%Iq1, iq2) = cea%G(cea%Iq1, iq2) + cea%Enn
       end if

       kmat = cea%Imat
       cea%Imat = cea%Imat - 1
    end if

    return
  end subroutine MATRIX



  subroutine NEWOF(cea)
    !***********************************************************************
    ! CALCULATE NEW VALUES OF B0 AND HSUB0 FOR NEW OF RATIO
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, j
    real(8):: assval, bigb, bratio, dbi, smalb, tem, v1, v2

    type(CEA_Point), pointer:: p !< current point

    cea%iOF = cea%iOF + 1
    cea%Npt = 1

    p => cea%points(cea%iOF, cea%Npt)

    if (.not. cea%Short) write(IOOUT, '(/" O/F = ", f10.6)') cea%Oxfl
    cea%Eqrat = 0
    tem = cea%Oxfl + 1
    v2 = (cea%Oxfl * cea%Vmin(1) + cea%Vmin(2)) / tem
    v1 = (cea%Oxfl * cea%Vpls(1) + cea%Vpls(2)) / tem
    if (v2 /= 0.) cea%Eqrat = abs(v1/v2)
    do i = 1, cea%Nlm
       p%B0(i) = (cea%Oxfl * p%B0p(i, 1) + p%B0p(i, 2)) / tem
       dbi = abs(p%B0(i))
       if (i == 1) then
          bigb = dbi
          smalb = dbi
       else if (dbi /= 0.) then
          if (dbi < smalb) smalb = dbi
          if (dbi > bigb) bigb = dbi
       end if
    end do
    cea%Bcheck = bigb * .000001d0

    ! CALCUALTE MOLECULAR WEIGHT OF TOTAL REACTANT, WMIX.
    if (all(cea%Am /= 0)) then
       cea%Wmix = (cea%Oxfl + 1) * cea%Am(1) * cea%Am(2) / (cea%Am(1) + cea%Oxfl * cea%Am(2))
    else
       cea%Wmix = cea%Am(2)
       if (cea%Am(2) == 0.0) cea%Wmix = cea%Am(1)
    end if

    ! IF ASSIGNED U OR H NOT GIVEN IN PROB DATA, INITIAL HSUB0 = 1.d30
    if (cea%Size == 0) assval = cea%Hsub0
    if (assval >= 1.d30) cea%Hsub0 = (cea%Oxfl * cea%Hpp(1) + cea%Hpp(2)) / tem

    ! NOTE THAT "BRATIO" IS "BRATIO" IN SEC 3.2 IN RP-1311.
    bratio = smalb / bigb
    cea%Size = 18.420681d0
    if (bratio < 1.d-5) cea%Size = log(1000/bratio)
    cea%Jsol = 0
    cea%Jliq = 0

    if (.not. cea%Short) then
       write(IOOUT, '(/, 23x, "EFFECTIVE FUEL", 5x, "EFFECTIVE OXIDANT", 8x, "MIXTURE")')
       if (cea%Vol) write(IOOUT, '(" INTERNAL ENERGY", 11x, "u(2)/R", 14x, "u(1)/R", 14x, "u0/R")')
       if (.not. cea%Vol) write(IOOUT, '(" ENTHALPY", 18x, "h(2)/R", 14x, "h(1)/R", 15x, "h0/R")')
       write(IOOUT, '(" (KG-MOL)(K)/KG", 4x, e18.8, 2e20.8)') cea%Hpp(2), cea%Hpp(1), cea%Hsub0
       write(IOOUT, '(/" KG-FORM.WT./KG", 13x, "bi(2)", 15x, "bi(1)", 15x, "b0i")')
    end if

    do i = 1, cea%Nlm
       j = cea%Jcm(i)
       if (.not. cea%Short) write(IOOUT, '(1x, a16, 3e20.8)') cea%Prod(j), p%B0p(i, 2), p%B0p(i, 1), p%B0(i)
    end do

    return
  end subroutine NEWOF


  subroutine ROCKET(cea)
    !***********************************************************************
    ! EXECUTIVE ROUTINE FOR ROCKET PROBLEMS.
    !***********************************************************************
    use mod_types
    use mod_legacy_io, only: RKTOUT
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, i01, i12, iof, iplt1, iplte, ipp, isub, isup1, isupsv, itnum, &
         itrot, nar, nipp, niter, nn, npr1, nptth
    logical:: done, seql, thi
    real(8):: a1l = -1.26505, b1 = 1.0257, c1 = -1.2318, pa = 1e5
    real(8):: acatsv, aeatl, appl, aratio, asq, check, cprf, dd, dh, &
         dlnp, dlnpe, dlt, dp, eln, mat, msq, p1, pcpa, pcplt, pinf, pinj, &
         pinjas, pjrat, ppa, pr, pracat, prat, pratsv, pvg, test, tmelt, usq
    integer:: ip, it

    type(CEA_Point), pointer:: p !< current point

    iplte = cea%Iplt
    isup1 = 1

    do concurrent (iof = 1:cea%Nof)
       cea%points(iof, 1)%App = 1
    end do

    cea%Iopt = 0
    cea%Npp = cea%Npp + 2
    nn = cea%Npp
    i01 = 0
    i12 = 1
    nipp = 1
    nptth = 2
    if (cea%Fac) then
       cea%Eql = .true.
       cea%Npp = cea%Npp + 1
       if (cea%Acat /= 0) then
          cea%Iopt = 1
       else if (cea%Ma /= 0) then
          cea%Iopt = 2
       else
          write(IOOUT, '(/" FATAL ERROR!! EITHER mdot OR ac/at MISSING FOR fac PROBLEM (ROCKET)")')
          cea%Tt = 0
          go to 1400
       end if
       i01 = 1
       i12 = 2
       nipp = 2
       nptth = 3
       do i = cea%Nsub, 1, - 1
          cea%Subar(i+1) = cea%Subar(i)
       end do
       cea%Nsub = cea%Nsub + 1
       if (cea%Iopt /= 1) then
          if (cea%Acat == 0) cea%Acat = 2
       end if
       cea%Subar(1) = cea%Acat
    else if (.not. cea%Eql .and. cea%Nfz > 1 .and. cea%Nsub > 0) then
       cea%Nsub = 0
       write(IOOUT, '(/" WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ", /, &
            & " RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)")')
    end if
    nn = nn + cea%Nsub + cea%Nsup
    if (cea%Nfz > 2 .and. nn > Ncol-2) then
       write(IOOUT, '(/" WARNING!!  nfz NOT ALLOWED TO BE > 2 IF THE TOTAL", /, &
            & " NUMBER OF POINTS IS >", i3, " (ROCKET)")') Ncol - 2
       cea%Nfz = 1
       cea%Froz = .false.
    end if
    seql = cea%Eql
    cea%Tt = cea%Tcest
    cea%Pp = cea%P(1)

    do concurrent (iof = 1:cea%Nof)
       cea%points(iof, i12)%App = 1
    end do

    ! LOOP FOR EACH O/F
    iof = 0
100 it = 1
    iof = iof + 1
    cea%Oxfl = cea%Oxf(iof)
    if (cea%T(1) /= 0.) then
       cea%Tp = .true.
    else
       cea%Hp = .true.
    end if
    cea%Sp = .false.
    call NEWOF(cea)
    if (cea%T(1) /= 0) cea%Tt = cea%T(1)
    ! LOOP FOR CHAMBER PRESSURES
200 do ip = 1, cea%Np
       itnum = 0
       cea%Area = .false.
       if (cea%T(1) == 0) cea%Hp = .true.
       if (cea%T(1) /= 0) cea%Tp = .true.
       cea%Sp = .false.
       cea%Eql = .true.
       isub = 1
       cea%Isup = 1
       cea%Pp = cea%P(ip)
       pinf = cea%Pp
       ipp = 1
       itrot = 3
       isupsv = 1
       niter = 1
       cea%Page1 = .true.
       iplt1 = iplte
       cea%Iplt = iplte
       done = .false.
       ! LOOP FOR OUTPUT COLUMNS
250    nar = cea%Npt
       if (cea%Eql) then
          call EQLBRM(cea)
          if (cea%Npt == cea%Nfz) cprf = cea%Cpsum
       else
          call FROZEN(cea)
       end if
       ! TT = 0 IF NO CONVERGENCE
       if (cea%Tt /= 0.) then
          ! TEST FOR FINITE AREA COMBUSTOR
          if (.not. cea%Fac) go to 400
          pinjas = cea%P(ip) * pa
          pinj = pinjas
          if (cea%Npt <= 2) then
             if (cea%Npt == 1 .and. cea%Trnspt) call TRANP(cea)
             if (cea%Npt == 2) pinf = cea%Ppp(2)
          end if
          if (cea%Npt /= 1) go to 400
          ! INITIAL ESTIMATE FOR PC (AND ACAT IF NOT ASSIGNED)
          do i = 1, 4
             prat = (b1 + c1 * cea%Acat) / (1 + a1l * cea%Acat)
             ppa = pinj * prat
             if (cea%Iopt == 1) go to 260
             cea%Acat = ppa / (cea%Ma * 2350)
             if (cea%Acat >= 1) then
                pratsv = prat
                if (cea%Debugf) then
                   if (i <= 1) write(IOOUT, '(/"  ITERATION", 9x, "PC", 7x, "CONTRACTION RATIO")')
                   write(IOOUT, '(5x, i2, 7x, f12.2, 3x, f12.6)') i, ppa, cea%Acat
                end if
             else
                write(IOOUT, '(/" INPUT VALUE OF mdot/a =", f12.3, " IS TOO LARGE."/ &
                     & " GIVES CONTRACTION RATIO ESTIMATE LESS THAN 1 (ROCKET)")') cea%Ma
                cea%Tt = 0
                go to 1400
             end if
          end do
          cea%Subar(1) = cea%Acat
260       p => cea%points(cea%iOF, 1)
          cea%Pp = ppa / pa
          p%App = cea%Pp / cea%Ppp(1)
          go to 1100
       else
          if (cea%Npt < 1) go to 1400
          if (.not. cea%Area) go to 600
          cea%Npt = nar - 1
          cea%Isup = cea%Nsup + 2
          cea%Isv = 0
          itnum = 0
          go to 950
       end if
300    cea%Hp = .true.
       cea%Sp = .false.
       niter = niter + 1
       cea%Isv = 0
       cea%Npt = 2
       ipp = 2
       call SETEN(cea)
       go to 250
350    done = .true.
       p => cea%points(cea%iOF, 1)
       p%App = cea%Ppp(2) / cea%Ppp(1)
       cea%Area = .false.
       if (cea%Nsub > 1) isub = 2
       cea%Isv = 4
       cea%Npt = 2
       ipp = min(4, cea%Npp)
       call SETEN(cea)
       p => cea%points(cea%iOF, 2)
       cea%Cpr(2) = cea%Cpr(4)
       cea%Dlvpt(2) = cea%Dlvpt(4)
       cea%Dlvtp(2) = cea%Dlvtp(4)
       cea%Gammas(2) = cea%Gammas(4)
       cea%Hsum(2) = cea%Hsum(4)
       cea%Ppp(2) = cea%Ppp(4)
       p%App = cea%Ppp(1)/pinf
       cea%Ssum(2) = cea%Ssum(4)
       cea%Totn(2) = cea%Totn(4)
       cea%Ttt(2) = cea%Ttt(4)
       cea%Vlm(2) = cea%Vlm(4)
       cea%Wm(2) = cea%Wm(4)
       if (.not. cea%Short) write(IOOUT, '(" END OF CHAMBER ITERATIONS")')
       go to 600
       ! INITIALIZE FOR THROAT
400    if (ipp > nipp) then
          usq = 2 * (cea%Hsum(1) - cea%Hsum(cea%Npt)) * R0
          if (ipp > nptth) go to 600
          ! THROAT
          if (.not. thi) then
             cea%Vv = cea%Vlm(nptth)
             pvg = cea%Pp * cea%Vv * cea%Gammas(nptth)
             if (pvg == 0) then
                write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                go to 550
             else
                msq = usq / pvg
                if (cea%Debug(1) .or. cea%Debug(2)) write(IOOUT, '(/" USQ=", e15.8, 5x, "PVG=", e15.8)') usq, pvg
                dh = abs(msq - 1)
                if (dh <= 0.4d-4) go to 550
                if (itrot > 0) then
                   p1 = cea%Pp
                   if (cea%Jsol /= 0) then
                      tmelt = cea%Tt
                      cea%Pp = cea%Pp * (1 + msq * cea%Gammas(nptth)) / (cea%Gammas(nptth) + 1)
                   else if (tmelt == 0) then
                      cea%Pp = cea%Pp * (1 + msq * cea%Gammas(nptth)) / (cea%Gammas(nptth) + 1)
                   else
                      write(IOOUT, '(/" WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)")')
                      p => cea%points(cea%iOF, nptth)
                      dlt = log(tmelt / cea%Tt)
                      dd = dlt * cea%Cpr(nptth) / (cea%Enn * cea%Dlvtp(nptth))
                      cea%Pp = cea%Pp * EXP(dd)
                      p%App = cea%P(ip) / cea%Pp
                      if (cea%Fac) p%App = pinf / cea%Pp
                      if (cea%Eql .and. .not. cea%Short) write(IOOUT, '(" Pinf/Pt =", F9.6)') p%App
                      thi = .true.
                      go to 250
                   end if
                   go to 500
                else if (itrot < 0) then
                   if (itrot < -19) then
                      write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                      go to 550
                   else
                      if (cea%Npr /= npr1) go to 550
                      cea%Pp = cea%Pp - dp
                      go to 500
                   end if
                else if (cea%Npr == npr1) then
                   write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                   go to 550
                else
                   dp = abs(cea%Pp - p1) / 20
                   cea%Pp = max(cea%Pp, p1)
                   write(IOOUT, '(/" WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)")')
                   cea%Pp = cea%Pp - dp
                   go to 500
                end if
             end if
          else
             cea%Gammas(nptth) = 0
             go to 550
          end if
       else
          if (.not. cea%Fac .and. cea%Trnspt) call TRANP(cea)
          if (cea%Npt == cea%Nfz) cea%Eql = seql
          cea%Tp = .false.
          cea%Hp = .false.
          cea%Sp = .true.
          cea%S0 = cea%Ssum(i12)
       end if
450    tmelt = 0
       itrot = 3
       thi = .false.
       p => cea%points(cea%iOF, nptth)
       p%App = ((cea%Gammas(i12) + 1) / 2)**(cea%Gammas(i12) / (cea%Gammas(i12) - 1))
       if (cea%Eql .and. .not. cea%Short) write(IOOUT, '(" Pinf/Pt =", f9.6)') p%App
       cea%Pp = Pinf / p%App
       cea%Isv = -i12
       go to 1200
500    npr1 = cea%Npr
       p => cea%points(cea%iOF, nptth)
       p%App = cea%P(ip) / cea%Pp
       if (cea%Fac) p%App = Pinf / cea%Pp
       if (cea%Eql .and. .not. cea%Short) write(IOOUT, '(" Pinf/Pt =", f9.6)') p%App
       itrot = itrot - 1
       go to 250
550    cea%Awt = cea%Enn * cea%Tt / (cea%Pp * sqrt(usq))
       p => cea%points(cea%iOF, nptth)
       pcplt = log(p%App)
600    cea%Isv = 0
       cea%points(cea%iOF, cea%Npt)%AeAt = cea%Enn * cea%Ttt(cea%Npt) / (cea%Pp * sqrt(usq) * cea%Awt)
       if (cea%Tt == 0) go to 1150
       if (cea%Area) go to 750
       if (cea%Trnspt .and. (.not. cea%Fac .or. done .or. cea%Npt > 2)) call TRANP(cea)
       if (cea%Npt == cea%Nfz) cea%Eql = seql
       if (cea%Fac) then
          if (cea%Npt == nptth) then
             cea%Area = .true.
             go to 750
          else if (cea%Npt == 2 .and. done) then
             cea%Npt = 3
             !  The following statement was corrected 1/30/2004.  Only fac parameters 
             !    after combustion were affected--generally extra or missing points.
             !  (remove) if (ipp <= cea%Npp) ipp = ipp - 1
             if (ipp < cea%Npp .or. cea%Npp == 4) ipp = ipp - 1
          end if
       end if
650    if (ipp < cea%Npp) go to 1100
700    if (cea%Nsub == i01 .and. cea%Nsup == 0) go to 1150
       cea%Area = .true.
       ! PCP ESTIMATES FOR AREA RATIOS
750    if (itnum == 0) then
          dlnp = 1
          itnum = 1
          aratio = cea%Subar(isub)
          if ((.not. cea%Fac .or. done) .and. cea%Nsub <= i01) aratio = cea%Supar(cea%Isup)
          if (.not. cea%Eql .and. cea%Nfz >= 3) then
             if (aratio <= cea%points(cea%iOF, cea%Nfz)%Aeat) then
                write(IOOUT, '(/, " WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED", &
                     & " WHERE THE ASSIGNED", /, " SUPERSONIC AREA RATIOS WERE ", &
                     & "LESS THAN THE VALUE AT POINT nfz =", i3, " (ROCKET)")') cea%Nfz
                go to 1050
             end if
          end if
          if (aratio  <  1) then
             write(IOOUT, '(/" AN ASSIGNED AREA RATIO IS < 1 (ROCKET)")')
             go to 1050
          end if
          eln = log(aratio)
          if (cea%Fac) then
             if (.not. done) go to 800
          end if
          if (cea%Nsub <= i01) then
             if (cea%Nfz == ipp) isupsv = cea%Isup
             if (cea%Supar(cea%Isup) < 2) then
                appl = sqrt(eln*(1.535d0 + 3.294d0 * eln)) + pcplt
                go to 1100
             else
                if (cea%Isup > isup1) then
                   if (cea%Supar(cea%Isup-1) >= 2) go to 850
                end if
                appl = cea%Gammas(nptth) + eln * 1.4
                go to 1100
             end if
          end if
          ! TEST FOR CONVERGENCE ON AREA RATIO.
       else if (cea%Gammas(cea%Npt) > 0) then
          check = 0.00004
          p => cea%points(cea%iOF, cea%Npt)
          if (cea%Debug(cea%Npt)) write(IOOUT, '(/" ITER=", i2, 2x, "ASSIGNED AE/AT=", f14.7, 3x, "AE/AT=", f14.7, &
               & /, 2x, "PC/P=", f14.7, 2x, "DELTA LN PCP=", f14.7)') &
               itnum, aratio, p%AeAt, p%App, dlnp
          if (abs(cea%points(cea%iOF, cea%Npt)%AeAt - Aratio) / Aratio <= check) go to 900
          if (abs(dlnp) < 0.00004) go to 900
          aeatl = log(cea%points(cea%iOF, cea%Npt)%AeAt)
          itnum = itnum + 1
          if (itnum > 10) then
             write(IOOUT, '(/" WARNING!!  DID NOT CONVERGE FOR AREA RATIO =", F10.5, " (ROCKET)")') aratio
             go to 900
          else
             ! IMPROVED PCP ESTIMATES.
             asq = cea%Gammas(cea%Npt) * cea%Enn * R0 * cea%Tt
             dlnpe = cea%Gammas(cea%Npt) * usq / (usq - asq)
             go to 850
          end if
       else
          write(IOOUT, '(/" WARNING!!  AREA RATIO CALCULATION CANNOT BE DONE ", &
               & "BECAUSE GAMMAs", /, " CALCULATION IMPOSSIBLE. (ROCKET)")')
          cea%Npt = cea%Npt - 1
          if (cea%Nsub <= 0) isup1 = 100
          if (cea%Nsub < 0.) cea%Nsup = cea%Isup - 1
          if (cea%Nsub > 0) cea%Nsub = isub - 1
          go to 1000
       end if
800    appl = pcplt / (cea%Subar(isub) + (10.587 * eln**2 + 9.454) * eln)
       if (Aratio < 1.09) appl = 0.9 * appl
       if (Aratio > 10) appl = appl / Aratio
       if (isub > 1 .or. cea%Npt == Ncol) go to 1100
       go to 1200
850    dlnp = dlnpe * eln - dlnpe * aeatl
       appl = appl + dlnp
       if (itnum == 1) go to 1100
       if (appl < 0.) appl = 0.000001
       p => cea%points(cea%iOF, cea%Npt)
       p%App = EXP(appl)
       cea%Pp = Pinf / p%App
       go to 250
       ! CONVERGENCE HAS BEEN REACHED FOR ASSIGNED AREA RATIO
900    cea%points(cea%iOF, cea%Npt)%AeAt = Aratio
       if (cea%Fac) then
          if (.not. done) then
             if (cea%Iopt == 1) then
                ! OPTION 1 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
                ! PRESSURE AND CONTRACTION RATIO. IMPROVED ESTIMATE FOR PC
                cea%Area = .false.
                itnum = 0
                ppa = cea%Ppp(cea%Npt) * pa
                pinj = ppa + 1.d05 * usq / cea%Vlm(cea%Npt)
                test = (pinj - pinjas) / pinjas
                pcpa = pinf * pa
                if (cea%Debugf) then
                   write(IOOUT, '(" ITER", 3x, "TEST", 3x, "ASSIGNED PINJ", 1x, "CALC PINJ", 5x, &
                        & "PC", 7x, "P AT ACAT", 3x, "PREV ACAT", 2x, "ACAT")')
                   write(IOOUT, '(i3, f10.6, 1x, 4f12.2, 2f9.5)') niter, test, pinjas, pinj, pcpa, ppa, acatsv, cea%Acat
                end if
                if (abs(test) < 0.00002) go to 350
                prat = pinjas / pinj
                cea%Pp = pinf * prat
                go to 300
             else if (cea%Iopt == 2) then
                ! OPTION 2 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
                ! PRESSURE AND MASS FLOW PER UNIT AREA. IMPROVED ESTIMATE FOR PC
                ! AND ACAT
                acatsv = cea%Acat
                pratsv = prat
                cea%Area = .false.
                itnum = 0
                ppa = cea%Ppp(4) * pa
                pinj = ppa + 1.d05 * usq / cea%Vlm(4)
                mat = pa / (cea%Awt * R0)
                cea%Acat = mat / cea%Ma
                prat = (b1 + c1 * cea%Acat) / (1 + a1l * cea%Acat)
                test = (pinj - pinjas) / pinjas
                pcpa = pinf * pa
                if (cea%Debugf) then
                   write(IOOUT, '(" ITER", 3x, "TEST", 3x, "ASSIGNED PINJ", 1x, "CALC PINJ", 5x, &
                        & "PC", 7x, "P AT ACAT", 3x, "PREV ACAT", 2x, "ACAT")')
                   write(IOOUT, '(i3, f10.6, 1x, 4f12.2, 2f9.5)') niter, test, pinjas, pinj, pcpa, ppa, acatsv, cea%Acat
                end if
                if (abs(test) < 0.00002) go to 350
                pjrat = pinj / pinjas
                cea%Pp = pinf
                do i = 1, 2
                   pracat = pratsv / prat
                   pr = pjrat * pracat
                   cea%Pp = cea%Pp / pr
                   pcpa = cea%Pp * pa
                   cea%Acat = cea%Acat / pr
                   cea%Subar(1) = cea%Acat
                   pratsv = prat
                   pjrat = 1
                   prat = (b1 + c1 * cea%Acat) / (1 + a1l * cea%Acat)
                   if (cea%Debugf) write(IOOUT, '(" NEW PC = ", f10.2, 2x, "NEW ACAT = ", f9.6, 2x, "PJRAT =", &
                        & f10.7, " PRACAT =", f10.7)') pcpa, cea%Acat, pjrat, pracat
                end do
                go to 300
             end if
          end if
       end if
950    if (cea%Trnspt) call TRANP(cea)
       if (cea%Npt == cea%Nfz) cea%Eql = seql
1000   itnum = 0
       if (cea%Nsub > i01) then
          isub = isub + 1
          if (isub <= cea%Nsub) go to 750
          isub = 1
          cea%Nsub = -cea%Nsub
          if (cea%Isup <= cea%Nsup) go to 750
          cea%Area = .false.
          go to 1150
       end if
1050   cea%Isup = cea%Isup + 1
       itnum = 0
       if (cea%Isup <= cea%Nsup) go to 750
       cea%Isup = isupsv
       cea%Area = .false.
       go to 1150
       ! TEST FOR OUTPUT -- SCHEDULES COMPLETE OR NPT=Ncol
1100   cea%Isv = cea%Npt
       if (cea%Npt /= Ncol) go to 1200
1150   if (.not. cea%Eql) then
          if (cea%Nfz <= 1) then
             cea%Cpr(cea%Nfz) = cprf
             cea%Gammas(cea%Nfz) = cprf / (cprf - 1 / cea%Wm(cea%Nfz))
          end if
       end if
       call RKTOUT(cea, it)
       cea%Iplt = cea%Iplt + cea%Npt
       if (.not. cea%Page1) then
          cea%Iplt = cea%Iplt - 2
          if (cea%Iopt /= 0) cea%Iplt = cea%Iplt - 1
          cea%Iplt = min(cea%Iplt, 500)
       else
          cea%Page1 = .false.
       end if
       iplte = max(iplte, cea%Iplt)
       dlnp = 1
       if (cea%Tt == 0) cea%Area = .false.
       if (.not. cea%Eql .and. cea%Tt == 0.) write(IOOUT, '(/" WARNING!!  CALCULATIONS WERE STOPPED BECAUSE NEXT ", &
            & "POINT IS MORE", /, " THAN 50 K BELOW THE TEMPERATURE", &
            & " RANGE OF A CONDENSED SPECIES (ROCKET)")')
       if (cea%Isv == 0) then
          ! PCP, SUBAR, AND SUPAR SCHEDULES COMPLETED
          if (cea%Nsub < 0) cea%Nsub = -cea%Nsub
          if (.not. cea%Froz .or. .not. cea%Eql) go to 1300
          ! SET UP FOR FROZEN.
          if (cea%Eql) cea%Iplt = iplt1
          cea%Eql = .false.
          cea%Page1 = .true.
          call SETEN(cea)
          cea%Tt = cea%Ttt(cea%Nfz)
          ipp = cea%Nfz
          if (cea%Nfz == cea%Npt) go to 1150
          cea%Npt = cea%Nfz
          cea%Enn = 1./cea%Wm(cea%Nfz)
          if (cea%Nfz == 1) go to 450
          if (cea%Nsub > 0) then
             cea%Nsub = -cea%Nsub
             write(IOOUT, '(/" WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ", /, &
                  & " RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)")')
          end if
          p => cea%points(cea%iOF, cea%Nfz)
          if (p%App < cea%points(cea%iOF, nptth)%App) then
             write(IOOUT, '(/" WARNING!!  FREEZING IS NOT ALLOWED AT A SUBSONIC ", &
                  & "PRESSURE RATIO FOR nfz GREATER"/" THAN 1. FROZEN ", &
                  & "PERFORMANCE CALCULATIONS WERE OMITTED (ROCKET)")')
          else
             if (cea%Nfz < cea%Npp) go to 1200
             go to 700
          end if
          go to 1300
       else
          if (cea%Eql) write(IOOUT, '(////)')
          cea%Npt = nptth
       end if
       ! SET INDICES AND ESTIMATES FOR NEXT POINT.
1200   cea%Npt = cea%Npt + 1
       if (cea%Eql .or. (cea%Isv == -i12 .and. .not. seql)) then
          ! THE FOLLOWING STATEMENT WAS ADDED TO TAKE CARE OF A SITUATION
          ! WHERE EQLBRM WENT SINGULAR WHEN STARTING FROM ESTIMATES WHERE
          ! BOTH SOLID AND LIQUID WERE INCLUDED.  JULY 27, 1990.
          if (cea%Jliq /= 0 .and. cea%Isv > 0) cea%Isv = 0
          call SETEN(cea)
       end if
1250   ipp = ipp + 1
       p => cea%points(cea%iOF, cea%Npt)
       if (cea%Npt > nptth) then
          if (cea%Area) then
             p%App = EXP(appl)
          else
             p%App = cea%Pcp(ipp - nptth)
             if (cea%Fac) p%App = p%App * Pinf / cea%Ppp(1)
             if (.not. cea%Eql .and. p%App < cea%points(cea%iOF, cea%Nfz)%App) then
                write(IOOUT, '(/, " WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED", &
                     & " WHERE THE ASSIGNED", /, &
                     & " PRESSURE RATIOS WERE LESS THAN ", &
                     & "THE VALUE AT POINT nfz =", i3, " (ROCKET)")') cea%Nfz
                go to 1250
             end if
          end if
          cea%Pp = pinf / p%App
          if (cea%Fac) then
             if (cea%Area) then
                if (isub <= cea%Nsub .and. isub > i01 .and. aratio >= cea%points(cea%iOF, 2)%AeAt) then
                   write(IOOUT, '(/" WARNING!!  ASSIGNED subae/at =", f10.5, " IS NOT ", &
                        & "PERMITTED TO BE GREATER"/" THAN ac/at =", f9.5, &
                        & ".  POINT OMITTED (ROCKET)")') aratio, cea%points(cea%iOF, 2)%AeAt
                   cea%Npt = cea%Npt - 1
                   go to 1000
                end if
             else if (cea%Npt > nptth .and. cea%Pcp(ipp-3) < cea%Ppp(1) / cea%Ppp(2)) then
                write(IOOUT, '(/" WARNING!!  ASSIGNED pip =", F10.5, &
                     & " IS NOT PERMITTED"/" TO BE LESS THAN  Pinj/Pc =", f9.5, &
                     & ". POINT OMITTED", " (ROCKET)")') cea%Pcp(ipp-3), cea%Ppp(1) / cea%Ppp(2)
                cea%Npt = cea%Npt - 1
                go to 650
             end if
          end if
       end if
       go to 250
1300   cea%Npt = 1
       ! CHECK FOR COMPLETED SCHEDULES -
       ! 1) CHAMBER PRESSURES(IP = NP)
       ! 2) CHAMBER TEMPERATURES(IT = NT)
       ! 3) O/F VALUES(IOF = NOF)
       if (ip == cea%Np .and. it == cea%Nt .and. iof == cea%Nof) go to 1400
       write(IOOUT, '(////)')
       call SETEN(cea)
       cea%Tt = cea%Ttt(i12)
    end do
    if (it < cea%Nt) then
       it = it + 1
       cea%Tt = cea%T(it)
       go to 200
    else if (iof < cea%Nof) then
       go to 100
    end if
1400 cea%Iplt = max(cea%Iplt, iplte)
    return
  end subroutine ROCKET


  subroutine read_libraries(cea)
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea
    integer:: i, j
    real(8):: xi, xln

    if (cea%invalid_case) return

    if (cea%Ions) then
       if (cea%Elmt(cea%Nlm) /= 'E') then
          cea%Nlm = cea%Nlm + 1
          cea%Elmt(cea%Nlm) = 'E'

          do concurrent (i = 1:maxMix, j = 1:Ncol)
             cea%points(i, j)%B0p(cea%Nlm, 1) = 0
             cea%points(i, j)%B0p(cea%Nlm, 2) = 0
          end do
       end if
    else if (cea%Elmt(cea%Nlm) == 'E') then
       cea%Nlm = cea%Nlm - 1
    end if

    cea%Jray(1:cea%Nreac) = 0

    call SEARCH(cea)

    if (cea%Ngc == 0) then
       cea%invalid_case = .true.
       return
    end if

    if (cea%Trnspt) call READTR(cea)

    ! INITIAL ESTIMATES
    cea%Npr = 0
    cea%Gonly = .true.
    cea%Enn = 0.1d0
    cea%Ennl = -2.3025851
    cea%Sumn = cea%Enn
    xi = cea%Ng
    if (xi == 0.) xi = 1
    xi = cea%Enn/xi
    xln = log(xi)

    cea%En(cea%Ng+1:cea%Ng+cea%Nc, 1) = 0
    cea%Enln(cea%Ng+1:cea%Ng+cea%Nc) = 0

    cea%En(1:cea%Ng, 1) = xi
    cea%Enln(1:cea%Ng) = xln

    if (cea%Nc /= 0 .and. cea%Nsert /= 0) then
       innerLoop: do i = 1, cea%Nsert
          do j = cea%Ngc, cea%Ngp1, - 1
             if (cea%Prod(j) == cea%ensert(i)) then
                cea%Npr = cea%Npr + 1
                cea%Jcond(cea%Npr) = j
                if (.not. cea%Short) write(cea%io_log, '(1X, A16, "INSERTED")') cea%Prod(j)
                cycle innerLoop
             end if
          end do
          write(cea%io_log, '(/" WARNING!!!", A16, "NOT FOUND FOR INSERTION")') cea%ensert(i)
       end do innerLoop
    end if

    return
  end subroutine read_libraries


  subroutine SEARCH(cea)
    !***********************************************************************
    ! SEARCH THERMO.LIB FOR THERMO DATA FOR SPECIES TO BE CONSIDERED.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(6):: date(maxNgc)
    character(2):: el(5)
    character(15):: sub
    integer:: i, j, k, ii
    integer:: i5, ifaz, itot, nall, ne, nint, ntgas, ntot
    real(8):: b(5), t1, t2, thermo(9, 3)
    integer:: io_thermo
    logical:: is_opened

    cea%Nc = 0
    ne = 0
    do i = 1, cea%Nlm
       cea%Jx(i) = 0
    end do
    do j = 1, maxNgc
       cea%S(j) = 0
       cea%H0(j) = 0
       cea%Deln(j) = 0
       do i = 1, cea%Nlm
          cea%A(i, j) = 0
       end do
    end do
    ! READ TEMPERATURE RANGES FOR COEFFICIENTS OF GASEOUS SPECIES.
    ! SOME DEFINITIONS:
    !   NTGAS = NUMBER OF GASEOUS SPECIES IN THERMO.LIB.
    !   NTOT =  NTGAS PLUS NUMBER OF TEMPERATURE INTERVALS FOR CONDENSED.
    !   NALL =  NTOT PLUS THE NUMBER OF REACTANT SPECIES IN THERMO.LIB.
    !   NG =    NUMBER OF GASES WITH STORED COEFFICIENTS.
    !   NC =    NUMBER OF CONDENSED INTERVALS WITH STORED COEFFICIENTS.
    !   NGC =    NG + NC
    !   THDATE = DATE READ FROM THERMO.INP FILE
    open(newunit = io_thermo, file = cea%filename_thermo_lib, status = 'old', form = 'unformatted', action = 'read')
    read(io_thermo) cea%Tg, ntgas, ntot, nall, cea%Thdate
    cea%Ngc = 1
    cea%Nc = 1

    ! BEGIN LOOP FOR READING SPECIES DATA FROM THERMO.LIB.
    outerLoop: do itot = 1, ntot
       if (itot > ntgas) then
          read(io_thermo) sub, nint, date(cea%Ngc), (el(j), b(j), j = 1, 5), cea%Ifz(cea%Nc), &
               cea%Temp(1, cea%Nc), cea%Temp(2, cea%Nc), cea%Mw(cea%Ngc), (cea%Cft(k, cea%Nc), k = 1, 9)
       else
          read(io_thermo) sub, nint, date(cea%Ngc), (el(j), b(j), j = 1, 5), ifaz, T1, T2, cea%Mw(cea%Ngc), thermo
       end if
       if (cea%Nonly /= 0) then
          i = 1
20        if (cea%Prod(i) /= sub .and. '*' // cea%Prod(i) /= sub) then
             i = i + 1
             if (i <= cea%Nonly) go to 20
             cycle outerLoop
          else
             if (sub == cea%Prod(cea%Ngc-1)) then
                cea%Nonly = cea%Nonly + 1
                do k = cea%Nonly, i+1, -1
                   cea%Prod(k) = cea%Prod(k-1)
                end do
             else
                cea%Prod(i) = cea%Prod(cea%Ngc)
             end if
             cea%Prod(cea%Ngc) = sub
          end if
       else if (cea%Nomit /= 0) then
          do i = 1, cea%Nomit
             if (cea%Omit(i) == sub .or. '*' // cea%Omit(i) == sub) cycle outerLoop
          end do
       end if
       innerLoop: do k = 1, 5
          if (b(k) == 0) exit
          do i = 1, cea%Nlm
             if (cea%Elmt(i) == el(k)) then
                cea%A(i, cea%Ngc) = b(k)
                cycle innerLoop
             end if
          end do
          do j = 1, cea%Nlm
             cea%A(j, cea%Ngc) = 0
          end do
          cycle outerLoop
       end do innerLoop
       cea%Prod(cea%Ngc) = sub
       if (itot > ntgas) then
          cea%Nc = cea%Nc + 1
          if (cea%Nc > maxNc) go to 400
       else
          cea%Ng = cea%Ngc
          if (cea%Ng > maxNg) go to 400
          do i = 1, 3
             do j = 1, 9
                cea%Coef(j, cea%Ng, i) = thermo(j, i)
             end do
          end do
          ! IF SPECIES IS AN ATOMIC GAS, STORE INDEX IN JX
          if (b(2) == 0 .and. b(1) == 1) then
             do i = 1, cea%Nlm
                if (cea%Elmt(i) == el(1)) then
                   ne = ne + 1
                   cea%Jx(i) = cea%Ngc
                   cea%Jcm(i) = cea%Ngc
                   go to 150
                end if
             end do
          end if
       end if
150    cea%Ngc = cea%Ngc + 1
       if (cea%Ngc > maxNgc) go to 400
    end do outerLoop

    close(io_thermo)

    ! FINISHED READING THERMO DATA FROM I/O UNIT io_thermo.
    cea%Ifz(cea%Nc) = 0
    cea%Nc = cea%Nc - 1
    cea%Ngc = cea%Ngc - 1
    cea%Ngp1 = cea%Ng + 1
    if (cea%Ngc < cea%Nonly) then
       do k = cea%Ngc+1, cea%Nonly
          write(cea%io_log, '(/" WARNING!!  ", a15, " NOT A PRODUCT IN thermo.lib FILE (SEARCH)")') cea%Prod(k)
       end do
    end if
    ! FIND MISSING ELEMENTS (IF ANY) FOR COMPONENTS
    cea%Nspx = cea%Ngc
    if (ne < cea%Nlm) then
       do i = 1, cea%Nlm
          if (cea%Nspx > maxNgc) go to 400
          if (cea%Jx(i) == 0) then
             cea%Nspx = cea%Nspx + 1
             do k = 1, cea%Nlm
                cea%A(k, cea%Nspx) = 0
             end do
             cea%A(i, cea%Nspx) = 1
             cea%Prod(cea%Nspx) = cea%Elmt(i)
             do k = 1, 100
                if (cea%Elmt(i) == atomic_symbol(k)) then
                   cea%Mw(cea%Nspx) = atomic_mass(k)
                   cea%Atwt(i) = atomic_mass(k)
                   cea%Cp(cea%Nspx) = 2.5d0
                   go to 210
                end if
             end do
210          cea%Jx(i) = cea%Nspx
             cea%Jcm(i) = cea%Nspx
          end if
       end do
    end if
    ! ARE ALL ELEMENTS IN PRODUCT SPECIES?
    outerLoop2: do i = 1, cea%Nlm
       do j = 1, cea%Ngc
          if (cea%A(i, j) /= 0) cycle outerLoop2
          ii = i
       end do
       write(cea%io_log, '(/" PRODUCT SPECIES CONTAINING THE ELEMENT", a3, " MISSING", &
            & //, 13x, "FATAL ERROR (SEARCH)")') cea%Elmt(ii)
       cea%Ngc = 0
       return
    end do outerLoop2
    ! WRITE POSSIBLE PRODUCT LIST
    if (.not. cea%Short) then
       write(cea%io_log, '(/2x, "SPECIES BEING CONSIDERED IN THIS SYSTEM", &
            & /" (CONDENSED PHASE MAY HAVE NAME LISTED SEVERAL TIMES)", &
            & /"  LAST thermo.inp UPDATE: ", a10, /)') cea%Thdate
       do i = 1, cea%Ngc, 3
          i5 = i + 2
          if (cea%Ngc < i5) i5 = cea%Ngc
          write(cea%io_log, '(3(2X, A6, 2X, A15))') (date(j), cea%Prod(j), j = i, i5)
       end do
    end if
    return

400 write(cea%io_log, '(/" INSUFFICIENT STORAGE FOR PRODUCTS-SEE RP-1311,", /"   PART 2, PAGE 39. (SEARCH)")')
    cea%Ngc = 0

    inquire(io_thermo, opened = is_opened)
    if (is_opened) close(io_thermo)

    return
  end subroutine SEARCH


  subroutine READTR(cea)
    ! SEARCH FOR TRANSPORT PROPERTIES FOR THIS CHEMICAL SYSTEM
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea

    character(16):: bin(2, 40), pure(6), spece(2)
    integer:: i, j, k, jj(2), jk, ir, lineb, npure, nrec
    real(8):: trdata(36)
    integer:: io_transport

    open(newunit = io_transport, file = trim(cea%filename_trans_lib), status = 'old', form = 'unformatted', action = 'read')

    cea%Ntape = 0
    npure = 0
    lineb = 1

    if (.not. cea%Short) write(cea%io_log, '(/" SPECIES WITH TRANSPORT PROPERTIES"//8X, "PURE SPECIES"/)')

    read(io_transport) nrec

    allocate(cea%transport_properties(nrec))

    do ir = 1, nrec
       read(io_transport) spece, trdata
       k = 1
       do j = 1, cea%Ng
          if (spece(k) == cea%Prod(j) .or. '*' // spece(k) == cea%Prod(j)) then
             jj(k) = j
             if (k == 2) then
                ! STORE NAMES FOR BINARIES IN BIN ARRAY.
                do k = 1, 2
                   bin(k, lineb) = spece(k)
                end do
                lineb = lineb + 1
                exit
             else
                jj(2) = j
                if (spece(2) == ' ') then
                   ! WRITE NAMES FOR PURE SPECIES.
                   npure = npure + 1
                   pure(npure) = spece(1)
                   exit
                else
                   k = 2
                   cycle
                end if
             end if
          end if
       end do

       if (j <= cea%Ng) then
          cea%Ntape = cea%Ntape + 1
          cea%transport_properties(cea%Ntape)%j = jj
          cea%transport_properties(cea%Ntape)%data = reshape(trdata, [6, 3, 2])
       end if

       if (npure /= 0 .and. (npure >= 6 .or. ir >= nrec)) then
          if (.not. cea%Short) write(cea%io_log, '(4(2x, A16))') (pure(jk), jk = 1, npure)
          npure = 0
       end if
    end do

    lineb = lineb - 1

    if (.not. cea%Short) then
       write(cea%io_log, '(/"     BINARY INTERACTIONS"/)')
       do j = 1, lineb
          write(cea%io_log, '(5X, 2A16)') (bin(i, j), i = 1, 2)
       end do
    end if

    write(cea%io_log, *)

    close(io_transport)

    return
  end subroutine READTR



  subroutine SETEN(cea)
    !***********************************************************************
    ! USE COMPOSITIONS FROM PREVIOUS POINT AS INITIAL ESTIMATES FOR
    ! CURRENT POINT NPT.  IF -
    !  ISV>0  USE COMPOSITIONS FROM POINT ISV.
    !  ISV<0  SAVE COMPOSITIONS FROM POINT -ISV FOR POSSIBLE LATER USE.
    !         ALSO USE COMPOSITIONS FROM POINT -ISV FOR NPT.
    !  ISV=0  USE COMPOSITIONS SAVED WHEN ISV<0.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: j

    if (cea%Isv > 0) then
       ! USE COMPOSITIONS FROM PREVIOUS POINT
       do j = 1, cea%Ngc
          cea%En(j, cea%Npt) = cea%En(j, cea%Isv)
       end do
    else if (cea%Isv < 0) then
       ! FIRST T--SAVE COMPOSITIONS FOR FUTURE POINTS WITH THIS T
       cea%Isv = -cea%Isv
       cea%Tsave = cea%Ttt(cea%Isv)
       cea%Ensave = cea%Enn
       cea%Enlsav = cea%Ennl
       cea%lsav = cea%lsave
       do j = 1, cea%Ng
          cea%Sln(j) = cea%Enln(j)
       end do
       do j = 1, cea%Ng
          cea%En(j, cea%Npt) = cea%En(j, cea%Isv)
       end do
       cea%Npr = 0
       do j = cea%Ngp1, cea%Ngc
          cea%Sln(j) = cea%En(j, cea%Isv)
          cea%En(j, cea%Npt) = cea%Sln(j)
          if (cea%Jliq == j) then
             cea%En(cea%Jsol, cea%Npt) = cea%En(cea%Jsol, cea%Isv) + cea%En(cea%Jliq, cea%Isv)
             cea%En(cea%Jliq, cea%Npt) = 0
             cea%Jsol = 0
             cea%Jliq = 0
             cea%Tsave = cea%Tsave - 5
             cea%Tt = cea%Tsave
             cea%Sln(j) = 0
          else if (cea%En(j, cea%Npt) > 0) then
             cea%Npr = cea%Npr + 1
             cea%Jcond(cea%Npr) = j
          end if
       end do
    else ! if (cea%Isv == 0) then
       ! NEXT POINT FIRST T IN SCHEDULE, USE PREVIOUS COMPOSITIONS FOR THIS T
       cea%Jsol = 0
       cea%Jliq = 0
       cea%Enn = cea%Ensave
       cea%Ennl = cea%Enlsav
       cea%lsave = cea%lsav
       cea%Npr = 0
       do j = cea%Ngp1, cea%Ngc
          cea%En(j, cea%Npt) = cea%Sln(j)
          if (cea%En(j, cea%Npt) > 0.d0) then
             cea%Npr = cea%Npr + 1
             cea%Jcond(cea%Npr) = j
          end if
       end do
       do j = 1, cea%Ng
          cea%En(j, cea%Npt) = 0
          cea%Enln(j) = cea%Sln(j)
          if (cea%Sln(j) /= 0) then
             if ((cea%Enln(j) - cea%Ennl + 18.5) > 0) cea%En(j, cea%Npt) = exp(cea%Enln(j))
          end if
       end do
       if (.not. cea%Tp) cea%Tt = cea%Tsave
       cea%Sumn = cea%Enn
    end if
  end subroutine SETEN



  subroutine SHCK(cea)
    !***********************************************************************
    ! PRIMARY ROUTINE FOR SHOCK PROBLEMS.
    !***********************************************************************
    use mod_types
    use mod_legacy_io
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(1):: cr12, cr52
    integer:: i, iof, it1, it2, itr, j, n
    logical:: refl, seql, srefl
    real(8):: ax, axx, b2, cormax, gg, hs, m2m1(Ncol), mis(13), mu12rt, p1, p21, &
         p21l, p2p1(Ncol), pmn, rho12, rho52, rrho(Ncol), sg(78), T1, T21, &
         t21l, t2t1(Ncol), ttmax, u1u2(Ncol), uis(13), utwo(Ncol), uu, wmx, ww

    type(CEA_Point), pointer:: p

    if (cea%Trace == 0) cea%Trace = 5.E-9
    cea%Tp = .true.
    cea%Cpmix = 0
    srefl = .false.
    if (.not. cea%Short) then
       write(IOOUT, '(/"   *** INPUT FOR SHOCK PROBLEMS ***")')
       write(IOOUT, '(/" INCDEQ =", L2, "   REFLEQ =", L2, "   INCDFZ =", L2, "    REFLFZ =", L2)') &
            cea%Incdeq, cea%Refleq, cea%Incdfz, cea%Reflfz
    end if
    if (cea%Refleq .or. cea%Reflfz) srefl = .true.
    seql = cea%Incdeq
    if (cea%T(1) == 0) cea%T(1) = cea%Rtemp(1)
    do i = 1, cea%Nsk
       uis(i) = cea%points(1, i)%U1
       mis(i) = cea%points(1, i)%Mach1
       if (cea%points(1, i)%Mach1 == 0 .and. cea%points(1, i)%U1 == 0) exit
    end do
    if (cea%Nsk > Ncol) then
       write(IOOUT, '(/" WARNING!!  ONLY ", I2, " u1 OR mach1 VALUES ALLOWED (SHCK)")') Ncol
       cea%Nsk = Ncol
    end if
    if (.not. cea%Short) then
       write(IOOUT, '(/1p, " U1 =   ", 5E13.6, /(8X, 5E13.6))') (cea%points(1, i)%U1, i = 1, cea%Nsk)
       write(IOOUT, '(/1p, " MACH1 =", 5E13.6, /(8X, 5E13.6))') (cea%points(1, i)%Mach1, i = 1, cea%Nsk)
    end if

    do iof = 1, cea%Nof
       cea%Oxfl = cea%Oxf(iof)
       call NEWOF(cea)
       cea%Incdeq = seql

       do
          refl = .false.
          it2 = 2
          it1 = 1
          cea%Pp = cea%P(1)
          cea%Tt = cea%T(1)

          if (.not. cea%Incdeq) then
             ! FROZEN
             do i = 1, cea%Nsk
                cea%Dlvtp(i) = 1
                cea%Dlvpt(i) = -1
             end do
          end if

          do i = 1, cea%Nsk
             p => cea%points(iof, i)

             cea%Ppp(i) = cea%P(i)
             cea%Ttt(i) = cea%T(i)

             if (i > 1) then
                if (cea%Ppp(i) == 0) cea%Ppp(i) = cea%Ppp(i-1)
                if (cea%Ttt(i) == 0) cea%Ttt(i) = cea%Ttt(i-1)
                cea%Ssum(i) = cea%Ssum(i-1)
                cea%Hsum(i) = cea%Hsum(i-1)
             end if

             if (i == 1 .or. cea%Ttt(i) /= cea%Tt .or. cea%Ppp(i) /= cea%Pp) then
                cea%Pp = cea%Ppp(i)
                cea%Tt = cea%Ttt(i)
                if (cea%Tt >= cea%Tg(1) * 0.8d0) then
                   cea%Npt = i
                   call HCALC(cea)
                   cea%Hsum(i) = cea%Hsub0
                else
                   write(IOOUT, '(/" TEMPERATURE=", E12.4, " IS OUT OF EXTENDED RANGE ", "FOR POINT", I5, " (SHCK)")') cea%Tt, i
                   return
                end if
             end if

             if (cea%Cpmix /= 0) cea%Gamma1 = cea%Cpmix / (cea%Cpmix - 1/cea%Wmix)
             cea%A1 = sqrt(R0 * cea%Gamma1 * cea%Tt / cea%Wmix)
             if (p%U1 == 0) p%U1 = cea%A1 * p%Mach1
             if (p%Mach1 == 0) p%Mach1 = p%U1 / cea%A1
             cea%Wm(i) = cea%Wmix
             cea%Cpr(i) = cea%Cpmix
             cea%Gammas(i) = cea%Gamma1
             cea%Vlm(i) = R0 * cea%Tt / (cea%Wmix * cea%Pp)
          end do

          ! OUTPUT--1ST CONDITION
          write(IOOUT, '(////25X, "SHOCK WAVE PARAMETERS ASSUMING")')
          if (.not. cea%Incdeq) then
             write(IOOUT, '(/, 17X, " FROZEN COMPOSITION FOR INCIDENT SHOCKED CONDITI1ONS"//)')
          else
             write(IOOUT, '(/, 16X, " EQUILIBRIUM COMPOSITION FOR INCIDENT SHOCKED CONDITIONS"//)')
          end if
          cea%Eql = .false.
          call OUT1(cea)
          write(IOOUT, '(/" INITIAL GAS (1)")')
          cea%fmt(4) = '13'
          cea%fmt(5) = ' '
          cea%fmt(7) = '4,'
          write(IOOUT, cea%fmt) 'MACH NUMBER1   ', (cea%points(iof, j)%Mach1, j = 1, cea%Nsk)
          cea%fmt(7) = '2,'
          write(IOOUT, cea%fmt) 'U1, M/SEC      ', (cea%points(iof, j)%U1, j = 1, cea%Nsk)
          call OUT2(cea, cea%Nsk)

          ! BEGIN CALCULATIONS FOR 2ND CONDITION
          if (cea%Incdeq) cea%Eql = .true.

          cea%ipt = 1
          do
             p => cea%points(iof, cea%ipt)

             cea%Gamma1 = cea%Gammas(cea%ipt)
             uu = p%U1
             wmx = cea%Wm(cea%ipt)
             p1 = cea%Ppp(cea%ipt)
             T1 = cea%Ttt(cea%ipt)
             hs = cea%Hsum(cea%ipt)
             if (refl) uu = u1u2(cea%ipt)
             mu12rt = wmx * uu**2 / (R0 * T1)
             if (refl) then
                ! REFLECTED--SUBSCRIPTS 2=1, 5=2, P52=P21
                T21 = 2
                b2 = (-1 - mu12rt - T21) / 2
                p21 = -b2 + sqrt(b2**2 - T21)
             else
                p21 = (2 * cea%Gamma1 * p%Mach1**2 - cea%Gamma1 + 1) / (cea%Gamma1 + 1)
                ! THE FOLLOWING IMPROVED FORMULATION FOR THE INITIAL ESTIMATE FOR THE
                ! 2ND CONDITION WAS MADE AND TESTED BY S. GORDON 7/10/89.
                if (.not. cea%Eql) then
                   T21 = p21 * (2 / p%Mach1**2 + cea%Gamma1 - 1) / (cea%Gamma1 + 1)
                else
                   cea%Pp = p21 * p1
                   cea%Tp = .false.
                   cea%Hp = .true.
                   cea%Hsub0 = hs + uu**2 / (2 * R0)

                   cea%Npt = cea%ipt
                   call EQLBRM(cea)
                   cea%ipt = cea%Npt

                   T21 = cea%Ttt(cea%ipt) / T1
                   cea%Hp = .false.
                   cea%Tp = .true.
                end if
             end if
             p21l = log(p21)
             ttmax = 1.05 * cea%Tg(4) / T1
             T21 = min(T21, ttmax)
             t21l = log(T21)

             ax = 1
             do itr = 1, 60
                if (cea%Shkdbg) write(IOOUT, '(/" ITR NO.=", I3, 3X, "P", I1, "/P", I1, " =", F9.4, 3X, "T", I1, &
                     & "/T", I1, " =", F9.4, "   RHO2/RHO1 =", F9.6)') itr, it2, it1, p21, it2, it1, T21, rho52
                cea%Tt = T21 * T1
                cea%Pp = p21 * p1
                if (.not. cea%Eql) then
                   ! FROZEN
                   cea%Tln = log(cea%Tt)
                   if (.not. cea%Incdeq) then
                      cea%Npt = cea%ipt
                      call HCALC(cea)
                      if (cea%Tt == 0) exit
                      cea%Hsum(cea%ipt) = cea%Hsub0
                      cea%Cpr(cea%ipt) = cea%Cpmix
                   else
                      call CPHS(cea)
                      cea%Cpr(cea%ipt) = cea%Cpsum
                      cea%Hsum(cea%ipt) = 0
                      do j = 1, cea%Ng
                         cea%Hsum(cea%ipt) = cea%Hsum(cea%ipt) + cea%H0(j) * cea%En(j, cea%ipt)
                      end do
                      cea%Hsum(cea%ipt) = cea%Hsum(cea%ipt) * cea%Tt
                   end if
                else
                   cea%Npt = cea%ipt
                   call EQLBRM(cea)
                   cea%ipt = cea%Npt
                   if (cea%Tt == 0) exit
                end if
                rho12 = wmx * T21 / (cea%Wm(cea%ipt) * p21)
                gg = rho12 * mu12rt
                rho52 = 1 / rho12
                if (refl) gg = -mu12rt * rho52 / (rho52 - 1)**2
                cea%G(1, 1) = -gg * cea%Dlvpt(cea%ipt) - p21
                cea%G(1, 2) = -gg * cea%Dlvtp(cea%ipt)
                cea%G(1, 3) = p21 - 1 + gg - mu12rt
                if (refl) cea%G(1, 3) = p21 - 1 + gg * (rho52 - 1)
                gg = gg * T1 / wmx
                if (.not. refl) gg = gg * rho12
                cea%G(2, 1) = -gg * cea%Dlvpt(cea%ipt) + cea%Tt * (cea%Dlvtp(cea%ipt) - 1) / cea%Wm(cea%ipt)
                cea%G(2, 2) = -gg * cea%Dlvtp(cea%ipt) - cea%Tt * cea%Cpr(cea%ipt)
                gg = 1 - rho12**2
                if (refl) gg = (rho52 + 1) / (rho52 - 1)
                cea%G(2, 3) = cea%Hsum(cea%ipt) - hs - uu**2 * gg / (2 * R0)
                cea%X(3) = cea%G(1, 1) * cea%G(2, 2) - cea%G(1, 2) * cea%G(2, 1)
                cea%X(1) = (cea%G(1, 3) * cea%G(2, 2) - cea%G(2, 3) * cea%G(1, 2)) / cea%X(3)
                cea%X(2) = (cea%G(1, 1) * cea%G(2, 3) - cea%G(2, 1) * cea%G(1, 3)) / cea%X(3)
                if (cea%Shkdbg) then
                   write(IOOUT, '(/" G(I,J)  ", 3E15.8)') cea%G(1, 1), cea%G(1, 2), cea%G(1, 3)
                   write(IOOUT, '(/" G(I,J)  ", 3E15.8)') cea%G(2, 1), cea%G(2, 2), cea%G(2, 3)
                   write(IOOUT, '(/" X       ", 2E15.8)') cea%X(1), cea%X(2)
                   write(IOOUT, '(/" HSUM HS UU U2 ", 4E15.8)') cea%Hsum(cea%ipt), hs, uu, uu * rho12
                end if
                ax = abs(cea%X(1))
                axx = abs(cea%X(2))
                if (axx > ax) ax = axx

                if (ax < 0.00005) exit

                cormax = 0.40546511
                if (itr > 4) cormax = 0.22314355
                if (itr > 12) cormax = 0.09531018
                if (itr > 20) cormax = 0.04879016
                ax = ax / cormax
                if (ax > 1) then
                   cea%X(1) = cea%X(1) / ax
                   cea%X(2) = cea%X(2) / ax
                end if
                p21l = p21l + cea%X(1)
                t21l = t21l + cea%X(2)
                p21 = exp(p21l)
                T21 = exp(t21l)

                if (cea%Shkdbg) write(IOOUT, '(/" MAX.COR.=", e13.6, " X(1)=", e13.6, " X(2)=", e13.6)') cormax, cea%X(1), cea%X(2)

                if (itr == 1 .and. T21 >= ttmax) then
                   cea%Tt = 0
                   cea%ipt = cea%ipt - 1
                end if
             end do

             p => cea%points(iof, cea%ipt)

             if (.not. cea%Eql .or. cea%Tt /= 0) then
                if (itr > 60) then
                   write(IOOUT, '(/6x, " WARNING!!  NO CONVERGENCE FOR u1=", F8.1, &
                        & /"  ANSWERS NOT RELIABLE, SOLUTION MAY NOT EXIST (SHCK)")') p%U1
                end if

                if ((itr > 60) .or. (ax < 0.00005) .or. ((.not. cea%Eql) .and. (.not. cea%Incdeq) .and. (cea%Tt == 0))) then
                   ! CONVERGED OR TOOK 60 ITERATIONS WITHOUT CONVERGING.
                   ! STORE RESULTS.
                   rrho(cea%ipt) = rho52
                   m2m1(cea%ipt) = cea%Wm(cea%ipt) / wmx
                   p2p1(cea%ipt) = p21
                   t2t1(cea%ipt) = T21
                   utwo(cea%ipt) = uu * rho12
                   u1u2(cea%ipt) = uu - utwo(cea%ipt)
                   if (cea%Tt >= cea%Tg(1) * 0.8d0 .and. cea%Tt <= cea%Tg(4) * 1.1d0) then
                      if (.not. cea%Eql) then
                         ! FROZEN
                         cea%Ppp(cea%ipt) = cea%Pp
                         cea%Ttt(cea%ipt) = cea%Tt
                         cea%Gammas(cea%ipt) = cea%Cpr(cea%ipt) / (cea%Cpr(cea%ipt) - 1 / wmx)
                         cea%Vlm(cea%ipt) = R0 * cea%Tt / (wmx * cea%Pp)
                         if (cea%Incdeq) then
                            cea%Ssum(cea%ipt) = 0
                            do j = 1, cea%Ngc
                               pmn = cea%Pp * wmx * cea%En(j, cea%ipt)
                               if (cea%En(j, cea%ipt) > 0) cea%Ssum(cea%ipt) = cea%Ssum(cea%ipt) + cea%En(j, cea%ipt) * (cea%S(j) - log(pmn))
                            end do
                         end if
                      end if
                      go to 900
                   end if
                end if

                write(IOOUT, '(/" TEMPERATURE=", E12.4, " IS OUT OF EXTENDED RANGE ", &
                     & "FOR POINT", I5, " (SHCK)")') cea%Tt, cea%ipt
                cea%Tt = 0
             end if

             if (cea%ipt < 1) return

             cea%Nsk = cea%ipt

900          if (cea%Trnspt) then
                cea%Npt = cea%ipt
                call TRANP(cea)
             end if

             cea%Isv = 0
             if (cea%ipt < cea%Nsk) cea%Isv = cea%ipt
             if (cea%ipt == 1) cea%Isv = -1
             cea%ipt = cea%ipt + 1

             if (cea%Eql) then
                cea%Npt = cea%ipt
                write(stderr, *) '[DEBUG] SHCK: call SETEN(', cea%Isv, cea%Npt, ')'
                call SETEN(cea)
             end if

             if (cea%ipt <= cea%Nsk) cycle
             cea%ipt = cea%Nsk
             if (refl) then
                if (.not. cea%Eql) write(IOOUT, '(/" SHOCKED GAS (5)--REFLECTED--FROZEN")')
                if (cea%Eql) write(IOOUT, '(/" SHOCKED GAS (5)--REFLECTED--EQUILIBRIUM")')
                cr12 = '2'
                cr52 = '5'
             else
                if (.not. cea%Eql) write(IOOUT, '(/" SHOCKED GAS (2)--INCIDENT--FROZEN")')
                if (cea%Eql) write(IOOUT, '(/" SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM")')
                cr12 = '1'
                cr52 = '2'
             end if
             cea%fmt(7) = '2,'
             write(IOOUT, cea%fmt) 'U' // cr52 // ', M/SEC      ', (utwo(j), j = 1, cea%ipt)

             call OUT2(cea, cea%ipt)

             if (cea%Trnspt) call OUT4(cea, cea%ipt)
             write(IOOUT, *)
             cea%fmt(7) = '3,'
             write(IOOUT, cea%fmt) 'P' // cr52 // '/P' // cr12 // '           ', (p2p1(j), j = 1, cea%ipt)
             write(IOOUT, cea%fmt) 'T' // cr52 // '/T' // cr12 // '           ', (t2t1(j), j = 1, cea%ipt)
             cea%fmt(7) = '4,'
             write(IOOUT, cea%fmt) 'M' // cr52 // '/M' // cr12 // '           ', (m2m1(j), j = 1, cea%ipt)
             write(IOOUT, cea%fmt) 'RHO' // cr52 // '/RHO' // cr12 // '       ', (rrho(j), j = 1, cea%ipt)
             cea%fmt(7) = '2,'
             if (.not. refl) write(IOOUT, cea%fmt) 'V2, M/SEC      ', (u1u2(j), j = 1, cea%ipt)
             if (refl) write(IOOUT, cea%fmt) 'U5+V2,M/SEC    ', (u1u2(j), j = 1, cea%ipt)
             if (.not. cea%Eql) then
                ! WRITE FROZEN MOLE (OR MASS) FRACTIONS
                cea%fmt(7) = '5,'
                if (.not. cea%Incdeq) then
                   if (cea%Massf) then
                      write(IOOUT, '(/1x, A4, " FRACTIONS"/)') 'MASS'
                   else
                      write(IOOUT, '(/1x, A4, " FRACTIONS"/)') 'MOLE'
                      ww = wmx
                   end if
                   do n = 1, cea%Nreac
                      j = cea%Jray(n)
                      if (cea%Massf) ww = cea%Mw(j)
                      write(IOOUT, '(" ", A16, F8.5, 12F9.5)') cea%Prod(j), (cea%En(j, i) * ww, i = 1, cea%ipt)
                   end do
                else
                   cea%Eql = .true.
                   call OUT3(cea, cea%ipt)
                   cea%Eql = .false.
                end if
             else
                call OUT3(cea, cea%ipt)
             end if

             cea%Iplt = min(cea%Iplt + cea%ipt, 500)

             if (srefl) then
                if (.not. refl) then
                   refl = .true.
                   it2 = 5
                   it1 = 2
                   cea%Eql = .true.
                   if (cea%Reflfz) then
                      cea%Eql = .false.
                      if (cea%Refleq) then
                         j = 1
                         do i = 1, cea%ipt
                            sg(j)   = u1u2(i)
                            sg(j+1) = cea%Wm(i)
                            sg(j+2) = cea%Ppp(i)
                            sg(j+3) = cea%Ttt(i)
                            sg(j+4) = cea%Hsum(i)
                            sg(j+5) = cea%Gammas(i)
                            j = j + 6
                         end do
                      end if
                   end if
                   cea%ipt = 1
                   cycle
                else if (.not. cea%Eql .and. cea%Refleq) then
                   j = 1
                   do i = 1, cea%ipt
                      u1u2(i) = sg(j)
                      cea%Wm(i) = sg(j+1)
                      cea%Ppp(i) = sg(j+2)
                      cea%Ttt(i) = sg(j+3)
                      cea%Hsum(i) = sg(j+4)
                      cea%Gammas(i) = sg(j+5)
                      j = j + 6
                   end do
                   cea%Eql = .true.
                   cea%ipt = 1
                   cycle
                end if
             end if

             exit
          end do

          if (.not. cea%Incdeq .or. .not. cea%Incdfz) exit

          cea%Incdeq = .false.
          cea%Eql = .false.
       end do

       do concurrent (i = 1:cea%Nsk)
          cea%points(iof, i)%U1 = uis(i)
          cea%points(iof, i)%Mach1 = mis(i)
       end do

#ifndef NDEBUG
       call print_debug(cea, iof)
#endif
    end do

    cea%Tp = .false.
    do n = 1, cea%Nreac
       cea%Rtemp(n) = cea%T(1)
    end do

    return
  end subroutine SHCK



  subroutine THERMP(cea)
    !***********************************************************************
    ! ASSIGNED THERMODYNAMIC STATES.  HP, SP, TP, UV, SV, AND TV PROBLEMS.
    !***********************************************************************
    use mod_types
    use mod_legacy_io
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: iof
    integer:: ip, it

#ifndef NDEBUG
    write(stderr, *) '[DEBUG] ', trim(cea%Case)
#endif

    cea%Eql = .true.
    iofLoop: do iof = 1, cea%Nof
       cea%Oxfl = cea%Oxf(iof)
       call NEWOF(cea)
       ! SET ASSIGNED P OR VOLUME
       do ip = 1, cea%Np
          cea%Pp = cea%P(ip)
          ! SET ASSIGNED T
          do it = 1, cea%Nt
             cea%ipt = (ip - 1) * cea%Nt + it
             cea%Npt = mod(cea%ipt - 1, Ncol) + 1

             if (iof > 1 .or. cea%ipt > 1) then
                if (.not. cea%Tp .and. cea%Tt /= 0) cea%T(1) = cea%Tt

                if (cea%ipt == 2) then
                   cea%Isv = -1
                else if (cea%Nt /= 1 .and. (it == 1 .or. cea%Tt == 0.)) then
                   cea%Isv = 0
                else
                   cea%Isv = mod(cea%ipt - 2, Ncol) + 1
                end if
                call SETEN(cea)
             end if

             cea%Vv = cea%V(ip)
             cea%Tt = cea%T(it)
             call EQLBRM(cea)

             if (cea%Npt == 0) then
#ifndef NDEBUG
                call print_debug(cea, iof)
#endif
                return
             end if

             if (cea%Trnspt .and. cea%Tt /= 0) call TRANP(cea)

             if (ip /= cea%Np .or. it /= cea%Nt .and. cea%Tt /= 0) then
                if (cea%Npt /= Ncol) cycle
             end if
             if (.not. cea%Hp) write(IOOUT, '(////15X, "THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED")')
             if (cea%Hp) write(IOOUT, '(////9X, "THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED")')
             if (.not. cea%Vol) then
                if (cea%Hp) write(IOOUT, '(/34X, " PRESSURES"/)')
                if (cea%Tp) write(IOOUT, '(/27X, "TEMPERATURE AND PRESSURE"/)')
                if (cea%Sp) write(IOOUT, '(/29X, "ENTROPY AND PRESSURE"/)')
             else
                if (cea%Hp) write(IOOUT, '(/36X, " VOLUME"/)')
                if (cea%Tp) write(IOOUT, '(/28X, "TEMPERATURE AND VOLUME"/)')
                if (cea%Sp) write(IOOUT, '(/30X, "ENTROPY AND VOLUME"/)')
             end if
             call OUT1(cea)
             write(IOOUT, '(/" THERMODYNAMIC PROPERTIES"/)')
             call OUT2(cea, cea%Npt)
             if (cea%Trnspt) call OUT4(cea, cea%Npt)
             call OUT3(cea, cea%Npt)
             cea%Iplt = min(cea%Iplt + cea%Npt, 500)

             if ((ip == cea%Np .and. it == cea%Nt .or. cea%Tt == 0) .and. iof == cea%Nof) then
#ifndef NDEBUG
                call print_debug(cea, iof)
#endif
                return
             end if

             write(IOOUT, '(////)')
          end do
       end do
    end do iofLoop

    return
  end subroutine THERMP



  subroutine TRANIN(cea)
    !***********************************************************************
    ! BRINGS IN AND SORTS OUT INPUT FOR TRANSPORT CALCULATIONS
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, ii, inds(maxTr), ir, j, jtape(2), k, k1, k2, kt, kvc, l, loop, m, nms
    logical:: change, elc1, elc2, ion1, ion2, setx
    real(8), save:: coeff, debye, ekt, enel, enmin, ionic, lambda, omega, prop, qc, ratio, &
         stcf(maxTr, maxTr), stcoef(maxTr), te, testen, testot, total, &
         trc(6, 3, 2), wmols(maxTr), wmred, xsel, xss(maxTr)

    if (.not. cea%Eql) then
       if (.not. cea%Shock) then
          if (.not. setx) then
             setx = .true.
             cea%Nm = nms
             do i = 1, cea%Nm
                cea%Xs(i) = xss(i)
                cea%Wmol(i) = wmols(i)
                cea%Ind(i) = inds(i)
             end do
          end if
          go to 300
       else if (.not. cea%Incdeq) then
          if (cea%Npt <= 1) then
             cea%Nm = cea%Nreac
             do i = 1, cea%Nm
                j = cea%Jray(i)
                cea%Ind(i) = j
                cea%Wmol(i) = cea%Mw(j)
                cea%Xs(i) = cea%En(j, 1) * cea%Wm(1)
             end do
          end if
          go to 300
       end if
    end if

    ! PICK OUT IMPORTANT SPECIES
    cea%Nm = 0
    total = 0
    enmin = 1.0d-11 / cea%Wm(cea%Npt)
    testot = 0.999999999d0 / cea%Wm(cea%Npt)
    do i = 1, cea%Lsave
       j = cea%Jcm(i)
       if (cea%En(j, cea%Npt) <= 0 .and. j <= cea%Ngc) then
          if ((cea%Enln(j) - cea%Ennl + 25.328436d0) > 0) cea%En(j, cea%Npt) = exp(cea%Enln(j))
       end if
       cea%Nm = cea%Nm + 1
       cea%Ind(cea%Nm) = j
       total = total + cea%En(j, cea%Npt)
       if (cea%Mw(j) < 1) enel = cea%En(j, cea%Npt)
       cea%En(j, cea%Npt) = -cea%En(j, cea%Npt)
    end do
    testen = 1 / (cea%Ng * cea%Wm(cea%Npt))

    if (total <= testot) then
       outerLoop1: do loop = 1, cea%Ng
          testen = testen / 10

          do j = 1, cea%Ng
             if (cea%En(j, cea%Npt) >= testen) then
                if (cea%Nm >= maxTr) then
                   write(IOOUT, '(/" WARNING!!  MAXIMUM ALLOWED NO. OF SPECIES", I3, " WAS USED IN ", &
                        & /" TRANSPORT PROPERTY CALCULATIONS FOR POINT", I3, "(TRANIN))")') cea%Nm, cea%Npt
                   exit outerLoop1
                else
                   total = total + cea%En(j, cea%Npt)
                   cea%Nm = cea%Nm + 1
                   cea%Ind(cea%Nm) = j
                   cea%En(j, cea%Npt) = -cea%En(j, cea%Npt)
                end if
             end if
          end do

          if (testen <= enmin) exit
       end do outerLoop1
    end if

    ! CALCULATE MOLE FRACTIONS FROM THE EN(J, NPT)
    do j = 1, cea%Ng
       cea%En(j, cea%Npt) = abs(cea%En(j, cea%Npt))
    end do
    do i = 1, cea%Nm
       j = cea%Ind(i)
       cea%Wmol(i) = cea%Mw(j)
       cea%Xs(i) = cea%En(j, cea%Npt) / total
    end do
    if (cea%Npt == cea%Nfz) then
       nms = cea%Nm
       do i = 1, cea%Nm
          xss(i) = cea%Xs(i)
          wmols(i) = cea%Wmol(i)
          inds(i) = cea%Ind(i)
       end do
       setx = .false.
    end if
    ! REWRITE REACTIONS TO ELIMINATE TRACE ELEMENTS
    cea%Nr = cea%Nm - cea%Lsave
    if (cea%Nr /= 0) then
       do k = 1, maxTr
          do m = 1, maxTr
             cea%Stc(k, m) = 0
          end do
       end do
       k = 1
       do i = cea%Lsave + 1, cea%Nm
          cea%Stc(k, i) = -1
          j = cea%Ind(i)
          do m = 1, cea%Lsave
             cea%Stc(k, m) = cea%A(m, j)
          end do
          k = k + 1
       end do
       do i = 1, cea%Nm
          if (cea%Xs(i) < 1.0d-10) then
             m = 1
             change = .false.

             do j = 1, cea%Nr
                coeff = cea%Stc(j, i)
                if (abs(coeff) > 1.0d-05) then
                   if (.not. change) then
                      change = .true.
                      do k = 1, cea%Nm
                         stcoef(k) = cea%Stc(j, k) / coeff
                      end do
                      cycle
                   else
                      do k = 1, cea%Nm
                         cea%Stc(j, k) = cea%Stc(j, k) / coeff - stcoef(k)
                      end do
                   end if
                end if
                do k = 1, cea%Nm
                   stcf(m, k) = cea%Stc(j, k)
                end do
                m = m + 1
             end do

             do ii = 1, cea%Nm
                do j = 1, cea%Nr
                   cea%Stc(j, ii) = stcf(j, ii)
                end do
             end do
             cea%Nr = m - 1
          end if
       end do
    end if

    ! FIND TRANSPORT DATA FOR IMPORTANT INTERACTIONS
300 do i = 1, cea%Nm
       cea%Con(i) = 0
       do j = 1, cea%Nm
          cea%Eta(i, j) = 0
       end do
    end do

    outerLoop2: do ir = 1, cea%Ntape
       jtape = cea%transport_properties(ir)%j
       trc = cea%transport_properties(ir)%data

       innerLoop: do k = 1, 2
          do i = 1, cea%Nm
             j = cea%Ind(i)
             if (j == jtape(k)) then
                l = i
                if (k == 2) then
                   kvc = 1

                   do
                      kt = 1
                      if (trc(2, 1, kvc) /= 0) then
                         if (trc(2, 2, kvc) /= 0) then
                            if (cea%Tt > trc(2, 1, kvc)) kt = 2
                            if (trc(2, 3, kvc) /= 0) then
                               if (cea%Tt > trc(2, 2, kvc)) kt = 3
                            end if
                         end if

                         prop = exp(trc(6, kt, kvc) + (trc(5, kt, kvc) / cea%Tt + trc(4, kt, kvc)) / cea%Tt &
                              + trc(3, kt, kvc) * cea%Tln)

                         if (kvc == 2) then
                            cea%Con(l) = prop
                            cycle outerLoop2
                         else
                            cea%Eta(l, m) = prop
                            if (l /= m) cea%Eta(m, l) = cea%Eta(l, m)
                         end if

                      else if (kvc == 2) then
                         cycle outerLoop2
                      end if

                      kvc = 2
                   end do

                else
                   m = i
                   cycle innerLoop
                end if
             end if
          end do

          exit
       end do innerLoop
    end do outerLoop2

    ! MAKE ESTIMATES FOR MISSING DATA
    !
    ! INCLUDES ION CROSS SECTION ESTIMATES
    ! ESTIMATES FOR  E-ION, ION-ION, E-NEUTRAL, ION-NEUTRAL
    ! DEBYE SHIELDING WITH IONIC CUTOFF DISTANCE
    if (cea%Ions) then
       te = cea%Tt / 1000
       ekt = 4.8032d0**2 / (Boltz * te)
       qc = 100 * ekt**2
       xsel = enel / total
       if (xsel < 1.0d-12) xsel = 1.0d-12
       debye = ((22.5d0 / pi) * (R0 / Avgdr * 100) * (te/xsel)) / ekt**3
       ionic = ((810 / (4*pi)) * (R0 / Avgdr * 100d0) * (te/xsel))**(2/3.) / ekt**2
       lambda = sqrt(debye + ionic)
       lambda = max(lambda, 2.71828183d0)
    end if
    do i = 1, cea%Nm
       k = cea%Ind(i)
       cea%Cprr(i) = cea%Cp(k)
       if (.not. (cea%Ions .and. (abs(cea%A(cea%Nlm, k)) == 1) .and. (cea%Eta(i, i) == 0))) then
          if (cea%Eta(i, i) == 0) then
             omega = log(50 * cea%Wmol(i)**4.6 / cea%Tt**1.4)
             omega = max(omega, 1d0)
             cea%Eta(i, i) = cea%Viscns * sqrt(cea%Wmol(i) * cea%Tt) / omega
          end if
          if (cea%Con(i) == 0) cea%Con(i) = cea%Eta(i, i) * R0 * (0.00375d0 + 0.00132d0 * (cea%Cprr(i) - 2.5d0)) / cea%Wmol(i)
       end if
    end do

    do i = 1, cea%Nm
       do j = i, cea%Nm
          ion1 = .false.
          ion2 = .false.
          elc1 = .false.
          elc2 = .false.
          omega = 0
          if (cea%Eta(i, j) == 0) cea%Eta(i, j) = cea%Eta(j, i)
          if (cea%Eta(j, i) == 0) cea%Eta(j, i) = cea%Eta(i, j)
          if (cea%Eta(i, j) == 0) then
             if (cea%Ions) then
                ! ESTIMATE FOR IONS
                k1 = cea%Ind(i)
                k2 = cea%Ind(j)
                if (abs(cea%A(cea%Nlm, k1)) == 1) ion1 = .true.
                if (abs(cea%A(cea%Nlm, k2)) == 1) ion2 = .true.
                if (cea%Wmol(i) < 1) elc1 = .true.
                if (cea%Wmol(j) < 1) elc2 = .true.
                if (ion1 .and. ion2) omega = 1.36d0 * qc * log(lambda)
                if ((ion1 .and. elc2) .or. (ion2 .and. elc1)) omega = 1.29d0 * qc * log(lambda)
                if ((ion1 .and. .not. ion2) .or. (ion2 .and. .not. ion1)) omega = exp(6.776 - 0.4 * cea%Tln)
                if (omega /= 0) then
                   wmred = sqrt(2 * cea%Tt * cea%Wmol(i) * cea%Wmol(j) / (cea%Wmol(i) + cea%Wmol(j)))
                   cea%Eta(i, j) = cea%Viscns * wmred * pi / omega
                   cea%Eta(j, i) = cea%Eta(i, j)
                   if (i == j) then
                      cea%Cprr(i) = cea%Cp(k1)
                      cea%Con(i) = cea%Eta(i, i) * R0 * (0.00375d0 + 0.00132d0 * (cea%Cprr(i) - 2.5d0)) / cea%Wmol(i)
                   end if
                   cycle
                end if
             end if
             ! ESTIMATE FOR UNLIKE INTERACTIONS FROM RIGID SPHERE ANALOGY
             ratio = sqrt(cea%Wmol(j) / cea%Wmol(i))
             cea%Eta(i, j) = 5.656854d0 * cea%Eta(i, i) * sqrt(cea%Wmol(j) / (cea%Wmol(i) + cea%Wmol(j)))
             cea%Eta(i, j) = cea%Eta(i, j) / (1 + sqrt(ratio * cea%Eta(i, i) / cea%Eta(j, j)))**2
             cea%Eta(j, i) = cea%Eta(i, j)
          end if
       end do
    end do

    return
  end subroutine TRANIN



  subroutine TRANP(cea)
    !***********************************************************************
    ! CALCULATES GAS TRANSPORT PROPERTIES
    !
    !   NUMBER OF GASEOUS SPECIES = NM   (MAXIMUM maxTr)
    !   NUMBER OF CHEMICAL REACTIONS = NR (NM - NLM)
    !   ARRAY OF STOICHIOMETRIC COEFFICIENTS = STC
    !***********************************************************************
    use mod_types
    use mod_general
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, i1, j, jj, k, m, mm, nlmm, nmm
    real(8), save:: cpreac, delh(maxTr), gmat(maxMat, maxMat+1), phi(maxTr, maxTr), &
         psi(maxTr, maxTr), reacon, rtpd(maxTr, maxTr), stx(maxTr), &
         stxij(maxTr, maxTr), sumc, sumv, wtmol, xskm(maxTr, maxTr)

    call TRANIN(cea)
    ! CALCULATE VISCOSITY AND FROZEN THERMAL CONDUCTIVITY
    nmm = cea%Nm - 1
    do i = 1, cea%Nm
       rtpd(i, i) = 0
       phi(i, i) = 1
       psi(i, i) = 1
    end do
    cea%Confro(cea%Npt) = 0
    cea%Vis(cea%Npt) = 0
    do i = 1, nmm
       i1 = i + 1
       !DIR$ IVDEP
       do j = i1, cea%Nm
          sumc = 2 / (cea%Eta(i, j) * (cea%Wmol(i) + cea%Wmol(j)))
          phi(i, j) = sumc * cea%Wmol(j) * cea%Eta(i, i)
          phi(j, i) = sumc * cea%Wmol(i) * cea%Eta(j, j)
          sumc = (cea%Wmol(i) + cea%Wmol(j))**2
          psi(i, j) = phi(i, j) * (1 + 2.41d0 * (cea%Wmol(i) - cea%Wmol(j)) * (cea%Wmol(i) - 0.142d0 * cea%Wmol(j)) / sumc)
          psi(j, i) = phi(j, i) * (1 + 2.41d0 * (cea%Wmol(j) - cea%Wmol(i)) * (cea%Wmol(j) - 0.142d0 * cea%Wmol(i)) / sumc)
       end do
    end do
    do i = 1, cea%Nm
       sumc = 0
       sumv = 0
       do j = 1, cea%Nm
          sumc = sumc + psi(i, j) * cea%Xs(j)
          sumv = sumv + phi(i, j) * cea%Xs(j)
       end do
       cea%Vis(cea%Npt) = cea%Vis(cea%Npt) + cea%Eta(i, i) * cea%Xs(i) / sumv
       cea%Confro(cea%Npt) = cea%Confro(cea%Npt) + cea%Con(i) * cea%Xs(i) / sumc
    end do
    if (cea%Eql .and. cea%Nr > 0) then
       ! CALCULATE REACTION HEAT CAPACITY AND THERMAL CONDUCTIVITY
       m = cea%Nr + 1
       do i = 1, cea%Nr
          delh(i) = 0
          do k = 1, cea%Lsave
             j = cea%Jcm(k)
             delh(i) = cea%Stc(i, k) * cea%H0(j) + delh(i)
          end do
          nlmm = cea%Lsave + 1
          do k = nlmm, cea%Nm
             j = cea%Ind(k)
             delh(i) = cea%Stc(i, k) * cea%H0(j) + delh(i)
          end do
          cea%G(i, m) = delh(i)
       end do
       do i = 1, maxTr
          do j = 1, maxTr
             if (abs(cea%Stc(i, j)) < 1.0d-6) cea%Stc(i, j) = 0
          end do
       end do
       jj = cea%Nm - 1
       do k = 1, jj
          mm = k + 1
          do m = mm, cea%Nm
             rtpd(k, m) = cea%Wmol(k) * cea%Wmol(m) / (1.1 * cea%Eta(k, m) * (cea%Wmol(k) + cea%Wmol(m)))
             xskm(k, m) = cea%Xs(k) * cea%Xs(m)
             xskm(m, k) = xskm(k, m)
             rtpd(m, k) = rtpd(k, m)
          end do
       end do
       do i = 1, cea%Nr
          do j = i, cea%Nr
             cea%G(i, j) = 0
             gmat(i, j) = 0
          end do
       end do
       do k = 1, jj
          mm = k + 1
          do m = mm, cea%Nm
             if (cea%Xs(k) >= 1.0d-10 .and. cea%Xs(m) >= 1.0d-10) then
                do j = 1, cea%Nr
                   if ((cea%Stc(j, k) == 0) .and. (cea%Stc(j, m) == 0)) stx(j) = 0
                   if ((cea%Stc(j, k) /= 0) .or. (cea%Stc(j, m) /= 0)) stx(j) = cea%Xs(m) * cea%Stc(j, k) - cea%Xs(k) * cea%Stc(j, m)
                end do
                do i = 1, cea%Nr
                   do j = i, cea%Nr
                      stxij(i, j) = stx(i) * stx(j) / xskm(k, m)
                      cea%G(i, j) = cea%G(i, j) + stxij(i, j)
                      gmat(i, j) = gmat(i, j) + stxij(i, j) * rtpd(k, m)
                   end do
                end do
             end if
          end do
       end do
       m = 1 + cea%Nr
       do i = 1, cea%Nr
          !DIR$ IVDEP
          do j = i, cea%Nr
             cea%G(j, i) = cea%G(i, j)
          end do
          cea%G(i, m) = delh(i)
       end do
       cea%Imat = cea%Nr

       call gauss_elimination(cea%G(1:cea%Imat, 1:cea%Imat+1), cea%X(1:cea%Imat))

       cpreac = 0
       do i = 1, cea%Nr
          cea%G(i, m) = delh(i)
          cpreac = cpreac + cea%R * delh(i) * cea%X(i)
          !DIR$ IVDEP
          do j = i, cea%Nr
             cea%G(i, j) = gmat(i, j)
             cea%G(j, i) = cea%G(i, j)
          end do
       end do

       call gauss_elimination(cea%G(1:cea%Imat, 1:cea%Imat+1), cea%X(1:cea%Imat))

       reacon = 0
       do i = 1, cea%Nr
          reacon = reacon + cea%R*delh(i)*cea%X(i)
       end do
       reacon = 0.6d0 * reacon
    else
       cpreac = 0
       reacon = 0
    end if
    ! CALCULATE OTHER ANSWERS
    cea%Cpfro(cea%Npt) = 0
    wtmol = 0
    do i = 1, cea%Nm
       cea%Cpfro(cea%Npt) = cea%Cpfro(cea%Npt) + cea%Xs(i) * cea%Cprr(i)
       wtmol = wtmol + cea%Xs(i) * cea%Wmol(i)
    end do
    cea%Cpfro(cea%Npt) = cea%Cpfro(cea%Npt) * cea%R / wtmol
    cea%Confro(cea%Npt) = cea%Confro(cea%Npt) / 1000
    if (.not. cea%SIunit) cea%Confro(cea%Npt) = cea%Confro(cea%Npt) / 4.184d0
    cea%Vis(cea%Npt) = cea%Vis(cea%Npt) / 1000
    cea%Prfro(cea%Npt) = cea%Vis(cea%Npt) * cea%Cpfro(cea%Npt) / cea%Confro(cea%Npt)
    if (cea%Eql) then
       cpreac = cpreac / wtmol
       reacon = reacon / 1000
       cea%Cpeql(cea%Npt) = cpreac + cea%Cpfro(cea%Npt)
       cea%Coneql(cea%Npt) = cea%Confro(cea%Npt) + reacon
       cea%Preql(cea%Npt) = cea%Vis(cea%Npt) * cea%Cpeql(cea%Npt) / cea%Coneql(cea%Npt)
    end if
  end subroutine TRANP


#ifndef NDEBUG
  subroutine print_debug(cea, iof)
    use mod_constants, only: stderr
    use mod_types, only: CEA_Problem
    implicit none

    type(CEA_Problem), intent(in):: cea
    integer, intent(in):: iof

    integer:: i

    write(stderr, *) "[DEBUG] Case: ", trim(cea%Case)
    write(stderr, *) "[DEBUG] iof = ", iof, ", Nof = ", cea%Nof
    write(stderr, *) "[DEBUG] ", (trim(cea%ensert(i)), i = 1, 20)
    write(stderr, *) "[DEBUG] ", cea%Enn, cea%Ennl, cea%Enlsav, cea%Ensave, cea%Sumn
!!$    write(stderr, *) "[DEBUG] ", cea%Deln(), cea%Enln(), cea%Sln()
!!$    write(stderr, *) "[DEBUG] ", cea%En()
    write(stderr, *) "[DEBUG] ", cea%Iplt, cea%Nc, cea%Ng, cea%Ngp1, cea%Nlm, cea%Nplt, cea%Nof
    write(stderr, *) "[DEBUG] ", cea%Nomit, cea%Nonly, cea%Np, cea%Npr, cea%Ngc, cea%Nsert, cea%Nspr, cea%Nspx, cea%Nt
    write(stderr, *) "[DEBUG] ", cea%Cpmix, cea%Wmix, cea%Bcheck
    write(stderr, *) "[DEBUG] ", cea%Imat, cea%Iq1, cea%Isv, cea%Jliq, cea%Jsol, cea%Lsave, cea%Msing
    write(stderr, *) "[DEBUG] ", cea%Convg, cea%Detdbg, cea%Detn, cea%Eql, cea%Gonly, cea%Hp, cea%Ions, cea%Massf, cea%Moles
    write(stderr, *) "[DEBUG] ", cea%Pderiv, cea%Shock, cea%Short, cea%SIunit, cea%Sp, cea%Tp, cea%Trnspt, cea%Vol
    write(stderr, *) "[DEBUG] ", cea%Eqrat, cea%Hsub0, cea%Oxfl, cea%Pp, cea%R, cea%Size, cea%S0
    write(stderr, *) "[DEBUG] ", cea%Tln, cea%Tm, cea%Trace, cea%Tt, cea%Viscns, cea%Vv
!!$    write(stderr, *) "[DEBUG] ", 

    return
  end subroutine print_debug
#endif

end module mod_cea
