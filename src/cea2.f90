!***********************************************************************
!                     P R O G R A M      C E A 2
!
!             CHEMICAL EQULIBRIUM WITH APPLICATIONS         5/21/04
!***********************************************************************
program main
  use mod_legacy_cea
  use mod_legacy_io
  implicit none

  integer:: icase = 0
  type(CEA_Problem):: cea(255)

  character(15):: ensert(20)
  character(200):: infile, ofile
  character(196):: prefix
  logical:: caseOK, ex, readOK
  integer:: i, inc, iof, j, ln, n
  real(8):: xi, xln

  write(*, '(//" ENTER INPUT FILE NAME WITHOUT .inp EXTENSION."/  &
       & "   THE OUTPUT FILES FOR LISTING AND PLOTTING WILL HAVE", / &
       & " THE SAME NAME WITH EXTENSIONS .out AND .plt RESPECTIVELY" &
       & //)')

  read(*, '(a)') prefix

  ln = index(prefix, ' ') - 1
  infile = prefix(1:ln) // '.inp'
  ofile  = prefix(1:ln) // '.out'
  Pfile  = prefix(1:ln) // '.plt'

  inquire(file=infile, exist=ex)
  if (.not. ex) then
     print *, infile, ' DOES NOT EXIST'
     error stop
  end if

  open(IOINP, file=infile, status='old', form='formatted')
  open(IOOUT, file=ofile, status='unknown', form='formatted')
  open(IOSCH, status='scratch', form='unformatted')
  open(IOTHM, file='thermo.lib', form='unformatted')
  open(IOTRN, file='trans.lib', form='unformatted')

  write(IOOUT, '(/" *******************************************************************************")')
  write(IOOUT, '(/, 9x, "NASA-GLENN CHEMICAL EQUILIBRIUM PROGRAM CEA2,", &
       & " MAY 21, 2004", /19x, "BY  BONNIE MCBRIDE", &
       & " AND SANFORD GORDON", /5x, &
       & " REFS: NASA RP-1311, PART I, 1994", &
       & " AND NASA RP-1311, PART II, 1996")')
  write(IOOUT, '(/" *******************************************************************************")')

  readOK = .true.
  Newr = .false.

  outerLoop: do while (readOK)
     icase = icase + 1

     Iplt = 0
     Nplt = 0

     call INPUT(readOK, caseOK, ensert)

     if ((.not. caseOK) .or. (.not. readOK)) cycle

     do iof = 1, Nof
        if (Oxf(iof) == 0. .and. B0p(1, 1) /= 0.) then
           do i = 1, Nlm
              if (B0p(i, 1) == 0. .or. B0p(i, 2) == 0.) then
                 write(IOOUT, '(/, "OXIDANT NOT PERMITTED WHEN SPECIFYING 100% FUEL(main)")')
                 cycle outerLoop
              end if
           end do
        end if
     end do

     if (Ions) then
        if (Elmt(Nlm) /= 'E') then
           Nlm = Nlm + 1
           Elmt(Nlm) = 'E'
           B0p(Nlm, 1) = 0
           B0p(Nlm, 2) = 0
        end if
     else if (Elmt(Nlm) == 'E') then
        Nlm = Nlm - 1
     end if

     Jray(1:Nreac) = 0

     call SEARCH

     if (Ngc == 0) exit outerLoop

     Newr = .false.

     if (Trnspt) call READTR

     ! INITIAL ESTIMATES
     Npr = 0
     Gonly = .true.
     Enn = 0.1d0
     Ennl = -2.3025851
     Sumn = Enn
     xi = Ng
     if (xi == 0.) xi = 1
     xi = Enn/xi
     xln = log(xi)

     En(Ng+1:Ng+Nc, 1) = 0
     Enln(Ng+1:Ng+Nc) = 0

     En(1:Ng, 1) = xi
     Enln(1:Ng) = xln

     if (Nc /= 0 .and. Nsert /= 0) then
        innerLoop: do i = 1, Nsert
           do j = Ngc, Ngp1, - 1
              if (Prod(j) == ensert(i)) then
                 Npr = Npr + 1
                 Jcond(Npr) = j
                 if (.not. Short) write(IOOUT, '(1X, A16, "INSERTED")') Prod(j)
                 cycle innerLoop
              end if
           end do
           write(IOOUT, '(/" WARNING!!!", A16, "NOT FOUND FOR INSERTION")') ensert(i)
        end do innerLoop
     end if

     if (Rkt) then
        call ROCKET(cea(icase))
     else if (Tp .or. Hp .or. Sp) then
        call THERMP(cea(icase))
     else if (Detn) then
        call DETON(cea(icase))
     else if (Shock) then
        call SHCK(cea(icase))
     end if

     if (Nplt > 0) then
        open(IOPLT, file = Pfile, form = 'formatted')
        call write_plt_file(IOPLT, Iplt, Nplt, Pltvar, Pltout)
     end if

  end do outerLoop

  close(IOINP)
  close(IOOUT)
  close(IOSCH)
  close(IOTHM)
  close(IOTRN)
  close(IOPLT)

  stop
end program main


subroutine CPHS
!***********************************************************************
! CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
!***********************************************************************
  use mod_legacy_cea
  implicit none
  integer:: i, ij, j, jj, k

  cx(:) = [0d0, 0d0, 1d0, 0.5d0, 0.6666666666666667d0, 0.75d0, 0.8d0]
  hcx(:) = [0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0]
  scx(:) = 0d0

  k = 1
  if (Tt > Tg(2)) k = 2
  if (Tt > Tg(3)) k = 3
  cx(2) = 1 / Tt
  cx(1) = cx(2)**2
  scx(3) = Tln
  scx(2) = -cx(2)
  hcx(2) = Tln * cx(2)
  hcx(1) = -cx(1)
  scx(1) = hcx(1) / 2

  forall(i = 4:7)
     hcx(i) = cx(i) * Tt
     scx(i) = cx(i-1) * Tt
  end forall

  H0(1:Ng) = 0
  S(1:Ng) = 0

  do i = 7, 4, -1
     forall(j = 1:Ng)
        S(j) = (S(j) + Coef(j, i, k)) * scx(i)
        H0(j) = (H0(j) + Coef(j, i, k)) * hcx(i)
     end forall
  end do

  do i = 1, 3
     forall(j = 1:Ng)
        S(j) = S(j) + Coef(j, i, k) * scx(i)
        H0(j) = H0(j) + Coef(j, i, k) * hcx(i)
     end forall
  end do

  forall(j = 1:Ng)
     S(j) = S(j) + Coef(j, 9, k)
     H0(j) = H0(j) + Coef(j, 8, k) * cx(2)
  end forall

  if (.not. Tp .or. Convg) then
     Cp(1:Ng) = 0

     do i = 7, 4, -1
        forall(j = 1:Ng) Cp(j) = (Cp(j) + Coef(j, i, k)) * Tt
     end do

     forall(j = 1:Ng) Cp(j) = Cp(j) + sum(Coef(j, 1:3, k) * cx(1:3))
  end if

  if (Npr /= 0 .and. k /= 3 .and. Ng /= Ngc) then
     do ij = 1, Npr
        j = Jcond(ij)
        jj = Jcond(ij) - Ng
        Cp(j) = 0
        H0(j) = 0
        S(j) = 0

        do i = 7, 4, -1
           S(j) = (S(j) + Cft(jj, i)) * scx(i)
           H0(j) = (H0(j) + Cft(jj, i)) * hcx(i)
           Cp(j) = (Cp(j) + Cft(jj, i)) * Tt
        end do

        do i = 1, 3
           S(j) = S(j) + Cft(jj, i) * scx(i)
           H0(j) = H0(j) + Cft(jj, i) * hcx(i)
           Cp(j) = Cp(j) + Cft(jj, i) * cx(i)
        end do

        S(j) = S(j) + Cft(jj, 9)
        H0(j) = H0(j) + Cft(jj, 8) * cx(2)
     end do
  end if

  return
end subroutine CPHS


subroutine ALLCON
!***********************************************************************
! CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
!***********************************************************************
  use mod_legacy_cea
  implicit none
  integer:: i, j, jj

  do jj = 1, Nc
     j = jj + Ng
     Cp(j) = 0
     H0(j) = 0
     S(j) = 0

     do i = 7, 4, -1
        S(j) = (S(j) + Cft(jj, i)) * scx(i)
        H0(j) = (H0(j) + Cft(jj, i)) * hcx(i)
        Cp(j) = (Cp(j) + Cft(jj, i)) * Tt
     end do

     do i = 1, 3
        S(j) = S(j) + Cft(jj, i) * scx(i)
        H0(j) = H0(j) + Cft(jj, i) * hcx(i)
        Cp(j) = Cp(j) + Cft(jj, i) * cx(i)
     end do

     S(j) = S(j) + Cft(jj, 9)
     H0(j) = H0(j) + Cft(jj, 8) * cx(2)
  end do

  return
end subroutine ALLCON



subroutine DETON(cea)
!***********************************************************************
! CHAPMAN-JOUGUET DETONATIONS.
!***********************************************************************
  use mod_legacy_cea
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
  real(8):: a11, a12, a21, a22, alam, alpha, amm, b1, b2, cpl(Ncol), d, gam, &
       gm1(Ncol), h1(Ncol), p1, pp1, pub(Ncol), rk, rr1, rrho(Ncol), T1, &
       tem, tt1, tub(Ncol), ud, x1, x2


  iof = 0
  Eql = .true.

  if (T(1) == 0) then
     T(1) = Rtemp(1)
     Nt = 1
  end if

100 Tt = T(1)

  iof = iof + 1
  Oxfl = Oxf(iof)

  call NEWOF

! BEGIN T LOOP.
  do It = 1, Nt
     T1 = T(It)

! BEGIN P LOOP.
     do Ip = 1, Np
        p1 = P(Ip)
        Tt = T1
        Pp = p1

        call HCALC

        if (Tt == 0) return
        if (Detdbg) call OUT1

        h1(Npt) = Hsub0 * R
        tub(Npt) = T1
        pub(Npt) = p1
        cpl(Npt) = Cpmix * R
        itr = 0
        Tt = 3800
        pp1 = 15
        Pp = pp1 * p1

! CALCULATE ENTHALPY FOR INITIAL ESTIMATE OF T2(TT AFTER EQLBRM)
        Hsub0 = h1(Npt) / R + 0.75 * T1 * pp1 / Wmix
        Tp = .false.
        Hp = .true.

        call EQLBRM

        Hsub0 = h1(Npt) / R
        Hp = .false.

        if (Tt /= 0) then
           gam = Gammas(Npt)
           tt1 = Tt / T1
           ii = 0
           tem = tt1 - 0.75 * pp1 / (Cpr(Npt) * Wmix)
           amm = Wm(Npt) / Wmix

           if (Detdbg) write(IOOUT, '(/" T EST.=", F8.2/11X, "P/P1", 17X, "T/T1")') Tt

! LOOP FOR IMPROVING T2/T1 AND P2/P1 INITIAL ESTIMATE.
           do ii = 1, 3
              alpha = amm / tt1
              pp1 = (1 + gam) * (1 + sqrt(1 - 4 * gam * alpha / (1 + gam)**2)) / (2 * gam * alpha)
              rk = pp1 * alpha
              tt1 = tem + 0.5 * pp1 * gam * (rk**2 - 1) / (Wmix * Cpr(Npt) * rk)
              if (Detdbg) write(IOOUT, '(i5, 2e20.8)') ii, pp1, tt1
           end do

           Tp = .true.
           Tt = T1 * tt1
           rr1 = pp1 * amm / tt1

! BEGIN MAIN ITERATION LOOP.
110        itr = itr + 1
           Pp = p1 * pp1

           call EQLBRM

           if (Npt == 0) go to 200

           if (Tt /= 0) then
              gam = Gammas(Npt)
              amm = Wm(Npt) / Wmix
              rr1 = pp1 * amm / tt1
              a11 = 1 / pp1 + gam * rr1 * Dlvpt(Npt)
              a12 = gam * rr1 * Dlvtp(Npt)
              a21 = 0.5 * gam * (rr1**2 - 1 - Dlvpt(Npt) * (1 + rr1**2)) + Dlvtp(Npt) - 1
              a22 = -0.5 * gam * Dlvtp(Npt) * (rr1**2 + 1) - Wm(Npt) * Cpr(Npt)
              b1 = 1 / pp1 - 1 + gam * (rr1 - 1)
              b2 = Wm(Npt) * (Hsum(Npt) - h1(Npt) / R) / Tt - 0.5 * gam * (rr1**2 - 1)
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

              Tt = T1 * tt1
              ud = rr1 * sqrt(Rr * gam * Tt/Wm(Npt))

              if (Detdbg) write(IOOUT, '(/" ITER =", i2, 5x, "P/P1 =", e15.8, /7x, "T/T1 =", e15.8, 5x, &
                   & "RHO/RHO1 =", e15.8, /7x, "DEL LN P/P1 =", e15.8, 5x, &
                   & "DEL LN T/T1 =", e15.8)') itr, pp1, tt1, rr1, x1, x2
! CONVERGENCE TEST
              if (itr < 8 .and. tem > 0.5E-04) go to 110

              if (itr < 8) then
                 rrho(Npt) = rr1
                 if (cpl(Npt) == 0) then
                    gm1(Npt) = 0
                    Vmoc(Npt) = 0
                 else
                    gm1(Npt) = cpl(Npt) / (cpl(Npt) - R / Wmix)
                    Vmoc(Npt) = ud / sqrt(Rr * gm1(Npt) * T1 / Wmix)
                 end if
              else
                 write(IOOUT, '(/" CONSERVATION EQNS NOT SATISFIED IN 8 ITERATIONS (DETON)")')
                 Npt = Npt - 1
                 Tt = 0
              end if

              if (Trnspt) call TRANP

              Isv = 0

              if (Ip /= Np .or. It /= Nt .and. Tt /= 0) then
                 Isv = Npt
                 if (Npt /= Ncol) go to 120
              end if
           end if

! OUTPUT
           write(IOOUT, '(//, 21X, "DETONATION PROPERTIES OF AN IDEAL REACTING GAS")')
           call OUT1

! SET MXX ARRAY FOR PLOTTING PARAMETERS
           mp    = 0
           mt    = 0
           mgam  = 0
           mh    = 0
           mdv   = 0
           mson  = 0
           mmach = 0

           do i = 1, Nplt
              if (index(Pltvar(i)(2:), '1') /= 0) then
                 if (Pltvar(i)(:3) == 'son') then
                    mson = i
                 else if (Pltvar(i)(:3) == 'gam') then
                    mgam = i
                 else if (Pltvar(i)(:1) == 'h') then
                    mh = i
                 else if (Pltvar(i)(:1) == 't') then
                    mt = i
                 else if (Pltvar(i)(:1) == 'p') then
                    mp = i
                 end if
              else if (index(Pltvar(i), 'vel') /= 0) then
                 mdv = i
              else if (index(Pltvar(i), 'mach') /= 0) then
                 mmach = i
              end if
           end do

           write(IOOUT, '(/" UNBURNED GAS"/)')

           cea%fmt(4) = '13'
           cea%fmt(5) = ' '
           cea%fmt(7) = '4,'

           do i = 1, Npt
              if (SIunit) then
                 V(i) = pub(i)
                 unit = 'BAR'
              else
                 V(i) = pub(i) / 1.01325d0
                 unit = 'ATM'
              end if
              if (mp > 0) Pltout(i+Iplt, mp) = V(i)
           end do

           write(IOOUT, cea%fmt) 'P1, ' // unit // '        ', (V(j), j = 1, Npt)

           cea%fmt(7) = '2,'
           write(IOOUT, cea%fmt) ft1, (tub(j), j = 1, Npt)

           if (.not. SIunit) write(IOOUT, cea%fmt) fh1, (h1(j), j = 1, Npt)
           if (SIunit) write(IOOUT, cea%fmt) fhs1, (h1(j), j = 1, Npt)

           forall(i = 1:Npt)
              V(i) = Wmix
              Sonvel(i) = sqrt(Rr * gm1(i) * tub(i) / Wmix)
           end forall

           cea%fmt(7) = '3,'
           write(IOOUT, cea%fmt) fm1, (V(j), j = 1, Npt)
           cea%fmt(7) = '4,'
           write(IOOUT, cea%fmt) fg1, (gm1(j), j = 1, Npt)
           cea%fmt(7) = '1,'
           write(IOOUT, cea%fmt) 'SON VEL1,M/SEC ', (Sonvel(j), j = 1, Npt)

           if (Nplt > 0) then
              do i = 1, Npt
                 if (mt > 0)   Pltout(i+Iplt, mt) = tub(i)
                 if (mgam > 0) Pltout(i+Iplt, mgam) = gm1(i)
                 if (mh > 0)   Pltout(i+Iplt, mh) = h1(i)
                 if (mson > 0) Pltout(i+Iplt, mson) = Sonvel(i)
              end do
           end if

           write(IOOUT, '(/" BURNED GAS"/)')

           cea%fmt(4) = cea%fmt(6)
           call OUT2(cea)

           if (Trnspt) call OUT4(cea)

           write(IOOUT, '(/" DETONATION PARAMETERS"/)')

           cea%fmt(7) = '3,'

           do i = 1, Npt
              V(i) = Ppp(i) / pub(i)
              Pcp(i) = Ttt(i) / tub(i)
              Sonvel(i) = Sonvel(i) * rrho(i)
              if (mmach > 0) Pltout(i+Iplt, mmach) = Vmoc(i)
              if (mdv > 0)   Pltout(i+Iplt, mdv) = Sonvel(i)
           end do

           write(IOOUT, cea%fmt) fpp1, (V(j), j = 1, Npt)
           write(IOOUT, cea%fmt) ftt1, (Pcp(j), j = 1, Npt)

           forall(i = 1:Npt) V(i) = Wm(i) / Wmix

           cea%fmt(7) = '4,'
           write(IOOUT, cea%fmt) fmm1, (V(j), j = 1, Npt)
           write(IOOUT, cea%fmt) frr1, (rrho(j), j = 1, Npt)
           write(IOOUT, cea%fmt) 'DET MACH NUMBER', (Vmoc(j), j = 1, Npt)

           cea%fmt(7) = '1,'
           write(IOOUT, cea%fmt) fdv, (Sonvel(j), j = 1, Npt)

           Eql = .true.

           call OUT3(cea)

           Iplt = min(Iplt+Npt, 500)

           if (Isv == 0 .and. iof == Nof) go to 200
           if (Np == 1 .and. Nt == 1) go to 100

           write(IOOUT, '(///)')

           Npt = 0
120        Npt = Npt + 1

           if (Isv == 1) Isv = -1

           call SETEN
        end if
     end do
  end do

  Iplt = min(Iplt + Npt - 1, 500)

  if (iof < Nof) go to 100

200 Tp = .false.

  return
end subroutine DETON



subroutine EFMT(Fone, Aa, Vx)
!***********************************************************************
! WRITE OUTPUT RECORD WITH NUMERICAL VALUES IN SPECIAL EXPONENT FORM.
!***********************************************************************
  use mod_legacy_cea
  implicit none
! DUMMY ARGUMENTS
  character(15):: Aa
  character(4):: Fone
  real(8):: Vx(maxMat)
! LOCAL VARIABLES
  character(4), parameter:: fmix(5) = [character(4):: 'I3,', '6.4,', 'I2,', '9X,', '5.3,']
  character(4):: frmt(8) = [character(4):: '(1H ', ',A15', ',', '9X,', '13(F', '6.4,', 'I2,', '1X))']
  integer:: i, j, j1, ne(Ncol)
  real(8):: ee, fe, w(Ncol)


  frmt(6) = fmix(2)
  frmt(7) = fmix(3)
  j1 = 1
  frmt(4) = '1x,'

  if (Fone == '9X,') then
     j1 = 2
     frmt(4) = fmix(4)
  end if

  do i = j1, Npt
     if (Vx(i) /= 0.) then
        ee = log10(abs(Vx(i)))
        ne(i) = int(ee)
        fe = ne(i)

        if (ee < -0.2181E-05 .and. fe /= ee) ne(i) = ne(i) - 1

        if (abs(ne(i)) >= 10) then
           frmt(6) = fmix(5)
           frmt(7) = fmix(1)
        end if
        w(i) = Vx(i) / 10.**ne(i)
     else
        w(i) = 0
        ne(i) = 0
     end if
  end do

  write(IOOUT, frmt) Aa, (w(j), ne(j), j = j1, Npt)

  return
end subroutine EFMT



subroutine EQLBRM
!***********************************************************************
! CALCULATE EQUILIBRIUM COMPOSITION AND PROPERTIES.
!***********************************************************************
  use mod_legacy_cea
  implicit none
! LOCAL VARIABLES
  character(12):: ae, cmp(maxEl)
  character(16):: amb
  logical, save:: cpcalc, i2many, newcom, reduce
  integer:: i, il, ilamb, ilamb1, inc, ipr, iq2, iter, ix, ixsing, iz, j, ja, jb, &
       jbx, jc, jcondi, jcons, jdelg, jex, jj, jkg, jneg, jsw, k, kc, kg, kk, &
       kmat, kneg, l, lc, lcs(maxEl), le, lelim, lk, ll, lncvg, ls, lsing, &
       lz, maxitn, ncvg, njc, nn, numb
  real(8), save:: aa, ambda, ambda1, bigen, bigneg, delg, dlnt, dpie, ensol, esize, &
       gap, gasfrc, pie, pisave(maxMat-2), siz9, sizeg, &
       sum0, sum1, szgj, tem, tmelt, tsize, ween, xi, xln, xsize, xx(maxMat)
  real(8):: smalno = 1e-6, smnol = -13.815511
  logical:: mask(Ng)

  ixsing = 0
  lsing = 0
  jsw = 0
  jdelg = 0
  maxitn = 50
  ncvg = 0
  lncvg = 3 * Nlm
  reduce = .false.
  siz9 = Size - 9.2103404d0
  tsize = Size
  xsize = Size + 6.90775528d0

  if (Trace /= 0.) then
     maxitn = maxitn + Ngc / 2
     xsize = -log(Trace)
     if (xsize < Size) xsize = Size + 0.1
  end if

  if (xsize > 80) xsize = 80

  esize = min(80., xsize + 6.90775528d0)
  jcons = 0
  pie = 0
  i2many = .false.
  Pderiv = .false.
  Convg = .false.
  numb = 0
  cpcalc = .true.

  if (Tp) cpcalc = .false.

  if (Tt /= 0) then
     if (Npr == 0 .or. (Tt /= T(1) .and. .not. Tp)) go to 400
     k = 1
  else
     Tt = 3800
     if (Npr == 0) go to 400
     k = 1
  end if

100 j = Jcond(k)
  jc = j - Ng
  kg = -Ifz(jc)

  do i = 1, 9
     kg = kg + 1
     kc = jc + kg

     if (Tt <= Temp(2, kc)) then
        if (kg /= 0) then
           Jcond(k) = j + kg
           En(j+kg, Npt) = En(j, Npt)
           En(j, Npt) = 0

           if (Prod(j) /= Prod(j+kg) .and. .not. Short) &
                & write(IOOUT, '(" PHASE CHANGE, REPLACE ", A16, "WITH ", A16)') Prod(j), Prod(j+kg)
        end if
        go to 300
     else if (kc >= Nc .or. Ifz(kc+1) <= Ifz(kc)) then
        exit
     end if
  end do

  if (.not. Tp) then
     Tt = Temp(2, kc) - 10
     k = 1
     go to 100
  end if

  write(IOOUT, '(" REMOVE ", A16)') Prod(j)

  En(j, Npt) = 0
  Enln(j) = 0
  Deln(j) = 0

  do i = k, Npr
     Jcond(i) = Jcond(i+1)
  end do

  Npr = Npr - 1

300 k = k + 1

  if (k <= Npr) go to 100

400 Tln = log(Tt)

  if (Vol) Pp = Rr * Enn * Tt / Vv

  call CPHS

  Tm = log(Pp / Enn)
  le = Nlm

  if (Lsave /= 0 .and. Nlm /= Lsave) then
     tem = exp(-tsize)
     do i = Lsave + 1, Nlm
        forall(j = 1:Ng, A(i, j) /= 0)
           En(j, Npt) = tem
           Enln(j) = -tsize
        end forall
     end do
  end if

  ls = Nlm
  lelim = 0
  lz = ls

  if (Ions) lz = ls - 1

  if (Npt == 1 .and. .not. Shock .and. .not. Short) then
     write(IOOUT, '(/" POINT ITN", 6X, "T", 10X, 4(A4, 8X)/(18X, 5(A4, 8X)))') (Elmt(i), i = 1, Nlm)
  end if

  if (Debug(Npt)) then
     do i = 1, Nlm
        cmp(i) = Elmt(i)
     end do
  end if

! BEGIN ITERATION
500 if (cpcalc) then
     Cpsum = sum(En(1:Ng, Npt) * Cp(1:Ng))

     if (Npr /= 0) then
        Cpsum = Cpsum + sum(En(Jcond(1:Npr), Npt) * Cp(Jcond(1:Npr)))
        cpcalc = .false.
     end if
  end if

  numb = numb + 1

  call MATRIX

  iq2 = Iq1 + 1

  if (Convg) Imat = Imat - 1

  if (Debug(Npt)) then
     if (.not. Convg) then
        write(IOOUT, '(/" ITERATION", I3, 6X, "MATRIX ")') numb
     else
        if (.not. Pderiv) write(IOOUT, '(/" T DERIV MATRIX")')
        if (Pderiv) write(IOOUT, '(/" P DERIV MATRIX")')
     end if

     kmat = Imat + 1

     do i = 1, Imat
        write(IOOUT, '(3X, 5E15.6)') (G(i, k), k = 1, kmat)
     end do
  end if

  Msing = 0

  call GAUSS

  if (Msing == 0) then
     if (Debug(Npt)) then
        write(IOOUT, '(/" SOLUTION VECTOR", /, 6x, 5A15/8X, 5A15)') (cmp(k), k = 1, le)
        write(IOOUT, '(3X, 5E15.6)') (X(i), i = 1, Imat)
     end if

     if (.not. Convg) then
! OBTAIN CORRECTIONS TO THE ESTIMATES
        if (Vol) X(iq2) = X(Iq1)
        if (Tp) X(iq2) = 0
        dlnt = X(iq2)
        sum0 = X(Iq1)

        if (Vol) then
           X(Iq1) = 0
           sum0 = -dlnt
        end if

        outerLoop0: do j = 1, Ng
           if (lelim /= 0) then
              Deln(j) = 0
              do i = lelim, ls
                 if (A(i, j) /= 0) cycle outerLoop0
              end do
           end if

           Deln(j) = -Mu(j) + H0(j) * dlnt + sum0


!!$           Deln(j) = Deln(j) + sum(A(1:Nlm, j) * X(1:Nlm)) ???
           do k = 1, Nlm
              Deln(j) = Deln(j) + A(k, j) * X(k)
           end do

           if (pie /= 0) Deln(j) = Deln(j) + A(ls, j) * pie
        end do outerLoop0

        if (Npr /= 0) then
           forall(k = 1:Npr) Deln(Jcond(k)) = X(Nlm+k)
        end if

! CALCULATE CONTROL FACTOR, AMBDA
        ambda = 1
        ambda1 = 1
        ilamb = 0
        ilamb1 = 0
        sum0 = 5 * max(abs(X(Iq1)), abs(dlnt))

        do j = 1, Ng
           if (Deln(j) > 0) then
              if ((Enln(j) - Ennl + Size) <= 0) then

                 sum1 = abs(Deln(j)-X(Iq1))

                 if (sum1 >= siz9) then
                    sum1 = abs(-9.2103404d0 - Enln(j) + Ennl) / sum1

                    if (sum1 < ambda1) then
                       ambda1 = sum1
                       ilamb1 = j
                    end if
                 end if

              else if (Deln(j) > sum0) then
                 sum0 = Deln(j)
                 ilamb = j
              end if
           end if
        end do

        if (sum0 > 2) ambda = 2 / sum0

        if (ambda1 <= ambda) then
           ambda = ambda1
           ilamb = ilamb1
        end if

        if (Debug(Npt)) then
! INTERMEDIATE OUTPUT
           write(IOOUT, '(/" T=", E15.8, " ENN=", E15.8, " ENNL=", E15.8, " PP=", E15.8, &
                & /" LN P/N=", E15.8, " AMBDA=", E15.8)') Tt, Enn, Ennl, Pp, Tm, ambda

           if (ambda /= 1) then
              amb = 'ENN'
              if (abs(X(iq2)) > abs(X(Iq1))) amb = 'TEMP'
              if (ilamb /= 0) amb = Prod(ilamb)
              write(IOOUT, '(/" AMBDA SET BY ", A16)') amb
           end if

           if (Vol) write(IOOUT, '(" VOLUME=", E15.8, "CC/G")') Vv * 0.001d0

           write(IOOUT, '(/24X, "Nj", 12X, "LN Nj", 8X, "DEL LN Nj", 6X, "H0j/RT", /, 41X, &
                & "S0j/R", 10X, " G0j/RT", 8X, " Gj/RT")')

           do j = 1, Ngc
              write(IOOUT, '(1X, A16, 4E15.6, /35x, 3E15.6)') &
                   Prod(j), En(j, Npt), Enln(j), Deln(j), H0(j), S(j), H0(j) - S(j), Mu(j)
           end do
        end if

! APPLY CORRECTIONS TO ESTIMATES
        Totn(Npt) = 0

        forall(j = 1:Ng) Enln(j) = Enln(j) + ambda * Deln(j)

        En(1:Ng, Npt) = 0

        forall(j = 1:Ng, (lelim == 0 .or. all(A(lelim:ls, j) == 0)) .and. (Enln(j) - Ennl + tsize) > 0)
           En(j, Npt) = exp(Enln(j))
        end forall

        Totn(Npt) = Totn(Npt) + sum(En(1:Ng, Npt))

        if (Ions .and. Elmt(Nlm) == 'E') then
           mask = .false.
           forall(j = 1:Ng, A(ls, j) /= 0 .and. En(j, Npt) == 0 .and. (Enln(j) - Ennl + esize) > 0)
              En(j, Npt) = exp(Enln(j))
              mask(j) = .true.
           end forall
           Totn(Npt) = Totn(Npt) + sum(En(1:Ng, Npt), mask=mask)
        end if

        Sumn = Totn(Npt)

        if (Npr /= 0) then
           forall(k = 1:Npr) En(Jcond(k), Npt) = En(Jcond(k), Npt) + ambda * Deln(Jcond(k))
           Totn(Npt) = Totn(Npt) + sum(En(Jcond(1:Npr), Npt))
        end if

        if (.not. Tp) then
           Tln = Tln + ambda * dlnt
           Tt = exp(Tln)
           cpcalc = .true.
           call CPHS
        end if

        if (Vol) then
           Enn = Sumn
           Ennl = log(Enn)
           if (Vol) Pp = Rr * Tt * Enn / Vv
        else
           Ennl = Ennl + ambda * X(Iq1)
           Enn = exp(Ennl)
        end if

        Tm = log(Pp / Enn)

        if (Elmt(Nlm) == 'E') then
! CHECK ON REMOVING IONS
           if (all(A(Nlm, 1:Ngc) == 0 .or. En(1:Ngc, Npt) <= 0)) then
              pie = X(Nlm)
              lelim = Nlm
              Nlm = Nlm - 1
              go to 500
           end if
        end if

! TEST FOR CONVERGENCE
        if (numb > maxitn) then
           write(IOOUT, '(/, I4, " ITERATIONS DID NOT SATISFY CONVERGENCE", /, 15x, &
                & " REQUIREMENTS FOR THE POINT", I5, " (EQLBRM)")') maxitn, Npt

           if (Nc == 0 .or. i2many) go to 1500

           i2many = .true.

           if (.not. Hp .or. Npt /= 1 .or. Tt > 100.) then
              if (Npr /= 1 .or. Enn > 1.E-4) go to 1500
! HIGH TEMPERATURE, INCLUDED CONDENSED CONDITION
              write(IOOUT, '(/" TRY REMOVING CONDENSED SPECIES (EQLBRM)")')

              Enn = 0.1
              Ennl = -2.3025851
              Sumn = Enn
              xi = Ng
              xi = Enn/xi
              xln = log(xi)

              forall(j = 1:Ng)
                 En(j, Npt) = xi
                 Enln(j) = xln
              end forall

              j = Jcond(1)
              k = 1
              go to 1000

           else
              write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", &
                   & "VE BEEN INSERTED,", &
                   & /" RESTART WITH insert DATASET (EQLBRM)")')
              go to 1500
           end if

        else
           if (abs(X(Iq1) * Enn / Totn(Npt)) > 0.5E-5 .or. &
                any(abs(Deln(1:Ng)) * En(1:Ng, Npt) / Totn(Npt) > 0.5d-5) .or. &
                abs(dlnt) > 1.d-04) go to 500

           if (Npr /= 0) then
              do k = 1, Npr
                 j = Jcond(k)
                 if (abs(Deln(j)/Totn(Npt)) > 0.5d-5) go to 500
                 if (En(j, Npt) < 0) go to 700
              end do
           end if

           le = Nlm

           if (any(abs(B0(1:Nlm)) >= 1.d-06 .and. &
                abs(B0(1:Nlm) - sum(spread(En(1:Ngc, Npt), 1, Nlm) * A(1:Nlm, 1:Ngc), DIM=2)) > Bcheck)) go to 500

           if (Trace /= 0.) then
              tsize = xsize
              tem = 1

              if (numb /= 1) then
                 lk = lz
                 if (Nlm < lz) lk = Nlm
                 do i = 1, lk
                    if (i /= lsing) then
                       tem = 0
                       if (X(i) /= 0) then
                          tem = abs((pisave(i) - X(i)) / X(i))
                          if (tem > 0.001) exit
                       end if
                    end if
                 end do
              end if

              forall(i = 1:Nlm) pisave(i) = X(i)

              if (tem > 0.001) go to 500

              if (Ions) then
! CHECK ON ELECTRON BALANCE
                 iter = 1

                 if (pie /= 0) then
                    X(Nlm+1) = pie
                 end if

566              sum1 = 0
                 sum0 = 0
                 pie = X(le)

                 do j = 1, Ng
                    if (A(ls, j) /= 0) then
                       En(j, Npt) = 0
                       tem = 0

                       if (Enln(j) > -87) tem = exp(Enln(j))

                       if ((Enln(j)-Ennl+tsize) > 0 .and. Elmt(Nlm) == 'E') then
                          pie = 0
                          En(j, Npt) = tem
                       end if

                       aa = A(ls, j) * tem
                       sum0 = sum0 + aa
                       sum1 = sum1 + aa * A(ls, j)
                    end if
                 end do

                 if (sum1 /= 0) then
                    dpie = -sum0 / sum1

                    forall(j = 1:Ng, A(ls, j) /= 0) Enln(j) = Enln(j) + A(ls, j) * dpie

                    if (Debug(Npt)) write(IOOUT, '(/" ELECTRON BALANCE ITER NO. =", i4, "  DELTA PI =", e14.7)') iter, dpie

                    if (abs(dpie) > 0.0001) then
                       X(le) = X(le) + dpie
                       iter = iter + 1

                       if (iter <= 80) go to 566
                       write(IOOUT, '(/" DID NOT CONVERGE ON ELECTRON BALANCE (EQLBRM)")')
                       go to 1500

                    else if (Elmt(Nlm) == 'E' .and. pie /= 0) then
                       Nlm = Nlm - 1
                       newcom = .true.
                    end if
                 end if
              end if
           end if
        end if

     else if (.not. Pderiv) then
! TEMPERATURE DERIVATIVES--CONVG=T, PDERIV=F
        Dlvtp(Npt) = 1. - X(Iq1)
        Cpr(Npt) = G(iq2, iq2)

        Cpr(Npt) = Cpr(Npt) - sum(G(iq2, 1:Iq1) * X(1:Iq1))

! PRESSURE DERIVATIVE--CONVG=T, PDERIV=T
        Pderiv = .true.
        go to 500

     else
        Dlvpt(Npt) = -1 + X(Iq1)
        if (Jliq == 0) then
           Gammas(Npt) = -1 / (Dlvpt(Npt) + (Dlvtp(Npt)**2) * Enn / Cpr(Npt))
        else
           En(Jsol, Npt) = ensol
           Hsum(Npt) = Hsum(Npt) + En(Jliq, Npt) * (H0(Jliq) - H0(Jsol))
           Gammas(Npt) = -1. / Dlvpt(Npt)
           Npr = Npr + 1
           Jcond(Npr) = Jliq
        end if
        go to 1400
     end if

! SINGULAR MATRIX
  else
     if (Convg) then
        write(IOOUT, '(/" DERIVATIVE MATRIX SINGULAR (EQLBRM)")')
        Dlvpt(Npt) = -1
        Dlvtp(Npt) = 1
        Cpr(Npt) = Cpsum
        Gammas(Npt) = -1 / (Dlvpt(Npt) + Dlvtp(Npt)**2 * Enn / Cpr(Npt))
        go to 1400

     else
        write(IOOUT, '(/" SINGULAR MATRIX, ITERATION", I3, "  VARIABLE", I3, "(EQLBRM)")') numb, Msing
        lsing = Msing
        ixsing = ixsing + 1
        if (ixsing <= 8) then
           xsize = 80
           tsize = xsize
           if (Msing > Nlm .and. numb < 1 .and. Npr > 1 .and. jdelg > 0) then
              ween = 1000
              j = 0

              do i = 1, Npr
                 jcondi = Jcond(i)
                 if (jcondi /= jdelg) then
                    do ll = 1, Nlm
                       if (A(ll, jdelg) /= 0 .and. A(ll, jcondi) /= 0) then
                          if (En(jcondi, Npt) <= ween) then
                             ween = En(jcondi, Npt)
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

           else if (.not. Hp .or. Npt /= 1 .or. Nc == 0 .or. Tt > 100) then
              if (ixsing >= 3) then
                 if (Msing < Iq1) then
                    if (reduce .and. Msing <= Nlm) then
                       if (Nlm < lelim) go to 1500
                       write(IOOUT, '(/" WARNING!! POINT", I3, &
                            & " USES A REDUCED SET OF COMPONENTS", / &
                            & " SPECIES CONTAINING THE ELIMINATED COMPONENT ARE OMITTED.", &
                            & / &
                            & " IT MAY BE NECESSARY TO RERUN WITH INSERTED CONDENSED SPECIES", &
                            & /" CONTAINING COMPONENT ", A8, "(EQLBRM)")') Npt, Elmt(Nlm)
                       Nlm = Nlm - 1
                       go to 500

                    else if (Msing <= Nlm) then
! FIND NEW COMPONENTS
                       if (.not. Ions) go to 1100
                       if (Elmt(Nlm) /= 'E') go to 1100

                       forall(j = 1:Ng, A(Nlm, j) /= 0) En(j, Npt) = 0

                       pie = X(Nlm)
                       Nlm = Nlm - 1
                       if (Msing > Nlm) go to 500
                       go to 1100
                    else
! REMOVE CONDENSED SPECIES TO CORRECT SINGULARITY
                       k = Msing - Nlm
                       j = Jcond(k)

                       if (j /= jcons) then
                          jcons = j
                          go to 1000
                       end if
                    end if
                 end if
              end if

              forall(j = 1:Ng, .not. (Ions .and. Elmt(Nlm) /= 'E' .and. A(ls, j) /= 0) .and. En(j, Npt) == 0)
                 En(j, Npt) = smalno
                 Enln(j) = smnol
              end forall

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
600 Ssum(Npt) = sum(En(1:Ng, Npt) * (S(1:Ng) - Enln(1:Ng) - Tm))

  if (Npr > 0) then
     Ssum(Npt) = Ssum(Npt) + sum(En(Jcond(1:Npr), Npt) * S(Jcond(1:Npr)))
  end if

  if (.not. Sp) then
     Convg = .true.
  else
     tem = Ssum(Npt) - S0
     if (abs(tem) > 0.0005) go to 500
     if (Debug(Npt)) write(IOOUT, '(/" DELTA S/R =", e15.8)') tem
     Convg = .true.
  end if

! CONVERGENCE TESTS ARE SATISFIED, TEST CONDENSED SPECIES.
700 ncvg = ncvg + 1

  if (ncvg > lncvg) then
! ERROR, SET TT=0
     write(IOOUT, '(/, I3, " CONVERGENCES FAILED TO ESTABLISH SET OF CONDENSED", " SPECIES (EQLBRM)")') lncvg
     go to 1500
  else
     if (.not. Shock) then
        forall(il = 1:le) xx(il) = X(il)

        if (.not. Short) then
           if (newcom) write(IOOUT, '(/" POINT ITN", 6x, "T", 10x, 4a12/(18x, 5a12))') (cmp(k), k = 1, le)
           write(IOOUT, '(i4, i5, 5f12.3, /(12x, 5f12.3))') Npt, numb, Tt, (xx(il), il = 1, le)
        end if

        if (.not. Tp .and. Npr == 0 .and. Tt <= Tg(1) * 0.2d0) then
           write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", "VE BEEN INSERTED,", &
                & /" RESTART WITH insert DATASET (EQLBRM)")')
           go to 1500
        end if

        newcom = .false.
     end if

     if (Npr /= 0) then
        bigneg = 0
        jneg = 0

        do k = 1, Npr
           j = Jcond(k)
           if (En(j, Npt) * Cp(j) <= bigneg) then
              bigneg = En(j, Npt) * Cp(j)
              jneg = j
              kneg = k
           end if
        end do

        if (jneg /= 0) then
           j = jneg
           k = kneg
           if (j == Jsol .or. j == Jliq) then
              Jsol = 0
              Jliq = 0
           end if
           go to 1000
        end if
     end if

     if (Ngc /= Ng .or. Tp) then
        Ng = Ngc

        call CPHS

        Ng = Ngp1 - 1
        cpcalc = .true.

        if (Ngc == Ng) go to 750

        call ALLCON

        if (Npr /= 0 .and. .not. Tp) then
           gap = 50

           outerLoop1: do ipr = 1, Npr
              j = Jcond(ipr)
              if (j /= Jsol .and. j /= Jliq) then
                 inc = j - Ng
                 kg = -Ifz(inc)

                 do iz = 1, 20
                    kg = kg + 1
                    kc = inc + kg

                    if (Tt <= Temp(2, kc)) then
                       if (kg /= 0) then
                          jkg = j + kg
                          if (abs(kg) > 1 .or. Prod(j) == Prod(jkg)) go to 740
                          if (jkg == jsw) go to 720
                          if (Tt < Temp(1, inc) - gap .or. Tt > Temp(2, inc) + gap) go to 740
                          go to 720
                       end if
                       cycle outerLoop1
                    else if (Ifz(kc+1) <= Ifz(kc)) then
                       cycle outerLoop1
                    end if
                 end do
                 if (Tt > Temp(2, kc) * 1.2d0) go to 1000
              end if
           end do outerLoop1
        end if

        sizeg = 0
        szgj = 0

        do inc = 1, Nc
           j = inc + Ng

           if (Debug(Npt)) write(IOOUT, '(/1x, a15, 2f10.3, 3x, e15.7)') Prod(j), Temp(1, inc), Temp(2, inc), En(j, Npt)

           if (En(j, Npt) <= 0) then
              if (Tt > Temp(1, inc) .or. Temp(1, inc) == Tg(1)) then
                 if (Tt <= Temp(2, inc)) then
                    delg = (H0(j) - S(j) - sum(A(1:Nlm, j) * X(1:Nlm))) / Mw(j)

                    if (delg < sizeg .and. delg < 0) then
                       if (j /= jcons) then
                          sizeg = delg
                          jdelg = j
                       else
                          szgj = delg
                       end if
                       ipr = ipr - 1
                    end if

                    if (Debug(Npt)) write(IOOUT, '(" [G0j-SUM(Aij*PIi)]/Mj =", E15.7, 9X, "MAX NEG DELTA G =", E15.7)') delg, sizeg
                 end if
              end if
           end if
        end do

        if (sizeg == 0 .and. szgj == 0) go to 750

        if (sizeg /= 0) then
           j = jdelg
           go to 800
        else
           write(IOOUT, '(/" REINSERTION OF ", A16, " LIKELY TO CAUSE SINGULARITY, ", "(EQLBRM)")') Prod(jcons)
           go to 1500
        end if

720     kk = max(0, kg)
        tmelt = Temp(kk+1, inc)
        Tt = tmelt
        Tln = log(Tt)
        Jsol = min(j, jkg)
        Jliq = Jsol + 1
        En(jkg, Npt) = 0.5d0 * En(j, Npt)
        En(j, Npt) = En(jkg, Npt)
        j = jkg
        go to 800

! WRONG PHASE INCLUDED FOR T INTERVAL, SWITCH EN
740     En(jkg, Npt) = En(j, Npt)
        Jcond(ipr) = jkg
        En(j, Npt) = 0
        jsw = j

        if (Prod(j) /= Prod(jkg) .and. .not. Short) &
             write(IOOUT, '(" PHASE CHANGE, REPLACE ", A16, "WITH ", A16)') Prod(j), Prod(jkg)

        j = jkg
        go to 900
     end if

! CONVERGED WITH NO CONDENSED CHANGES.  IF BOTH SOLID & LIQ PRESENT, 
! TEMPORARILY REMOVE LIQ TO PREVENT SINGULAR DERIVATIVE MATRIX.
750  Sumn = Enn
     if (Jsol /= 0) then
        ensol = En(Jsol, Npt)
        En(Jsol, Npt) = En(Jsol, Npt) + En(Jliq, Npt)
        Dlvtp(Npt) = 0
        Cpr(Npt) = 0
        Gammas(Npt) = 0
        Pderiv = .true.

        do k = 1, Npr
           if (Jcond(k) == Jliq) exit
        end do

        Jcond(k:Npr) = Jcond(k+1:Npr+1)

        Npr = Npr - 1
     end if

     go to 500
  end if

! ADD CONDENSED SPECIES
800 Npr = Npr + 1

  Jcond(2:Npr) = Jcond(1:Npr-1)
  Jcond(1) = j

  if (.not. Short) write(IOOUT, '(" ADD ", a16)') Prod(j)

900 inc = j - Ng
  Convg = .false.
  if (Tp) cpcalc = .false.
  numb = -1
  go to 500

! REMOVE CONDENSED SPECIES
1000 En(j, Npt) = 0
  Deln(j) = 0
  Enln(j) = 0

  Jcond(k:Npr) = Jcond(k+1:Npr+1)

  if (.not. Short) write(IOOUT, '(" REMOVE ", A16)') Prod(j)

  Npr = Npr - 1
  do i = 1, Nlm
     if (cmp(i) == Prod(j)) then
        numb = -1
        Convg = .false.
        if (Tp) cpcalc = .false.
        go to 1100
     end if
  end do

  go to 900

1100 newcom = .false.
  nn = Nlm

  if (Elmt(Nlm) == 'E') nn = Nlm - 1

! FIND ORDER OF SPECIES FOR COMPONENTS - BIGGEST TO SMALLEST
  njc = 0
  lcs(1:nn) = 0

1200 bigen = -1d-35

  do j = 1, Ng
     if (En(j, Npt) > bigen) then
        if (.not. Ions .or. A(ls, j) == 0) then
           bigen = En(j, Npt)
           jbx = j
        end if
     end if
  end do

  if (bigen > 0.) then
     do 1250 lc = 1, nn
        if (jbx == 0) jbx = Jx(lc)

        if (A(lc, jbx) > smalno) then
           if (njc /= 0) then
              do i = 1, njc
                 l = lcs(i)
                 if (l == lc) go to 1250
                 if (l == 0) exit
                 j = Jcm(l)
                 if (all(A(1:nn, jbx) == A(1:nn, j))) go to 1250
              end do
           end if

           do i = 1, nn
              if (i /= lc .and. abs(A(lc, jbx) * A(i, jx(i)) - A(lc, jx(i)) * A(i, jbx)) <= smalno) go to 1250
           end do

           njc = njc + 1
           if (jbx /= Jcm(lc)) newcom = .true.
           Jcm(lc) = jbx
           lcs(njc) = lc
           go to 1300
        end if
1250 continue

1300 En(jbx, Npt) = -En(jbx, Npt)
     if (njc < nn) go to 1200
  end if

  forall(j = 1:Ng) En(j, Npt) = abs(En(j, Npt))

  if (newcom) then
! SWITCH COMPONENTS
     do lc = 1, nn
        jb = Jcm(lc)

        if (A(lc, jb) == 0) then
           jb = Jx(lc)
           Jcm(lc) = jb
        end if

        tem = A(lc, jb)

        if (tem /= 0) then
           pisave(lc) = H0(jb) - S(jb)

           if (jb <= Ng) pisave(lc) = pisave(lc) + Enln(jb) + Tm
           cmp(lc) = trim(Prod(jb))

! CALCULATE NEW COEFFICIENTS
           if (tem /= 1) then
              B0(lc) = B0(lc) / tem
              B0p(lc, 1) = B0p(lc, 1) / tem
              B0p(lc, 2) = B0p(lc, 2) / tem

              forall(j = 1:Nspx) A(lc, j) = A(lc, j) / tem
           end if

           do i = 1, nn
              if (A(i, jb) /= 0. .and. i /= lc) then
                 tem = A(i, jb)

                 forall(j = 1:Nspx) A(i, j) = A(i, j) - A(lc, j) * tem
                 forall(j = 1:Nspx, abs(A(i, j)) < 1.E-5) A(i, j) = 0

                 B0(i) = B0(i) - B0(lc) * tem
                 B0p(i, 1) = B0p(i, 1) - B0p(lc, 1) * tem
                 B0p(i, 2) = B0p(i, 2) - B0p(lc, 2) * tem
              end if
           end do
        end if
     end do

     if (Debug(Npt)) then
        write(IOOUT, '(/" NEW COMPONENTS")')
        write(IOOUT, '(/2x, 6A12)') (cmp(k), k = 1, nn)
     end if
  end if

  if (Msing /= 0) then
! SWITCH ORDER OF MSING AND NLM COMPONENTS
     reduce = .true.
     lelim = Nlm
     lsing = Nlm

     if (Msing /= Nlm) then

        do j = 1, Nspx
           aa = A(Msing, j)
           A(Msing, j) = A(Nlm, j)
           A(Nlm, j) = aa
        end do

        ja = Jcm(Msing)
        Jcm(Msing) = Jcm(Nlm)
        Jcm(Nlm) = ja
        ae = cmp(Msing)
        cmp(Msing) = cmp(Nlm)
        cmp(Nlm) = ae
        ae = Elmt(Msing)
        Elmt(Msing) = Elmt(Nlm)
        Elmt(Nlm) = trim(ae)
        ja = Jx(Msing)
        Jx(Msing) = Jx(Nlm)
        Jx(Nlm) = ja
        aa = Atwt(Msing)
        Atwt(Msing) = Atwt(Nlm)
        Atwt(Nlm) = aa
        aa = B0(Msing)
        B0(Msing) = B0(Nlm)
        B0(Nlm) = aa
        aa = pisave(Msing)
        pisave(Msing) = pisave(Nlm)
        pisave(Nlm) = aa

        do i = 1, 2
           aa = B0p(Msing, i)
           B0p(Msing, i) = B0p(Nlm, i)
           B0p(Nlm, i) = aa
        end do
     end if
  else if (.not. newcom .and. Trace == 0.) then
     go to 600
  end if

  Msing = 0
  tsize = xsize
  go to 500

1400 Ttt(Npt) = Tt
  Ppp(Npt) = Pp
  Vlm(Npt) = Rr * Enn * Tt / Pp
  Hsum(Npt) = Hsum(Npt) * Tt
  Wm(Npt) = 1. / Enn
  gasfrc = Enn/Totn(Npt)

  if (gasfrc < 0.0001) write(IOOUT, '(/" WARNING!  RESULTS MAY BE WRONG FOR POINT", i3, " DUE TO", &
       & /" LOW MOLE FRACTION OF GASES (", e15.8, ") (EQLBRM)")') Npt, gasfrc

  if (Trace /= 0) then
     forall(j = 1:Ng, lelim == 0 .or. all(A(lelim:ls, j) == 0) .and. Enln(j) > -87) En(j, Npt) = exp(Enln(j))
  end if

  if (Debug(Npt)) write(IOOUT, '(/" POINT=", i3, 3x, "P=", e13.6, 3x, "T=", e13.6, /3x, "H/R=", &
       & e13.6, 3x, "S/R=", e13.6, /3x, "M=", e13.6, 3x, "CP/R=", e13.6, 3x, &
       & "DLVPT=", e13.6, /3x, "DLVTP=", e13.6, 3x, "GAMMA(S)=", e13.6, 3x, "V=", e13.6)') &
       Npt, Pp, Tt, Hsum(Npt), Ssum(Npt), Wm(Npt), Cpr(Npt), Dlvpt(Npt), Dlvtp(Npt), Gammas(Npt), Vlm(Npt)

  if (Tt >= Tg(1) .and. Tt <= Tg(4)) go to 1600

  if (Shock) go to 1600

  write(IOOUT, '(" THE TEMPERATURE=", e12.4, " IS OUT OF RANGE FOR POINT", i5, "(EQLBRM)")') Tt, Npt

  if (Tt >= Tg(1) * 0.8d0 .and. Tt <= Tg(4) * 1.1d0) go to 1600

  Npt = Npt + 1

1500 Tt = 0
  Npt = Npt - 1
  write(IOOUT, '(/" CALCULATIONS STOPPED AFTER POINT", I3, "(EQLBRM)")') Npt

1600 Lsave = Nlm
  Nlm = ls

  if (Npr > 0) Gonly = .false.

  return
end subroutine EQLBRM



subroutine FROZEN
!***********************************************************************
! CALCULATE PROPERTIES WITH FROZEN COMPOSITION AT ASSIGNED ENTROPY
! AND PRESSURE.  CALLED FROM ROCKET.
!***********************************************************************
  use mod_legacy_cea
  implicit none
! LOCAL VARIABLES
  integer:: i, inc, iter, j, k, nnn
  real(8):: dlnt, dlpm

  Convg = .false.
  Tln = log(Tt)
  dlpm = log(Pp * Wm(Nfz))
  nnn = Npt
  Npt = Nfz

  forall(j = 1:Ng, En(j, Nfz) /= 0) Deln(j) = -(log(En(j, Nfz)) + dlpm)

  do iter = 1, 8
     call CPHS

     Cpsum = sum(En(1:Ng, Nfz) * Cp(1:Ng))
     Ssum(nnn) = sum(En(1:Ng, Nfz) * (S(1:Ng) + Deln(1:Ng)))

     if (Npr /= 0) then
        Cpsum = Cpsum + sum(En(Jcond(1:Npr), Nfz) * Cp(Jcond(1:Npr)))
        Ssum(nnn) = Ssum(nnn) + sum(En(Jcond(1:Npr), Nfz) * S(Jcond(1:Npr)))
     end if

     if (Convg) then
        Npt = nnn

        Hsum(Npt) = sum(En(1:Ngc, Nfz) * H0(1:Ngc)) * Tt

        Ttt(Npt) = Tt
        Gammas(Npt) = Cpsum / (Cpsum - 1 / Wm(Nfz))
        Vlm(Npt) = Rr * Tt / (Wm(Nfz) * Pp)
        Wm(Npt) = Wm(Nfz)
        Dlvpt(Npt) = -1
        Dlvtp(Npt) = 1
        Totn(Npt) = Totn(Nfz)
        Ppp(Npt) = Pp
        Cpr(Npt) = Cpsum

        if (Tt >= Tg(1) * 0.8d0) then
           if (all(En(Ngp1:Ngc, Nfz) == 0)) return
           if (all(Temp(1, Ngp1-Ng:Ngc-Ng) - 50 <= Tt .and. Tt <= Temp(2, Ngp1-Ng:Ngc-Ng) + 50)) return
        end if

        Tt = 0
        Npt = Npt - 1
        return

     else
        dlnt = (Ssum(Nfz) - Ssum(nnn)) / Cpsum
        Tln = Tln + dlnt
        if (abs(dlnt) < 0.5d-4) Convg = .true.
        Tt = exp(Tln)
     end if
  end do

  write(IOOUT, '(/" FROZEN DID NOT CONVERGE IN 8 ITERATIONS (FROZEN)")')

  Tt = 0
  Npt = Npt - 1

  return
end subroutine FROZEN



subroutine GAUSS
!***********************************************************************
! SOLVE ANY LINEAR SET OF UP TO maxMat EQUATIONS
! NUMBER OF EQUATIONS = IMAT
!***********************************************************************
  use mod_legacy_cea
  implicit none

! LOCAL VARIABLES
  integer:: i, imatp1, j, k, nn, nnp1
  real(8), parameter:: bigNo = 1e25
  real(8):: coefx(50)
  real(8):: tmp(Imat+1)

! BEGIN ELIMINATION OF NNTH VARIABLE
  imatp1 = Imat + 1

  do nn = 1, Imat
     if (nn /= Imat) then
! SEARCH FOR MAXIMUM COEFFICIENT IN EACH ROW
        nnp1 = nn + 1

        do i = nn, Imat
           coefx(i) = bigNo

           if (G(i, nn) /= 0) then
              coefx(i) = maxval(abs(G(i, nnp1:imatp1)))

              if (coefx(i) < bigNo * abs(G(i, nn))) then
                 coefx(i) = coefx(i) / abs(G(i, nn))
              else
                 coefx(i) = bigNo
              end if
           end if
        end do

        if (all(coefx(nn:Imat) == bigNo)) then
           Msing = nn
           return
        end if

! LOCATE ROW WITH SMALLEST MAXIMUM COEFFICIENT
        tmp(1:1) = minloc(coefx(nn:Imat))
        i = tmp(1) + nn - 1

! INDEX I LOCATES EQUATION TO BE USED FOR ELIMINATING THE NTH
! VARIABLE FROM THE REMAINING EQUATIONS
! INTERCHANGE EQUATIONS I AND NN
        if (i /= nn) then
           tmp(nn:imatp1) = G(i, nn:imatp1)
           G(i, nn:imatp1) = G(nn, nn:imatp1)
           G(nn, nn:imatp1) = tmp(nn:imatp1)
        end if

     else if (G(nn, nn) == 0) then
        Msing = nn
        return
     end if

! DIVIDE NTH ROW BY NTH DIAGONAL ELEMENT AND ELIMINATE THE NTH
! VARIABLE FROM THE REMAINING EQUATIONS
     if (G(nn, nn) == 0) then
        Msing = nn
        return
     else
        forall(j = nn+1:imatp1) G(nn, j) = G(nn, j) / G(nn, nn)

        if (nn + 1 /= imatp1) then
!DIR$ IVDEP
           forall(i = nn+1:Imat, j = nn+1:imatp1) G(i, j) = G(i, j) - G(i, nn) * G(nn, j)
        end if
     end if
  end do

  ! BACKSOLVE FOR THE VARIABLES
  do k = Imat, 1, -1
     X(k) = G(k, imatp1)
     if (Imat >= k+1) then
        X(k) = X(k) - sum(G(k, k+1:Imat) * X(k+1:Imat))
     end if
  end do

  return
end subroutine GAUSS



subroutine HCALC
!***********************************************************************
! CALCULATE PROPERTIES FOR TOTAL REACTANT USING THERMO DATA FOR
! ONE OR MORE REACTANTS. USED ONLY FOR SHOCK AND DETON PROBLEMS.
!***********************************************************************
  use mod_legacy_cea
  implicit none

! LOCAL VARIABLES
  character(6):: date(maxNgc)
  character(2):: el(5)
  character(15):: sub
  integer, save:: i, icf, ifaz, itot, j, k, l, m, n, nall, nint, ntgas, ntot
  real(8):: bb(5), enj, er, sj, t1, t2, tem, thermo(9, 3), Tsave

  Tsave = Tt

  Tm = 0
  if (Pp > 0) Tm = log(Pp * Wmix)

  Ssum(Npt) = 0
  Hpp(1) = 0
  Hpp(2) = 0
  Hsub0 = 0
  Cpmix = 0
  tem = (1 + Oxfl)

! LOOP ON REACTANTS.
! if oxidant, k = 1
! if fuel,    k = 2
  Nspr = Nspx
  do n = 1, Nreac
     k = 2
     if (Fox(n)(:1) == 'O' .or. Fox(n)(:1) == 'o') k = 1
     if (Tt == 0) Tt = Rtemp(n)

     j = Jray(n)

     if (j == 0) then
! SEARCH FOR REACTANT IN STORED THERMO SPECIES. STORE INDEX IN JRAY(N).
        ifaz = 0
        do j = 1, Ngc
           if (Rname(n) == Prod(j) .or. '*' // Rname(n) == Prod(j)) then
              Jray(n) = j
              if (j > Ng) then
                 write(IOOUT, '(/" REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)")')
                 Tt = 0
                 Cpmix = 0
                 return
              end if
              go to 50
           end if
        end do

! SEARCH THERMO.LIB FOR SPECIES.
        rewind IOTHM

        read(IOTHM) Tg, ntgas, ntot, nall

        Nspr = Nspr + 1
        do itot = 1, nall
           if (itot <= ntot) then
              icf = 3
              if (itot > ntgas) icf = 1
              read(IOTHM) sub, nint, date(Nspr), (el(j), bb(j), j = 1, 5), ifaz, &
                   T1, T2, Mw(Nspr), ((thermo(l, m), l = 1, 9), m = 1, icf)
           else
              read(IOTHM) sub, nint, date(Nspr), (el(j), bb(j), j = 1, 5), ifaz, T1, T2, Mw(Nspr), er
              if (nint /= 0) then
                 read(IOTHM) ((thermo(i, j), i = 1, 9), j = 1, nint)
                 icf = nint
              end if
           end if

           if (sub == Rname(n) .or. sub == '*' // Rname(n)) then
              if (ifaz <= 0 .and. nint > 0) then
                 do j = 1, 5
                    if (bb(j) == 0) exit
                    Nfla(n) = j
                    Ratom(n, j) = el(j)
                    Rnum(n, j) = bb(j)
                 end do

                 Jray(n) = Nspr
                 j = Nspr

                 forall(l = 1:icf, m = 1:9) Coef(j, m, l) = thermo(m, l)
                 go to 50

              else
                 if (ifaz > 0) write(IOOUT, '(/" REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)")')
                 if (nint == 0) write(IOOUT, '(/" COEFFICIENTS FOR ", a15, " ARE NOT AVAILABLE (HCALC)")') Rname(n)
                 Tt = 0
                 Cpmix = 0
                 return
              end if
           end if
        end do

        Nspr = Nspr - 1
        write(IOOUT, '(/" ERROR IN DATA FOR ", a15, " CHECK NAME AND TEMPERATURE", &
             & " RANGE IN", /, " thermo.inp (HCALC)")') Rname(n)

        Energy(n) = ' '
        Tt = 0
        Cpmix = 0

        return
     end if

! CALCULATE EN FOR REACTANT AND CALCULATE PROPERTIES.
50   if (Moles) then
        enj = Pecwt(n) / Wp(k)
     else
        enj = Pecwt(n) / Rmw(n)
     end if
     enj = enj / tem
     if (k == 1) enj = enj * Oxfl

     Tln = log(Tt)
     En(j, Npt) = enj
     l = 1

     if (ifaz <= 0) then
        if (Tt > Tg(2)) l = 2
        if (Tt > Tg(3) .and. ifaz < 0) l = 3
     end if

     S(j) = Coef(j, 7, l) / 4 * Tt**4 + Coef(j, 6, l) / 3 * Tt**3 + Coef(j, 5, l) / 2 * Tt**2 + Coef(j, 4, l) * Tt &
          - Coef(j, 1, l) * 0.5d0 / Tt**2 - Coef(j, 2, l) / Tt + Coef(j, 3, l) * Tln + Coef(j, 9, l)

     H0(j) = Coef(j, 7, l) / 5 * Tt**4 + Coef(j, 6, l) / 4 * Tt**3 + Coef(j, 5, l) / 3 * Tt**2 + Coef(j, 4, l) / 2 * Tt &
          - Coef(j, 1, l) / Tt**2 + (Coef(j, 2, l) * Tln + Coef(j, 8, l)) / Tt + Coef(j, 3, l)

     Cp(j) = Coef(j, 7, l) * Tt**4 + Coef(j, 6, l) * Tt**3 + Coef(j, 5, l) * Tt**2 + Coef(j, 4, l) * Tt &
          + Coef(j, 1, l) / Tt**2 + Coef(j, 2, l) / Tt + Coef(j, 3, l)

     if (abs(H0(j)) < 0.01) H0(j) = 0

! ADD CONTRIBUTION TO CP, H, AND S OF TOTAL REACTANT.
     Cpmix = Cpmix + Cp(j) * enj

! FOR CONDENSED SPECIES:  SJ = S(J)
     sj = S(j) - log(enj) - Tm
     Ssum(Npt) = Ssum(Npt) + enj * sj
     er = H0(j) * enj * Tt
     Hsub0 = Hsub0 + er
     Hpp(k) = Hpp(k) + er
  end do

  if (Tsave /= 0) Tt = Tsave

  return
end subroutine HCALC


subroutine MATRIX
!***********************************************************************
! SET UP ITERATION OR DERIVATIVE MATRIX.
!***********************************************************************
  use mod_legacy_cea
  implicit none

! LOCAL VARIABLES
  integer:: i, iq, iq2, iq3, isym, j, k, kk, kmat
  real(8):: energyl, f, h, ss, sss, term, term1

  iq = Nlm + Npr
  Iq1 = iq + 1
  iq2 = Iq1 + 1
  iq3 = iq2 + 1
  kmat = iq3
  if (.not. Convg .and. Tp) kmat = iq2
  imat = kmat - 1
  G = 0
  sss = 0
  Hsum(Npt) = 0

! BEGIN SET-UP OF ITERATION OR DERIVATIVE MATRIX
  do j = 1, Ng
     Mu(j) = H0(j) - S(j) + Enln(j) + Tm
     if (En(j, Npt) /= 0) then
        h = H0(j) * En(j, Npt)
        f = Mu(j) * En(j, Npt)
        ss = h - f
        term1 = h
        if (kmat == iq2) term1 = f

        do i = 1, Nlm
           if (A(i, j) /= 0) then
              term = A(i, j) * En(j, Npt)
              forall(k = i:Nlm) G(i, k) = G(i, k) + A(k, j) * term
              G(i, Iq1) = G(i, Iq1) + term
              G(i, iq2) = G(i, iq2) + A(i, j) * term1
              if (.not. (Convg .or. Tp)) then
                 G(i, iq3) = G(i, iq3) + A(i, j) * f
                 if (Sp) G(iq2, i) = G(iq2, i) + A(i, j) * ss
              end if
           end if
        end do

        if (kmat /= iq2) then
           if (Convg .or. Hp) then
              G(iq2, iq2) = G(iq2, iq2) + H0(j) * h
              if (.not. Convg) then
                 G(iq2, iq3) = G(iq2, iq3) + H0(j) * f
                 G(Iq1, iq3) = G(Iq1, iq3) + f
              end if
           else
              G(iq2, Iq1) = G(iq2, Iq1) + ss
              G(iq2, iq2) = G(iq2, iq2) + H0(j) * ss
              G(iq2, iq3) = G(iq2, iq3) + Mu(j) * ss
              G(Iq1, iq3) = G(Iq1, iq3) + f
           end if
        end if
        G(Iq1, iq2) = G(Iq1, iq2) + term1
     end if
  end do

! CONDENSED SPECIES
  if (Npr /= 0) then
     do k = 1, Npr
        j = Jcond(k)
        kk = Nlm + k
        Mu(j) = H0(j) - S(j)

        forall(i = 1:Nlm)
           G(i, kk) = A(i, j)
           G(i, kmat) = G(i, kmat) - A(i, j) * En(j, Npt)
        end forall

        G(kk, iq2) = H0(j)
        G(kk, kmat) = Mu(j)
        Hsum(Npt) = Hsum(Npt) + H0(j) * En(j, Npt)

        if (Sp) then
           sss = sss + S(j) * En(j, Npt)
           G(iq2, kk) = S(j)
        end if
     end do
  end if

  sss = sss + G(iq2, Iq1)
  Hsum(Npt) = Hsum(Npt) + G(Iq1, iq2)
  G(Iq1, Iq1) = Sumn - Enn

! REFLECT SYMMETRIC PORTIONS OF THE MATRIX
  isym = Iq1
  if (Hp .or. Convg) isym = iq2

  do i = 1, isym
!DIR$ IVDEP
     forall(j = i:isym) G(j, i) = G(i, j)
  end do

! COMPLETE THE RIGHT HAND SIDE
  if (.not. Convg) then
     forall(i = 1:Nlm) G(i, kmat) = G(i, kmat) + B0(i) - G(i, Iq1)
     G(Iq1, kmat) = G(Iq1, kmat) + Enn - Sumn

! COMPLETE ENERGY ROW AND TEMPERATURE COLUMN
     if (kmat /= iq2) then
        if (Sp) energyl = S0 + Enn - Sumn - sss
        if (Hp) energyl = Hsub0/Tt - Hsum(Npt)
        G(iq2, iq3) = G(iq2, iq3) + energyl
        G(iq2, iq2) = G(iq2, iq2) + Cpsum
     end if

  else
     if (Pderiv) then
! PDERIV = .true.-- SET UP MATRIX TO SOLVE FOR DLVPT
        G(Iq1, iq2) = Enn
        forall(i = 1:iq) G(i, iq2) = G(i, Iq1)
     end if
     G(iq2, iq2) = G(iq2, iq2) + Cpsum
  end if

  if (Vol .and. .not. Convg) then
! CONSTANT VOLUME MATRIX
     if (kmat == iq2) then
        forall(i = 1:iq) G(i, Iq1) = G(i, iq2)
     else

!DIR$ IVDEP
        forall(i = 1:iq)
           G(Iq1, i) = G(iq2, i) - G(Iq1, i)
           G(i, Iq1) = G(i, iq2) - G(i, Iq1)
           G(i, iq2) = G(i, iq3)
        end forall

        G(Iq1, Iq1) = G(iq2, iq2) - G(Iq1, iq2) - G(iq2, Iq1)
        G(Iq1, iq2) = G(iq2, iq3) - G(Iq1, iq3)
        if (Hp) G(Iq1, iq2) = G(Iq1, iq2) + Enn
     end if

     kmat = imat
     imat = imat - 1
  end if

  return
end subroutine MATRIX



subroutine NEWOF
!***********************************************************************
! CALCULATE NEW VALUES OF B0 AND HSUB0 FOR NEW OF RATIO
!***********************************************************************
  use mod_legacy_cea
  implicit none

! LOCAL VARIABLES
  integer:: i, j
  real(8):: assval, bigb, bratio, dbi, smalb, tem, v1, v2

  if (.not. Short) write(IOOUT, '(/" O/F = ", f10.6)') Oxfl
  Eqrat = 0
  tem = Oxfl + 1
  v2 = (Oxfl * Vmin(1) + Vmin(2)) / tem
  v1 = (Oxfl * Vpls(1) + Vpls(2)) / tem
  if (v2 /= 0.) Eqrat = abs(v1/v2)
  do i = 1, Nlm
     B0(i) = (Oxfl * B0p(i, 1) + B0p(i, 2)) / tem
     dbi = abs(B0(i))
     if (i == 1) then
        bigb = dbi
        smalb = dbi
     else if (dbi /= 0.) then
        if (dbi < smalb) smalb = dbi
        if (dbi > bigb) bigb = dbi
     end if
  end do
  Bcheck = bigb * .000001d0

! CALCUALTE MOLECULAR WEIGHT OF TOTAL REACTANT, WMIX.
  if (all(Am /= 0)) then
     Wmix = (Oxfl + 1) * Am(1) * Am(2) / (Am(1) + Oxfl * Am(2))
  else
     Wmix = Am(2)
     if (Am(2) == 0.0) Wmix = Am(1)
  end if
  Npt = 1

! IF ASSIGNED U OR H NOT GIVEN IN PROB DATA, INITIAL HSUB0 = 1.d30
  if (Size == 0) assval = Hsub0
  if (assval >= 1.d30) Hsub0 = (Oxfl * Hpp(1) + Hpp(2)) / tem

! NOTE THAT "BRATIO" IS "BRATIO" IN SEC 3.2 IN RP-1311.
  bratio = smalb / bigb
  Size = 18.420681d0
  if (bratio < 1.d-5) Size = log(1000/bratio)
  Jsol = 0
  Jliq = 0

  if (.not. Short) then
     write(IOOUT, '(/, 23x, "EFFECTIVE FUEL", 5x, "EFFECTIVE OXIDANT", 8x, "MIXTURE")')
     if (Vol) write(IOOUT, '(" INTERNAL ENERGY", 11x, "u(2)/R", 14x, "u(1)/R", 14x, "u0/R")')
     if (.not. Vol) write(IOOUT, '(" ENTHALPY", 18x, "h(2)/R", 14x, "h(1)/R", 15x, "h0/R")')
     write(IOOUT, '(" (KG-MOL)(K)/KG", 4x, e18.8, 2e20.8)') Hpp(2), Hpp(1), Hsub0
     write(IOOUT, '(/" KG-FORM.WT./KG", 13x, "bi(2)", 15x, "bi(1)", 15x, "b0i")')
  end if

  do i = 1, Nlm
     j = Jcm(i)
     if (.not. Short) write(IOOUT, '(1x, a16, 3e20.8)') Prod(j), B0p(i, 2), B0p(i, 1), B0(i)
  end do

  return
end subroutine NEWOF



subroutine REACT
!***********************************************************************
! READ AND PROCESS REACTANT RECORDS.  CALLED FROM subroutine INPUT.
!***********************************************************************
  use mod_cea
  use mod_legacy_cea
  implicit none
! LOCAL VARIABLES
  character(6):: date
  character(2):: el(5)
  character(15):: sub
  integer:: i, icf, ifaz, ifrmla, itot, j, jj, k, kk, kr, l, n, nall, nint, nj, ntgas, ntot
  logical:: fuel, rcoefs, wdone(2)
  logical:: hOK
  real(8):: bb(5), dat(35), dift, eform, pcwt, rcf(9, 3), rm, T1, T2
  real(8):: T1save, T2save

  wdone = .false.
  Wp = 0
  Hpp = 0
  Vpls = 0
  Vmin = 0
  Am = 0
  Rh = 0
  Elmt = ' '
  B0p = 0
  dat = 0

! IF OXIDANT, KR = 1
! IF FUEL, KR = 2
  do n = 1, Nreac
     hOK = .false.
     T1save = 20000
     T2save = 0
     rcoefs = .true.

     if (Energy(n) == 'lib' .or. Rnum(n, 1) == 0) then
        Tt = Rtemp(n)
        rewind IOTHM
        read(IOTHM) Tg, ntgas, ntot, nall

        do itot = 1, nall
           if (itot <= ntot) then
              icf = 3
              if (itot > ntgas) icf = 1
              read(IOTHM) sub, nint, date, (el(j), bb(j), j = 1, 5), ifaz, T1, T2, rm, ((rcf(i, j), i = 1, 9), j = 1, icf)
           else
              read(IOTHM) sub, nint, date, (el(j), bb(j), j = 1, 5), ifaz, T1, T2, rm, eform
              if (nint > 0) read(IOTHM) ((rcf(i, j), i = 1, 9), j = 1, nint)
           end if

           if (sub == Rname(n) .or. sub == '*' // Rname(n)) then
              if (nint == 0) then
                 rcoefs = .false.
                 hOK = .true.
                 Enth(n) = eform * 1000 / Rr

                 if (Tt == 0) then
                    Tt = T1
                    Rtemp(n) = T1
                 else
                    dift = abs(Tt - T1)
                    if (dift > 1) then
                       if (dift > 10) then
                          write(IOOUT, '(/" REACTANT ", a15, "HAS BEEN DEFINED FOR THE TEMPERATURE", &
                               & f8.2, "K ONLY."/" YOUR TEMPERATURE ASSIGNMENT", f8.2, &
                               & " IS MORE THAN 10 K FROM THIS VALUE. (REACT)")') Rname(n), T1, Tt
                          Nlm = 0
                          hOK = .false.
                          return
                       else
                          write(IOOUT, '(/" NOTE! REACTANT ", a15, "HAS BEEN DEFINED FOR ", &
                               & "TEMPERATURE", f8.2, "K ONLY."/" YOUR TEMPERATURE ASSIGNMENT", &
                               & f8.2, " IS NOT = BUT <10 K FROM THIS VALUE. (REACT)")') Rname(n), T1, Tt
                          Tt = T1
                          Rtemp(n) = T1
                       end if
                    end if
                 end if
              else
                 if (ifaz <= 0) then
                    T1save = min(T1save, 0.8d0 * Tg(1))
                    T2save = max(T2save, 1.2d0 * T2)
                 else
                    T1save = min(T1save, T1 - 0.001d0)
                    T2save = max(T2save, T2 + 0.001d0)
                 endif
                 if (T1save < Tt .and. Tt < T2save) hOK = .true.
              end if

              do j = 1, 5
                 if (bb(j) == 0) exit
                 Nfla(n) = j
                 Ratom(n, j) = el(j)
                 Rnum(n, j) = bb(j)
              end do

              if (Tt == 0) then
                 if (.not. Hp) go to 50
                 write(IOOUT, '(/" TEMPERATURE MISSING FOR REACTANT NO.", i2, "(REACT)")') n
                 Nlm = 0
                 return
              end if

              if (rcoefs .and. hOK) then
                 Tln = log(Tt)
                 l = 1
                 if (ifaz <= 0) then
                    if (Tt > Tg(2)) l = 2
                    if (Tt > Tg(3)) l = 3
                 end if

                 Enth(n) = (((((rcf(7, l)/5) * Tt + rcf(6, l) / 4) * Tt + rcf(5, l) / 3) * Tt &
                      + rcf(4, l) / 2) * Tt + rcf(3, l)) * Tt - rcf(1, l) / Tt + rcf(2, l) * Tln + rcf(8, l)

                 if (Vol .and. ifaz <= 0) Enth(n) = Enth(n) - Tt
              end if

              if (hOK) go to 50
           end if
        end do

        if (.not. hOK) then
           write(IOOUT, '(/" YOUR ASSIGNED TEMPERATURE", f8.2, "K FOR ", a15, /, &
                & "IS OUTSIDE ITS TEMPERATURE RANGE", f8.2, " TO", f9.2, "K (REACT)")') Tt, Rname(n), T1save, T2save
           Energy(n) = ' '
           Nlm = 0
           return
        endif
     end if

50   ifrmla = Nfla(n)
     if (Fox(n)(:1) == 'f') then
        fuel = .true.
        kr = 2
        Fox(n) = 'FUEL'
     else if (Fox(n)(:4) == 'name') then
        fuel = .true.
        kr = 2
        Fox(n) = 'NAME'
     else
        kr = 1
        Fox(n) = 'OXIDANT'
     end if

     dat = 0

! STORE ATOMIC SYMBOLS IN ELMT ARRAY.
! CALCULATE MOLECULAR WEIGHT.
! TEMPORARILY STORE ATOMIC VALENCE IN X.
     rm = 0
     do jj = 1, ifrmla
        do j = 1, maxEl
           nj = j
           if (Elmt(j) == ' ') exit
           if (Ratom(n, jj) == Elmt(j)) go to 80
        end do

        Nlm = nj
        Elmt(j) = Ratom(n, jj)

80      do kk = 1, 100
           if (atomic_symbol(kk) == Ratom(n, jj)) then
              rm = rm + Rnum(n, jj) * atomic_mass(kk)
              Atwt(j) = atomic_mass(kk)
              X(j) = atomic_valence(kk)
              dat(j) = dat(j) + Rnum(n, jj)
              go to 100
           end if
        end do

        write(IOOUT, '(/1x, a2, " NOT FOUND IN BLOCKDATA (REACT)")') Ratom(n, jj)
        Nlm = 0
        return
100  end do

     if (Pecwt(n) < 0) then
        Pecwt(n) = 0
        if (.not. Moles .and. .not. wdone(kr)) then
           wdone(kr) = .true.
           Pecwt(n) = 100.
           write(IOOUT, '(/" WARNING!!  AMOUNT MISSING FOR REACTANT", i3, ".", &
                & /" PROGRAM SETS WEIGHT PERCENT = 100. (REACT)")') n
        else
           write(IOOUT, '(/" AMOUNT MISSING FOR REACTANT NO.", i2, "(REACT)")') n
           Nlm = 0
           return
        end if
     end if

! ADD CONTRIBUTIONS TO WP(K), HPP(K), AM(K), AND B0P(I, K)
     if (Pecwt(n) > 0) wdone(kr) = .true.
     pcwt = Pecwt(n)
     if (Moles) pcwt = pcwt * rm
     Wp(kr) = Wp(kr) + pcwt

     if (rm <= 0) then
        Nlm = 0
        return
     else
        Hpp(kr) = Hpp(kr) + Enth(n) * pcwt / rm
        Am(kr) = Am(kr) + pcwt / rm
        if (Dens(n) /= 0) then
           Rh(kr) = Rh(kr) + pcwt / Dens(n)
        else
           Rh = 0
        end if

        forall(j = 1:Nlm)
           B0p(j, kr) = dat(j) * pcwt / rm + B0p(j, kr)
        end forall

        Rmw(n) = rm
     end if

  end do

  if (.not. fuel) then
! 100 PERCENT OXIDANT, SWITCH INDICES
     Fox = ' '
     Wp(2) = Wp(1)
     Wp(1) = 0
     Hpp(2) = Hpp(1)
     Am(2) = Am(1)
     Am(1) = 0
     B0p(:, 2) = B0p(:, 1)
  end if

  if (Nlm /= 0) then
! NORMALIZE HPP(KKR), AM(KR), B0P(I, KR), AND PECWT(N).
! CALCULATE V+(KR), AND V-(KR)
     do kr = 1, 2
        if (Wp(kr) /= 0) then
           Hpp(kr) = Hpp(kr) / Wp(kr)
           Am(kr) = Wp(kr) / Am(kr)
           if (Rh(kr) /= 0) Rh(kr) = Wp(kr) / Rh(kr)
           do j = 1, Nlm
              B0p(j, kr) = B0p(j, kr) / Wp(kr)
              if (X(j) < 0) Vmin(kr) = Vmin(kr) + B0p(j, kr) * X(j)
              if (X(j) > 0) Vpls(kr) = Vpls(kr) + B0p(j, kr) * X(j)
           end do
        end if
     end do

     if (.not. Moles) then
        do n = 1, Nreac
           if (Fox(n)(:1) == 'O') then
              Pecwt(n) = Pecwt(n) / Wp(1)
           else
              Pecwt(n) = Pecwt(n) / Wp(2)
           end if
        end do
     end if

     if (.not. Short) then
        if (Moles) then
           write(IOOUT, '(/4x, "REACTANT", 10x, a7, 3x, "(ENERGY/R),K", 3x, &
                & "TEMP,K  DENSITY"/, 8x, "EXPLODED FORMULA")') ' MOLES '
        else
           write(IOOUT, '(/4x, "REACTANT", 10x, a7, 3x, "(ENERGY/R),K", 3x, &
                & "TEMP,K  DENSITY"/, 8x, "EXPLODED FORMULA")') 'WT.FRAC'
        end if
        do n = 1, Nreac
           write(IOOUT, '(1x, a1, ": ", a15, f10.6, e15.6, f9.2, f8.4, /8x, 5(2x, a2, f8.5))') &
                Fox(n), Rname(n), Pecwt(n), Enth(n), Rtemp(n), Dens(n), (Ratom(n, i), Rnum(n, i), i = 1, Nfla(n))
        end do
     end if

  end if

  return
end subroutine REACT



subroutine RKTOUT(cea)
!***********************************************************************
! SPECIAL OUTPUT FOR ROCKET PROBLEMS.
!***********************************************************************
  use mod_legacy_cea
  use mod_legacy_io
  implicit none

  type(CEA_Problem), intent(inout):: cea

! LOCAL VARIABLES
  character(4):: exit(11) = 'EXIT'
  character(15):: fi, fiv, fr, z(4)
  integer, save:: i, i23, i46, i57, i68, i79, ione, ixfr, ixfz, j, k, line, ln, mae, mcf, misp, mivac, mmach, mppf, mppj, nex
  real(8):: agv, aw, gc, tem, tra, vaci(Ncol), ww


  if (.not. Eql) then
     write(IOOUT, '(/////10x, " THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION")')
     if (Nfz > 1) write(IOOUT, '(33x, "AFTER POINT", i2)') Nfz
  else
     write(IOOUT, '(/////13x, " THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM")')
     if (Iopt /= 0) then
        write(IOOUT, '(/11x, " COMPOSITION DURING EXPANSION FROM FINITE AREA COMBUSTOR")')
     else
        write(IOOUT, '(/10x, " COMPOSITION DURING EXPANSION FROM INFINITE AREA COMBUSTOR")')
     end if
  end if

  if (Ttt(1) == T(It)) write(IOOUT, '(25X, "AT AN ASSIGNED TEMPERATURE  ")')

  tem = Ppp(1) * 14.696006d0 / 1.01325d0
  write(IOOUT, '(/1x, a3, " =", f8.1, " PSIA")') 'Pin', tem

  i23 = 2
  if (Iopt > 0) then
     if (Iopt == 1) write(IOOUT, '(" Ac/At =", f8.4, 6x, "Pinj/Pinf =", f10.6)') Subar(1), App(2)
     if (Iopt == 2) write(IOOUT, '(" MDOT/Ac =", f10.3, " (KG/S)/M**2", 6x, "Pinj/Pinf =", f10.6)') Ma, App(2)
     i23 = 3
  end if

  call OUT1

  cea%fmt(4) = cea%fmt(6)
  nex = Npt - 2
  if (Page1) then
     ione = 0
     i46 = 4
     i57 = 5
     i68 = 6
     i79 = 7
  else
     ione = i23
  end if

! PRESSURE RATIOS
  if (Iopt == 0) then
     write(IOOUT, '(/17X, "CHAMBER   THROAT", 11(5X, A4))') (exit(i), i = 1, nex)
     call VARFMT(cea, App)
     write(IOOUT, cea%fmt) 'Pinf/P         ', (App(j), j = 1, Npt)
  else
     nex = nex - 1
     write(IOOUT, '(/, 17X, "INJECTOR  COMB END  THROAT", 10(5X, A4))') (exit(i), i = 1, nex)
     X(1) = 1.d0
     do i = 2, Npt
        X(i) = Ppp(1) / Ppp(i)
     end do
     call VARFMT(cea, X)
     write(IOOUT, cea%fmt) 'Pinj/P         ', (X(i), i = 1, Npt)
  end if

  call OUT2(cea)

  mppf  = 0
  mppj  = 0
  mmach = 0
  mae   = 0
  mcf   = 0
  mivac = 0
  misp  = 0

  do i = 1, Nplt
     ixfz = index(Pltvar(i)(2:), 'fz')
     ixfr = index(Pltvar(i)(2:), 'fr')

     if (ixfz /= 0 .or. ixfr /= 0) then
        if (Eql) cycle
     else if (.not. Eql) then
        cycle
     end if

     if (Pltvar(i)(:4) == 'pi/p' .or. Pltvar(i)(:3) == 'pip') then
        if (Iopt == 0) mppf = i
        if (Iopt /= 0) mppj = i
     else if (Pltvar(i)(:4) == 'mach') then
        mmach = i
     else if (Pltvar(i)(:2) == 'ae') then
        mae = i
     else if (Pltvar(i)(:2) == 'cf') then
        mcf = i
     else if (Pltvar(i)(:4) == 'ivac') then
        mivac = i
     else if (Pltvar(i)(:3) == 'isp') then
        misp = i
     end if
  end do

  if (SIunit) then
     agv = 1
     gc = 1
     fr = 'CSTAR, M/SEC'
     fiv = 'Ivac, M/SEC'
     fi = 'Isp, M/SEC'
  else
     gc = 32.174
     agv = 9.80665
     fr = 'CSTAR, FT/SEC'
     fiv = 'Ivac,LB-SEC/LB'
     fi = 'Isp, LB-SEC/LB'
  end if

  do k = 2, Npt
     Spim(k) = sqrt(2 * Rr * (Hsum(1) - Hsum(k))) / agv
! AW IS THE LEFT SIDE OF EQ.(6.12) IN RP-1311, PT I.
     aw = Rr * Ttt(k) / (Ppp(k) * Wm(k) * Spim(k) * agv**2)
     if (k == i23) then
        if (Iopt == 0) Cstr = gc * Ppp(1) * aw
        if (Iopt /= 0) Cstr = gc * Ppp(1) / App(2) * aw
     end if
     vaci(k) = Spim(k) + Ppp(k) * aw
     Vmoc(k) = 0
     if (Sonvel(k) /= 0) Vmoc(k) = Spim(k) * agv / Sonvel(k)
  end do

! MACH NUMBER
  Vmoc(1) = 0
  if (Gammas(i23) == 0) Vmoc(i23) = 0
  cea%fmt(7) = '3,'
  write(IOOUT, cea%fmt) 'MACH NUMBER    ', (Vmoc(j), j = 1, Npt)
  if (Trnspt) call OUT4(cea)
  write(IOOUT, '(/" PERFORMANCE PARAMETERS"/)')

! AREA RATIO
  cea%fmt(4) = '9x,'
  cea%fmt(i46) = '9x,'
  call VARFMT(cea, Aeat)
  cea%fmt(5) = ' '
  cea%fmt(i57) = ' '
  write(IOOUT, cea%fmt) 'Ae/At          ', (Aeat(j), j = 2, Npt)

! C*
  cea%fmt(i57) = '13'
  cea%fmt(i68) = cea%fmt(i68 + 2)
  cea%fmt(i79) = '1,'
  write(IOOUT, cea%fmt) fr, (Cstr, j = 2, Npt)

! CF - THRUST COEFICIENT
  cea%fmt(i79) = '4,'
  do i = 2, Npt
     X(i) = gc * Spim(i) / Cstr
  end do
  write(IOOUT, cea%fmt) 'CF             ', (X(j), j = 2, Npt)

! VACUUM IMPULSE
  cea%fmt(i57) = '13'
  cea%fmt(i79) = '1,'
  write(IOOUT, cea%fmt) fiv, (vaci(j), j = 2, Npt)

! SPECIFIC IMPULSE
  write(IOOUT, cea%fmt) fi, (Spim(j), j = 2, Npt)

  if (Nplt > 0) then
     Spim(1) = 0
     Aeat(1) = 0
     Vmoc(1) = 0
     vaci(1) = 0
     X(1) = 0
     Spim(1) = 0
     do i = ione + 1, Npt
        if (mppj > 0)  Pltout(i+Iplt-ione, mppj)  = Ppp(1) / Ppp(i)
        if (mppf > 0)  Pltout(i+Iplt-ione, mppf)  = App(i)
        if (mmach > 0) Pltout(i+Iplt-ione, mmach) = Vmoc(i)
        if (mae > 0)   Pltout(i+Iplt-ione, mae)   = Aeat(i)
        if (mcf > 0)   Pltout(i+Iplt-ione, mcf)   = X(i)
        if (mivac > 0) Pltout(i+Iplt-ione, mivac) = vaci(i)
        if (misp > 0)  Pltout(i+Iplt-ione, misp)  = Spim(i)
     end do
  end if

  write(IOOUT, *)
  cea%fmt(4) = ' '
  cea%fmt(5) = '13'
  cea%fmt(7) = '5,'

  if (Iopt /= 0) then
     cea%fmt(i46) = cea%fmt(8)
     cea%fmt(i57) = cea%fmt(9)
  end if

  if (.not. Eql) then
     if (Massf) then
        write(IOOUT, '(1x, A4, " FRACTIONS"/)') 'MASS'
     else
        write(IOOUT, '(1x, A4, " FRACTIONS"/)') 'MOLE'
        ww = 1 / Totn(Nfz)
     end if

! MOLE (OR MASS) FRACTIONS - FROZEN
     tra = 5.E-6
     if (Trace /= 0) tra = Trace
     line = 0

     do k = 1, Ngc
        if (Massf) ww = Mw(k)
        X(line+1) = En(k, Nfz) * ww

        if (X(line+1) >= tra) then
           line = line + 1
           z(line) = Prod(k)
        end if

        if (line == 3 .or. k == Ngc) then
           if (line == 0) then
              call OUT3(cea)
              return
           end if
           write(IOOUT, '(1x, 3(a15, f8.5, 3x))') (z(ln), X(ln), ln = 1, line)
           line = 0
        end if
     end do
  end if

  call OUT3(cea)
  return
end subroutine RKTOUT



subroutine ROCKET(cea)
!***********************************************************************
! EXECUTIVE ROUTINE FOR ROCKET PROBLEMS.
!***********************************************************************
  use mod_cea
  use mod_legacy_cea
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

  iplte = Iplt
  isup1 = 1
  App(1) = 1
  Iopt = 0
  Npp = Npp + 2
  nn = Npp
  i01 = 0
  i12 = 1
  nipp = 1
  nptth = 2
  if (Fac) then
     Eql = .true.
     Npp = Npp + 1
     if (Acat /= 0) then
        Iopt = 1
     else if (Ma /= 0) then
        Iopt = 2
     else
        write(IOOUT, '(/" FATAL ERROR!! EITHER mdot OR ac/at MISSING FOR fac PROBLEM (ROCKET)")')
        Tt = 0
        go to 1400
     end if
     i01 = 1
     i12 = 2
     nipp = 2
     nptth = 3
     do i = Nsub, 1, - 1
        Subar(i+1) = Subar(i)
     end do
     Nsub = Nsub + 1
     if (Iopt /= 1) then
        if (Acat == 0) Acat = 2
     end if
     Subar(1) = Acat
  else if (.not. Eql .and. Nfz > 1 .and. Nsub > 0) then
     Nsub = 0
     write(IOOUT, '(/" WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ", /, &
          & " RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)")')
  end if
  nn = nn + Nsub + Nsup
  if (Nfz > 2 .and. nn > Ncol-2) then
     write(IOOUT, '(/" WARNING!!  nfz NOT ALLOWED TO BE > 2 IF THE TOTAL", /, &
          & " NUMBER OF POINTS IS >", i3, " (ROCKET)")') Ncol - 2
     Nfz = 1
     Froz = .false.
  end if
  seql = Eql
  iof = 0
  Tt = Tcest
  Pp = P(1)
  App(i12) = 1
! LOOP FOR EACH O/F
100 It = 1
  iof = iof + 1
  Oxfl = Oxf(iof)
  if (T(1) /= 0.) then
     Tp = .true.
  else
     Hp = .true.
  end if
  Sp = .false.
  call NEWOF
  if (T(1) /= 0) Tt = T(1)
! LOOP FOR CHAMBER PRESSURES
200 do Ip = 1, Np
     itnum = 0
     Area = .false.
     if (T(1) == 0) Hp = .true.
     if (T(1) /= 0) Tp = .true.
     Sp = .false.
     Eql = .true.
     isub = 1
     Isup = 1
     Pp = P(Ip)
     pinf = Pp
     ipp = 1
     itrot = 3
     isupsv = 1
     niter = 1
     Page1 = .true.
     iplt1 = iplte
     Iplt = iplte
     done = .false.
! LOOP FOR OUTPUT COLUMNS
250  nar = Npt
     if (Eql) then
        call EQLBRM
        if (Npt == Nfz) cprf = Cpsum
     else
        call FROZEN
     end if
! TT = 0 IF NO CONVERGENCE
     if (Tt /= 0.) then
! TEST FOR FINITE AREA COMBUSTOR
        if (.not. Fac) go to 400
        pinjas = P(Ip) * pa
        pinj = pinjas
        if (Npt <= 2) then
           if (Npt == 1 .and. Trnspt) call TRANP
           if (Npt == 2) pinf = Ppp(2)
        end if
        if (Npt /= 1) go to 400
! INITIAL ESTIMATE FOR PC (AND ACAT IF NOT ASSIGNED)
        do i = 1, 4
           prat = (b1 + c1 * Acat) / (1 + a1l * Acat)
           ppa = pinj * prat
           if (Iopt == 1) go to 260
           Acat = ppa / (Ma * 2350)
           if (Acat >= 1) then
              pratsv = prat
              if (Debugf) then
                 if (i <= 1) write(IOOUT, '(/"  ITERATION", 9x, "PC", 7x, "CONTRACTION RATIO")')
                 write(IOOUT, '(5x, i2, 7x, f12.2, 3x, f12.6)') i, ppa, Acat
              end if
           else
              write(IOOUT, '(/" INPUT VALUE OF mdot/a =", f12.3, " IS TOO LARGE."/ &
                   & " GIVES CONTRACTION RATIO ESTIMATE LESS THAN 1 (ROCKET)")') Ma
              Tt = 0
              go to 1400
           end if
        end do
        Subar(1) = Acat
260     Pp = ppa / pa
        App(1) = Pp / Ppp(1)
        go to 1100
     else
        if (Npt < 1) go to 1400
        if (.not. Area) go to 600
        Npt = nar - 1
        Isup = Nsup + 2
        Isv = 0
        itnum = 0
        go to 950
     end if
300  Hp = .true.
     Sp = .false.
     niter = niter + 1
     Isv = 0
     Npt = 2
     ipp = 2
     call SETEN
     go to 250
350  done = .true.
     App(1) = Ppp(2) / Ppp(1)
     Area = .false.
     if (Nsub > 1) isub = 2
     Isv = 4
     Npt = 2
     ipp = min(4, Npp)
     call SETEN
     Cpr(2) = Cpr(4)
     Dlvpt(2) = Dlvpt(4)
     Dlvtp(2) = Dlvtp(4)
     Gammas(2) = Gammas(4)
     Hsum(2) = Hsum(4)
     Ppp(2) = Ppp(4)
     App(2) = Ppp(1)/pinf
     Ssum(2) = Ssum(4)
     Totn(2) = Totn(4)
     Ttt(2) = Ttt(4)
     Vlm(2) = Vlm(4)
     Wm(2) = Wm(4)
     if (.not. Short) write(IOOUT, '(" END OF CHAMBER ITERATIONS")')
     go to 600
! INITIALIZE FOR THROAT
400  if (ipp > nipp) then
        usq = 2 * (Hsum(1) - Hsum(Npt)) * Rr
        if (ipp > nptth) go to 600
! THROAT
        if (.not. thi) then
           Vv = Vlm(nptth)
           pvg = Pp * Vv * Gammas(nptth)
           if (pvg == 0) then
              write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
              go to 550
           else
              msq = usq / pvg
              if (Debug(1) .or. Debug(2)) write(IOOUT, '(/" USQ=", e15.8, 5x, "PVG=", e15.8)') usq, pvg
              dh = abs(msq - 1)
              if (dh <= 0.4d-4) go to 550
              if (itrot > 0) then
                 p1 = Pp
                 if (Jsol /= 0) then
                    tmelt = Tt
                    Pp = Pp * (1 + msq * Gammas(nptth)) / (Gammas(nptth) + 1)
                 else if (tmelt == 0) then
                    Pp = Pp * (1 + msq * Gammas(nptth)) / (Gammas(nptth) + 1)
                 else
                    write(IOOUT, '(/" WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)")')
                    dlt = log(tmelt / Tt)
                    dd = dlt * Cpr(nptth) / (Enn * Dlvtp(nptth))
                    Pp = Pp * EXP(dd)
                    App(nptth) = P(Ip) / Pp
                    if (Fac) App(nptth) = pinf / Pp
                    if (Eql .and. .not. Short) write(IOOUT, '(" Pinf/Pt =", F9.6)') App(nptth)
                    thi = .true.
                    go to 250
                 end if
                 go to 500
              else if (itrot < 0) then
                 if (itrot < -19) then
                    write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                    go to 550
                 else
                    if (Npr /= npr1) go to 550
                    Pp = Pp - dp
                    go to 500
                 end if
              else if (Npr == npr1) then
                 write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                 go to 550
              else
                 dp = abs(Pp - p1) / 20
                 Pp = max(Pp, p1)
                 write(IOOUT, '(/" WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)")')
                 Pp = Pp - dp
                 go to 500
              end if
           end if
        else
           Gammas(nptth) = 0
           go to 550
        end if
     else
        if (.not. Fac .and. Trnspt) call TRANP
        if (Npt == Nfz) Eql = seql
        Tp = .false.
        Hp = .false.
        Sp = .true.
        S0 = Ssum(i12)
     end if
450  tmelt = 0
     itrot = 3
     thi = .false.
     App(nptth) = ((Gammas(i12) + 1) / 2)**(Gammas(i12) / (Gammas(i12) - 1))
     if (Eql .and. .not. Short) write(IOOUT, '(" Pinf/Pt =", f9.6)') App(nptth)
     Pp = Pinf / App(nptth)
     Isv = -i12
     go to 1200
500  npr1 = Npr
     App(nptth) = P(Ip) / Pp
     if (Fac) App(nptth) = Pinf / Pp
     if (Eql .and. .not. Short) write(IOOUT, '(" Pinf/Pt =", f9.6)') App(nptth)
     itrot = itrot - 1
     go to 250
550  Awt = Enn * Tt / (Pp * sqrt(usq))
     pcplt = log(App(nptth))
600  Isv = 0
     Aeat(Npt) = Enn * Ttt(Npt) / (Pp * sqrt(usq) * Awt)
     if (Tt == 0) go to 1150
     if (Area) go to 750
     if (Trnspt .and. (.not. Fac .or. done .or. Npt > 2)) call TRANP
     if (Npt == Nfz) Eql = seql
     if (Fac) then
        if (Npt == nptth) then
           Area = .true.
           go to 750
        else if (Npt == 2 .and. done) then
           Npt = 3
!  The following statement was corrected 1/30/2004.  Only fac parameters 
!    after combustion were affected--generally extra or missing points.
!  (remove) if (ipp <= Npp) ipp = ipp - 1
           if (ipp < Npp .or. npp == 4) ipp = ipp - 1
        end if
     end if
650  if (ipp < Npp) go to 1100
700  if (Nsub == i01 .and. Nsup == 0) go to 1150
     Area = .true.
! PCP ESTIMATES FOR AREA RATIOS
750  if (itnum == 0) then
        dlnp = 1
        itnum = 1
        aratio = Subar(isub)
        if ((.not. Fac .or. done) .and. Nsub <= i01) aratio = Supar(Isup)
        if (.not. Eql .and. Nfz >= 3) then
           if (aratio <= Aeat(Nfz)) then
              write(IOOUT, '(/, " WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED", &
                   & " WHERE THE ASSIGNED", /, " SUPERSONIC AREA RATIOS WERE ", &
                   & "LESS THAN THE VALUE AT POINT nfz =", i3, " (ROCKET)")') Nfz
              go to 1050
           end if
        end if
        if (aratio  <  1) then
           write(IOOUT, '(/" AN ASSIGNED AREA RATIO IS < 1 (ROCKET)")')
           go to 1050
        end if
        eln = log(aratio)
        if (Fac) then
           if (.not. done) go to 800
        end if
        if (Nsub <= i01) then
           if (Nfz == ipp) isupsv = Isup
           if (Supar(Isup) < 2) then
              appl = sqrt(eln*(1.535d0 + 3.294d0 * eln)) + pcplt
              go to 1100
           else
              if (Isup > isup1 .and. Supar(Isup-1) >= 2) go to 850
              appl = Gammas(nptth) + eln * 1.4
              go to 1100
           end if
        end if
! TEST FOR CONVERGENCE ON AREA RATIO.
     else if (Gammas(Npt) > 0) then
        check = 0.00004
        if (Debug(Npt)) write(IOOUT, '(/" ITER=", i2, 2x, "ASSIGNED AE/AT=", f14.7, 3x, "AE/AT=", f14.7, &
             & /, 2x, "PC/P=", f14.7, 2x, "DELTA LN PCP=", f14.7)') itnum, aratio, Aeat(Npt), App(Npt), dlnp
        if (abs(Aeat(Npt) - Aratio) / Aratio <= check) go to 900
        if (abs(dlnp) < 0.00004) go to 900
        aeatl = log(Aeat(Npt))
        itnum = itnum + 1
        if (itnum > 10) then
           write(IOOUT, '(/" WARNING!!  DID NOT CONVERGE FOR AREA RATIO =", F10.5, " (ROCKET)")') aratio
           go to 900
        else
! IMPROVED PCP ESTIMATES.
           asq = Gammas(Npt) * Enn * Rr * Tt
           dlnpe = Gammas(Npt) * usq / (usq - asq)
           go to 850
        end if
     else
        write(IOOUT, '(/" WARNING!!  AREA RATIO CALCULATION CANNOT BE DONE ", &
             & "BECAUSE GAMMAs", /, " CALCULATION IMPOSSIBLE. (ROCKET)")')
        Npt = Npt - 1
        if (Nsub <= 0) isup1 = 100
        if (Nsub < 0.) Nsup = Isup - 1
        if (Nsub > 0) Nsub = isub - 1
        go to 1000
     end if
800  appl = pcplt / (Subar(isub) + (10.587 * eln**2 + 9.454) * eln)
     if (Aratio < 1.09) appl = 0.9 * appl
     if (Aratio > 10) appl = appl / Aratio
     if (isub > 1 .or. Npt == Ncol) go to 1100
     go to 1200
850  dlnp = dlnpe * eln - dlnpe * aeatl
     appl = appl + dlnp
     if (itnum == 1) go to 1100
     if (appl < 0.) appl = 0.000001
     App(Npt) = EXP(appl)
     Pp = Pinf / App(Npt)
     go to 250
! CONVERGENCE HAS BEEN REACHED FOR ASSIGNED AREA RATIO
900  Aeat(Npt) = Aratio
     if (Fac) then
        if (.not. done) then
           if (Iopt == 1) then
! OPTION 1 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
! PRESSURE AND CONTRACTION RATIO. IMPROVED ESTIMATE FOR PC
              Area = .false.
              itnum = 0
              ppa = Ppp(Npt) * pa
              pinj = ppa + 1.d05 * usq / Vlm(Npt)
              test = (pinj - pinjas) / pinjas
              pcpa = pinf * pa
              if (Debugf) then
                 write(IOOUT, '(" ITER", 3x, "TEST", 3x, "ASSIGNED PINJ", 1x, "CALC PINJ", 5x, &
                      & "PC", 7x, "P AT ACAT", 3x, "PREV ACAT", 2x, "ACAT")')
                 write(IOOUT, '(i3, f10.6, 1x, 4f12.2, 2f9.5)') niter, test, pinjas, pinj, pcpa, ppa, acatsv, Acat
              end if
              if (abs(test) < 0.00002) go to 350
              prat = pinjas / pinj
              Pp = pinf * prat
              go to 300
           else if (Iopt == 2) then
! OPTION 2 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
! PRESSURE AND MASS FLOW PER UNIT AREA. IMPROVED ESTIMATE FOR PC
! AND ACAT
              acatsv = Acat
              pratsv = prat
              Area = .false.
              itnum = 0
              ppa = Ppp(4) * pa
              pinj = ppa + 1.d05 * usq / Vlm(4)
              mat = pa / (Awt * Rr)
              Acat = mat / Ma
              prat = (b1 + c1 * Acat) / (1 + a1l * Acat)
              test = (pinj - pinjas) / pinjas
              pcpa = pinf * pa
              if (Debugf) then
                 write(IOOUT, '(" ITER", 3x, "TEST", 3x, "ASSIGNED PINJ", 1x, "CALC PINJ", 5x, &
                      & "PC", 7x, "P AT ACAT", 3x, "PREV ACAT", 2x, "ACAT")')
                 write(IOOUT, '(i3, f10.6, 1x, 4f12.2, 2f9.5)') niter, test, pinjas, pinj, pcpa, ppa, acatsv, Acat
              end if
              if (abs(test) < 0.00002) go to 350
              pjrat = pinj / pinjas
              Pp = pinf
              do i = 1, 2
                 pracat = pratsv / prat
                 pr = pjrat * pracat
                 Pp = Pp / pr
                 pcpa = Pp * pa
                 Acat = Acat / pr
                 Subar(1) = Acat
                 pratsv = prat
                 pjrat = 1
                 prat = (b1 + c1 * Acat) / (1 + a1l * Acat)
                 if (Debugf) write(IOOUT, '(" NEW PC = ", f10.2, 2x, "NEW ACAT = ", f9.6, 2x, "PJRAT =", &
                      & f10.7, " PRACAT =", f10.7)') pcpa, Acat, pjrat, pracat
              end do
              go to 300
           end if
        end if
     end if
950  if (Trnspt) call TRANP
     if (Npt == Nfz) Eql = seql
1000 itnum = 0
     if (Nsub > i01) then
        isub = isub + 1
        if (isub <= Nsub) go to 750
        isub = 1
        Nsub = -Nsub
        if (Isup <= Nsup) go to 750
        Area = .false.
        go to 1150
     end if
1050 Isup = Isup + 1
     itnum = 0
     if (Isup <= Nsup) go to 750
     Isup = isupsv
     Area = .false.
     go to 1150
! TEST FOR OUTPUT -- SCHEDULES COMPLETE OR NPT=Ncol
1100 Isv = Npt
     if (Npt /= Ncol) go to 1200
1150 if (.not. Eql) then
        if (Nfz <= 1) then
           Cpr(Nfz) = cprf
           Gammas(Nfz) = cprf / (cprf - 1 / Wm(Nfz))
        end if
     end if
     call RKTOUT(cea)
     Iplt = Iplt + Npt
     if (.not. Page1) then
        Iplt = Iplt - 2
        if (Iopt /= 0) Iplt = Iplt - 1
        Iplt = min(Iplt, 500)
     else
        Page1 = .false.
     end if
     iplte = max(iplte, Iplt)
     dlnp = 1
     if (Tt == 0) Area = .false.
     if (.not. Eql .and. Tt == 0.) write(IOOUT, '(/" WARNING!!  CALCULATIONS WERE STOPPED BECAUSE NEXT ", &
          & "POINT IS MORE", /, " THAN 50 K BELOW THE TEMPERATURE", &
          & " RANGE OF A CONDENSED SPECIES (ROCKET)")')
     if (Isv == 0) then
! PCP, SUBAR, AND SUPAR SCHEDULES COMPLETED
        if (Nsub < 0) Nsub = -Nsub
        if (.not. Froz .or. .not. Eql) go to 1300
! SET UP FOR FROZEN.
        if (Eql) Iplt = iplt1
        Eql = .false.
        Page1 = .true.
        call SETEN
        Tt = Ttt(Nfz)
        ipp = Nfz
        if (Nfz == Npt) go to 1150
        Npt = Nfz
        Enn = 1./Wm(Nfz)
        if (Nfz == 1) go to 450
        if (Nsub > 0) then
           Nsub = -Nsub
           write(IOOUT, '(/" WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ", /, &
                & " RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)")')
        end if
        if (App(Nfz) < App(nptth)) then
           write(IOOUT, '(/" WARNING!!  FREEZING IS NOT ALLOWED AT A SUBSONIC ", &
                & "PRESSURE RATIO FOR nfz GREATER"/" THAN 1. FROZEN ", &
                & "PERFORMANCE CALCULATIONS WERE OMITTED (ROCKET)")')
        else
           if (Nfz < Npp) go to 1200
           go to 700
        end if
        go to 1300
     else
        if (Eql) write(IOOUT, '(////)')
        Npt = nptth
     end if
! SET INDICES AND ESTIMATES FOR NEXT POINT.
1200 Npt = Npt + 1
     if (Eql .or. (Isv == -i12 .and. .not. seql)) then
! THE FOLLOWING STATEMENT WAS ADDED TO TAKE CARE OF A SITUATION
! WHERE EQLBRM WENT SINGULAR WHEN STARTING FROM ESTIMATES WHERE
! BOTH SOLID AND LIQUID WERE INCLUDED.  JULY 27, 1990.
        if (Jliq /= 0 .and. Isv > 0) Isv = 0
        call SETEN
     end if
1250 ipp = ipp + 1
     if (Npt > nptth) then
        if (Area) then
           App(Npt) = EXP(appl)
        else
           App(Npt) = Pcp(ipp - nptth)
           if (Fac) App(Npt) = App(Npt) * Pinf / Ppp(1)
           if (.not. Eql .and. App(Npt) < App(Nfz)) then
              write(IOOUT, '(/, " WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED", &
                   & " WHERE THE ASSIGNED", /, &
                   & " PRESSURE RATIOS WERE LESS THAN ", &
                   & "THE VALUE AT POINT nfz =", i3, " (ROCKET)")') Nfz
              go to 1250
           end if
        end if
        Pp = pinf / App(Npt)
        if (Fac) then
           if (Area) then
              if (isub <= Nsub .and. isub > i01 .and. aratio >= Aeat(2)) then
                 write(IOOUT, '(/" WARNING!!  ASSIGNED subae/at =", f10.5, " IS NOT ", &
                      & "PERMITTED TO BE GREATER"/" THAN ac/at =", f9.5, &
                      & ".  POINT OMITTED (ROCKET)")') aratio, Aeat(2)
                 Npt = Npt - 1
                 go to 1000
              end if
           else if (Npt > nptth .and. Pcp(ipp-3) < Ppp(1) / Ppp(2)) then
              write(IOOUT, '(/" WARNING!!  ASSIGNED pip =", F10.5, &
                   & " IS NOT PERMITTED"/" TO BE LESS THAN  Pinj/Pc =", f9.5, &
                   & ". POINT OMITTED", " (ROCKET)")') Pcp(ipp-3), Ppp(1) / Ppp(2)
              Npt = Npt - 1
              go to 650
           end if
        end if
     end if
     go to 250
1300 Npt = 1
! CHECK FOR COMPLETED SCHEDULES -
! 1) CHAMBER PRESSURES(IP = NP)
! 2) CHAMBER TEMPERATURES(IT = NT)
! 3) O/F VALUES(IOF = NOF)
     if (Ip == Np .and. It == Nt .and. iof == Nof) go to 1400
     write(IOOUT, '(////)')
     call SETEN
     Tt = Ttt(i12)
  end do
  if (It < Nt) then
     It = It + 1
     Tt = T(It)
     go to 200
  else if (iof < Nof) then
     go to 100
  end if
1400 Iplt = max(Iplt, iplte)
  return
end subroutine ROCKET



subroutine SEARCH
!***********************************************************************
! SEARCH THERMO.LIB FOR THERMO DATA FOR SPECIES TO BE CONSIDERED.
!***********************************************************************
  use mod_cea
  use mod_legacy_cea
  implicit none

  ! LOCAL VARIABLES
  character(6):: date(maxNgc)
  character(2):: el(5)
  character(15):: sub
  integer:: i, j, k, ii
  integer:: i5, ifaz, itot, nall, ne, nint, ntgas, ntot
  real(8):: b(5), t1, t2, thermo(9, 3)

  Nc = 0
  ne = 0
  do i = 1, Nlm
     Jx(i) = 0
  end do
  do j = 1, maxNgc
     S(j) = 0
     H0(j) = 0
     Deln(j) = 0
     do i = 1, Nlm
        A(i, j) = 0
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
  rewind IOTHM
  read(IOTHM) Tg, ntgas, ntot, nall, Thdate
  Ngc = 1
  Nc = 1
! BEGIN LOOP FOR READING SPECIES DATA FROM THERMO.LIB.
  do 200 itot = 1, ntot
     if (itot > ntgas) then
        read(IOTHM) sub, nint, date(Ngc), (el(j), b(j), j = 1, 5), Ifz(Nc), &
             Temp(1, Nc), Temp(2, Nc), Mw(Ngc), (Cft(Nc, k), k = 1, 9)
     else
        read(IOTHM) sub, nint, date(Ngc), (el(j), b(j), j = 1, 5), ifaz, T1, T2, Mw(Ngc), thermo
     end if
     if (Nonly /= 0) then
        i = 1
20      if (Prod(i) /= sub .and. '*' // Prod(i) /= sub) then
           i = i + 1
           if (i <= Nonly) go to 20
           go to 200
        else
           if (sub == Prod(Ngc-1)) then
              Nonly = Nonly + 1
              do k = Nonly, i+1, -1
                 Prod(k) = Prod(k-1)
              end do
           else
              Prod(i) = Prod(Ngc)
           end if
           Prod(Ngc) = sub
        end if
     else if (Nomit /= 0) then
        do i = 1, Nomit
           if (Omit(i) == sub .or. '*' // Omit(i) == sub) go to 200
        end do
     end if
     do 50 k = 1, 5
        if (b(k) == 0) go to 100
        do i = 1, Nlm
           if (Elmt(i) == el(k)) then
              A(i, Ngc) = b(k)
              go to 50
           end if
        end do
        do j = 1, Nlm
           A(j, Ngc) = 0
        end do
        go to 200
50   continue
100  Prod(Ngc) = sub
     if (itot > ntgas) then
        Nc = Nc + 1
        if (Nc > maxNc) go to 400
     else
        Ng = Ngc
        if (Ng > maxNg) go to 400
        do i = 1, 3
           do j = 1, 9
              Coef(Ng, j, i) = thermo(j, i)
           end do
        end do
! IF SPECIES IS AN ATOMIC GAS, STORE INDEX IN JX
        if (b(2) == 0 .and. b(1) == 1) then
           do i = 1, Nlm
              if (Elmt(i) == el(1)) then
                 ne = ne + 1
                 Jx(i) = Ngc
                 Jcm(i) = Ngc
                 go to 150
              end if
           end do
        end if
     end if
150  Ngc = Ngc + 1
     if (Ngc > maxNgc) go to 400
200 continue
! FINISHED READING THERMO DATA FROM I/O UNIT IOTHM.
  Ifz(Nc) = 0
  Nc = Nc - 1
  Ngc = Ngc - 1
  Ngp1 = Ng + 1
  if (Ngc < Nonly) then
     do k = Ngc+1, Nonly
        write(IOOUT, '(/" WARNING!!  ", a15, " NOT A PRODUCT IN thermo.lib FILE (SEARCH)")') Prod(k)
     end do
  end if
! FIND MISSING ELEMENTS (IF ANY) FOR COMPONENTS
  Nspx = Ngc
  if (ne < Nlm) then
     do i = 1, Nlm
        if (Nspx > maxNgc) go to 400
        if (Jx(i) == 0) then
           Nspx = Nspx + 1
           do k = 1, Nlm
              A(k, Nspx) = 0
           end do
           A(i, Nspx) = 1
           Prod(Nspx) = Elmt(i)
           do k = 1, 100
              if (Elmt(i) == atomic_symbol(k)) then
                 Mw(Nspx) = atomic_mass(k)
                 Atwt(i) = atomic_mass(k)
                 Cp(Nspx) = 2.5d0
                 go to 210
              end if
           end do
210        Jx(i) = Nspx
           Jcm(i) = Nspx
        end if
     end do
  end if
! ARE ALL ELEMENTS IN PRODUCT SPECIES?
  do 300 i = 1, Nlm
     do j = 1, Ngc
        if (A(i, j) /= 0) go to 300
        ii = i
     end do
     write(IOOUT, '(/" PRODUCT SPECIES CONTAINING THE ELEMENT", a3, " MISSING", &
          & //, 13x, "FATAL ERROR (SEARCH)")') Elmt(ii)
     Ngc = 0
     return
300 continue
! WRITE POSSIBLE PRODUCT LIST
  if (.not. Short) then
     write(IOOUT, '(/2x, "SPECIES BEING CONSIDERED IN THIS SYSTEM", &
          & /" (CONDENSED PHASE MAY HAVE NAME LISTED SEVERAL TIMES)", &
          & /"  LAST thermo.inp UPDATE: ", a10, /)') Thdate
     do i = 1, Ngc, 3
        i5 = i + 2
        if (Ngc < i5) i5 = Ngc
        write(IOOUT, '(3(2X, A6, 2X, A15))') (date(j), Prod(j), j = i, i5)
     end do
  end if
  return

400 write(IOOUT, '(/" INSUFFICIENT STORAGE FOR PRODUCTS-SEE RP-1311,", /"   PART 2, PAGE 39. (SEARCH)")')
  Ngc = 0
  return
end subroutine


subroutine READTR
  ! SEARCH FOR TRANSPORT PROPERTIES FOR THIS CHEMICAL SYSTEM
  use mod_legacy_cea
  implicit none

  character(16):: bin(2, 40), pure(6), spece(2)
  integer:: i, j, k, jj(2), jk, ir, lineb, npure, nrec
  real(8):: trdata(36)

  rewind IOTRN
  rewind IOSCH
  Ntape = 0
  npure = 0
  lineb = 1
  if (.not. Short) write(IOOUT, '(/" SPECIES WITH TRANSPORT PROPERTIES"//8X, "PURE SPECIES"/)')
  read(IOTRN) nrec
  do ir = 1, nrec
     read(IOTRN) spece, trdata
     k = 1
450  do j = 1, Ng
        if (spece(k) == Prod(j) .or. '*' // spece(k) == Prod(j)) then
           jj(k) = j
           if (k == 2) then
! STORE NAMES FOR BINARIES IN BIN ARRAY.
              do k = 1, 2
                 bin(k, lineb) = spece(k)
              end do
              lineb = lineb + 1
              go to 500
           else
              jj(2) = j
              if (spece(2) == ' ') then
! WRITE NAMES FOR PURE SPECIES.
                 npure = npure + 1
                 pure(npure) = spece(1)
                 go to 500
              else
                 k = 2
                 go to 450
              end if
           end if
        end if
     end do
     go to 550
500  write(IOSCH) jj, trdata
     Ntape = Ntape + 1
550  if (npure /= 0 .and. (npure >= 6 .or. ir >= nrec)) then
        if (.not. Short) write(IOOUT, '(4(2x, A16))') (pure(jk), jk = 1, npure)
        npure = 0
     end if
  end do
  lineb = lineb - 1
  if (.not. Short) then
     write(IOOUT, '(/"     BINARY INTERACTIONS"/)')
     do j = 1, lineb
        write(IOOUT, '(5X, 2A16)') (bin(i, j), i = 1, 2)
     end do
  end if
  write(IOOUT, *)
  return
end subroutine



subroutine SETEN
!***********************************************************************
! USE COMPOSITIONS FROM PREVIOUS POINT AS INITIAL ESTIMATES FOR
! CURRENT POINT NPT.  IF -
!  ISV>0  USE COMPOSITIONS FROM POINT ISV.
!  ISV<0  SAVE COMPOSITIONS FROM POINT -ISV FOR POSSIBLE LATER USE.
!         ALSO USE COMPOSITIONS FROM POINT -ISV FOR NPT.
!  ISV=0  USE COMPOSITIONS SAVED WHEN ISV<0.
!***********************************************************************
  use mod_legacy_cea
  implicit none
! LOCAL VARIABLES
  integer, save:: j, lsav
  real(8), save:: Tsave

  if (Isv < 0) then
! FIRST T--SAVE COMPOSITIONS FOR FUTURE POINTS WITH THIS T
     Isv = -Isv
     Tsave = Ttt(Isv)
     Ensave = Enn
     Enlsav = Ennl
     lsav = Lsave
     do j = 1, Ng
        Sln(j) = Enln(j)
     end do
     do j = 1, Ng
        En(j, Npt) = En(j, Isv)
     end do
     Npr = 0
     do j = Ngp1, Ngc
        Sln(j) = En(j, Isv)
        En(j, Npt) = Sln(j)
        if (Jliq == j) then
           En(Jsol, Npt) = En(Jsol, Isv) + En(Jliq, Isv)
           En(Jliq, Npt) = 0
           Jsol = 0
           Jliq = 0
           Tsave = Tsave - 5
           Tt = Tsave
           Sln(j) = 0
        else if (En(j, Npt) > 0) then
           Npr = Npr + 1
           Jcond(Npr) = j
        end if
     end do
  else if (Isv == 0) then
! NEXT POINT FIRST T IN SCHEDULE, USE PREVIOUS COMPOSITIONS FOR THIS T
     Jsol = 0
     Jliq = 0
     Enn = Ensave
     Ennl = Enlsav
     Lsave = lsav
     Npr = 0
     do j = Ngp1, Ngc
        En(j, Npt) = Sln(j)
        if (En(j, Npt) > 0.d0) then
           Npr = Npr + 1
           Jcond(Npr) = j
        end if
     end do
     do j = 1, Ng
        En(j, Npt) = 0
        Enln(j) = Sln(j)
        if (Sln(j) /= 0) then
           if ((Enln(j) - Ennl + 18.5) > 0) En(j, Npt) = exp(Enln(j))
        end if
     end do
     if (.not. Tp) Tt = Tsave
     Sumn = Enn
  else if (Isv > 0) then
! USE COMPOSITIONS FROM PREVIOUS POINT
     do j = 1, Ngc
        En(j, Npt) = En(j, Isv)
     end do
  end if
end subroutine SETEN



subroutine SHCK(cea)
!***********************************************************************
! PRIMARY ROUTINE FOR SHOCK PROBLEMS.
!***********************************************************************
  use mod_legacy_cea
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

  if (Trace == 0) Trace = 5.E-9
  Tp = .true.
  Cpmix = 0
  srefl = .false.
  if (.not. Short) then
     write(IOOUT, '(/"   *** INPUT FOR SHOCK PROBLEMS ***")')
     write(IOOUT, '(/" INCDEQ =", L2, "   REFLEQ =", L2, "   INCDFZ =", L2, "    REFLFZ =", L2)') &
          Incdeq, Refleq, Incdfz, Reflfz
  end if
  if (Refleq .or. Reflfz) srefl = .true.
  seql = Incdeq
  if (T(1) == 0) T(1) = Rtemp(1)
  do i = 1, Nsk
     uis(i) = U1(i)
     mis(i) = Mach1(i)
     if (Mach1(i) == 0 .and. U1(i) == 0) exit
  end do
  if (Nsk > Ncol) then
     write(IOOUT, '(/" WARNING!!  ONLY ", I2, " u1 OR mach1 VALUES ALLOWED (SHCK)")') Ncol
     Nsk = Ncol
  end if
  if (.not. Short) then
     write(IOOUT, '(/1p, " U1 =   ", 5E13.6, /(8X, 5E13.6))') (U1(i), i = 1, Nsk)
     write(IOOUT, '(/1p, " MACH1 =", 5E13.6, /(8X, 5E13.6))') (Mach1(i), i = 1, Nsk)
  end if
  iof = 0
200 iof = iof + 1
  Oxfl = Oxf(iof)
  call NEWOF
  Incdeq = seql
300 refl = .false.
  it2 = 2
  it1 = 1
  Pp = P(1)
  Tt = T(1)
  if (.not. Incdeq) then
! FROZEN
     do n = 1, Nsk
        Dlvtp(n) = 1
        Dlvpt(n) = -1
     end do
  end if
  do Npt = 1, Nsk
     Ppp(Npt) = P(Npt)
     Ttt(Npt) = T(Npt)
     if (Npt > 1) then
        if (Ppp(Npt) == 0) Ppp(Npt) = Ppp(Npt-1)
        if (Ttt(Npt) == 0) Ttt(Npt) = Ttt(Npt-1)
        Ssum(Npt) = Ssum(Npt-1)
        Hsum(Npt) = Hsum(Npt-1)
        if (Ttt(Npt) == Tt .and. Ppp(Npt) == Pp) go to 350
     end if
     Pp = Ppp(Npt)
     Tt = Ttt(Npt)
     if (Tt >= Tg(1) * 0.8d0) then
        call HCALC
        Hsum(Npt) = Hsub0
     else
        write(IOOUT, '(/" TEMPERATURE=", E12.4, " IS OUT OF EXTENDED RANGE ", "FOR POINT", I5, " (SHCK)")') Tt, Npt
        go to 1000
     end if
350  if (Cpmix /= 0) Gamma1 = Cpmix / (Cpmix - 1/Wmix)
     A1 = sqrt(Rr * Gamma1 * Tt / Wmix)
     if (U1(Npt) == 0) U1(Npt) = A1 * Mach1(Npt)
     if (Mach1(Npt) == 0.) Mach1(Npt) = U1(Npt) / A1
     Wm(Npt) = Wmix
     Cpr(Npt) = Cpmix
     Gammas(Npt) = Gamma1
     Vlm(Npt) = Rr * Tt / (Wmix * Pp)
  end do
  Npt = Nsk
! OUTPUT--1ST CONDITION
  write(IOOUT, '(////25X, "SHOCK WAVE PARAMETERS ASSUMING")')
  if (.not. Incdeq) then
     write(IOOUT, '(/, 17X, " FROZEN COMPOSITION FOR INCIDENT SHOCKED CONDITI1ONS"//)')
  else
     write(IOOUT, '(/, 16X, " EQUILIBRIUM COMPOSITION FOR INCIDENT SHOCKED CONDITIONS"//)')
  end if
  Eql = .false.
  call OUT1
  write(IOOUT, '(/" INITIAL GAS (1)")')
  cea%fmt(4) = '13'
  cea%fmt(5) = ' '
  cea%fmt(7) = '4,'
  write(IOOUT, cea%fmt) 'MACH NUMBER1   ', (Mach1(j), j = 1, Npt)
  cea%fmt(7) = '2,'
  write(IOOUT, cea%fmt) 'U1, M/SEC      ', (U1(j), j = 1, Npt)
  call OUT2(cea)
! BEGIN CALCULATIONS FOR 2ND CONDITION
  if (Incdeq) Eql = .true.
  Npt = 1
400 Gamma1 = Gammas(Npt)
  uu = U1(Npt)
  wmx = Wm(Npt)
  p1 = Ppp(Npt)
  T1 = Ttt(Npt)
  hs = Hsum(Npt)
  if (refl) uu = u1u2(Npt)
  mu12rt = wmx * uu**2 / (Rr * T1)
  if (refl) then
! REFLECTED--SUBSCRIPTS 2=1, 5=2, P52=P21
     T21 = 2
     b2 = (-1 - mu12rt - T21) / 2
     p21 = -b2 + sqrt(b2**2 - T21)
  else
     p21 = (2 * Gamma1 * Mach1(Npt)**2 - Gamma1 + 1) / (Gamma1 + 1)
! THE FOLLOWING IMPROVED FORMULATION FOR THE INITIAL ESTIMATE FOR THE
! 2ND CONDITION WAS MADE AND TESTED BY S. GORDON 7/10/89.
     if (.not. Eql) then
        T21 = p21 * (2 / Mach1(Npt)**2 + Gamma1 - 1) / (Gamma1 + 1)
     else
        Pp = p21 * p1
        Tp = .false.
        Hp = .true.
        Hsub0 = hs + uu**2 / (2 * Rr)
        call EQLBRM
        T21 = Ttt(Npt) / T1
        Hp = .false.
        Tp = .true.
     end if
  end if
  p21l = log(p21)
  ttmax = 1.05 * Tg(4) / T1
  T21 = min(T21, ttmax)
  t21l = log(T21)
  itr = 1
500 if (Shkdbg) write(IOOUT, '(/" ITR NO.=", I3, 3X, "P", I1, "/P", I1, " =", F9.4, 3X, "T", I1, &
       & "/T", I1, " =", F9.4, "   RHO2/RHO1 =", F9.6)') itr, it2, it1, p21, it2, it1, T21, rho52
  Tt = T21 * T1
  Pp = p21 * p1
  if (.not. Eql) then
! FROZEN
     Tln = log(Tt)
     if (.not. Incdeq) then
        call HCALC
        if (Tt == 0) go to 600
        Hsum(Npt) = Hsub0
        Cpr(Npt) = Cpmix
     else
        call CPHS
        Cpr(Npt) = Cpsum
        Hsum(Npt) = 0
        do j = 1, Ng
           Hsum(Npt) = Hsum(Npt) + H0(j) * En(j, Npt)
        end do
        Hsum(Npt) = Hsum(Npt) * Tt
     end if
  else
     call EQLBRM
     if (Tt == 0) go to 800
  end if
  rho12 = wmx * T21 / (Wm(Npt) * p21)
  gg = rho12 * mu12rt
  rho52 = 1 / rho12
  if (refl) gg = -mu12rt * rho52 / (rho52 - 1)**2
  G(1, 1) = -gg * Dlvpt(Npt) - p21
  G(1, 2) = -gg * Dlvtp(Npt)
  G(1, 3) = p21 - 1 + gg - mu12rt
  if (refl) G(1, 3) = p21 - 1 + gg * (rho52 - 1)
  gg = gg * T1 / wmx
  if (.not. refl) gg = gg * rho12
  G(2, 1) = -gg * Dlvpt(Npt) + Tt * (Dlvtp(Npt) - 1) / Wm(Npt)
  G(2, 2) = -gg * Dlvtp(Npt) - Tt * Cpr(Npt)
  gg = 1 - rho12**2
  if (refl) gg = (rho52 + 1) / (rho52 - 1)
  G(2, 3) = Hsum(Npt) - hs - uu**2 * gg / (2 * Rr)
  X(3) = G(1, 1) * G(2, 2) - G(1, 2) * G(2, 1)
  X(1) = (G(1, 3) * G(2, 2) - G(2, 3) * G(1, 2)) / X(3)
  X(2) = (G(1, 1) * G(2, 3) - G(2, 1) * G(1, 3)) / X(3)
  if (Shkdbg) then
     write(IOOUT, '(/" G(I,J)  ", 3E15.8)') G(1, 1), G(1, 2), G(1, 3)
     write(IOOUT, '(/" G(I,J)  ", 3E15.8)') G(2, 1), G(2, 2), G(2, 3)
     write(IOOUT, '(/" X       ", 2E15.8)') X(1), X(2)
     write(IOOUT, '(/" HSUM HS UU U2 ", 4E15.8)') Hsum(Npt), hs, uu, uu * rho12
  end if
  ax = abs(X(1))
  axx = abs(X(2))
  if (axx > ax) ax = axx
  if (ax >= 0.00005) then
     cormax = 0.40546511
     if (itr > 4) cormax = 0.22314355
     if (itr > 12) cormax = 0.09531018
     if (itr > 20) cormax = 0.04879016
     ax = ax / cormax
     if (ax > 1) then
        X(1) = X(1) / ax
        X(2) = X(2) / ax
     end if
     p21l = p21l + X(1)
     t21l = t21l + X(2)
     p21 = exp(p21l)
     T21 = exp(t21l)
     if (Shkdbg) write(IOOUT, '(/" MAX.COR.=", e13.6, " X(1)=", e13.6, " X(2)=", e13.6)') cormax, X(1), X(2)
     if (itr /= 1 .or. T21 < ttmax) then
        itr = itr + 1
        if (itr < 61) go to 500
        write(IOOUT, '(/6x, " WARNING!!  NO CONVERGENCE FOR u1=", F8.1, &
             & /"  ANSWERS NOT RELIABLE, SOLUTION MAY NOT EXIST (SHCK)")') U1(Npt)
     else
        Tt = 0
        Npt = Npt - 1
        go to 700
     end if
  end if
! CONVERGED OR TOOK 60 ITERATIONS WITHOUT CONVERGING.
! STORE RESULTS.
600 rrho(Npt) = rho52
  m2m1(Npt) = Wm(Npt) / wmx
  p2p1(Npt) = p21
  t2t1(Npt) = T21
  utwo(Npt) = uu * rho12
  u1u2(Npt) = uu - utwo(Npt)
  if (Tt >= Tg(1) * 0.8d0 .and. Tt <= Tg(4) * 1.1d0) then
     if (.not. Eql) then
! FROZEN
        Ppp(Npt) = Pp
        Ttt(Npt) = Tt
        Gammas(Npt) = Cpr(Npt) / (Cpr(Npt) - 1 / wmx)
        Vlm(Npt) = Rr * Tt / (wmx * Pp)
        if (Incdeq) then
           Ssum(Npt) = 0
           do j = 1, Ngc
              pmn = Pp * wmx * En(j, Npt)
              if (En(j, Npt) > 0) Ssum(Npt) = Ssum(Npt) + En(j, Npt) * (S(j) - log(pmn))
           end do
        end if
     end if
     go to 900
  end if
700 write(IOOUT, '(/" TEMPERATURE=", E12.4, " IS OUT OF EXTENDED RANGE ", &
       & "FOR POINT", I5, " (SHCK)")') Tt, Npt
  Tt = 0
800 if (Npt < 1) go to 1000
  Nsk = Npt
900 if (Trnspt) call TRANP
  Isv = 0
  if (Npt < Nsk) Isv = Npt
  if (Npt == 1) Isv = -1
  Npt = Npt + 1
  if (Eql) call SETEN
  if (Npt <= Nsk) go to 400
  Npt = Nsk
  if (refl) then
     if (.not. Eql) write(IOOUT, '(/" SHOCKED GAS (5)--REFLECTED--FROZEN")')
     if (Eql) write(IOOUT, '(/" SHOCKED GAS (5)--REFLECTED--EQUILIBRIUM")')
     cr12 = '2'
     cr52 = '5'
  else
     if (.not. Eql) write(IOOUT, '(/" SHOCKED GAS (2)--INCIDENT--FROZEN")')
     if (Eql) write(IOOUT, '(/" SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM")')
     cr12 = '1'
     cr52 = '2'
  end if
  cea%fmt(7) = '2,'
  write(IOOUT, cea%fmt) 'U' // cr52 // ', M/SEC      ', (utwo(j), j = 1, Npt)
  call OUT2(cea)
  if (Trnspt) call OUT4(cea)
  write(IOOUT, *)
  cea%fmt(7) = '3,'
  write(IOOUT, cea%fmt) 'P' // cr52 // '/P' // cr12 // '           ', (p2p1(j), j = 1, Npt)
  write(IOOUT, cea%fmt) 'T' // cr52 // '/T' // cr12 // '           ', (t2t1(j), j = 1, Npt)
  cea%fmt(7) = '4,'
  write(IOOUT, cea%fmt) 'M' // cr52 // '/M' // cr12 // '           ', (m2m1(j), j = 1, Npt)
  write(IOOUT, cea%fmt) 'RHO' // cr52 // '/RHO' // cr12 // '       ', (rrho(j), j = 1, Npt)
  cea%fmt(7) = '2,'
  if (.not. refl) write(IOOUT, cea%fmt) 'V2, M/SEC      ', (u1u2(j), j = 1, Npt)
  if (refl) write(IOOUT, cea%fmt) 'U5+V2,M/SEC    ', (u1u2(j), j = 1, Npt)
  if (.not. Eql) then
! WRITE FROZEN MOLE (OR MASS) FRACTIONS
     cea%fmt(7) = '5,'
     if (.not. Incdeq) then
        if (Massf) then
           write(IOOUT, '(/1x, A4, " FRACTIONS"/)') 'MASS'
        else
           write(IOOUT, '(/1x, A4, " FRACTIONS"/)') 'MOLE'
           ww = wmx
        end if
        do n = 1, Nreac
           j = Jray(n)
           if (Massf) ww = Mw(j)
           write(IOOUT, '(" ", A16, F8.5, 12F9.5)') Prod(j), (En(j, i) * ww, i = 1, Npt)
        end do
     else
        Eql = .true.
        call OUT3(cea)
        Eql = .false.
     end if
  else
     call OUT3(cea)
  end if
  Iplt = min(Iplt + Npt, 500)
  if (srefl) then
     if (.not. refl) then
        refl = .true.
        it2 = 5
        it1 = 2
        Eql = .true.
        if (Reflfz) then
           Eql = .false.
           if (Refleq) then
              j = 0
              do i = 1, Npt
                 j = j + 1
                 sg(j) = u1u2(i)
                 j = j + 1
                 sg(j) = Wm(i)
                 j = j + 1
                 sg(j) = Ppp(i)
                 j = j + 1
                 sg(j) = Ttt(i)
                 j = j + 1
                 sg(j) = Hsum(i)
                 j = j + 1
                 sg(j) = Gammas(i)
              end do
           end if
        end if
        Npt = 1
        go to 400
     else if (.not. Eql .and. Refleq) then
        j = 1
        do i = 1, Npt
           u1u2(i) = sg(j)
           Wm(i) = sg(j+1)
           Ppp(i) = sg(j+2)
           Ttt(i) = sg(j+3)
           Hsum(i) = sg(j+4)
           Gammas(i) = sg(j+5)
           j = j + 6
        end do
        Eql = .true.
        Npt = 1
        go to 400
     end if
  end if
  if (Incdeq .and. Incdfz) then
     Incdeq = .false.
     Eql = .false.
     go to 300
  else if (iof >= Nof) then
     Tp = .false.
     do n = 1, Nreac
        Rtemp(n) = T(1)
     end do
  else
     do i = 1, Nsk
        U1(i) = uis(i)
        Mach1(i) = mis(i)
     end do
     go to 200
  end if
1000 return
end subroutine SHCK



subroutine THERMP(cea)
!***********************************************************************
! ASSIGNED THERMODYNAMIC STATES.  HP, SP, TP, UV, SV, AND TV PROBLEMS.
!***********************************************************************
  use mod_legacy_cea
  use mod_legacy_io
  implicit none

  type(CEA_Problem), intent(inout):: cea

! LOCAL VARIABLES
  integer:: iof
  logical:: Uv, Tv, Sv

  Uv = transfer(Hp, Uv)
  Tv = transfer(Tp, Tv)
  Sv = transfer(Sp, Sv)
  Eql = .true.
  outerLoop: do iof = 1, Nof
     Oxfl = Oxf(iof)
     call NEWOF
! SET ASSIGNED P OR VOLUME
     do Ip = 1, Np
        Pp = P(Ip)
! SET ASSIGNED T
        do It = 1, Nt
           Vv = V(Ip)
           Tt = T(It)
           call EQLBRM
           if (Npt == 0) return
           if (Trnspt .and. Tt /= 0) call TRANP
           Isv = 0
           if (Ip /= Np .or. It /= Nt .and. Tt /= 0) then
              Isv = Npt
              if (Npt /= Ncol) go to 10
           end if
           if (.not. Hp) write(IOOUT, '(////15X, "THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED")')
           if (Hp) write(IOOUT, '(////9X, "THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED")')
           if (.not. Vol) then
              if (Hp) write(IOOUT, '(/34X, " PRESSURES"/)')
              if (Tp) write(IOOUT, '(/27X, "TEMPERATURE AND PRESSURE"/)')
              if (Sp) write(IOOUT, '(/29X, "ENTROPY AND PRESSURE"/)')
           else
              if (Uv) write(IOOUT, '(/36X, " VOLUME"/)')
              if (Tv) write(IOOUT, '(/28X, "TEMPERATURE AND VOLUME"/)')
              if (Sv) write(IOOUT, '(/30X, "ENTROPY AND VOLUME"/)')
           end if
           call OUT1
           write(IOOUT, '(/" THERMODYNAMIC PROPERTIES"/)')
           call OUT2(cea)
           if (Trnspt) call OUT4(cea)
           call OUT3(cea)
           Iplt = min(Iplt + Npt, 500)
           if (Isv == 0 .and. iof == Nof) return
           write(IOOUT, '(////)')
           Npt = 0
10         Npt = Npt + 1
           if (.not. Tp .and. Tt /= 0) T(1) = Tt
           if (Nt == 1 .and. Np == 1) cycle outerLoop
           if (Ip == 1 .and. It == 1) Isv = -Isv
           if (Nt /= 1) then
              if (It == Nt .or. Tt == 0.) Isv = 0
           end if
           call SETEN
        end do
     end do
  end do outerLoop

     return
end subroutine THERMP



subroutine TRANIN
!***********************************************************************
! BRINGS IN AND SORTS OUT INPUT FOR TRANSPORT CALCULATIONS
!***********************************************************************
  use mod_legacy_cea
  implicit none
! LOCAL VARIABLES
  integer:: i, ii, inds(maxTr), ir, j, jtape(2), k, k1, k2, kt, kvc, l, loop, m, nms
  logical:: change, elc1, elc2, ion1, ion2, setx
  real(8):: coeff, debye, ekt, enel, enmin, ionic, lambda, omega, prop, qc, ratio, &
       stcf(maxTr, maxTr), stcoef(maxTr), te, testen, testot, total, &
       trc(6, 3, 2), wmols(maxTr), wmred, xsel, xss(maxTr)

  if (.not. Eql) then
     if (.not. Shock) then
        if (.not. setx) then
           setx = .true.
           Nm = nms
           do i = 1, Nm
              Xs(i) = xss(i)
              Wmol(i) = wmols(i)
              Ind(i) = inds(i)
           end do
        end if
        go to 300
     else if (.not. Incdeq) then
        if (Npt <= 1) then
           Nm = Nreac
           do i = 1, Nm
              j = Jray(i)
              Ind(i) = j
              Wmol(i) = Mw(j)
              Xs(i) = En(j, 1)*Wm(1)
           end do
        end if
        go to 300
     end if
  end if

! PICK OUT IMPORTANT SPECIES
  Nm = 0
  total = 0
  enmin = 1.0d-11 / Wm(Npt)
  testot = 0.999999999d0 / Wm(Npt)
  do i = 1, Lsave
     j = Jcm(i)
     if (En(j, Npt) <= 0 .and. j <= Ngc) then
        if ((Enln(j) - Ennl + 25.328436d0) > 0) En(j, Npt) = exp(Enln(j))
     end if
     Nm = Nm + 1
     Ind(Nm) = j
     total = total + En(j, Npt)
     if (Mw(j) < 1) enel = En(j, Npt)
     En(j, Npt) = -En(j, Npt)
  end do
  testen = 1 / (Ng * Wm(Npt))

  if (total <= testot) then
     outerLoop1: do loop = 1, Ng
        testen = testen / 10

        do j = 1, Ng
           if (En(j, Npt) >= testen) then
              if (Nm >= maxTr) then
                 write(IOOUT, '(/" WARNING!!  MAXIMUM ALLOWED NO. OF SPECIES", I3, " WAS USED IN ", &
                      & /" TRANSPORT PROPERTY CALCULATIONS FOR POINT", I3, "(TRANIN))")') Nm, Npt
                 exit outerLoop1
              else
                 total = total + En(j, Npt)
                 Nm = Nm + 1
                 Ind(Nm) = j
                 En(j, Npt) = -En(j, Npt)
              end if
           end if
        end do

        if (testen <= enmin) exit
     end do outerLoop1
  end if

! CALCULATE MOLE FRACTIONS FROM THE EN(J, NPT)
  do j = 1, Ng
     En(j, Npt) = abs(En(j, Npt))
  end do
  do i = 1, Nm
     j = Ind(i)
     Wmol(i) = Mw(j)
     Xs(i) = En(j, Npt) / total
  end do
  if (Npt == Nfz) then
     nms = Nm
     do i = 1, Nm
        xss(i) = Xs(i)
        wmols(i) = Wmol(i)
        inds(i) = Ind(i)
     end do
     setx = .false.
  end if
! REWRITE REACTIONS TO ELIMINATE TRACE ELEMENTS
  Nr = Nm - Lsave
  if (Nr /= 0) then
     do k = 1, maxTr
        do m = 1, maxTr
           Stc(k, m) = 0
        end do
     end do
     k = 1
     do i = Lsave + 1, Nm
        Stc(k, i) = -1
        j = Ind(i)
        do m = 1, Lsave
           Stc(k, m) = A(m, j)
        end do
        k = k + 1
     end do
     do i = 1, Nm
        if (Xs(i) < 1.0d-10) then
           m = 1
           change = .false.

           do j = 1, Nr
              coeff = Stc(j, i)
              if (abs(coeff) > 1.0d-05) then
                 if (.not. change) then
                    change = .true.
                    do k = 1, Nm
                       stcoef(k) = Stc(j, k) / coeff
                    end do
                    cycle
                 else
                    do k = 1, Nm
                       Stc(j, k) = Stc(j, k) / coeff - stcoef(k)
                    end do
                 end if
              end if
              do k = 1, Nm
                 stcf(m, k) = Stc(j, k)
              end do
              m = m + 1
           end do

           do ii = 1, Nm
              do j = 1, Nr
                 Stc(j, ii) = stcf(j, ii)
              end do
           end do
           Nr = m - 1
        end if
     end do
  end if

! FIND TRANSPORT DATA FOR IMPORTANT INTERACTIONS
300 do i = 1, Nm
     Con(i) = 0
     do j = 1, Nm
        Eta(i, j) = 0
     end do
  end do

  rewind IOSCH

  outerLoop2: do ir = 1, Ntape
     read(IOSCH) jtape, trc

     innerLoop: do k = 1, 2
        do i = 1, Nm
           j = Ind(i)
           if (j == jtape(k)) then
              l = i
              if (k == 2) then
                 kvc = 1

                 do
                    kt = 1
                    if (trc(2, 1, kvc) /= 0) then
                       if (trc(2, 2, kvc) /= 0) then
                          if (Tt > trc(2, 1, kvc)) kt = 2
                          if (trc(2, 3, kvc) /= 0) then
                             if (Tt > trc(2, 2, kvc)) kt = 3
                          end if
                       end if

                       prop = exp(trc(6, kt, kvc) + (trc(5, kt, kvc) / Tt + trc(4, kt, kvc)) / Tt + trc(3, kt, kvc) * Tln)

                       if (kvc == 2) then
                          Con(l) = prop
                          cycle outerLoop2
                       else
                          Eta(l, m) = prop
                          if (l /= m) Eta(m, l) = Eta(l, m)
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
  if (Ions) then
     te = Tt / 1000
     ekt = 4.8032d0**2 / (Boltz * te)
     qc = 100 * ekt**2
     xsel = enel / total
     if (xsel < 1.0d-12) xsel = 1.0d-12
     debye = ((22.5d0 / Pi) * (Rr / Avgdr * 100) * (te/xsel)) / ekt**3
     ionic = ((810 / (4*Pi)) * (Rr / Avgdr * 100d0) * (te/xsel))**(2/3.) / ekt**2
     lambda = sqrt(debye + ionic)
     lambda = max(lambda, 2.71828183d0)
  end if
  do i = 1, Nm
     k = Ind(i)
     Cprr(i) = Cp(k)
     if (.not. (Ions .and. (abs(A(Nlm, k)) == 1) .and. (Eta(i, i) == 0))) then
        if (Eta(i, i) == 0) then
           omega = log(50 * Wmol(i)**4.6 / Tt**1.4)
           omega = max(omega, 1.)
           Eta(i, i) = Viscns * sqrt(Wmol(i) * Tt) / omega
        end if
        if (Con(i) == 0) Con(i) = Eta(i, i) * Rr * (0.00375d0 + 0.00132d0 * (Cprr(i) - 2.5d0)) / Wmol(i)
     end if
  end do

  do i = 1, Nm
     do j = i, Nm
        ion1 = .false.
        ion2 = .false.
        elc1 = .false.
        elc2 = .false.
        omega = 0
        if (Eta(i, j) == 0) Eta(i, j) = Eta(j, i)
        if (Eta(j, i) == 0) Eta(j, i) = Eta(i, j)
        if (Eta(i, j) == 0) then
           if (Ions) then
! ESTIMATE FOR IONS
              k1 = Ind(i)
              k2 = Ind(j)
              if (abs(A(Nlm, k1)) == 1) ion1 = .true.
              if (abs(A(Nlm, k2)) == 1) ion2 = .true.
              if (Wmol(i) < 1) elc1 = .true.
              if (Wmol(j) < 1) elc2 = .true.
              if (ion1 .and. ion2) omega = 1.36d0 * qc * log(lambda)
              if ((ion1 .and. elc2) .or. (ion2 .and. elc1)) omega = 1.29d0 * qc * log(lambda)
              if ((ion1 .and. .not. ion2) .or. (ion2 .and. .not. ion1)) omega = exp(6.776 - 0.4 * Tln)
              if (omega /= 0) then
                 wmred = sqrt(2 * Tt * Wmol(i) * Wmol(j) / (Wmol(i) + Wmol(j)))
                 Eta(i, j) = Viscns * wmred * Pi / omega
                 Eta(j, i) = Eta(i, j)
                 if (i == j) then
                    Cprr(i) = Cp(k1)
                    Con(i) = Eta(i, i) * Rr * (0.00375d0 + 0.00132d0 * (Cprr(i) - 2.5d0)) / Wmol(i)
                 end if
                 cycle
              end if
           end if
! ESTIMATE FOR UNLIKE INTERACTIONS FROM RIGID SPHERE ANALOGY
           ratio = sqrt(Wmol(j) / Wmol(i))
           Eta(i, j) = 5.656854d0 * Eta(i, i) * sqrt(Wmol(j) / (Wmol(i) + Wmol(j)))
           Eta(i, j) = Eta(i, j) / (1 + sqrt(ratio * Eta(i, i) / Eta(j, j)))**2
           Eta(j, i) = Eta(i, j)
        end if
     end do
  end do

  return
end subroutine TRANIN



subroutine TRANP
!***********************************************************************
! CALCULATES GAS TRANSPORT PROPERTIES
!
!   NUMBER OF GASEOUS SPECIES = NM   (MAXIMUM maxTr)
!   NUMBER OF CHEMICAL REACTIONS = NR (NM - NLM)
!   ARRAY OF STOICHIOMETRIC COEFFICIENTS = STC
!***********************************************************************
  use mod_legacy_cea
  implicit none
! LOCAL VARIABLES
  integer:: i, i1, j, jj, k, m, mm, nlmm, nmm
  real(8):: cpreac, delh(maxTr), gmat(maxMat, maxMat+1), phi(maxTr, maxTr), &
       psi(maxTr, maxTr), reacon, rtpd(maxTr, maxTr), stx(maxTr), &
       stxij(maxTr, maxTr), sumc, sumv, wtmol, xskm(maxTr, maxTr)

  call TRANIN
! CALCULATE VISCOSITY AND FROZEN THERMAL CONDUCTIVITY
  nmm = Nm - 1
  do i = 1, Nm
     rtpd(i, i) = 0
     phi(i, i) = 1
     psi(i, i) = 1
  end do
  Confro(Npt) = 0
  Vis(Npt) = 0
  do i = 1, nmm
     i1 = i + 1
!DIR$ IVDEP
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
     Vis(Npt) = Vis(Npt) + Eta(i, i) * Xs(i) / sumv
     Confro(Npt) = Confro(Npt) + Con(i) * Xs(i) / sumc
  end do
  if (Eql .and. Nr > 0) then
! CALCULATE REACTION HEAT CAPACITY AND THERMAL CONDUCTIVITY
     m = Nr + 1
     do i = 1, Nr
        delh(i) = 0
        do k = 1, Lsave
           j = Jcm(k)
           delh(i) = Stc(i, k) * H0(j) + delh(i)
        end do
        nlmm = Lsave + 1
        do k = nlmm, Nm
           j = Ind(k)
           delh(i) = Stc(i, k) * H0(j) + delh(i)
        end do
        G(i, m) = delh(i)
     end do
     do i = 1, maxTr
        do j = 1, maxTr
           if (abs(Stc(i, j)) < 1.0d-6) Stc(i, j) = 0
        end do
     end do
     jj = Nm - 1
     do k = 1, jj
        mm = k + 1
        do m = mm, Nm
           rtpd(k, m) = Wmol(k) * Wmol(m) / (1.1 * Eta(k, m) * (Wmol(k) + Wmol(m)))
           xskm(k, m) = Xs(k) * Xs(m)
           xskm(m, k) = xskm(k, m)
           rtpd(m, k) = rtpd(k, m)
        end do
     end do
     do i = 1, Nr
        do j = i, Nr
           G(i, j) = 0
           gmat(i, j) = 0
        end do
     end do
     do k = 1, jj
        mm = k + 1
        do m = mm, Nm
           if (Xs(k) >= 1.0d-10 .and. Xs(m) >= 1.0d-10) then
              do j = 1, Nr
                 if ((Stc(j, k) == 0) .and. (Stc(j, m) == 0)) stx(j) = 0
                 if ((Stc(j, k) /= 0) .or. (Stc(j, m) /= 0)) stx(j) = Xs(m) * Stc(j, k) - Xs(k) * Stc(j, m)
              end do
              do i = 1, Nr
                 do j = i, Nr
                    stxij(i, j) = stx(i) * stx(j) / xskm(k, m)
                    G(i, j) = G(i, j) + stxij(i, j)
                    gmat(i, j) = gmat(i, j) + stxij(i, j) * rtpd(k, m)
                 end do
              end do
           end if
        end do
     end do
     m = 1 + Nr
     do i = 1, Nr
!DIR$ IVDEP
        do j = i, Nr
           G(j, i) = G(i, j)
        end do
        G(i, m) = delh(i)
     end do
     imat = Nr
     call GAUSS
     cpreac = 0
     do i = 1, Nr
        G(i, m) = delh(i)
        cpreac = cpreac + R * delh(i) * X(i)
!DIR$ IVDEP
        do j = i, Nr
           G(i, j) = gmat(i, j)
           G(j, i) = G(i, j)
        end do
     end do
     call GAUSS
     reacon = 0
     do i = 1, Nr
        reacon = reacon + R*delh(i)*X(i)
     end do
     reacon = 0.6d0 * reacon
  else
     cpreac = 0
     reacon = 0
  end if
! CALCULATE OTHER ANSWERS
  Cpfro(Npt) = 0
  wtmol = 0
  do i = 1, Nm
     Cpfro(Npt) = Cpfro(Npt) + Xs(i) * Cprr(i)
     wtmol = wtmol + Xs(i) * Wmol(i)
  end do
  Cpfro(Npt) = Cpfro(Npt) * R / wtmol
  Confro(Npt) = Confro(Npt) / 1000
  if (.not. SIunit) Confro(Npt) = Confro(Npt) / 4.184d0
  Vis(Npt) = Vis(Npt) / 1000
  Prfro(Npt) = Vis(Npt) * Cpfro(Npt) / Confro(Npt)
  if (Eql) then
     cpreac = cpreac / wtmol
     reacon = reacon / 1000
     Cpeql(Npt) = cpreac + Cpfro(Npt)
     Coneql(Npt) = Confro(Npt) + reacon
     Preql(Npt) = Vis(Npt) * Cpeql(Npt) / Coneql(Npt)
  end if
end subroutine TRANP



subroutine UTHERM(readOK)
!***********************************************************************
! READ THERMO DATA FROM I/O UNIT 7 IN RECORD format AND WRITE
! UNFORMATTED ON I/O UNIT IOTHM.  DATA ARE REORDERED GASES FIRST.
!
! USES SCRATCH I/O UNIT IOSCH.
!
! UTHERM IS CALLED FROM subroutine INPUT.
!
! NOTE:  THIS ROUTINE MAY BE CALLED DIRECTLY AND USED BY ITSELF TO
! PROCESS THE THERMO DATA.
!
! GASEOUS SPECIES:
! THE STANDARD TEMPERATURE RANGES TGL ARE GIVEN ON THE FIRST
! RECORD, FOLLOWED BY THE DATE OF THE LAST DATA CHANGE THDATE.
!
! WHEN COEFFICIENTS ARE NOT GIVEN FOR THE THIRD TEMPERATURE
! INTERVAL, A STRAIGHT LINE FOR CP/R IS USED.  FOR HIGH TEMPS,
! THE EXTRAPOLATION GOES BETWEEN THE LAST POINT GIVEN AND THE
! FOLLOWING VALUES AA AT TINF = 1.d06 K:
!      MONATOMICS  2.5
!      DIATOMICS   4.5
!      POLYATOMICS 3*N-1.75  (AVERAGE 1.5 AND 2)
!
! THE FOLLOWING EXTRAPOLATION IS NOT CURRENTLY PROGRAMED (12/9/98):
!   FOR LOW TEMPS, THE EXTRAPOLATION GOES BETWEEN THE FIRST VALUE
!   DOWN TO THE FOLLOWING VALUES AA AT 0 K:
!      MONATOMICS  2.5
!      DIATOMICS   3.5
!      POLYATOMICS 3.75 (AVERAGE 3.5 AND 4.0)
!
! IF DATA ARE AVAILABLE FOR THE THIRD T INTERVAL, IFAZ (SEE
! DEFINITION) IS SET TO -1 AND THE NAME IS ALTERED TO START WITH *.
!
! CONDENSED SPECIES:
! NO EXTRAPOLATIONS ARE DONE.  TEMP INTERVALS VARY.
!
! SOME DEFINITIONS:
! TGL(I)  - TEMPERATURE INTERVALS FOR GASES (I.E. 200, 1000, 6000, 20000).
! FILL(I) - IF TRUE, DATA MISSING FOR INTERVAL.  CURRENTLY ONLY 3RD
!           INTERVAL CHECKED.
! NGL     - NUMBER OF GASEOUS PRODUCTS.
! NS      - NGL + NUMBER OF CONDENSED PRODUCT PHASES.
! NALL    - NS + NUMBER OF REACTANT SPECIES.
! IFAZ    - PHASE INDICATOR. GASES ARE 0, CONDENSED PHASES ARE NUMBERED
!           STARTING WITH 1 FOR THE LOWEST T RANGE, 2 FOR THE NEXT
!           CONTIGUOUS PHASE, ETC.
! NTL     - NUMBER OF T INTERVALS FOR A SPECIES SET.
!***********************************************************************
  use mod_legacy_cea
  implicit none
! DUMMY ARGUMENTS
  logical:: readOK
! LOCAL VARIABLES
  character(15):: name
  character(16):: namee
  character(65):: notes
  character(2):: sym(5)
  character(6):: date
  integer:: i, ifaz, ifzm1, inew, int, j, k, kk, l, nall, ncoef, ngl, ns, ntl
  logical:: fill(3)
  real(8):: aa, atms, cpfix, dlt, expn(8), fno(5), hform, hh, mwt, templ(9), tex, &
       tgl(4), thermo(9, 3), tinf, tl(2), ttl, tx

  ngl = 0
  ns = 0
  nall = 0
  ifzm1 = 0
  inew = 0
  tinf = 1.d06
  rewind IOSCH
  read(IOINP, '(4f10.3, a10)') tgl, Thdate

  do
     do i = 1, 3
        fill(i) = .true.
        do j = 1, 9
           thermo(j, i) = 0
        end do
     end do
     hform = 0
     tl(1) = 0
     tl(2) = 0
     read(IOINP, '(a15, a65)', END=300, ERR=400) name, notes
     if (name(:3) == 'END' .or. name(:3) == 'end') then
        if (index(name, 'ROD') == 0 .and. index(name, 'rod') == 0) exit
        ns = nall
        cycle
     end if
     read(IOINP, '(i2, 1x, a6, 1x, 5(a2, f6.2), i2, f13.5, f15.3)', ERR=400) &
          ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, mwt, hform
     write(IOOUT, '(" ", a15, 2x, a6, e15.6, 2x, a65)') name, date, hform, notes
     ! IF NTL=0, REACTANT WITHOUT COEFFICIENTS
     if (ntl == 0) then
        if (ns == 0) exit
        nall = nall + 1
        read(IOINP, '(2F11.3, i1, 8F5.1, 2x, f15.3)', ERR=400) tl, ncoef, expn, hh
        thermo(1, 1) = hform
        write(IOSCH) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
        cycle
     else if (name == 'Air') then
        sym(1) = 'N'
        fno(1) = 1.56168d0
        sym(2) = 'O'
        fno(2) = .419590d0
        sym(3) = 'AR'
        fno(3) = .009365d0
        sym(4) = 'C'
        fno(4) = .000319d0
     else if (name == 'e-') then
        mwt = 5.48579903d-04
     end if

     do i = 1, ntl
        read(IOINP, '(2F11.3, i1, 8F5.1, 2x, f15.3)', ERR=400) tl, ncoef, expn, hh
        read(IOINP, '(5d16.8/2d16.8, 16x, 2d16.8)', ERR=400) templ

        if (ifaz == 0 .and. i > 3) go to 400

        if (ifaz <= 0) then
           if (tl(2) > tgl(4)-.01d0) then
              ifaz = -1
              namee = '*' // name
              name = namee(:15)
           end if
           if (tl(1) >= tgl(i+1)) cycle
           int = i
           fill(i) = .false.
        else
           int = 1
           if (i > 1) then
              do k = 1, 7
                 thermo(k, 1) = 0
              end do
           end if
        end if

        do l = 1, ncoef
           do k = 1, 7
              if (expn(l) == real(k-3)) then
                 thermo(k, int) = templ(l)
                 exit
              end if
           end do
        end do

        thermo(8, int) = templ(8)
        thermo(9, int) = templ(9)
        if (ifaz > 0) then
           nall = nall + 1
           if (ifaz > ifzm1) then
              inew = inew + 1
           else
              inew = i
           end if
           write(IOSCH) name, ntl, date, (sym(j), fno(j), j = 1, 5), inew, tl, mwt, thermo
        end if
     end do

     ifzm1 = ifaz
     if (ifaz <= 0) then
        inew = 0
        nall = nall + 1
        if (ifaz <= 0 .and. ns == 0) then
           ngl = ngl + 1
           if (fill(3)) then
              atms = 0

              do i = 1, 5
                 if (sym(i) == ' ' .or. sym(i) == 'E') exit
                 atms = atms + fno(i)
              end do

              ! FOR GASES WITH NO COEFFICIENTS FOR TGL(3)-TGL(4) INTERVAL,
              ! CALCULATE ESTIMATED COEFFICIENTS. (STRAIGHT LINE FOR CP/R)
              aa = 2.5d0
              if (atms > 1.9) aa = 4.5d0
              if (atms > 2.1) aa = 3 * atms - 1.75d0
              ttl = tl(2)
              tx = ttl - tinf
              cpfix = 0
              templ(8) = 0
              templ(9) = 0
              dlt = log(ttl)
              do k = 7, 1, - 1
                 kk = k - 3
                 if (kk == 0) then
                    cpfix = cpfix + thermo(k, 2)
                    templ(8) = templ(8) + thermo(k, 2)
                    templ(9) = templ(9) + thermo(k, 2) * dlt
                 else
                    tex = ttl**kk
                    cpfix = cpfix + thermo(k, 2) * tex
                    templ(9) = templ(9) + thermo(k, 2) * tex / kk
                    if (kk == -1) then
                       templ(8) = templ(8) + thermo(k, 2)*dlt / ttl
                    else
                       templ(8) = templ(8) + thermo(k, 2)*tex / (kk + 1)
                    end if
                 end if
              end do
              templ(2) = (cpfix - aa) / tx
              thermo(4, 3) = templ(2)
              templ(1) = cpfix - ttl * templ(2)
              thermo(3, 3) = templ(1)
              thermo(8, 3) = thermo(8, 2) + ttl * (templ(8) - templ(1) - 0.5 * templ(2) * ttl)
              thermo(9, 3) = -templ(1) * dlt + thermo(9, 2) + templ(9) - templ(2) * ttl
           end if
        end if
        ! WRITE COEFFICIENTS ON SCRATCH I/O UNIT IOSCH
        write(IOSCH) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
     end if
  end do

! END OF DATA. COPY CONDENSED & REACTANT DATA FROM IOSCH & ADD TO IOTHM.
300 rewind IOSCH

  if (ns == 0) ns = nall
  write(IOTHM) tgl, ngl, ns, nall, Thdate
! WRITE GASEOUS PRODUCTS ON IOTHM
  if (ngl /= 0) then
     do i = 1, ns
        read(IOSCH) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
        if (ifaz <= 0) write(IOTHM) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
     end do
  end if
  if (ngl /= nall) then
! WRITE CONDENSED PRODUCTS AND REACTANTS ON IOTHM
     rewind IOSCH
     do i = 1, nall
        read(IOSCH) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
        if (i > ns) then
           write(IOTHM) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo(1, 1)
           if (ntl > 0) write(IOTHM) thermo
        else if (ifaz > 0) then
           write(IOTHM) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, (thermo(k, 1), k = 1, 9)
        end if
     end do
  end if

  return

400 write(IOOUT, '(/" ERROR IN PROCESSING thermo.inp AT OR NEAR ", A15, " (UTHERM)")') name
  readOK = .false.

  return
end subroutine UTHERM



subroutine UTRAN(readOK)
!***********************************************************************
! READ TRANSPORT PROPERTIES FORM I/O UNIT 7 IN RECORD format AND WRITE
! UNFORMATTED ON I/O UNIT IOTRN.  USES SCRATCH I/O UNIT IOSCH.
!
! UTRAN IS CALLED FROM subroutine INPUT AFTER A RECORD WITH 'tran'
! IN COLUMNS 1-4 HAS BEEN READ.
!
! NOTE:  THIS ROUTINE MAY BE CALLED DIRECTLY  AND USED BY ITSELF TO
! PROCESS THE TRANSPORT PROPERTY DATA.
!***********************************************************************
  use mod_legacy_cea
  implicit none
! DUMMY ARGUMENTS
  logical, intent(out):: readOK
! LOCAL VARIABLES
  character(16):: tname(2)
  character(1):: vorc
  integer:: i, ic, in, iv, j, k, ncc, nn, ns, nv
  real(8):: cc, tcin(6), trcoef(6, 3, 2), vvl


  ns = 0
  rewind IOSCH
  outerLoop: do
     trcoef(:, :, :) = 0
     read(IOINP, '(2A16, 2X, A1, I1, A1, I1)') tname, vvl, nv, cc, ncc
     if (tname(1) == 'end' .or. tname(1) == 'LAST') then
        write(IOTRN) ns
        rewind IOSCH
        do i = 1, ns
           read(IOSCH, ERR=200) tname, trcoef
           write(IOTRN) tname, trcoef
        end do
        return
     else
        ic = 0
        iv = 0
        nn = nv + ncc
        if (nv <= 3 .and. ncc <= 3) then
           do in = 1, nn
              read(IOINP, '(1X, A1, 2F9.2, 4E15.8)') vorc, tcin
              if (vorc == 'C') then
                 k = 2
                 ic = ic + 1
                 j = ic
              else
                 k = 1
                 iv = iv + 1
                 j = iv
              end if
              if (j > 3) exit outerLoop
              do i = 1, 6
                 trcoef(i, j, k) = tcin(i)
              end do
           end do
           ns = ns + 1
           write(IOSCH) tname, trcoef
           cycle
        end if
     end if
     exit
  end do outerLoop
200 write(IOOUT, '(/" ERROR IN PROCESSING trans.inp AT OR NEAR (UTRAN)", /1X, 2A16)') tname
  readOK = .false.
  return
end subroutine UTRAN
