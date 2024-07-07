module mod_legacy_io
  implicit none

  integer, parameter:: MAX_CHARS = 132 !< Maximum number of characters in a line (record).

  character(15), private:: fc
  integer, private:: ione, mcond, mcondf, mpn, mpnf, mvis
  real(8), private:: pfuel, phi, tem

contains

  subroutine count_cases(inp_filename, num_cases)
    character(*), intent(in):: inp_filename
    integer, intent(out):: num_cases
    integer:: io_inp
    character(MAX_CHARS):: line
    character(4):: head
    logical:: case_ended

    num_cases = 0
    case_ended = .true.

    open(newunit = io_inp, file = trim(inp_filename), status = 'old', form = 'formatted', action = 'read')

    do
       read(io_inp, '(a)', end = 999) line
       line = adjustl(line)
       head = line(1:4)

       if (head == 'ther' .or. head == 'tran' .or. head == 'prob' .or.  &
            head == 'reac' .or. head == 'outp' .or. head == 'omit' .or.  &
            head == 'only' .or. head == 'inse') then
          case_ended = .false.
       else if (head(1:3) == 'end' .and. .not. case_ended) then
          num_cases = num_cases + 1
          case_ended = .true.
       end if
    end do

999 close(io_inp)

    if (.not. case_ended) then
       num_cases = num_cases + 1
    end if

    return
  end subroutine count_cases


  subroutine INPUT(cea, readOK, caseOK, Ensert)
    !***********************************************************************
    ! DECIPHER KEYWORDS, LITERAL VARIABLES, & NUMERICAL VARIABLES IN INPUT.
    !***********************************************************************
    use mod_cea
    implicit none

    type(CEA_Problem), intent(inout):: cea
    logical, intent(out):: caseOK, readOK
    character(15), intent(out):: Ensert(20)

    ! LOCAL VARIABLES
    character(26), parameter:: lc = 'abcdefghijklmnopqrstuvwxyz', uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(15):: cin(maxNgc), cx15
    character(4):: code, cx4
    character(1):: cx1
    character(2):: cx2
    character(3):: cx3
    logical:: eqrats, incd, phi, pltdat, reacts, refl
    integer:: i, ifrmla, ii, in, iv, ix, j, jj, k, lcin(maxNgc), ncin, nmix
    real(8):: denmtr, dpin(maxNgc), eratio, hr, mix(maxNgc), ur, xyz


    open(newunit = cea%io_log, status = 'scratch', form = 'formatted')

    write(cea%io_log, '(/, /)')

    caseOK = .true.
    cea%Nonly = 0
    cea%Nomit = 0
    cea%Nsert = 0
    reacts = .false.
    cea%Trace = 0
    cea%Short = .false.
    cea%Massf = .false.
    cea%Debug(1:Ncol) = .false.
    cea%Nplt = 0
    cea%SIunit = .true.
    pltdat = .false.

    do
       do
          ! CALL INFREE TO READ DATASET
          call INFREE(readOK, cin, ncin, lcin, dpin, cea%io_log)

          if (.not. readOK) return

          code = trim(cin(1))

          if (code == '    ') then
             cycle
          end if

          ! STORE PRODUCT NAMES FROM 'ONLY' DATASET
          if (code == 'only') then
             cea%Nonly = min(maxNgc, ncin-1)
             do concurrent (i = 1:cea%Nonly)
                cea%Prod(i) = cin(i+1)
             end do

             ! STORE CONDENSED PRODUCT NAMES FROM 'INSERT' DATASET
          else if (code == 'inse') then
             cea%Nsert = min(20, ncin-1)
             do concurrent (i = 1:cea%Nsert)
                Ensert(i) = cin(i+1)
             end do

             ! STORE PRODUCT NAMES FROM 'OMIT' DATASET
          else if (code == 'omit') then
             ! CHECK OMIT DATASET
             cea%Nomit = min(maxNgc, ncin-1)
             do concurrent (i = 1:cea%Nomit)
                cea%Omit(i) = cin(i+1)
             end do

             ! KEYWORD 'THER' READ
             ! CALL UTHERM TO CONVERT FORMATTED THERMODYNAMIC DATA
          else if (code == 'ther') then
             call UTHERM(cea, readOK)
             if (.not. readOK) then
                write(cea%io_log, '(/" FATAL ERROR IN DATASET (INPUT)")')
                return
             end if

             ! KEYWORD 'TRAN' READ
             ! CALL UTRAN TO CONVERT FORMATTED TRANSPORT PROPERTIES
          else if (code == 'tran') then
             call UTRAN(cea%filename_trans_lib, IOINP, cea%io_log, readOK)
             if (.not. readOK) then
                write(cea%io_log, '(/" FATAL ERROR IN DATASET (INPUT)")')
                return
             end if

             ! PROCESS 'OUTP' DATASET.
          else if (code == 'outp') then
             do i = 2, ncin
                if (lcin(i) < 0) then
                   cx2 = cin(i)(1:2)
                   cx3 = cin(i)(1:3)
                   cx4 = cin(i)(1:4)

                   if (cx3 == 'cal') then
                      cea%SIunit = .false.

                   else if (cx4 == 'tran' .or. cx3 == 'trn') then
                      cea%Trnspt = .true.

                   else if (cx4 == 'trac') then
                      cea%Trace = dpin(i+1)

                   else if (cin(i)(:5) == 'short') then
                      cea%Short = .true.

                   else if (cin(i)(:5) == 'massf') then
                      cea%Massf = .true.

                   else if (cx3 == 'deb' .or. cx3 == 'dbg') then
                      do concurrent (j = i+1:ncin, lcin(j) == i .and. int(dpin(j)) <= Ncol)
                         cea%Debug(int(dpin(j))) = .true.
                      end do
                      do concurrent (j = i+1:ncin, lcin(j) == i)
                         lcin(j) = 0
                      end do

                   else if (cx2 == 'si') then
                      cea%SIunit = .true.

                   else if (pltdat .and. cea%Nplt < 20) then
                      cea%Nplt = cea%Nplt + 1
                      cea%Pltvar(cea%Nplt) = cin(i)

                   else if (cx2 == 'pl') then
                      pltdat = .true.

                   else
                      write(cea%io_log, '("  WARNING!!  DID NOT RECOGNIZE ", a15, " (INPUT)"/)') cin(i)
                   end if
                end if
             end do

             ! SORT AND STORE DATA FROM 'REAC' DATASET.
          else if (code == 'reac') then
             reacts = .true.
             cea%Moles = .false.
             cea%Nreac = 0
             cea%Pecwt(1:maxR) = -1

             i = 1
             do while (i < ncin)
                i = i + 1

                if (lcin(i) == 0) then
                   cycle
                end if

                if (lcin(i) > 0) then
                   write(cea%io_log, '(/" WARNING!!  LITERAL EXPECTED FOR ", a15, "(INPUT)")') cin(i)
                   cycle
                end if

                cx15 = cin(i)
                cx1 = cx15(:1)
                cx2 = cx15(:2)
                cx3 = cx15(:3)
                cx4 = cx15(:4)

                ! NEW REACTANT
                if (cx2 == 'na' .or. cx2 == 'ox' .or. cx2 == 'fu') then
                   cea%Nreac = min(cea%Nreac+1, maxR)
                   cea%Fox(cea%Nreac) = trim(cx15)
                   i = i + 1
                   if (lcin(i) < 0) cea%Rname(cea%Nreac) = cin(i)
                   ifrmla = 0
                   cea%Nfla(cea%Nreac) = 0
                   cea%Energy(cea%Nreac) = 'lib'
                   cea%Enth(cea%Nreac) = 0
                   cea%Jray(cea%Nreac) = 0
                   cea%Pecwt(cea%Nreac) = -1
                   cea%Rnum(cea%Nreac, 1) = 0
                   cea%Rmw(cea%Nreac) = 0
                   cea%Rtemp(cea%Nreac) = 0

                else
                   ! LOOK FOR PERCENTS
                   if (cx1 == 'm' .or. cx1 == 'w') then
                      if (lcin(i+1) > 0) then
                         i = i + 1
                         cea%Pecwt(cea%Nreac) = dpin(i)
                      else
                         caseOK = .false.
                         write(cea%io_log, '(/" REACTANT AMOUNT MISSING (INPUT)")')
                      end if

                      if (cx1 == 'm' .and. cea%Nreac == 1) cea%Moles = .true.

                      if (cx1 == 'm' .and. .not. cea%Moles .or. cx1 == 'w' .and. cea%Moles) then
                         caseOK = .false.
                         write(cea%io_log, '(/" MOLES AND WEIGHT PERCENTS SHOULD NOT BE MIXED (INPUT)")')
                      end if

                      ! LOOK FOR TEMPERATURES
                   else if (cx1 == 't') then
                      if (lcin(i+1) > 0) then
                         i = i + 1
                         cea%Rtemp(cea%Nreac) = dpin(i)

                         if (lcin(i-1) < 1) then
                            if (index(cx15, 'r') > 0) cea%Rtemp(cea%Nreac) = cea%Rtemp(cea%Nreac) / 1.8d0
                            if (index(cx15, 'c') > 0) cea%Rtemp(cea%Nreac) = cea%Rtemp(cea%Nreac) + 273.15d0
                            if (index(cx15, 'f') > 0) cea%Rtemp(cea%Nreac) = (cea%Rtemp(cea%Nreac) - 32) / 1.8d0 + 273.15d0
                         end if

                      else
                         write(cea%io_log, '(/" REACTANT TEMPERATURE MISSING (INPUT) ")')
                         caseOK = .false.
                      end if

                      ! LOOK FOR ENTHALPY
                   else if (cx1 == 'h' .or. cx1 == 'u') then
                      cea%Energy(cea%Nreac) = cx15

                      if (lcin(i+1) > 0) then
                         i = i + 1
                         cea%Enth(cea%Nreac) = dpin(i) * 1000 / R0

                         if (index(cin(i-1), 'c') > 0) cea%Enth(cea%Nreac) = cea%Enth(cea%Nreac) * 4.184d0
                         if (index(cin(i-1), 'k') > 0) cea%Enth(cea%Nreac) = cea%Enth(cea%Nreac) * 1000
                      end if

                      ! LOOK FOR DENSITY
                   else if (cx3 == 'rho' .or. cx3 == 'den') then
                      if (lcin(i+1) > 0) then
                         i = i + 1
                         cea%Dens(cea%Nreac) = dpin(i)
                         if (index(cx15, 'kg') > 0) cea%Dens(cea%Nreac) = cea%Dens(cea%Nreac) / 1000
                      end if

                      ! CHECK FOR CHEMICAL SYMBOLS IN EXPLODED FORMULA
                   else if ((lcin(i) == -1 .or. lcin(i) == -2) .and. index(uc, cx1) > 0) then
                      cea%Energy(cea%Nreac) = ' '
                      ifrmla = ifrmla + 1
                      cea%Nfla(cea%Nreac) = ifrmla

                      if (lcin(i) == -2) then
                         ix = index(lc, cx2(2:2))
                         if (ix > 0) cx2(2:2) = uc(ix:ix)
                      end if

                      cea%Ratom(cea%Nreac, ifrmla) = cx2
                      if (lcin(i+1) == i) then
                         cea%Rnum(cea%Nreac, ifrmla) = dpin(i+1)
                      else
                         cea%Rnum(cea%Nreac, ifrmla) = 1
                      end if

                      i = i + 1

                   else
                      write(cea%io_log, '(/" WARNING!! ", a15, " NOT RECOGNIZED (INPUT)")') cin(i)
                   end if
                end if
             end do

             ! SORT AND STORE INPUT FROM 'PROB' DATASET
          else if (code == 'prob') then
             cea%Case = ' '
             cea%P(1:maxPv) = 0
             cea%V(1:maxPv) = 0
             cea%T(1:maxT) = 0
             cea%P(1) = 1
             cea%Trace = 0
             cea%Lsave = 0
             cea%R = R0 / 4184
             cea%S0 = 0
             hr = 1.d30
             ur = 1.d30
             cea%Tp = .false.
             cea%Hp = .false.
             cea%Sp = .false.
             cea%Rkt = .false.
             cea%Shock = .false.
             cea%Detn = .false.
             cea%Vol = .false.
             cea%Ions = .false.
             cea%Eql = .false.
             cea%Froz = .false.
             cea%Fac = .false.
             cea%Debugf = .false.
             cea%Acat = 0
             cea%Ma = 0
             cea%Nfz = 1
             cea%Nsub = 0
             cea%Nsup = 0
             cea%Npp = 0
             cea%Tcest = 3800
             do concurrent (i = 1:Ncol)
                cea%Pcp(i) = 0
                cea%Pcp(i+Ncol) = 0
                cea%Supar(i) = 0
                cea%Subar(i) = 0
                cea%Mach1(i) = 0
                cea%U1(i) = 0
             end do
             cea%Gamma1 = 0
             phi = .false.
             eqrats = .false.
             incd = .false.
             refl = .false.
             cea%Shkdbg = .false.
             cea%Incdeq = .false.
             cea%Incdfz = .false.
             cea%Refleq = .false.
             cea%Reflfz = .false.
             cea%Np = 0
             cea%Nt = 1
             cea%Trnspt = .false.

             ! PROCESS LITERAL VARIABLES IN 'PROB' DATASET THAT DO NOT HAVE
             ! ASSOCIATED NUMERICAL DATA.
             do i = 2, ncin
                if (lcin(i) < 0) then
                   if (any(lcin(i+1:ncin) == i)) cycle

                   cx15 = cin(i)
                   cx2 = cx15(:2)
                   cx3 = cx15(:3)
                   cx4 = cx15(:4)

                   if (cx4 == 'case') then
                      cea%Case = cin(i+1)
                      lcin(i+1) = 0

                   else if (cx2 == 'tp' .or. cx2 == 'pt') then
                      cea%Tp = .true.

                   else if (cx2 == 'hp' .or. cx2 == 'ph') then
                      cea%Hp = .true.

                   else if (cx2 == 'sp' .or. cx2 == 'ps') then
                      cea%Sp = .true.

                   else if (cx2 == 'sv' .or. cx2 == 'vs') then
                      cea%Sp = .true.
                      cea%Vol = .true.

                   else if (cx2 == 'uv' .or. cx2 == 'vu') then
                      cea%Hp = .true.
                      cea%Vol = .true.

                   else if (cx2 == 'tv' .or. cx2 == 'vt') then
                      cea%Tp = .true.
                      cea%Vol = .true.

                   else if (cx2 == 'ro' .or. cx3 == 'rkt') then
                      cea%Rkt = .true.

                   else if (cx3 == 'dbg' .or. cx3 == 'deb') then
                      cea%Debugf = .true.
                      cea%Shkdbg = .true.
                      cea%Detdbg = .true.

                   else if (cx3 == 'fac') then
                      cea%Rkt = .true.
                      cea%Eql = .true.
                      cea%Fac = .true.
                      cea%Froz = .false.

                   else if (cx2 == 'eq') then
                      cea%Eql = .true.

                   else if (cx2 == 'fr' .or. cx2 == 'fz') then
                      cea%Froz = .true.

                   else if (cx2 == 'sh') then
                      cea%Shock = .true.

                   else if (cx3 == 'inc') then
                      cea%Shock = .true.
                      incd = .true.
                      if (index(cx15, 'eq') > 0) cea%Eql = .true.
                      if (index(cx15, 'fr') > 0) cea%Froz = .true.
                      if (index(cx15, 'fz') > 0) cea%Froz = .true.

                   else if (cx3 == 'ref') then
                      cea%Shock = .true.
                      refl = .true.
                      if (index(cx15, 'eq') > 0) cea%Eql = .true.
                      if (index(cx15, 'fr') > 0) cea%Froz = .true.
                      if (index(cx15, 'fz') > 0) cea%Froz = .true.

                   else if (cx3 == 'det') then
                      cea%Detn = .true.

                   else if (cx4 == 'ions') then
                      cea%Ions = .true.

                   else
                      write(cea%io_log, '("  WARNING!!  DID NOT RECOGNIZE ", a15, " (INPUT)"/)') cx15

                   end if
                   lcin(i) = 0
                end if
             end do

             iv = 2
             cea%Nof = 0
             exit

          else if (code(1:3) == 'end') then
             if (cea%Shock) then
                if (incd .and. cea%Froz) cea%Incdfz = .true.
                if (incd .and. cea%Eql) cea%Incdeq = .true.
                if (refl .and. cea%Froz) cea%Reflfz = .true.
                if (refl .and. cea%Eql) cea%Refleq = .true.
             end if

             cea%Hsub0 = min(hr, ur)
             cea%Size = 0

             if (hr > 0.9d30) hr = 0
             if (ur > 0.9d30) ur = 0
             if (cea%Trnspt) cea%Viscns = 0.3125 * sqrt(1.e5 * Boltz / (pi * Avgdr))
             if (cea%SIunit) cea%R = R0 / 1000

             if (.not. cea%Short) then
                write(cea%io_log, '(/" OPTIONS: TP=", l1, "  HP=", l1, "  SP=", l1, "  TV=", l1, &
                     & "  UV=", l1, "  SV=", l1, "  DETN=", l1, "  SHOCK=", l1, &
                     & "  REFL=", l1, "  INCD=", l1, /" RKT=", l1, "  FROZ=", l1, &
                     & "  EQL=", l1, "  IONS=", l1, "  SIUNIT=", l1, "  DEBUGF=", l1, &
                     & "  SHKDBG=", l1, "  DETDBG=", l1, "  TRNSPT=", l1)') &
                     cea%Tp, (cea%Hp .and. .not. cea%Vol), cea%Sp, (cea%Tp .and. cea%Vol), (cea%Hp .and. cea%Vol), &
                     (cea%Sp .and. cea%Vol), cea%Detn, cea%Shock, refl, incd, cea%Rkt, cea%Froz, cea%Eql, cea%Ions, &
                     cea%SIunit, cea%Debugf, cea%Shkdbg, cea%Detdbg, cea%Trnspt

                if (cea%T(1) > 0) then
                   write(cea%io_log, '(/" T,K =", 7f11.4)') (cea%T(jj), jj = 1, cea%Nt)
                end if

                write(cea%io_log, '(/1p, " TRACE=", e9.2, "  S/R=", e13.6, "  H/R=", e13.6, "  U/R=",  e13.6)') cea%Trace, cea%S0, hr, ur

                if (cea%Np > 0 .and. cea%Vol) then
                   write(cea%io_log, '(/" SPECIFIC VOLUME,M**3/KG =", 1p, (4e14.7))') (cea%V(jj)*1.d-05, jj = 1, cea%Np)
                end if
             end if

             if (cea%Rkt) then
                if (cea%Nt == 0) cea%Hp = .true.
                if (.not. cea%Short) then
                   write(cea%io_log, '(/" Pc,BAR =", 7f13.6)') (cea%P(jj), jj = 1, cea%Np)
                   write(cea%io_log, '(/" Pc/P =", 9f11.4)') (cea%Pcp(jj), jj = 1, cea%Npp)
                   write(cea%io_log, '(/" SUBSONIC AREA RATIOS =", (5f11.4))') (cea%Subar(i), i = 1, cea%Nsub)
                   write(cea%io_log, '(/" SUPERSONIC AREA RATIOS =", (5f11.4))') (cea%Supar(i), i = 1, cea%Nsup)
                   write(cea%io_log, '(/" NFZ=", i3, 1p, "  Mdot/Ac=", e13.6, "  Ac/At=", e13.6)') cea%Nfz, cea%Ma, cea%Acat
                end if
             else
                if (.not. cea%Vol .and. .not. cea%Short) write(cea%io_log, '(/" P,BAR =", 7f13.6)') (cea%P(jj), jj = 1, cea%Np)
             end if

             if (reacts) call REACT(cea)
             if (cea%Nreac == 0 .or. cea%Nlm <= 0) then
                write(cea%io_log, '(/" ERROR IN REACTANTS DATASET (INPUT)")')
                caseOK = .false.
                write(cea%io_log, '(/" FATAL ERROR IN DATASET (INPUT)")')
                return
             end if

             if (cea%Nof == 0) then
                cea%Nof = 1
                cea%Oxf(1) = 0
                if (cea%Wp(2) > 0) then
                   cea%Oxf(1) = cea%Wp(1) / cea%Wp(2)
                else
                   caseOK = .false.
                   write(cea%io_log, '(/" REACTANT AMOUNT MISSING (INPUT)")')
                   write(cea%io_log, '(/" FATAL ERROR IN DATASET (INPUT)")')
                   return
                end if

             else if (phi .or. eqrats) then
                do i = 1, cea%Nof
                   eratio = cea%Oxf(i)
                   if (eqrats) then
                      xyz = -eratio * cea%Vmin(2) - cea%Vpls(2)
                      denmtr = eratio * cea%Vmin(1) + cea%Vpls(1)
                   else
                      xyz = -cea%Vmin(2) - cea%Vpls(2)
                      denmtr = eratio * (cea%Vmin(1) + cea%Vpls(1))
                   end if

                   if (abs(denmtr) < 1.d-30) then
                      caseOK = .false.
                      write(cea%io_log, '(/" UNABLE TO PROCESS EQUIVALENCE RATIO =", E11.4, "(INPUT)")') eratio
                      write(cea%io_log, '(/" FATAL ERROR IN DATASET (INPUT)")')
                      return
                   end if
                   cea%Oxf(i) = xyz / denmtr
                end do
             end if

             if (.not. (cea%Sp .or. cea%Tp .or. cea%Hp .or. cea%Rkt .or. cea%Detn .or. cea%Shock)) then
                caseOK = .false.
                write(cea%io_log, '(/" TYPE OF PROBLEM NOT SPECIFIED (INPUT)")')

             else if (cea%Tp .and. cea%T(1) <= 0) then
                caseOK = .false.
                write(cea%io_log, '(/" ASSIGNED VALUES OF TEMPERATURE ARE MISSING IN prob", " DATASET (INPUT)")')

             else if (cea%Np <= 0) then
                caseOK = .false.
                write(cea%io_log, '(/" ASSIGNED PRESSURE (OR DENSITY) MISSING IN prob", " DATASET (INPUT)")')
             end if

             if (.not. (caseOK .and. cea%Nlm > 0)) write(cea%io_log, '(/" FATAL ERROR IN DATASET (INPUT)")')
             return

          else
             write(cea%io_log, '(/" WARNING!!  A KEYWORD IS MISSING (INPUT)")')
          end if

       end do

       ! PROCESS NUMERICAL DATA FOLLOWING 'PROB' LITERALS
       do
          in = 0
          nmix = 0
          ii = iv
          do i = ii, ncin
             iv = i
             if (lcin(i) /= 0) then
                if (lcin(i) < 0) then
                   if (in > 0) exit
                   in = i
                else
                   if (lcin(i) /= in) exit
                   nmix = nmix + 1
                   mix(nmix) = dpin(i)
                   lcin(i) = 0
                end if
             end if
          end do

          if (nmix <= 0) then
             if (iv < ncin) cycle
             exit
          end if

          cx15 = cin(in)
          cx1 = cx15(:1)
          cx2 = cx15(:2)
          cx3 = cx15(:3)
          cx4 = cx15(:4)

          if (cx1 == 't') then
             cea%Nt = nmix
             if (nmix > maxMix) then
                cea%Nt = maxMix
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 't', cea%Nt
             end if

             do i = 1, cea%Nt
                if (cx4 /= 'tces') then
                   cea%T(i) = mix(i)
                   if (lcin(in) < -1) then
                      if (index(cx15, 'r') > 0) cea%T(i) = cea%T(i) / 1.8d0
                      if (index(cx15, 'c') > 0) cea%T(i) = cea%T(i) + 273.15d0
                      if (index(cx15, 'f') > 0) cea%T(i) = (cea%T(i) - 32) / 1.8d0 + 273.15d0
                   end if
                end if
             end do

          else if ((cx2 == 'pc' .or. cx2 == 'pi') .and. index(cx15(3:15), 'p') > 0 .and. index(cx15, 'psi') == 0) then
             cea%Npp = nmix
             if (nmix > 2*Ncol) then
                cea%Npp = 2 * Ncol
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'pcp', cea%Npp
             end if
             cea%Pcp(1:cea%Npp) = mix(1:cea%Npp)

          else if (cx1 == 'p' .and. cx3 /= 'phi') then
             cea%Np = nmix
             if (nmix > maxPv) then
                cea%Np = maxPv
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'p', cea%Np
             end if

             do i = 1, cea%Np
                cea%P(i) = mix(i)
                if (index(cx15, 'psi') /= 0) then
                   cea%P(i) = cea%P(i) / 14.696006d0

                else if (index(cx15, 'mmh') /= 0) then
                   cea%P(i) = cea%P(i) / 760

                else if (index(cx15, 'atm') == 0) then
                   cycle
                end if

                cea%P(i) = cea%P(i) * 1.01325d0
             end do

          else if (cx3 == 'rho') then
             xyz = 1.d02
             if (index(cx15, 'kg') /= 0) xyz = 1.d05
             cea%Np = nmix
             if (nmix > maxPv) then
                cea%Np = maxPv
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'rho', cea%Np
             end if
             cea%V(1:cea%Np) = xyz / mix(1:cea%Np)

          else if (cx1 == 'v') then
             xyz = 1.d02
             if (index(cx15, 'kg') /= 0) xyz = 1.d05
             cea%Np = nmix
             if (nmix > maxPv) then
                cea%Np = maxPv
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'v', cea%Np
             end if
             cea%V(1:cea%Np) = mix(1:cea%Np) * xyz

          else if (cx3 == 'nfz' .or. cx3 == 'nfr') then
             cea%Nfz = int(mix(1))
             cea%Froz = .true.

          else if (cx4 == 'tces') then
             cea%Tcest = mix(1)

          else if (cx4 == 'trac') then
             cea%Trace = mix(1)

          else if (cx3 == 's/r') then
             cea%S0 = mix(1)

          else if (cx3 == 'u/r' .or. cx2 == 'ur') then
             ur = mix(1)

          else if (cx3 == 'h/r' .or. cx2 == 'hr') then
             hr = mix(1)

          else if (cx2 == 'u1') then
             cea%Nsk = nmix
             if (nmix > Ncol) then
                cea%Nsk = Ncol
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'u1', cea%Nsk
             end if
             cea%U1(1:cea%Nsk) = mix(1:cea%Nsk)

          else if (cx4 == 'mach') then
             cea%Nsk = nmix
             if (nmix > Ncol) then
                cea%Nsk = Ncol
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'mach1', cea%Nsk
             end if
             cea%Mach1(1:cea%Nsk) = mix(1:cea%Nsk)

          else if (cx3 == 'sub') then
             cea%Nsub = nmix
             if (nmix > 13) then
                cea%Nsub = 13
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'subar', cea%Nsub
             end if
             cea%Subar(1:cea%Nsub) = mix(1:cea%Nsub)

          else if (cx3 == 'sup') then
             cea%Nsup = nmix
             if (nmix > 13) then
                cea%Nsup = 13
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'supar', cea%Nsup
             end if
             cea%Supar(1:cea%Nsup) = mix(1:cea%Nsup)

          else if (cx2 == 'ac') then
             cea%Acat = mix(1)

          else if (cx4 == 'mdot' .or. cx2 == 'ma') then
             cea%Ma = mix(1)

          else if (cx4 == 'case') then
             cea%Case = cin(in+1)
             lcin(in+1) = 0

          else if (cea%Nof == 0 .and. (cx3 == 'phi' .or. cx3 == 'o/f' .or. cx3 == 'f/a' .or. cx2 == '%f' .or. cx1 == 'r')) then
             cea%Nof = nmix
             if (nmix > maxMix) then
                cea%Nof = maxMix
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'o/f', cea%Nof
             end if
             cea%Oxf(1:cea%Nof) = mix(1:cea%Nof)

             if (cx3 == 'phi') then
                phi = .true.

             else if (cx1 == 'r') then
                eqrats = .true.

             else if (cx3 == 'f/a') then
                do concurrent (k = 1:cea%Nof, cea%Oxf(k) > 0)
                   cea%Oxf(k) = 1/cea%Oxf(k)
                end do

             else if (cx4 == '%fue') then
                do concurrent (k = 1:cea%Nof, cea%Oxf(k) > 0)
                   cea%Oxf(k) = (100 - cea%Oxf(k)) / cea%Oxf(k)
                end do
             end if

          else
             write(cea%io_log, '("  WARNING!!  DID NOT RECOGNIZE ", a15, " (INPUT)"/)') cx15
          end if

          if (iv >= ncin) exit
       end do
    end do

    return
  end subroutine INPUT


  subroutine INFREE(readOK, Cin, Ncin, Lcin, Dpin, io_log)
    !***********************************************************************
    ! FREE-FORM READ FOR CEA.  READS AND DECIPHERS DATA FOR ONE DATASET.
    !
    ! DEFINITIONS:
    !   Ch1  - INDIVIDUAL CHARACTERS IN RECORD, MAXIMUM 132.
    !   Nch1 - COLUMN NUMBER FOR THE LAST NON-BLANK CHARACTER IN RECORD.
    !   Ncin - NUMBER OF VARIABLES IN DATASET.
    !   Cin  - CHARACTER STRINGS IN DATASET. MAXIMUM 15 CHARACTERS.
    !   LCIN - NEG. LENGTH OF LITERALS.  FOR NUMERICS, INDEX OF PREVIOUS
    !          LITERAL.  ZERO FOR UNACCEPTIBLE VARIABLES.  VARIABLE
    !          FOLLOWING "CASE" IS ALWAYS ASSUMED TO BE LITERAL.
    !   Nb   - NUMBER OF DELIMITERS IN STRING.
    !   Nx   - NUMBER OF CHARACTERS IN STRING.
    !   Dpin - NUMERICS AS DOUBLE PRECISION VARIABLE.
    !   Cnum - CHARACTER STRING REPRESENTING DATASET NUMBERS. MAXIMUM 24
    !          CHARACTERS.
    !***********************************************************************
    use mod_cea
    implicit none

    ! DUMMY ARGUMENTS
    character(15), intent(out):: Cin(maxNgc)
    integer, intent(out):: Ncin
    integer, intent(out):: Lcin(maxNgc)
    logical, intent(out):: readOK
    real(8), intent(out):: Dpin(maxNgc)
    integer, intent(in):: io_log

    ! LOCAL VARIABLES
    character(1), parameter:: nums(13) = ['+', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.']
    character(2), parameter:: numg(24) = &
         [character(2):: '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', &
         '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24']
    character(1):: ch1(132), cx_
    character(24):: cnum
    character(3):: fmtl(3)
    character(4):: w1
    integer:: i, ich1, j, kcin, nb, nch1, nx

    fmtl = [character(3):: '(g', '16', '.0)']

    Ncin = 1
    Lcin(1) = 0
    kcin = 0
    Dpin(1) = 0

    do
       nb = 1
       nx = 0
       cnum = ' '
       Cin(Ncin) = ' '
       ch1(1) = ' '
       nch1 = 1

       ! READ CHARACTERS, ONE AT A TIME
       read(IOINP, '(132a1)', END = 500, ERR = 500) ch1

       ! FIND FIRST AND LAST NON-BLANK CHARACTER
       do i = 132, 1, - 1
          nch1 = i
          if (ch1(i) /= ' ' .and. ch1(i) /= '	') exit
       end do

       do i = 1, nch1
          ich1 = i
          if (ch1(i) /= ' ' .and. ch1(i) /= '	') exit
       end do

       if (nch1 == 1 .or. ch1(ich1) == '#' .or. ch1(ich1) == '!') then
          write(io_log, '(1x, 80a1)') (ch1(i), i = 1, nch1)
          cycle
       end if

       w1 = ch1(ich1) // ch1(ich1+1) // ch1(ich1+2) // ch1(ich1+3)

       ! IS STRING A KEYWORD SIGNALLING START OR END OF DATASET?
       if (w1 == 'ther' .or. w1 == 'tran' .or. w1 == 'prob' .or.  &
            w1 == 'reac' .or. w1 == 'outp' .or. w1 == 'omit' .or.  &
            w1 == 'only' .or. w1 == 'inse' .or. w1(1:3) == 'end') then
          if (Ncin == 1) then
             Cin(Ncin) = w1
             if (w1(1:3) == 'end' .or. w1 == 'ther' .or. w1 == 'tran') then
                write(io_log, '(1x, 80a1)') (ch1(i), i = 1, nch1)
                return
             end if
             ich1 = ich1 + 4
             nx = 4
             Lcin(1) = -4
          else

             ! KEYWORD READ FOR NEXT DATASET. END PROCESSING
             backspace IOINP
             if (nx == 0) Ncin = Ncin - 1
             return
          end if

       else if (Ncin == 1) then
          write(io_log, '(/" FATAL ERROR IN INPUT format (INFREE)")')
          go to 500
       end if

       write(io_log, '(1x, 80a1)') (ch1(i), i = 1, nch1)

       do i = ich1, nch1
          cx_ = ch1(i)

          ! LOOK FOR DELIMITER STRINGS
          if (cx_ == ',' .and. (Lcin(Ncin) > 0 .or. nx == 0)) cx_ = ' '
          if (cx_ == '=' .and. (Lcin(Ncin) < 0 .or. nx == 0)) cx_ = ' '
          if (cx_ /= ' ' .and. cx_ /= '	') then

             ! LOOK FOR CHARACTER STRINGS
             nx = nx + 1
             if (Ncin > 1) then
                cnum(nx:nx) = cx_
                if (nx <= 15) Cin(Ncin) = trim(cnum)
                if (nx == 1) then

                   ! IS THIS A NUMERIC?
                   do j = 1, 13
                      if (ch1(i) == nums(j)) then
                         Lcin(Ncin) = kcin
                         go to 310
                      end if
                   end do

                   Lcin(Ncin) = -1
                   kcin = Ncin
                else if (Lcin(Ncin) < 0) then
                   Lcin(Ncin) = -nx
                end if

310             nb = 1
             end if

             if (i < nch1 .or. Lcin(Ncin) < 0) cycle
          end if

          if (nb == 1. .and. nx > 0) then
             if (Ncin > 0 .and. Lcin(Ncin) > 0) then

                ! CONVERT NUMERIC CHARACTER STRINGS TO real(8) VARIABLES (DPIN)
                fmtl(2) = numg(min(24, nx))

                ! INTERNAL READ TO CONVERT TO NUMERIC
                read(cnum, fmtl, ERR=320) Dpin(Ncin)
             end if

             go to 340

320          if (Cin(Ncin-1)(:4) /= 'case') write(io_log, '(/" WARNING!!  UNACCEPTABLE NUMBER ", a15, " (INFREE)")') Cin(i)
             Lcin(Ncin) = 0

340          Ncin = Ncin + 1
             Cin(Ncin) = ' '
             Lcin(Ncin) = 0
             Dpin(Ncin) = 0
             nx = 0
             cnum = ' '
          end if

          nb = nb + 1

       end do

       if (nx > 0) then
          Ncin = Ncin + 1
          Lcin(Ncin) = 0
          Dpin(Ncin) = 0
       end if

    end do

500 readOK = .false.

    return
  end subroutine INFREE


  subroutine OUT0(cea)
    use mod_cea
    implicit none

    type(CEA_Problem), intent(inout):: cea
    character(MAX_CHARS):: line

    rewind cea%io_log

    do
       read(cea%io_log, '(a)', end = 999) line
!!$       write(0, '(a)') trim(line)
       write(IOOUT, '(a)') trim(line)
    end do

999 close(cea%io_log)

    return
  end subroutine OUT0


  subroutine OUT1(cea)
    !***********************************************************************
    ! OUT1 WRITES REACTANT AND FUEL-OXIDANT RATIO INformatION.
    !
    ! NOTE - ROCKET, SHOCK, AND DETON PROBLEMS HAVE ADDITIONAL OUTPUT.
    !***********************************************************************
    use mod_cea
    implicit none

    type(CEA_Problem), intent(inout):: cea

    integer:: n
    real(8):: rho

    write(IOOUT, '(" CASE = ", a15)') cea%Case

    if (cea%Moles) then
       write(IOOUT, '(/13X, "REACTANT", 20x, a11, "      ENERGY", 6x, "TEMP")') '   MOLES   '
       if (.not. cea%SIunit) write(IOOUT, '(57X, " CAL/MOL ", 6x, "K")')
       if (cea%SIunit) write(IOOUT, '(57X, "KJ/KG-MOL", 6x, "K")')
    else
       write(IOOUT, '(/13X, "REACTANT", 20x, a11, "      ENERGY", 6x, "TEMP")') 'WT FRACTION'
       if (.not. cea%SIunit) write(IOOUT, '(42X, "(SEE NOTE)      CAL/MOL       K  ")')
       if (cea%SIunit) write(IOOUT, '(42X, "(SEE NOTE)     KJ/KG-MOL      K  ")')
    end if

    do n = 1, cea%Nreac
       write(IOOUT, '(1x, a8, 4x, a15, 11x, f12.7, f14.3, f11.3)') &
            cea%Fox(n), cea%Rname(n), cea%Pecwt(n), cea%Enth(n)*cea%R, cea%Rtemp(n)
    end do

    phi = 0
    tem = (cea%Vpls(1) + cea%Vmin(1)) * cea%Oxfl
    if (ABS(tem) >= 1.d-3) phi = -(cea%Vmin(2) + cea%Vpls(2)) / tem

    if (cea%Fox(1) == 'NAME') then
       pfuel = 0
    else
       pfuel = 100 / (1 + cea%Oxfl)
    end if

    if (any(cea%Rh /= 0)) then
       if (any(cea%Rh == 0)) then
          rho = max(cea%Rh(1), cea%Rh(2))
       else
          rho = (cea%Oxfl + 1) * cea%Rh(1) * cea%Rh(2) / (cea%Rh(1) + cea%Oxfl * cea%Rh(2))
       end if
       if (cea%SIunit) then
          rho = rho * 1000
          write(IOOUT, '(/" REACTANT DENSITY=", F8.2, " KG/CU M")') rho
       else
          write(IOOUT, '(/" REACTANT DENSITY=", F8.4, " G/CC")') rho
       end if
    end if

    write(IOOUT, '(/" O/F=", F11.5, 2X, "%FUEL=", F10.6, 2X, "R,EQ.RATIO=", F9.6, 2X, "PHI,EQ.RATIO=", F9.6)') &
         cea%Oxfl, pfuel, cea%Eqrat, phi
    return
  end subroutine OUT1


  subroutine OUT2(cea)
    !***********************************************************************
    ! OUT2 WRITES THERMODYNAMIC PROPERTIES.
    !***********************************************************************
    use mod_cea
    implicit none

    type(CEA_Problem), intent(inout):: cea

    character(15):: fgi, fh, fp, frh, fs, fu
    integer:: i, j
    integer, save:: mcp, mdvp, mdvt, meq, mfa, mg, mgam, mh, mie, mm, mmw, mof, mp, mpf, mph, mrho, ms, mson, mt
    real(8):: pfactor, vnum

    ione = 0
    if (cea%Rkt .and. .not. cea%Page1) then
       ione = 2
       if (cea%Iopt /= 0) ione = 3
    end if

    ! SET MXX ARRAY FOR PLOTTING PARAMETERS
    mp     = 0
    mt     = 0
    mrho   = 0
    mh     = 0
    mie    = 0
    mg     = 0
    ms     = 0
    mm     = 0
    mcp    = 0
    mgam   = 0
    mson   = 0
    mcond  = 0
    mvis   = 0
    mpn    = 0
    mpf    = 0
    mof    = 0
    mph    = 0
    meq    = 0
    mfa    = 0
    mmw    = 0
    mdvt   = 0
    mdvp   = 0
    mcondf = 0
    mpnf   = 0

    do i = 1, cea%Nplt
       if (index(cea%Pltvar(i)(2:), '1') == 0) then
          if (index(cea%Pltvar(i)(1:), 'dlnt') /= 0) then
             mdvt = i
          else if (index(cea%Pltvar(i)(1:), 'dlnp') /= 0) then
             mdvp = i
          else if (cea%Pltvar(i)(:4) == 'pran') then
             if (index(cea%Pltvar(i)(3:), 'fz') /= 0 .or. index(cea%Pltvar(i)(3:), 'fr') /= 0) then
                mpnf = i
             else
                mpn = i
             end if
          else if (cea%Pltvar(i)(:4) == 'cond') then
             if (index(cea%Pltvar(i)(3:), 'fz') /= 0 .or. index(cea%Pltvar(i)(3:), 'fr') /= 0) then
                mcondf = i
             else
                mcond = i
             end if
          else if (cea%Pltvar(i)(:3) == 'phi') then
             mph = i
          else if (cea%Pltvar(i)(:2) == 'p ') then
             mp = i
          else if (cea%Pltvar(i)(:1) == 't') then
             mt = i
          else if (cea%Pltvar(i)(:3) == 'rho') then
             mrho = i
          else if (cea%Pltvar(i)(:1) == 'h') then
             mh = i
          else if (cea%Pltvar(i)(:1) == 'u') then
             mie = i
          else if (cea%Pltvar(i)(:3) == 'gam') then
             mgam = i
          else if (cea%Pltvar(i)(:3) == 'son') then
             mson = i
          else if (cea%Pltvar(i)(:2) == 'g ') then
             mg = i
          else if (cea%Pltvar(i)(:2) == 's ') then
             ms = i
          else if (cea%Pltvar(i)(:1) == 'm' .and. cea%Pltvar(i)(:2) /= 'ma') then
             if (.not. cea%Gonly .and. cea%Pltvar(i)(:2) == 'mw') then
                mmw = i
             else
                mm = i
             end if
          else if (cea%Pltvar(i)(:2) == 'cp') then
             mcp = i
          else if (cea%Pltvar(i)(:3) == 'vis') then
             mvis = i
          else if (cea%Pltvar(i)(:3) == 'o/f') then
             mof = i
          else if (cea%Pltvar(i)(:2) == '%f') then
             mpf = i
          else if (cea%Pltvar(i)(:3) == 'f/a') then
             mfa = i
          else if (cea%Pltvar(i)(:1) == 'r') then
             meq = i
          end if
       end if
    end do
    do i = cea%Iplt + 1, cea%Iplt + cea%Npt
       if (mof > 0) cea%Pltout(i, mof) = cea%Oxfl
       if (mpf > 0) cea%Pltout(i, mpf) = pfuel
       if (mph > 0) cea%Pltout(i, mph) = phi
       if (mfa > 0) cea%Pltout(i, mfa) = 1.d0/cea%Oxfl
       if (meq > 0) cea%Pltout(i, meq) = cea%Eqrat
    end do

    if (cea%SIunit) then
       pfactor = 1
       fp = 'P, BAR'
       vnum = 1.d05
       frh = 'RHO, KG/CU M'
       fh = 'H, KJ/KG'
       fu = 'U, KJ/KG'
       fgi = 'G, KJ/KG'
       fs = 'S, KJ/(KG)(K)'
       fc = 'Cp, KJ/(KG)(K)'
    else
       pfactor = 1 / 1.01325d0
       fp = 'P, ATM'
       vnum = 100
       frh = 'RHO, G/CC'
       fh = 'H, CAL/G'
       fu = 'U, CAL/G'
       fgi = 'G, CAL/G'
       fs = 'S, CAL/(G)(K)'
       fc = 'Cp, CAL/(G)(K)'
    end if

    cea%fmt(4) = cea%fmt(6)

    ! PRESSURE
    call VARFMT(cea, cea%Ppp)

    do i = 1, cea%Npt
       cea%X(i) = cea%Ppp(i) * pfactor
       if (cea%Nplt /= 0 .and. i > ione) then
          if (mp > 0) cea%Pltout(i + cea%Iplt - ione, mp) = cea%X(i)
          if (mt > 0) cea%Pltout(i + cea%Iplt - ione, mt) = cea%Ttt(i)
       end if
    end do

    write(IOOUT, cea%fmt) fp, (cea%X(j), j = 1, cea%Npt)

    ! TEMPERATURE
    cea%fmt(4) = '13'
    cea%fmt(5) = ' '
    cea%fmt(7) = '2,'
    write(IOOUT, cea%fmt) 'T, K            ', (cea%Ttt(j), j = 1, cea%Npt)

    ! DENSITY
    do i = 1, cea%Npt
       if (cea%Vlm(i) /= 0) cea%X(i) = vnum / cea%Vlm(i)
       if (cea%Nplt /= 0 .and. i > ione .and. mrho > 0) cea%Pltout(i+cea%Iplt-ione, mrho) = cea%X(i)
    end do
    call EFMT(cea%fmt(4), frh, cea%X, cea%Npt)

    ! ENTHALPY
    do i = 1, cea%Npt
       cea%X(i) = cea%Hsum(i) * cea%R
       if (cea%Nplt /= 0 .and. i > ione .and. mh > 0) cea%Pltout(i+cea%Iplt-ione, mh) = cea%X(i)
    end do
    cea%fmt(4) = cea%fmt(6)
    call VARFMT(cea, cea%X)
    write(IOOUT, cea%fmt) fh, (cea%X(j), j = 1, cea%Npt)

    ! INTERNAL ENERGY
    do i = 1, cea%Npt
       cea%X(i) = (cea%Hsum(i) - cea%Ppp(i) * cea%Vlm(i) / R0) * cea%R
       if (cea%Nplt /= 0 .and. i > ione .and. mie > 0) cea%Pltout(i+cea%Iplt-ione, mie) = cea%X(i)
    end do
    call VARFMT(cea, cea%X)
    write(IOOUT, cea%fmt) fu, (cea%X(j), j = 1, cea%Npt)

    ! GIBBS ENERGY
    do i = 1, cea%Npt
       cea%X(i) = (cea%Hsum(i) - cea%Ttt(i) * cea%Ssum(i)) * cea%R
       if (cea%Nplt /= 0 .and. i > ione) then
          if (mg > 0) cea%Pltout(i+cea%Iplt-ione, mg) = cea%X(i)
          if (mm > 0) cea%Pltout(i+cea%Iplt-ione, mm) = cea%Wm(i)
          if (mmw > 0) cea%Pltout(i+cea%Iplt-ione, mmw) = 1 / cea%Totn(i)
          if (ms > 0) cea%Pltout(i+cea%Iplt-ione, ms) = cea%Ssum(i) * cea%R
          if (mcp > 0) cea%Pltout(i+cea%Iplt-ione, mcp) = cea%Cpr(i) * cea%R
          if (mgam > 0) cea%Pltout(i+cea%Iplt-ione, mgam) = cea%Gammas(i)
          if (mdvt > 0) cea%Pltout(i+cea%Iplt-ione, mdvt) = cea%Dlvtp(i)
          if (mdvp > 0) cea%Pltout(i+cea%Iplt-ione, mdvp) = cea%Dlvpt(i)
       end if
    end do
    call VARFMT(cea, cea%X)
    write(IOOUT, cea%fmt) fgi, (cea%X(j), j = 1, cea%Npt)

    ! ENTROPY
    cea%fmt(4) = '13'
    cea%fmt(5) = ' '
    cea%fmt(7) = '4,'
    write(IOOUT, cea%fmt) fs, (cea%Ssum(j) * cea%R, j = 1, cea%Npt)
    write(IOOUT, *)

    ! MOLECULAR WEIGHT
    cea%fmt(7) = '3,'
    write(IOOUT, cea%fmt) 'M, (1/n)        ', (cea%Wm(j), j = 1, cea%Npt)
    if (.not. cea%Gonly) write(IOOUT, cea%fmt) 'MW, MOL WT      ', (1/cea%Totn(j), j = 1, cea%Npt)

    ! (DLV/DLP)T
    cea%fmt(7) = '5,'
    if (cea%Eql) write(IOOUT, cea%fmt) '(dLV/dLP)t      ', (cea%Dlvpt(j), j = 1, cea%Npt)

    ! (DLV/DLT)P
    cea%fmt(7) = '4,'
    if (cea%Eql) write(IOOUT, cea%fmt) '(dLV/dLT)p      ', (cea%Dlvtp(j), j = 1, cea%Npt)

    ! HEAT CAPACITY
    write(IOOUT, cea%fmt) fc, (cea%Cpr(j) * cea%R, j = 1, cea%Npt)

    ! GAMMA(S)
    cea%fmt(7) = '4,'
    write(IOOUT, cea%fmt) 'GAMMAs          ', (cea%Gammas(j), j = 1, cea%Npt)

    ! SONIC VELOCITY
    cea%fmt(7) = '1,'
    do i = 1, cea%Npt
       cea%Sonvel(i) = sqrt(R0 * cea%Gammas(i) * cea%Ttt(i) / cea%Wm(i))
       if (cea%Nplt /= 0 .and. i > ione .and. mson > 0) cea%Pltout(i+cea%Iplt-ione, mson) = cea%Sonvel(i)
    end do
    write(IOOUT, cea%fmt) 'SON VEL,M/SEC   ', (cea%Sonvel(j), j = 1, cea%Npt)
    return
  end subroutine OUT2


  subroutine OUT3(cea)
    !***********************************************************************
    ! OUT3 WRITES MOLE FRACTIONS.
    !***********************************************************************
    use mod_cea
    implicit none

    type(CEA_Problem), intent(inout):: cea

    character(4):: mamo
    logical:: kOK
    integer:: i, k, m, im, kin
    integer, save:: notuse
    real(8):: tra

    tra = 5.d-6
    if (cea%Trace /= 0.) tra = cea%Trace

    ! MASS OR MOLE FRACTIONS
    if (cea%Massf) then
       mamo = 'MASS'
    else
       mamo = 'MOLE'
    end if

    if (cea%Eql) then
       write(IOOUT, '(/1x, A4, " FRACTIONS"/)') mamo
       notuse = 0

       do k = 1, cea%Ngc
          kOK = .true.

          if (k > cea%Ng .and. k < cea%Ngc .and. cea%Prod(k) == cea%Prod(k+1)) then
             kOK = .false.
             im = 0
          else
             do m = 1, cea%Nplt
                im = 0
                if (cea%Pltvar(m) == cea%Prod(k) .or. '*' // cea%Pltvar(m) == cea%Prod(k)) then
                   im = m
                   exit
                end if
             end do
          end if

          kin = 0
          do i = 1, cea%Npt
             if (cea%Massf) then
                tem = cea%Mw(k)
             else
                tem = 1 / cea%Totn(i)
             end if
             if (k <= cea%Ng) then
                cea%X(i) = cea%En(k, i) * tem
             else
                if (cea%Prod(k) /= cea%Prod(k-1)) cea%X(i) = 0
                if (cea%En(k, i) > 0) cea%X(i) = cea%En(k, i) * tem
             end if
             if (cea%Nplt /= 0 .and. i > ione .and. im > 0) cea%Pltout(i+cea%Iplt-ione, im) = cea%X(i)
             if (kOK .and. cea%X(i) >= tra) kin = 1
          end do

          if (kin == 1) then
             if (cea%Trace == 0) then
                write(IOOUT, '(1x, A15, F9.5, 12F9.5)') cea%Prod(k), (cea%X(i), i = 1, cea%Npt)
             else
                call EFMT(cea%fmt(4), cea%Prod(k), cea%X, cea%Npt)
             end if
             if (cea%Prod(k) == cea%Omit(notuse)) notuse = notuse - 1
          else if (cea%Prod(k) /= cea%Prod(k-1)) then
             notuse = notuse + 1
             cea%Omit(notuse) = cea%Prod(k)
          end if
       end do
    end if

    write(IOOUT, '(/"  * THERMODYNAMIC PROPERTIES FITTED TO", f7.0, "K")') cea%Tg(4)
    if (.not. cea%Short) then
       write(IOOUT, '(/"    PRODUCTS WHICH WERE CONSIDERED BUT WHOSE ", a4, &
            & " FRACTIONS", /"    WERE LESS THAN", 1pe13.6, &
            & " FOR ALL ASSIGNED CONDITIONS"/)') mamo, tra
       write(IOOUT, '(5(1x, a15))') (cea%Omit(i), i = 1, notuse)
    end if

    if (.not. cea%Moles) write(IOOUT, '(/" NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF", &
         & " OXIDANT IN TOTAL OXIDANTS")')
    return
  end subroutine OUT3


  subroutine OUT4(cea)
    !***********************************************************************
    ! OUT4 WRITES TRANSPORT PROPERTIES.
    !***********************************************************************
    use mod_cea
    implicit none

    type(CEA_Problem), intent(inout):: cea

    integer:: i, j

    write(IOOUT, *)
    write(IOOUT, '(" TRANSPORT PROPERTIES (GASES ONLY)")')
    if (cea%SIunit) then
       write(IOOUT, '("   CONDUCTIVITY IN UNITS OF MILLIWATTS/(CM)(K)"/)')
    else
       write(IOOUT, '("   CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)"/)')
    end if

    ! TRANSPORT PROPERTIES
    cea%fmt(4) = cea%fmt(6)
    if (cea%Nplt > 0) then
       do i = 1, cea%Npt
          if (i > ione) then
             if (mvis > 0) cea%Pltout(i+cea%Iplt-ione, mvis) = cea%Vis(i)
             if (mcond > 0) cea%Pltout(i+cea%Iplt-ione, mcond) = cea%Coneql(i)
             if (mpn > 0) cea%Pltout(i+cea%Iplt-ione, mpn) = cea%Preql(i)
             if (mcondf > 0) cea%Pltout(i+cea%Iplt-ione, mcondf) = cea%Confro(i)
             if (mpnf > 0) cea%Pltout(i+cea%Iplt-ione, mpnf) = cea%Prfro(i)
          end if
       end do
    end if

    call VARFMT(cea, cea%Vis)

    write(IOOUT, cea%fmt) 'VISC,MILLIPOISE', (cea%Vis(j), j = 1, cea%Npt)

    cea%fmt(4) = '13'
    cea%fmt(5) = ' '
    cea%fmt(7) = '4,'

    if (cea%Eql) then
       write(IOOUT, '(/"  WITH EQUILIBRIUM REACTIONS"/)')
       ! SPECIFIC HEAT
       write(IOOUT, cea%fmt) fc, (cea%Cpeql(j), j = 1, cea%Npt)
       ! CONDUCTIVITY
       write(IOOUT, cea%fmt) 'CONDUCTIVITY    ', (cea%Coneql(j), j = 1, cea%Npt)
       ! PRANDTL NUMBER
       write(IOOUT, cea%fmt) 'PRANDTL NUMBER  ', (cea%Preql(j), j = 1, cea%Npt)
    end if

    write(IOOUT, '(/"  WITH FROZEN REACTIONS"/)')
    ! SPECIFIC HEAT
    write(IOOUT, cea%fmt) fc, (cea%Cpfro(j), j = 1, cea%Npt)
    ! CONDUCTIVITY
    write(IOOUT, cea%fmt) 'CONDUCTIVITY    ', (cea%Confro(j), j = 1, cea%Npt)
    ! PRANDTL NUMBER
    write(IOOUT, cea%fmt) 'PRANDTL NUMBER  ', (cea%Prfro(j), j = 1, cea%Npt)

    return
  end subroutine OUT4


  subroutine EFMT(Fone, Aa, Vx, Npt)
    !***********************************************************************
    ! WRITE OUTPUT RECORD WITH NUMERICAL VALUES IN SPECIAL EXPONENT FORM.
    !***********************************************************************
    use mod_cea
    implicit none

    ! DUMMY ARGUMENTS
    character(4), intent(in):: Fone
    character(15), intent(in):: Aa
    real(8), intent(in):: Vx(:)
    integer, intent(in):: Npt

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


  subroutine VARFMT(cea, Vx)
    !***********************************************************************
    ! SET DECIMAL PLACES ACCORDING TO NUMBER SIZE FOR F-format IN
    ! VARIABLE format FMT.
    !***********************************************************************
    use mod_cea
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! DUMMY ARGUMENTS
    real(8), intent(in):: Vx(Ncol)
    ! LOCAL VARIABLES
    integer:: i, k
    real(8):: vi

    do i = 1, cea%Npt
       vi = abs(Vx(i))
       k = 2*i + 3
       cea%fmt(k) = '5,'
       if (vi >= 0.99995d0)  cea%fmt(k) = '4,'
       if (vi >= 9.99950d0)  cea%fmt(k) = '3,'
       if (vi >= 99.9950d0)  cea%fmt(k) = '2,'
       if (vi >= 9999.95d0)  cea%fmt(k) = '1,'
       if (vi >= 999999.5d0) cea%fmt(k) = '0,'
    end do
    cea%fmt(29)(2:) = ' '
  end subroutine VARFMT


  subroutine write_plt_file(cea, filename)
    use mod_cea
    implicit none

    type(CEA_Problem), intent(in):: cea(:)
    character(*), intent(in):: filename

    integer:: num_cases, IOPLT
    integer:: i, j, icase

    num_cases = size(cea)

    open(newunit = IOPLT, file = filename, form = 'formatted')

    do icase = 1, num_cases
       if (cea(icase)%Nplt > 0) then
          write(IOPLT, '("#", 2x, 20A12)') (cea(icase)%Pltvar(j), j = 1, cea(icase)%Nplt)
          do i = 1, cea(icase)%Iplt
             write(IOPLT, '(1x, 1p, 20E12.4)') (cea(icase)%Pltout(i, j), j = 1, cea(icase)%Nplt)
          end do
          write(IOPLT, '("#", 2x, 20A12)') (cea(icase)%Pltvar(j), j = 1, cea(icase)%Nplt)
       end if
    end do

    close(IOPLT)

    return
  end subroutine write_plt_file

end module mod_legacy_io
!!$
!!$
!!$  subroutine read_legacy_input(filename, problems)
!!$    character(*), intent(in):: filename
!!$    type(CEA_Problem), dimension(:), intent(out):: problems
!!$
!!$  end subroutine read_legacy_input
!!$
!!$  subroutine write_legacy_output(filename, problems)
!!$    character(*), intent(in):: filename
!!$    type(CEA_Problem), dimension(:), intent(in):: problems
!!$
!!$  end subroutine write_legacy_output
