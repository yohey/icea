module mod_legacy_io
  implicit none

  integer, parameter:: MAX_CHARS = 132 !< Maximum number of characters in a line (record).

  private:: INPUT, INFREE, UTHERM, UTRAN, EFMT, VARFMT

contains

  subroutine count_cases(inp_filename, num_cases)
    character(*), intent(in):: inp_filename
    integer, intent(out):: num_cases
    integer:: io_inp, iostat
    character(MAX_CHARS):: line
    character(4):: head
    logical:: case_ended

    num_cases = 0
    case_ended = .true.

    open(newunit = io_inp, file = trim(inp_filename), status = 'old', form = 'formatted', action = 'read')

    do
       read(io_inp, '(a)', iostat = iostat) line
       if (is_iostat_end(iostat)) exit

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

    close(io_inp)

    if (.not. case_ended) then
       num_cases = num_cases + 1
    end if

    return
  end subroutine count_cases


  subroutine read_legacy_input(cea, filename)
    use mod_types, only: CEA_Problem, IOINP, init_case
    implicit none

    type(CEA_Problem), allocatable, intent(out):: cea(:)
    character(*), intent(in):: filename

    integer:: num_cases
    integer:: icase
    logical:: file_exists

    inquire(file = filename, exist = file_exists)
    if (.not. file_exists) then
       print *, filename, ' DOES NOT EXIST'
       error stop
    end if

    call count_cases(filename, num_cases)

    allocate(cea(num_cases))

    open(newunit = IOINP, file = filename, status = 'old', form = 'formatted', action = 'read')

    do icase = 1, num_cases
       call init_case(cea(icase))

       !! TEMPORARY WORK AROUND TO REPRODUCE KNOWN BUG !!
       if (icase >= 2) then
          cea(icase)%Dens(:) = cea(icase-1)%Dens(:)
       end if
       !!!!!!!!!!!!!!!!! TO BE DELETED !!!!!!!!!!!!!!!!!!

       call read_legacy_case(cea(icase))
    end do

    close(IOINP)

    return
  end subroutine read_legacy_input


  subroutine read_legacy_case(cea)
    use mod_types, only: CEA_Problem
    implicit none

    type(CEA_Problem), intent(inout):: cea

    integer:: iof, i
    logical:: caseOK, readOK

    caseOK = .true.
    readOK = .true.

    cea%Iplt = 0
    cea%Nplt = 0

    cea%legacy_mode = .true.

    call INPUT(cea, readOK, caseOK)

    do iof = 1, cea%Nof
       if (cea%Oxf(iof) == 0. .and. cea%points(1, 1)%B0p(1, 1) /= 0.) then
          do i = 1, cea%Nlm
             if (cea%points(1, 1)%B0p(i, 1) == 0. .or. cea%points(1, 1)%B0p(i, 2) == 0.) then
                write(cea%io_log, '(/, "OXIDANT NOT PERMITTED WHEN SPECIFYING 100% FUEL(main)")')
                caseOK = .false.
             end if
          end do
       end if
    end do

    cea%invalid_case = (.not. caseOK) .or. (.not. readOK)

    return
  end subroutine read_legacy_case


  subroutine INPUT(cea, readOK, caseOK)
    !***********************************************************************
    ! DECIPHER KEYWORDS, LITERAL VARIABLES, & NUMERICAL VARIABLES IN INPUT.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea
    logical, intent(out):: caseOK, readOK

    ! LOCAL VARIABLES
    character(26), parameter:: lc = 'abcdefghijklmnopqrstuvwxyz', uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(15):: cin(maxNgc), cx15
    character(4):: code, cx4
    character(1):: cx1
    character(2):: cx2
    character(3):: cx3
    logical:: eqrats, incd, phi, pltdat, reacts, refl
    integer:: i, ifrmla, ii, in, iv, ix, j, jj, k, lcin(maxNgc), ncin, nmix, iof
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
    cea%Nplt = 0
    cea%SIunit = .true.
    pltdat = .false.

    do
       do
          ! CALL INFREE TO READ DATASET
          call INFREE(readOK, cin, ncin, lcin, dpin, ioinp, cea%io_log)

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
                cea%ensert(i) = cin(i+1)
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
                      do concurrent (iof = 1:maxMix, j = i+1:ncin, lcin(j) == i .and. int(dpin(j)) <= Ncol)
                         cea%points(iof, int(dpin(j)))%Debug = .true.
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

                         if (index(cin(i-1), 'c') > 0) cea%Enth(cea%Nreac) = cea%Enth(cea%Nreac) * cal_to_J
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
             cea%R = R0 / (cal_to_J * 1d3)
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
                   write(cea%io_log, '(/" Pc,BAR =", 7f13.6)') cea%P(1:cea%Np)
                   write(cea%io_log, '(/" Pc/P =", 9f11.4)') cea%Pcp(1:cea%Npp)
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
             do concurrent (i = 1:maxMix, j = 1:cea%Nsk)
                cea%points(1, j)%U1 = mix(j)
             end do

          else if (cx4 == 'mach') then
             cea%Nsk = nmix
             if (nmix > Ncol) then
                cea%Nsk = Ncol
                write(cea%io_log, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'mach1', cea%Nsk
             end if
             do concurrent (i = 1:maxMix, j = 1:cea%Nsk)
                cea%points(1, j)%Mach1 = mix(j)
             end do

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


  subroutine INFREE(readOK, Cin, Ncin, Lcin, Dpin, io_input, io_log)
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
    implicit none

    ! DUMMY ARGUMENTS
    character(15), intent(out):: Cin(:)
    integer, intent(out):: Ncin
    integer, intent(out):: Lcin(:)
    logical, intent(out):: readOK
    real(8), intent(out):: Dpin(:)
    integer, intent(in):: io_input
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
    integer:: i, ich1, j, kcin, nb, nch1, nx, iostat

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
       read(io_input, '(132a1)', iostat = iostat) ch1

       if (iostat /= 0) exit

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
             backspace io_input
             if (nx == 0) Ncin = Ncin - 1
             return
          end if

       else if (Ncin == 1) then
          write(io_log, '(/" FATAL ERROR IN INPUT format (INFREE)")')
          exit
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
                read(cnum, fmtl, iostat = iostat) Dpin(Ncin)
             end if

             if (iostat > 0) then
                if (Cin(Ncin-1)(:4) /= 'case') write(io_log, '(/" WARNING!!  UNACCEPTABLE NUMBER ", a15, " (INFREE)")') Cin(i)
                Lcin(Ncin) = 0
             end if

             Ncin = Ncin + 1
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

    readOK = .false.

    return
  end subroutine INFREE


  subroutine open_legacy_output(io_out, filename)
    use mod_types, only: MAX_FILENAME
    implicit none

    integer, intent(inout):: io_out
    character(*), intent(in), optional:: filename
    integer:: i
    logical:: is_opened

    if (present(filename)) then
       inquire(file = filename, number = i)

       if (i > 0) then
          io_out = i
       else
          open(newunit = io_out, file = filename, status = 'unknown', form = 'formatted')
       end if
    else
       inquire(io_out, opened = is_opened)

       if (.not. is_opened) then
          open(newunit = io_out, status = 'scratch', form = 'formatted')
       end if
    end if

    write(io_out, '(/" *******************************************************************************")')
    write(io_out, '(/, 9x, "NASA-GLENN CHEMICAL EQUILIBRIUM PROGRAM CEA2,", &
                   & " MAY 21, 2004", /19x, "BY  BONNIE MCBRIDE", &
                   & " AND SANFORD GORDON", /5x, &
                   & " REFS: NASA RP-1311, PART I, 1994", &
                   & " AND NASA RP-1311, PART II, 1996")')
    write(io_out, '(/" *******************************************************************************")')

    return
  end subroutine open_legacy_output


  subroutine write_input_log(io_log, io_out)
    implicit none

    integer, intent(in):: io_log, io_out
    character(MAX_CHARS):: line
    integer:: iostat

    rewind io_log

    do
       read(io_log, '(a)', iostat = iostat) line
       if (is_iostat_end(iostat)) exit

       write(io_out, '(a)') trim(line)
    end do

    close(io_log)

    return
  end subroutine write_input_log


  subroutine OUT1(cea, io_out)
    !***********************************************************************
    ! OUT1 WRITES REACTANT AND FUEL-OXIDANT RATIO INFORMATION.
    !
    ! NOTE - ROCKET, SHOCK, AND DETON PROBLEMS HAVE ADDITIONAL OUTPUT.
    !***********************************************************************
    use mod_types, only: CEA_Problem
    implicit none

    type(CEA_Problem), intent(in):: cea
    integer, intent(in):: io_out

    integer:: n
    real(8):: rho, tem, pfuel, phi

    write(io_out, '(" CASE = ", a15)') cea%Case

    if (cea%Moles) then
       write(io_out, '(/13X, "REACTANT", 20x, a11, "      ENERGY", 6x, "TEMP")') '   MOLES   '
       if (.not. cea%SIunit) write(io_out, '(57X, " CAL/MOL ", 6x, "K")')
       if (cea%SIunit) write(io_out, '(57X, "KJ/KG-MOL", 6x, "K")')
    else
       write(io_out, '(/13X, "REACTANT", 20x, a11, "      ENERGY", 6x, "TEMP")') 'WT FRACTION'
       if (.not. cea%SIunit) write(io_out, '(42X, "(SEE NOTE)      CAL/MOL       K  ")')
       if (cea%SIunit) write(io_out, '(42X, "(SEE NOTE)     KJ/KG-MOL      K  ")')
    end if

    do n = 1, cea%Nreac
       write(io_out, '(1x, a8, 4x, a15, 11x, f12.7, f14.3, f11.3)') &
            cea%Fox(n), cea%Rname(n), cea%Pecwt(n), cea%Enth(n) * cea%R, cea%Rtemp(n)
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
          write(io_out, '(/" REACTANT DENSITY=", F8.2, " KG/CU M")') rho
       else
          write(io_out, '(/" REACTANT DENSITY=", F8.4, " G/CC")') rho
       end if
    end if

    write(io_out, '(/" O/F=", F11.5, 2X, "%FUEL=", F10.6, 2X, "R,EQ.RATIO=", F9.6, 2X, "PHI,EQ.RATIO=", F9.6)') &
         cea%Oxfl, pfuel, cea%Eqrat, phi
    return
  end subroutine OUT1


  subroutine OUT2(cea, i_col_end, io_out)
    !***********************************************************************
    ! OUT2 WRITES THERMODYNAMIC PROPERTIES.
    !***********************************************************************
    use mod_constants, only: R0
    use mod_types, only: CEA_Problem, CEA_Point, Ncol
    implicit none

    type(CEA_Problem), intent(in):: cea
    integer, intent(in):: i_col_end
    integer, intent(in):: io_out

    character(4), dimension(30):: fmt
    character(15):: fgi, fh, fp, frh, fs, fu, fc
    integer:: i
    real(8):: pfactor, vnum

    type(CEA_Point), pointer:: p !< current point

    integer:: n_cols_print
    integer, allocatable:: i_cols_print(:)

    real(8), allocatable:: out_cp(:), out_gamma(:), out_Dlvpt(:), out_Dlvtp(:)
    real(8), allocatable:: out_Ppp(:), out_Ssum(:), out_Totn(:), out_Ttt(:), out_Wm(:)
    real(8), allocatable:: out_Sonvel(:)
    real(8), allocatable:: X(:)

    fmt = cea%fmt

    call get_print_cols(cea, i_col_end, i_cols_print)
    n_cols_print = size(i_cols_print)

    allocate(out_cp(n_cols_print))
    allocate(out_gamma(n_cols_print))
    allocate(out_Dlvpt(n_cols_print))
    allocate(out_Dlvtp(n_cols_print))
    allocate(out_Ppp(n_cols_print))
    allocate(out_Ssum(n_cols_print))
    allocate(out_Totn(n_cols_print))
    allocate(out_Ttt(n_cols_print))
    allocate(out_Wm(n_cols_print))
    allocate(out_Sonvel(n_cols_print))
    allocate(X(n_cols_print))

    do i = 1, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       out_cp(i) = p%Cpr * cea%R
       out_gamma(i) = p%Gammas
       out_Dlvpt(i) = p%Dlvpt
       out_Dlvtp(i) = p%Dlvtp
       out_Ppp(i) = p%Ppp
       out_Ssum(i) = p%Ssum
       out_Totn(i) = p%Totn
       out_Ttt(i) = p%Ttt
       out_Wm(i) = p%Wm
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

    fmt(4) = fmt(6)

    ! PRESSURE
    fmt = VARFMT(fmt, out_Ppp, n_cols_print)

    write(io_out, fmt) fp, out_Ppp(:) * pfactor

    ! TEMPERATURE
    fmt(4) = '13'
    fmt(5) = ' '
    fmt(7) = '2,'
    write(io_out, fmt) 'T, K            ', out_Ttt(:)

    ! DENSITY
    do i = 1, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       if (p%Vlm /= 0) X(i) = vnum / p%Vlm
    end do
    call EFMT(fmt(4), frh, X, n_cols_print)

    ! ENTHALPY
    do i = 1, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       X(i) = p%Hsum * cea%R
    end do
    fmt(4) = fmt(6)
    fmt = VARFMT(fmt, X, n_cols_print)
    write(io_out, fmt) fh, X(1:n_cols_print)

    ! INTERNAL ENERGY
    do i = 1, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       X(i) = (p%Hsum - p%Ppp * p%Vlm / R0) * cea%R
    end do
    fmt = VARFMT(fmt, X, n_cols_print)
    write(io_out, fmt) fu, X(1:n_cols_print)

    ! GIBBS ENERGY
    do i = 1, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       X(i) = (p%Hsum - p%Ttt * p%Ssum) * cea%R
    end do
    fmt = VARFMT(fmt, X, n_cols_print)
    write(io_out, fmt) fgi, X(1:n_cols_print)

    ! ENTROPY
    fmt(4) = '13'
    fmt(5) = ' '
    fmt(7) = '4,'
    write(io_out, fmt) fs, out_Ssum(:) * cea%R
    write(io_out, *)

    ! MOLECULAR WEIGHT
    fmt(7) = '3,'
    write(io_out, fmt) 'M, (1/n)        ', out_Wm(:)
    if (.not. cea%Gonly) write(io_out, fmt) 'MW, MOL WT      ', 1 / out_Totn(:)

    ! (DLV/DLP)T
    fmt(7) = '5,'
    if (cea%Eql) write(io_out, fmt) '(dLV/dLP)t      ', out_Dlvpt(:)

    ! (DLV/DLT)P
    fmt(7) = '4,'
    if (cea%Eql) write(io_out, fmt) '(dLV/dLT)p      ', out_Dlvtp(:)

    ! HEAT CAPACITY
    write(io_out, fmt) fc, out_cp(:)

    ! GAMMA(S)
    fmt(7) = '4,'
    write(io_out, fmt) 'GAMMAs          ', out_gamma(:)

    ! SONIC VELOCITY
    fmt(7) = '1,'
    do i = 1, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       p%Sonvel = sqrt(R0 * p%Gammas * p%Ttt / p%Wm)
       out_Sonvel(i) = p%Sonvel
    end do
    write(io_out, fmt) 'SON VEL,M/SEC   ', out_Sonvel(:)

    deallocate(out_cp)
    deallocate(out_gamma)
    deallocate(out_Dlvpt)
    deallocate(out_Dlvtp)
    deallocate(out_Ppp)
    deallocate(out_Ssum)
    deallocate(out_Totn)
    deallocate(out_Ttt)
    deallocate(out_Wm)

    return
  end subroutine OUT2


  subroutine OUT3(cea, i_col_end, io_out)
    !***********************************************************************
    ! OUT3 WRITES MOLE FRACTIONS.
    !***********************************************************************
    use mod_types, only: CEA_Problem, CEA_Point, Ncol
    implicit none

    type(CEA_Problem), intent(inout):: cea
    integer, intent(in):: i_col_end
    integer, intent(in):: io_out

    character(4):: mamo
    logical:: kOK
    integer:: i, k, m, im, kin
    integer, save:: notuse
    real(8):: tra, tem

    type(CEA_Point), pointer:: p

    integer:: n_cols_print
    integer, allocatable:: i_cols_print(:)
    real(8), allocatable:: X(:)

    call get_print_cols(cea, i_col_end, i_cols_print)
    n_cols_print = size(i_cols_print)
    allocate(X(n_cols_print))

    tra = 5.d-6
    if (cea%Trace /= 0.) tra = cea%Trace

    ! MASS OR MOLE FRACTIONS
    if (cea%Massf) then
       mamo = 'MASS'
    else
       mamo = 'MOLE'
    end if

    if (cea%Eql) then
       write(io_out, '(/1x, A4, " FRACTIONS"/)') mamo
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
          do i = 1, n_cols_print
             p => cea%points(cea%iOF, i_cols_print(i))

             if (cea%Massf) then
                tem = cea%Mw(k)
             else
                tem = 1 / p%Totn
             end if
             if (k <= cea%Ng) then
                X(i) = p%En(k) * tem
             else
                if (cea%Prod(k) /= cea%Prod(k-1)) X(i) = 0
                if (p%En(k) > 0) X(i) = p%En(k) * tem
             end if
             if (kOK .and. X(i) >= tra) kin = 1
          end do

          if (kin == 1) then
             if (cea%Trace == 0) then
                write(io_out, '(1x, A15, F9.5, 12F9.5)') cea%Prod(k), X(1:n_cols_print)
             else
                call EFMT(cea%fmt(4), cea%Prod(k), X, n_cols_print)
             end if
             if (cea%Prod(k) == cea%Omit(notuse)) notuse = notuse - 1
          else if (cea%Prod(k) /= cea%Prod(k-1)) then
             notuse = notuse + 1
             cea%Omit(notuse) = cea%Prod(k)
          end if
       end do

!!$       write(0, '(a, i5)') '[DEBUG] OUT3: notuse   = ', notuse
!!$       write(0, '(a, *(a, ", "))') '[DEBUG] OUT3: cea%Omit = ', (trim(cea%Omit(i)), i = 1, notuse)
    end if

    write(io_out, '(/"  * THERMODYNAMIC PROPERTIES FITTED TO", f7.0, "K")') cea%Tg(4)
    if (.not. cea%Short) then
       write(io_out, '(/"    PRODUCTS WHICH WERE CONSIDERED BUT WHOSE ", a4, &
            & " FRACTIONS", /"    WERE LESS THAN", 1pe13.6, &
            & " FOR ALL ASSIGNED CONDITIONS"/)') mamo, tra
       write(io_out, '(5(1x, a15))') (cea%Omit(i), i = 1, notuse)
    end if

    if (.not. cea%Moles) write(io_out, '(/" NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF", &
         & " OXIDANT IN TOTAL OXIDANTS")')

    deallocate(X)

    return
  end subroutine OUT3


  subroutine OUT4(cea, i_col_end, io_out)
    !***********************************************************************
    ! OUT4 WRITES TRANSPORT PROPERTIES.
    !***********************************************************************
    use mod_types, only: CEA_Problem, CEA_Point, Ncol
    implicit none

    type(CEA_Problem), intent(in):: cea
    integer, intent(in):: i_col_end
    integer, intent(in):: io_out

    integer:: i
    character(4), dimension(30):: fmt
    character(15):: fc

    type(CEA_Point), pointer:: p !< current point

    integer:: n_cols_print
    integer, allocatable:: i_cols_print(:)

    real(8), allocatable:: out_Coneql(:), out_Confro(:), out_Cpeql(:), out_Cpfro(:), out_Preql(:), out_Prfro(:), out_Vis(:)

    fmt = cea%fmt

    if (cea%SIunit) then
       fc = 'Cp, KJ/(KG)(K)'
    else
       fc = 'Cp, CAL/(G)(K)'
    end if

    call get_print_cols(cea, i_col_end, i_cols_print)
    n_cols_print = size(i_cols_print)

    allocate(out_Coneql(n_cols_print))
    allocate(out_Confro(n_cols_print))
    allocate(out_Cpeql(n_cols_print))
    allocate(out_Cpfro(n_cols_print))
    allocate(out_Preql(n_cols_print))
    allocate(out_Prfro(n_cols_print))
    allocate(out_Vis(n_cols_print))

    do i = 1, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       out_Coneql(i) = p%Coneql
       out_Confro(i) = p%Confro
       out_Cpeql(i) = p%Cpeql
       out_Cpfro(i) = p%Cpfro
       out_Preql(i) = p%Preql
       out_Prfro(i) = p%Prfro
       out_Vis(i) = p%Vis
    end do

    write(io_out, *)
    write(io_out, '(" TRANSPORT PROPERTIES (GASES ONLY)")')
    if (cea%SIunit) then
       write(io_out, '("   CONDUCTIVITY IN UNITS OF MILLIWATTS/(CM)(K)"/)')
    else
       write(io_out, '("   CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)"/)')
    end if

    ! TRANSPORT PROPERTIES
    fmt(4) = fmt(6)

    fmt = VARFMT(fmt, out_Vis, n_cols_print)

    write(io_out, fmt) 'VISC,MILLIPOISE', out_Vis(:)

    fmt(4) = '13'
    fmt(5) = ' '
    fmt(7) = '4,'

    if (cea%Eql) then
       write(io_out, '(/"  WITH EQUILIBRIUM REACTIONS"/)')
       ! SPECIFIC HEAT
       write(io_out, fmt) fc, out_Cpeql(:)
       ! CONDUCTIVITY
       write(io_out, fmt) 'CONDUCTIVITY    ', out_Coneql(:)
       ! PRANDTL NUMBER
       write(io_out, fmt) 'PRANDTL NUMBER  ', out_Preql(:)
    end if

    write(io_out, '(/"  WITH FROZEN REACTIONS"/)')
    ! SPECIFIC HEAT
    write(io_out, fmt) fc, out_Cpfro(:)
    ! CONDUCTIVITY
    write(io_out, fmt) 'CONDUCTIVITY    ', out_Confro(:)
    ! PRANDTL NUMBER
    write(io_out, fmt) 'PRANDTL NUMBER  ', out_Prfro(:)

    deallocate(out_Coneql)
    deallocate(out_Confro)
    deallocate(out_Cpeql)
    deallocate(out_Cpfro)
    deallocate(out_Preql)
    deallocate(out_Prfro)
    deallocate(out_Vis)

    return
  end subroutine OUT4


  subroutine get_print_cols(cea, i_col_end, i_cols_print)
    use mod_types, only: CEA_Problem, Ncol
    implicit none

    type(CEA_Problem), intent(in):: cea
    integer, intent(in):: i_col_end
    integer, allocatable, intent(out):: i_cols_print(:)
    integer:: i, i_col_start, n_cols_print, n_cols_fixed

    if (cea%Rkt) then
       if (cea%Fac) then
          n_cols_fixed = 3
       else
          n_cols_fixed = 2
       end if
    else
       n_cols_fixed = 0
    end if

    i_col_start = ((i_col_end - n_cols_fixed - 1) / (Ncol - n_cols_fixed)) * (Ncol - n_cols_fixed) + n_cols_fixed + 1
    n_cols_print = i_col_end - i_col_start + n_cols_fixed + 1

    allocate(i_cols_print(n_cols_print))

    i_cols_print(:) = -1
    do concurrent (i = 1:n_cols_fixed)
       i_cols_print(i) = i
    end do
    do concurrent (i = i_col_start:i_col_end)
       i_cols_print(i - i_col_start + n_cols_fixed + 1) = i
    end do

    return
  end subroutine get_print_cols


  subroutine REACT(cea)
    !***********************************************************************
    ! READ AND PROCESS REACTANT RECORDS.  CALLED FROM SUBROUTINE INPUT.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(6):: date
    character(2):: el(5)
    character(15):: sub
    integer:: i, icf, ifaz, ifrmla, itot, j, jj, kk, kr, l, n, nall, nint, nj, ntgas, ntot
    logical:: fuel, rcoefs, wdone(2)
    logical:: hOK
    real(8):: bb(5), dat(35), dift, eform, pcwt, rcf(9, 3), rm, T1, T2
    real(8):: T1save, T2save
    integer:: io_thermo
    logical:: is_opened

    real(8):: B0p(maxEl, 2)

    wdone = .false.
    cea%Wp = 0
    cea%Hpp = 0
    cea%Vpls = 0
    cea%Vmin = 0
    cea%Am = 0
    cea%Rh = 0
    cea%Elmt = ' '
    dat = 0

    B0p = 0

    ! IF OXIDANT, KR = 1
    ! IF FUEL, KR = 2
    do n = 1, cea%Nreac
       hOK = .false.
       T1save = 20000
       T2save = 0
       rcoefs = .true.

       if (cea%Energy(n) == 'lib' .or. cea%Rnum(n, 1) == 0) then
          cea%Tt = cea%Rtemp(n)
          open(newunit = io_thermo, file = cea%filename_thermo_lib, status = 'old', form = 'unformatted', action = 'read')
          read(io_thermo) cea%Tg, ntgas, ntot, nall

          do itot = 1, nall
             if (itot <= ntot) then
                icf = 3
                if (itot > ntgas) icf = 1
                read(io_thermo) sub, nint, date, (el(j), bb(j), j = 1, 5), ifaz, T1, T2, rm, ((rcf(i, j), i = 1, 9), j = 1, icf)
             else
                read(io_thermo) sub, nint, date, (el(j), bb(j), j = 1, 5), ifaz, T1, T2, rm, eform
                if (nint > 0) read(io_thermo) ((rcf(i, j), i = 1, 9), j = 1, nint)
             end if

             if (sub == cea%Rname(n) .or. sub == '*' // cea%Rname(n)) then
                if (nint == 0) then
                   rcoefs = .false.
                   hOK = .true.
                   cea%Enth(n) = eform * 1000 / R0

                   if (cea%Tt == 0) then
                      cea%Tt = T1
                      cea%Rtemp(n) = T1
                   else
                      dift = abs(cea%Tt - T1)
                      if (dift > 1) then
                         if (dift > 10) then
                            if (cea%legacy_mode) write(cea%io_log, '(/" REACTANT ", a15, "HAS BEEN DEFINED FOR THE TEMPERATURE", &
                                 & f8.2, "K ONLY."/" YOUR TEMPERATURE ASSIGNMENT", f8.2, &
                                 & " IS MORE THAN 10 K FROM THIS VALUE. (REACT)")') cea%Rname(n), T1, cea%Tt
                            cea%Nlm = 0
                            hOK = .false.
                            close(io_thermo)
                            return
                         else
                            if (cea%legacy_mode) write(cea%io_log, '(/" NOTE! REACTANT ", a15, "HAS BEEN DEFINED FOR ", &
                                 & "TEMPERATURE", f8.2, "K ONLY."/" YOUR TEMPERATURE ASSIGNMENT", &
                                 & f8.2, " IS NOT = BUT <10 K FROM THIS VALUE. (REACT)")') cea%Rname(n), T1, cea%Tt
                            cea%Tt = T1
                            cea%Rtemp(n) = T1
                         end if
                      end if
                   end if
                else
                   if (ifaz <= 0) then
                      T1save = min(T1save, 0.8d0 * cea%Tg(1))
                      T2save = max(T2save, 1.2d0 * T2)
                   else
                      T1save = min(T1save, T1 - 0.001d0)
                      T2save = max(T2save, T2 + 0.001d0)
                   endif
                   if (T1save < cea%Tt .and. cea%Tt < T2save) hOK = .true.
                end if

                do j = 1, 5
                   if (bb(j) == 0) exit
                   cea%Nfla(n) = j
                   cea%Ratom(n, j) = el(j)
                   cea%Rnum(n, j) = bb(j)
                end do

                if (cea%Tt == 0) then
                   if (.not. cea%Hp) go to 50
                   if (cea%legacy_mode) write(cea%io_log, '(/" TEMPERATURE MISSING FOR REACTANT NO.", i2, "(REACT)")') n
                   cea%Nlm = 0
                   close(io_thermo)
                   return
                end if

                if (rcoefs .and. hOK) then
                   cea%Tln = log(cea%Tt)
                   l = 1
                   if (ifaz <= 0) then
                      if (cea%Tt > cea%Tg(2)) l = 2
                      if (cea%Tt > cea%Tg(3)) l = 3
                   end if

                   cea%Enth(n) = (((((rcf(7, l)/5) * cea%Tt + rcf(6, l) / 4) * cea%Tt + rcf(5, l) / 3) * cea%Tt &
                        + rcf(4, l) / 2) * cea%Tt + rcf(3, l)) * cea%Tt - rcf(1, l) / cea%Tt + rcf(2, l) * cea%Tln + rcf(8, l)

                   if (cea%Vol .and. ifaz <= 0) cea%Enth(n) = cea%Enth(n) - cea%Tt
                end if

                if (hOK) go to 50
             end if
          end do

          close(io_thermo)

          if (.not. hOK) then
             if (cea%legacy_mode) write(cea%io_log, '(/" YOUR ASSIGNED TEMPERATURE", f8.2, "K FOR ", a15, /, &
                  & "IS OUTSIDE ITS TEMPERATURE RANGE", f8.2, " TO", f9.2, "K (REACT)")') cea%Tt, cea%Rname(n), T1save, T2save
             cea%Energy(n) = ' '
             cea%Nlm = 0
             return
          endif
       end if

50     continue

       inquire(io_thermo, opened = is_opened)
       if (is_opened) close(io_thermo)

       ifrmla = cea%Nfla(n)
       if (cea%Fox(n)(:1) == 'f') then
          fuel = .true.
          kr = 2
          cea%Fox(n) = 'FUEL'
       else if (cea%Fox(n)(:4) == 'name') then
          fuel = .true.
          kr = 2
          cea%Fox(n) = 'NAME'
       else
          kr = 1
          cea%Fox(n) = 'OXIDANT'
       end if

       dat = 0

       ! STORE ATOMIC SYMBOLS IN ELMT ARRAY.
       ! CALCULATE MOLECULAR WEIGHT.
       ! TEMPORARILY STORE ATOMIC VALENCE IN X.
       rm = 0
       do jj = 1, ifrmla
          do j = 1, maxEl
             nj = j
             if (cea%Elmt(j) == ' ') exit
             if (cea%Ratom(n, jj) == cea%Elmt(j)) go to 80
          end do

          cea%Nlm = nj
          cea%Elmt(j) = cea%Ratom(n, jj)

80        do kk = 1, 100
             if (atomic_symbol(kk) == cea%Ratom(n, jj)) then
                rm = rm + cea%Rnum(n, jj) * atomic_mass(kk)
                cea%Atwt(j) = atomic_mass(kk)
                cea%X(j) = atomic_valence(kk)
                dat(j) = dat(j) + cea%Rnum(n, jj)
                go to 100
             end if
          end do

          if (cea%legacy_mode) write(cea%io_log, '(/1x, a2, " NOT FOUND IN BLOCKDATA (REACT)")') cea%Ratom(n, jj)
          cea%Nlm = 0
          return
100    end do

       if (cea%Pecwt(n) < 0) then
          cea%Pecwt(n) = 0
          if (.not. cea%Moles .and. .not. wdone(kr)) then
             wdone(kr) = .true.
             cea%Pecwt(n) = 100.
             if (cea%legacy_mode) write(cea%io_log, '(/" WARNING!!  AMOUNT MISSING FOR REACTANT", i3, ".", &
                  & /" PROGRAM SETS WEIGHT PERCENT = 100. (REACT)")') n
          else
             if (cea%legacy_mode) write(cea%io_log, '(/" AMOUNT MISSING FOR REACTANT NO.", i2, "(REACT)")') n
             cea%Nlm = 0
             return
          end if
       end if

       ! ADD CONTRIBUTIONS TO WP(K), HPP(K), AM(K), AND B0P(I, K)
       if (cea%Pecwt(n) > 0) wdone(kr) = .true.
       pcwt = cea%Pecwt(n)
       if (cea%Moles) pcwt = pcwt * rm
       cea%Wp(kr) = cea%Wp(kr) + pcwt

       if (rm <= 0) then
          cea%Nlm = 0
          return
       else
          cea%Hpp(kr) = cea%Hpp(kr) + cea%Enth(n) * pcwt / rm
          cea%Am(kr) = cea%Am(kr) + pcwt / rm
          if (cea%Dens(n) /= 0) then
             cea%Rh(kr) = cea%Rh(kr) + pcwt / cea%Dens(n)
          else
             cea%Rh = 0
          end if

          do concurrent (j = 1:cea%Nlm)
             B0p(j, kr) = dat(j) * pcwt / rm + B0p(j, kr)
          end do

          cea%Rmw(n) = rm
       end if

    end do

    if (.not. fuel) then
       ! 100 PERCENT OXIDANT, SWITCH INDICES
       cea%Fox = ' '
       cea%Wp(2) = cea%Wp(1)
       cea%Wp(1) = 0
       cea%Hpp(2) = cea%Hpp(1)
       cea%Am(2) = cea%Am(1)
       cea%Am(1) = 0
       B0p(:, 2) = B0p(:, 1)
    end if

    if (cea%Nlm /= 0) then
       ! NORMALIZE HPP(KKR), AM(KR), B0P(I, KR), AND PECWT(N).
       ! CALCULATE V+(KR), AND V-(KR)
       do kr = 1, 2
          if (cea%Wp(kr) /= 0) then
             cea%Hpp(kr) = cea%Hpp(kr) / cea%Wp(kr)
             cea%Am(kr) = cea%Wp(kr) / cea%Am(kr)
             if (cea%Rh(kr) /= 0) cea%Rh(kr) = cea%Wp(kr) / cea%Rh(kr)
             do j = 1, cea%Nlm
                B0p(j, kr) = B0p(j, kr) / cea%Wp(kr)
                if (cea%X(j) < 0) cea%Vmin(kr) = cea%Vmin(kr) + B0p(j, kr) * cea%X(j)
                if (cea%X(j) > 0) cea%Vpls(kr) = cea%Vpls(kr) + B0p(j, kr) * cea%X(j)
             end do
          end if
       end do

       if (.not. cea%Moles) then
          do n = 1, cea%Nreac
             if (cea%Fox(n)(:1) == 'O') then
                cea%Pecwt(n) = cea%Pecwt(n) / cea%Wp(1)
             else
                cea%Pecwt(n) = cea%Pecwt(n) / cea%Wp(2)
             end if
          end do
       end if

       if (cea%legacy_mode .and. .not. cea%Short) then
          if (cea%Moles) then
             write(cea%io_log, '(/4x, "REACTANT", 10x, a7, 3x, "(ENERGY/R),K", 3x, &
                  & "TEMP,K  DENSITY"/, 8x, "EXPLODED FORMULA")') ' MOLES '
          else
             write(cea%io_log, '(/4x, "REACTANT", 10x, a7, 3x, "(ENERGY/R),K", 3x, &
                  & "TEMP,K  DENSITY"/, 8x, "EXPLODED FORMULA")') 'WT.FRAC'
          end if
          do n = 1, cea%Nreac
             write(cea%io_log, '(1x, a1, ": ", a15, f10.6, e15.6, f9.2, f8.4, /8x, 5(2x, a2, f8.5))') &
                  cea%Fox(n), cea%Rname(n), cea%Pecwt(n), cea%Enth(n), cea%Rtemp(n), cea%Dens(n), &
                  (cea%Ratom(n, i), cea%Rnum(n, i), i = 1, cea%Nfla(n))
          end do
       end if

    end if

    do concurrent (i = 1:maxMix, j = 1:Ncol)
       cea%points(i, j)%B0p(:, :) = B0p(:, :)
    end do

    return
  end subroutine REACT


  subroutine RKTOUT(cea, it)
    !***********************************************************************
    ! SPECIAL OUTPUT FOR ROCKET PROBLEMS.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Problem), intent(inout):: cea
    integer, intent(in):: it

    ! LOCAL VARIABLES
    character(4):: exit(11) = 'EXIT'
    character(15):: fi, fiv, fr, z(4)
    integer, save:: i, i23, i46, i57, i68, i79, ione, ixfr, ixfz, j, k, line, ln, mae, mcf, misp, mivac, mmach, mppf, mppj, nex
    real(8):: agv, aw, gc, tem, tra, vaci(Ncol), ww

    type(CEA_Point), pointer:: p !< current point
    type(CEA_Point), pointer:: pfz
    type(CEA_Point), pointer:: p1, p2, p23

    integer:: i_col_end
    integer:: n_cols_print
    integer, allocatable:: i_cols_print(:)

    real(8), allocatable:: out_AeAt(:), out_App(:), X(:)

    i_col_end = cea%ipt
    call get_print_cols(cea, i_col_end, i_cols_print)
    n_cols_print = size(i_cols_print)

    allocate(out_AeAt(n_cols_print))
    allocate(out_App(n_cols_print))
    allocate(X(n_cols_print))

    do i = 1, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       out_AeAt(i) = p%AeAt
       out_App(i) = p%App
    end do

    if (.not. cea%Eql) then
       write(IOOUT, '(/////10x, " THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION")')
       if (cea%Nfz > 1) write(IOOUT, '(33x, "AFTER POINT", i2)') cea%Nfz
    else
       write(IOOUT, '(/////13x, " THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM")')
       if (cea%Iopt /= 0) then
          write(IOOUT, '(/11x, " COMPOSITION DURING EXPANSION FROM FINITE AREA COMBUSTOR")')
       else
          write(IOOUT, '(/10x, " COMPOSITION DURING EXPANSION FROM INFINITE AREA COMBUSTOR")')
       end if
    end if

    pfz => cea%points(cea%iOF, cea%Nfz)
    p1 => cea%points(cea%iOF, 1)
    p2 => cea%points(cea%iOF, 2)

    if (p1%Ttt == cea%T(it)) write(IOOUT, '(25X, "AT AN ASSIGNED TEMPERATURE  ")')

    tem = p1%Ppp * 14.696006d0 / 1.01325d0
    write(IOOUT, '(/1x, a3, " =", f8.1, " PSIA")') 'Pin', tem

    i23 = 2
    if (cea%Iopt > 0) then
       if (cea%Iopt == 1) write(IOOUT, '(" Ac/At =", f8.4, 6x, "Pinj/Pinf =", f10.6)') cea%Subar(1), p2%App
       if (cea%Iopt == 2) write(IOOUT, '(" MDOT/Ac =", f10.3, " (KG/S)/M**2", 6x, "Pinj/Pinf =", f10.6)') cea%Ma, p2%App
       i23 = 3
    end if

    call OUT1(cea, IOOUT)

    cea%fmt(4) = cea%fmt(6)
    nex = n_cols_print - 2
    if (cea%Page1) then
       ione = 0
       i46 = 4
       i57 = 5
       i68 = 6
       i79 = 7
    else
       ione = i23
    end if

    ! PRESSURE RATIOS
    if (cea%Iopt == 0) then
       write(IOOUT, '(/17X, "CHAMBER   THROAT", 11(5X, A4))') (exit(i), i = 1, nex)
       cea%fmt = VARFMT(cea%fmt, out_App, n_cols_print)
       write(IOOUT, cea%fmt) 'Pinf/P         ', out_App(:)
    else
       nex = nex - 1
       write(IOOUT, '(/, 17X, "INJECTOR  COMB END  THROAT", 10(5X, A4))') (exit(i), i = 1, nex)
       X(1) = 1.d0
       do i = 2, n_cols_print
          p => cea%points(cea%iOF, i_cols_print(i))
          X(i) = p1%Ppp / p%Ppp
       end do
       cea%fmt = VARFMT(cea%fmt, X, n_cols_print)
       write(IOOUT, cea%fmt) 'Pinj/P         ', X(:)
    end if

    call OUT2(cea, cea%ipt, IOOUT)

    mppf  = 0
    mppj  = 0
    mmach = 0
    mae   = 0
    mcf   = 0
    mivac = 0
    misp  = 0

    do i = 1, cea%Nplt
       ixfz = index(cea%Pltvar(i)(2:), 'fz')
       ixfr = index(cea%Pltvar(i)(2:), 'fr')

       if (ixfz /= 0 .or. ixfr /= 0) then
          if (cea%Eql) cycle
       else if (.not. cea%Eql) then
          cycle
       end if

       if (cea%Pltvar(i)(:4) == 'pi/p' .or. cea%Pltvar(i)(:3) == 'pip') then
          if (cea%Iopt == 0) mppf = i
          if (cea%Iopt /= 0) mppj = i
       else if (cea%Pltvar(i)(:4) == 'mach') then
          mmach = i
       else if (cea%Pltvar(i)(:2) == 'ae') then
          mae = i
       else if (cea%Pltvar(i)(:2) == 'cf') then
          mcf = i
       else if (cea%Pltvar(i)(:4) == 'ivac') then
          mivac = i
       else if (cea%Pltvar(i)(:3) == 'isp') then
          misp = i
       end if
    end do

    if (cea%SIunit) then
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

    do k = 2, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(k))
       p%Spim = sqrt(2 * R0 * (p1%Hsum - p%Hsum)) / agv
       ! AW IS THE LEFT SIDE OF EQ.(6.12) IN RP-1311, PT I.
       aw = R0 * p%Ttt / (p%Ppp * p%Wm * p%Spim * agv**2)
       if (k == i23) then
          if (cea%Iopt == 0) cea%Cstr = gc * p1%Ppp * aw
          if (cea%Iopt /= 0) cea%Cstr = gc * p1%Ppp / cea%points(cea%iOF, 2)%App * aw
       end if
       vaci(k) = p%Spim + p%Ppp * aw
       p%Vmoc = 0
       if (p%Sonvel /= 0) p%Vmoc = p%Spim * agv / p%Sonvel
    end do

    ! MACH NUMBER
    p1%Vmoc = 0
    p23 => cea%points(cea%iOF, i_cols_print(i23))
    if (p23%Gammas == 0) p23%Vmoc = 0
    cea%fmt(4) = '13'
    cea%fmt(5) = ''
    cea%fmt(7) = '3,'
    write(IOOUT, cea%fmt) 'MACH NUMBER    ', (cea%points(cea%iOF, i_cols_print(i))%Vmoc, i = 1, n_cols_print)
    if (cea%Trnspt) call OUT4(cea, cea%ipt, IOOUT)
    write(IOOUT, '(/" PERFORMANCE PARAMETERS"/)')

    ! AREA RATIO
    cea%fmt(4) = '9x,'
    cea%fmt(i46) = '9x,'
    cea%fmt = VARFMT(cea%fmt, out_AeAt, n_cols_print)
    cea%fmt(5) = ' '
    cea%fmt(i57) = ' '
    write(IOOUT, cea%fmt) 'Ae/At          ', out_AeAt(2:)

    ! C*
    cea%fmt(i57) = '13'
    cea%fmt(i68) = cea%fmt(i68 + 2)
    cea%fmt(i79) = '1,'
    write(IOOUT, cea%fmt) fr, (cea%Cstr, j = 2, n_cols_print)

    ! CF - THRUST COEFICIENT
    cea%fmt(i79) = '4,'
    do i = 2, n_cols_print
       p => cea%points(cea%iOF, i_cols_print(i))
       X(i) = gc * p%Spim / cea%Cstr
    end do
    write(IOOUT, cea%fmt) 'CF             ', X(2:n_cols_print)

    ! VACUUM IMPULSE
    cea%fmt(i57) = '13'
    cea%fmt(i79) = '1,'
    write(IOOUT, cea%fmt) fiv, vaci(2:n_cols_print)

    ! SPECIFIC IMPULSE
    write(IOOUT, cea%fmt) fi, (cea%points(cea%iOF, i_cols_print(i))%Spim, i = 2, n_cols_print)

    if (cea%Nplt > 0) then
       p1 => cea%points(cea%iOF, 1)
       p1%Spim = 0
       p1%AeAt = 0
       p1%Vmoc = 0
       vaci(1) = 0
       X(1) = 0
       do i = ione + 1, n_cols_print
          p => cea%points(cea%iOF, i_cols_print(i))
          if (mppj > 0)  cea%Pltout(i+cea%Iplt-ione, mppj)  = p1%Ppp / p%Ppp
          if (mppf > 0)  cea%Pltout(i+cea%Iplt-ione, mppf)  = p%App
          if (mmach > 0) cea%Pltout(i+cea%Iplt-ione, mmach) = p%Vmoc
          if (mae > 0)   cea%Pltout(i+cea%Iplt-ione, mae)   = p%AeAt
          if (mcf > 0)   cea%Pltout(i+cea%Iplt-ione, mcf)   = X(i)
          if (mivac > 0) cea%Pltout(i+cea%Iplt-ione, mivac) = vaci(i)
          if (misp > 0)  cea%Pltout(i+cea%Iplt-ione, misp)  = p%Spim
       end do
    end if

    write(IOOUT, *)
    cea%fmt(4) = ' '
    cea%fmt(5) = '13'
    cea%fmt(7) = '5,'

    if (cea%Iopt /= 0) then
       cea%fmt(i46) = cea%fmt(8)
       cea%fmt(i57) = cea%fmt(9)
    end if

    if (.not. cea%Eql) then
       if (cea%Massf) then
          write(IOOUT, '(1x, A4, " FRACTIONS"/)') 'MASS'
       else
          write(IOOUT, '(1x, A4, " FRACTIONS"/)') 'MOLE'
          ww = 1 / pfz%Totn
       end if

       ! MOLE (OR MASS) FRACTIONS - FROZEN
       tra = 5.E-6
       if (cea%Trace /= 0) tra = cea%Trace
       line = 0

       do k = 1, cea%Ngc
          if (cea%Massf) ww = cea%Mw(k)
          X(line+1) = pfz%En(k) * ww

          if (X(line+1) >= tra) then
             line = line + 1
             z(line) = cea%Prod(k)
          end if

          if (line == 3 .or. k == cea%Ngc) then
             if (line == 0) then
                call OUT3(cea, cea%ipt, IOOUT)
                deallocate(out_AeAt)
                deallocate(out_App)
                return
             end if
             write(IOOUT, '(1x, 3(a15, f8.5, 3x))') (z(ln), X(ln), ln = 1, line)
             line = 0
          end if
       end do
    end if

    call OUT3(cea, cea%ipt, IOOUT)

    deallocate(out_AeAt)
    deallocate(out_App)

    return
  end subroutine RKTOUT


  subroutine UTHERM(cea, readOK)
    !***********************************************************************
    ! READ THERMO DATA FROM I/O UNIT 7 IN RECORD FORMAT AND WRITE
    ! UNFORMATTED ON NEW FILE.  DATA ARE REORDERED GASES FIRST.
    !
    ! UTHERM IS CALLED FROM SUBROUTINE INPUT.
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
    use mod_types
    implicit none

    ! DUMMY ARGUMENTS
    type(CEA_Problem), intent(inout):: cea
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
    integer:: io_scratch, io_thermo, iostat

    ngl = 0
    ns = 0
    nall = 0
    ifzm1 = 0
    inew = 0
    tinf = 1.d06

    open(newunit = io_scratch, status = 'scratch', form = 'unformatted')

    read(IOINP, '(4f10.3, a10)') tgl, cea%Thdate

    outerLoop: do
       do i = 1, 3
          fill(i) = .true.
          do j = 1, 9
             thermo(j, i) = 0
          end do
       end do
       hform = 0
       tl(1) = 0
       tl(2) = 0

       read(IOINP, '(a15, a65)', iostat = iostat) name, notes
       if (iostat /= 0) exit

       if (name(:3) == 'END' .or. name(:3) == 'end') then
          if (index(name, 'ROD') == 0 .and. index(name, 'rod') == 0) exit
          ns = nall
          cycle
       end if

       read(IOINP, '(i2, 1x, a6, 1x, 5(a2, f6.2), i2, f13.5, f15.3)', iostat = iostat) &
            ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, mwt, hform
       if (iostat /= 0) exit

       write(cea%io_log, '(" ", a15, 2x, a6, e15.6, 2x, a65)') name, date, hform, notes

       ! IF NTL=0, REACTANT WITHOUT COEFFICIENTS
       if (ntl == 0) then
          if (ns == 0) exit
          nall = nall + 1

          read(IOINP, '(2F11.3, i1, 8F5.1, 2x, f15.3)', iostat = iostat) tl, ncoef, expn, hh
          if (iostat /= 0) exit

          thermo(1, 1) = hform
          write(io_scratch) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
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
          read(IOINP, '(2F11.3, i1, 8F5.1, 2x, f15.3)', iostat = iostat) tl, ncoef, expn, hh
          if (iostat /= 0) exit outerLoop
          read(IOINP, '(5d16.8/2d16.8, 16x, 2d16.8)', iostat = iostat) templ
          if (iostat /= 0) exit outerLoop

          if (ifaz == 0 .and. i > 3) then
             iostat = 1
             exit outerLoop
          end if

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
             write(io_scratch) name, ntl, date, (sym(j), fno(j), j = 1, 5), inew, tl, mwt, thermo
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
          ! WRITE COEFFICIENTS ON SCRATCH I/O UNIT
          write(io_scratch) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
       end if
    end do outerLoop

    ! END OF DATA. COPY CONDENSED & REACTANT DATA FROM IO_SCRATCH & ADD TO io_thermo.
    if (iostat <= 0) then
       rewind io_scratch

       open(newunit = io_thermo, file = cea%filename_thermo_lib, status = 'new', form = 'unformatted', action = 'write')

       if (ns == 0) ns = nall
       write(io_thermo) tgl, ngl, ns, nall, cea%Thdate
       ! WRITE GASEOUS PRODUCTS ON io_thermo
       if (ngl /= 0) then
          do i = 1, ns
             read(io_scratch) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
             if (ifaz <= 0) write(io_thermo) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
          end do
       end if
       if (ngl /= nall) then
          ! WRITE CONDENSED PRODUCTS AND REACTANTS ON io_thermo
          rewind io_scratch
          do i = 1, nall
             read(io_scratch) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo
             if (i > ns) then
                write(io_thermo) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, thermo(1, 1)
                if (ntl > 0) write(io_thermo) thermo
             else if (ifaz > 0) then
                write(io_thermo) name, ntl, date, (sym(j), fno(j), j = 1, 5), ifaz, tl, mwt, (thermo(k, 1), k = 1, 9)
             end if
          end do
       end if

       close(io_thermo)

    else
       write(cea%io_log, '(/" ERROR IN PROCESSING thermo.inp AT OR NEAR ", A15, " (UTHERM)")') name
       readOK = .false.

    end if

    close(io_scratch)

    return
  end subroutine UTHERM



  subroutine UTRAN(filename, io_input, io_log, readOK)
    !***********************************************************************
    ! READ TRANSPORT PROPERTIES FORM I/O UNIT 7 IN RECORD FORMAT AND WRITE
    ! UNFORMATTED ON NEW FILE.
    !
    ! UTRAN IS CALLED FROM SUBROUTINE INPUT AFTER A RECORD WITH 'tran'
    ! IN COLUMNS 1-4 HAS BEEN READ.
    !
    ! NOTE:  THIS ROUTINE MAY BE CALLED DIRECTLY  AND USED BY ITSELF TO
    ! PROCESS THE TRANSPORT PROPERTY DATA.
    !***********************************************************************
    implicit none
    ! DUMMY ARGUMENTS
    character(*), intent(in):: filename
    integer, intent(in):: io_input
    integer, intent(in):: io_log
    logical, intent(out):: readOK
    ! LOCAL VARIABLES
    character(16):: tname(2)
    character(1):: vorc
    integer:: i, ic, in, iv, j, k, ncc, nn, ns, nv
    real(8):: cc, tcin(6), trcoef(6, 3, 2), vvl
    integer:: io_scratch, io_transport, iostat


    open(newunit = io_transport, file = trim(filename), status = 'new', form = 'unformatted', action = 'write')
    open(newunit = io_scratch, status = 'scratch', form = 'unformatted')

    ns = 0
    outerLoop: do
       trcoef(:, :, :) = 0

       read(io_input, '(2A16, 2X, A1, I1, A1, I1)') tname, vvl, nv, cc, ncc

       if (tname(1) == 'end' .or. tname(1) == 'LAST') then
          write(io_transport) ns

          rewind io_scratch

          do i = 1, ns
             read(io_scratch, iostat = iostat) tname, trcoef
             if (iostat > 0) exit outerLoop
             write(io_transport) tname, trcoef
          end do

          close(io_transport)
          close(io_scratch)
          return

       else
          ic = 0
          iv = 0
          nn = nv + ncc
          if (nv <= 3 .and. ncc <= 3) then
             do in = 1, nn
                read(io_input, '(1X, A1, 2F9.2, 4E15.8)') vorc, tcin
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

             write(io_scratch) tname, trcoef

             cycle
          end if
       end if

       exit
    end do outerLoop

    write(io_log, '(/" ERROR IN PROCESSING trans.inp AT OR NEAR (UTRAN)", /1X, 2A16)') tname

    readOK = .false.

    close(io_transport)
    close(io_scratch)

    return
  end subroutine UTRAN


  subroutine EFMT(Fone, Aa, Vx, Npt)
    !***********************************************************************
    ! WRITE OUTPUT RECORD WITH NUMERICAL VALUES IN SPECIAL EXPONENT FORM.
    !***********************************************************************
    use mod_types
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


  function VARFMT(fmt_old, Vx, Npt)
    !***********************************************************************
    ! SET DECIMAL PLACES ACCORDING TO NUMBER SIZE FOR F-FORMAT IN
    ! VARIABLE FORMAT FMT.
    !***********************************************************************
    implicit none

    character(4), dimension(30):: VARFMT

    ! DUMMY ARGUMENTS
    character(4), dimension(30), intent(in):: fmt_old
    real(8), intent(in):: Vx(*)
    integer, intent(in):: Npt

    ! LOCAL VARIABLES
    integer:: i, k
    real(8):: vi

    VARFMT = fmt_old

    do i = 1, Npt
       vi = abs(Vx(i))
       k = 2*i + 3
       VARFMT(k) = '5,'
       if (vi >= 0.99995d0)  VARFMT(k) = '4,'
       if (vi >= 9.99950d0)  VARFMT(k) = '3,'
       if (vi >= 99.9950d0)  VARFMT(k) = '2,'
       if (vi >= 9999.95d0)  VARFMT(k) = '1,'
       if (vi >= 999999.5d0) VARFMT(k) = '0,'
    end do

    VARFMT(29)(2:) = ' '

    return
  end function VARFMT


  subroutine write_plt_file(cea, filename)
    use mod_types
    implicit none

    type(CEA_Problem), intent(in):: cea(:)
    character(*), intent(in):: filename

    integer:: num_cases, IOPLT
    integer:: i, j, k, icase, iof, ipt, iplt, max_points
    integer:: mcond, mcondf, mpn, mpnf, mvis
    integer:: mp, mt, mrho, mh, mie, mg, ms, mm, mcp, mgam, mson, mpf, mof, mph, meq, mfa, mmw, mdvt, mdvp
    real(8):: pfactor, vnum

    integer, allocatable:: im(:)
    real(8), allocatable:: Pltout(:, :, :)
    type(CEA_Point), pointer:: p

    num_cases = size(cea)

    allocate(Pltout(num_cases, 500, 20))
    Pltout(:, :, :) = 0

    do icase = 1, num_cases
       if (cea(icase)%Shock) then
          max_points = cea(icase)%Nsk
       else
          max_points = cea(icase)%Nt * cea(icase)%Np
          if (cea(icase)%Rkt) then
             max_points = max_points * (cea(icase)%Npp + cea(icase)%Nsub + cea(icase)%Nsup)
          end if
       end if

       Pltout(icase, :, :) = cea(icase)%Pltout(:, :)

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

       do i = 1, cea(icase)%Nplt
          if (index(cea(icase)%Pltvar(i)(2:), '1') == 0) then
             if (index(cea(icase)%Pltvar(i)(1:), 'dlnt') /= 0) then
                mdvt = i
             else if (index(cea(icase)%Pltvar(i)(1:), 'dlnp') /= 0) then
                mdvp = i
             else if (cea(icase)%Pltvar(i)(:4) == 'pran') then
                if (index(cea(icase)%Pltvar(i)(3:), 'fz') /= 0 .or. index(cea(icase)%Pltvar(i)(3:), 'fr') /= 0) then
                   mpnf = i
                else
                   mpn = i
                end if
             else if (cea(icase)%Pltvar(i)(:4) == 'cond') then
                if (index(cea(icase)%Pltvar(i)(3:), 'fz') /= 0 .or. index(cea(icase)%Pltvar(i)(3:), 'fr') /= 0) then
                   mcondf = i
                else
                   mcond = i
                end if
             else if (cea(icase)%Pltvar(i)(:3) == 'phi') then
                mph = i
             else if (cea(icase)%Pltvar(i)(:2) == 'p ') then
                mp = i
             else if (cea(icase)%Pltvar(i)(:1) == 't') then
                mt = i
             else if (cea(icase)%Pltvar(i)(:3) == 'rho') then
                mrho = i
             else if (cea(icase)%Pltvar(i)(:1) == 'h') then
                mh = i
             else if (cea(icase)%Pltvar(i)(:1) == 'u') then
                mie = i
             else if (cea(icase)%Pltvar(i)(:3) == 'gam') then
                mgam = i
             else if (cea(icase)%Pltvar(i)(:3) == 'son') then
                mson = i
             else if (cea(icase)%Pltvar(i)(:2) == 'g ') then
                mg = i
             else if (cea(icase)%Pltvar(i)(:2) == 's ') then
                ms = i
             else if (cea(icase)%Pltvar(i)(:1) == 'm' .and. cea(icase)%Pltvar(i)(:2) /= 'ma') then
                if (.not. cea(icase)%Gonly .and. cea(icase)%Pltvar(i)(:2) == 'mw') then
                   mmw = i
                else
                   mm = i
                end if
             else if (cea(icase)%Pltvar(i)(:2) == 'cp') then
                mcp = i
             else if (cea(icase)%Pltvar(i)(:3) == 'vis') then
                mvis = i
             else if (cea(icase)%Pltvar(i)(:3) == 'o/f') then
                mof = i
             else if (cea(icase)%Pltvar(i)(:2) == '%f') then
                mpf = i
             else if (cea(icase)%Pltvar(i)(:3) == 'f/a') then
                mfa = i
             else if (cea(icase)%Pltvar(i)(:1) == 'r') then
                meq = i
             end if
          end if
       end do

       allocate(im(cea(icase)%Ngc))
       im(:) = 0

       do k = 1, cea(icase)%Ngc
          do i = 1, cea(icase)%Nplt
             if (cea(icase)%Pltvar(i) == cea(icase)%Prod(k) .or. '*' // cea(icase)%Pltvar(i) == cea(icase)%Prod(k)) then
                im(k) = i
                exit
             end if
          end do
       end do

       if (cea(icase)%Nplt > 0) then
          if (cea(icase)%SIunit) then
             pfactor = 1
             vnum = 1.d05
          else
             pfactor = 1 / 1.01325d0
             vnum = 100
          end if

          do iof = 1, cea(icase)%Nof
             do ipt = 1, max_points
                p => cea(icase)%points(iof, ipt)
                iplt = ipt + (iof - 1) * max_points

                if (mof > 0) Pltout(icase, iplt, mof) = cea(icase)%Oxf(iof)
                if (mpf > 0) Pltout(icase, iplt, mpf) = 100 / (1 + cea(icase)%Oxf(iof))
                if (mph > 0) Pltout(icase, iplt, mph) = -(cea(icase)%Vmin(2) + cea(icase)%Vpls(2)) &
                     / ((cea(icase)%Vpls(1) + cea(icase)%Vmin(1)) * cea(icase)%Oxf(iof))
                if (mfa > 0) Pltout(icase, iplt, mfa) = 1 / cea(icase)%Oxf(iof)
                if (meq > 0) Pltout(icase, iplt, meq) = cea(icase)%Eqrat

                if (mp > 0) Pltout(icase, iplt, mp) = p%Ppp * pfactor
                if (mt > 0) Pltout(icase, iplt, mt) = p%Ttt
                if (mrho > 0 .and. p%Vlm /= 0) Pltout(icase, iplt, mrho) = vnum / p%Vlm
                if (mh > 0) Pltout(icase, iplt, mh) = p%Hsum * cea(icase)%R
                if (mie > 0) Pltout(icase, iplt, mie) = (p%Hsum - p%Ppp * p%Vlm / R0) * cea(icase)%R

                if (mg > 0) Pltout(icase, iplt, mg) = (p%Hsum - p%Ttt * p%Ssum) * cea(icase)%R
                if (mm > 0) Pltout(icase, iplt, mm) = p%Wm
                if (mmw > 0) Pltout(icase, iplt, mmw) = 1 / p%Totn
                if (ms > 0) Pltout(icase, iplt, ms) = p%Ssum * cea(icase)%R
                if (mcp > 0) Pltout(icase, iplt, mcp) = p%Cpr * cea(icase)%R
                if (mgam > 0) Pltout(icase, iplt, mgam) = p%Gammas
                if (mdvt > 0) Pltout(icase, iplt, mdvt) = p%Dlvtp
                if (mdvp > 0) Pltout(icase, iplt, mdvp) = p%Dlvpt

                if (mson > 0) Pltout(icase, iplt, mson) = p%Sonvel

                do k = 1, cea(icase)%Ngc
                   if (k <= cea(icase)%Ng .or. p%En(k) > 0) then
                      if (cea(icase)%Massf) then
                         if (im(k) > 0) Pltout(icase, iplt, im(k)) = p%En(k) * cea(icase)%Mw(k)
                      else
                         if (im(k) > 0) Pltout(icase, iplt, im(k)) = p%En(k) / p%Totn
                      end if
                   end if
                end do

                ! TRANSPORT PROPERTIES
                if (cea(icase)%Trnspt) then
                   if (mvis > 0) Pltout(icase, iplt, mvis) = p%Vis
                   if (mcond > 0) Pltout(icase, iplt, mcond) = p%Coneql
                   if (mpn > 0) Pltout(icase, iplt, mpn) = p%Preql
                   if (mcondf > 0) Pltout(icase, iplt, mcondf) = p%Confro
                   if (mpnf > 0) Pltout(icase, iplt, mpnf) = p%Prfro
                end if
             end do
          end do
       end if

       deallocate(im)
    end do


    open(newunit = IOPLT, file = filename, form = 'formatted')

    do icase = 1, num_cases
       if (cea(icase)%Nplt > 0) then
          write(IOPLT, '("#", 2x, 20A12)') (cea(icase)%Pltvar(j), j = 1, cea(icase)%Nplt)
          do i = 1, cea(icase)%Iplt
             write(IOPLT, '(1x, 1p, 20E12.4)') (Pltout(icase, i, j), j = 1, cea(icase)%Nplt)
          end do
          write(IOPLT, '("#", 2x, 20A12)') (cea(icase)%Pltvar(j), j = 1, cea(icase)%Nplt)
       end if
    end do

    close(IOPLT)

    deallocate(Pltout)

    return
  end subroutine write_plt_file

end module mod_legacy_io
