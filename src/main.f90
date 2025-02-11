program main
  use mod_cea_core
  use mod_cea_types
  use mod_legacy_io
  use mod_io
  implicit none

  type(CEA_Problem), allocatable:: cea(:)

  character(MAX_FILENAME):: inp_filename = ''
  character(MAX_FILENAME):: out_filename = ''
  character(MAX_FILENAME):: plt_filename = ''
  character(MAX_FILENAME):: filename_thermo_lib = ''
  character(MAX_FILENAME):: filename_trans_lib = ''
  logical:: legacy_mode = .false.

  call parse_arguments()

  if (legacy_mode) then
     call read_legacy_input(cea, inp_filename, filename_thermo_lib, filename_trans_lib)

  else
     write(stderr, *) '[ERROR] Not implemented yet.'
     error stop

  end if

  call run_all_cases(cea, out_filename, plt_filename)

  deallocate(cea)

  stop

contains

  subroutine parse_arguments()
    use mod_constants, only: MAX_FILENAME, stderr
    implicit none

    character(MAX_FILENAME-4):: basename
    character(MAX_FILENAME), allocatable:: args(:)
    integer:: i, nargs
    logical:: exists

    nargs = command_argument_count()
    allocate(args(0:nargs))

    do i = 0, nargs
       call get_command_argument(i, args(i))

       if (trim(args(i)) == '-l' .or. trim(args(i)) == '--legacy-mode') legacy_mode = .true.
    end do

    i = 1
    do while (i < nargs)
       select case (args(i))
       case ('-o', '--output')
          if (legacy_mode) then
             write(stderr, *) '[WARNING] Ignoring option in legacy mode:', trim(args(i))
          else
             out_filename = trim(adjustl(args(i+1)))
          end if
          i = i + 1

       case ('-p', '--plt-file')
          if (legacy_mode) then
             write(stderr, *) '[WARNING] Ignoring option in legacy mode:', trim(args(i))
          else
             plt_filename = trim(adjustl(args(i+1)))
          end if
          i = i + 1

       case ('--thermo-lib')
          filename_thermo_lib = trim(adjustl(args(i+1)))
          inquire(file = filename_thermo_lib, exist = exists)
          if (.not. exists) then
             write(stderr, *) '[ERROR] Not exist: ', trim(filename_thermo_lib)
             error stop
          end if
          i = i + 1

       case ('--trans-lib')
          filename_trans_lib = trim(adjustl(args(i+1)))
          inquire(file = filename_trans_lib, exist = exists)
          if (.not. exists) then
             write(stderr, *) '[ERROR] Not exist: ', trim(filename_trans_lib)
             error stop
          end if
          i = i + 1

       case default
          if (legacy_mode) then
             write(stderr, *) '[WARNING] Ignoring argument in legacy mode:', trim(args(i))
          else if (len_trim(inp_filename) == 0) then
             inp_filename = trim(adjustl(args(i)))

             inquire(file = inp_filename, exist = exists)
             if (.not. exists) then
                write(stderr, *) '[ERROR] Input file does not exist: ', trim(inp_filename)
                error stop
             end if
          else
             write(stderr, *) '[WARNING] Multiple input file names. Using only first one.'
          end if

       end select

       i = i + 1
    end do

    if (legacy_mode) then
       write(*, '(//" ENTER INPUT FILE NAME WITHOUT .inp EXTENSION."/  &
            & "   THE OUTPUT FILES FOR LISTING AND PLOTTING WILL HAVE", / &
            & " THE SAME NAME WITH EXTENSIONS .out AND .plt RESPECTIVELY" &
            & //)')

       read(*, '(a)') basename

       inp_filename = trim(basename) // '.inp'
       out_filename = trim(basename) // '.out'
       plt_filename = trim(basename) // '.plt'
    end if

    return
  end subroutine parse_arguments

end program
