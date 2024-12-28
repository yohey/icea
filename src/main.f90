program main
  use mod_cea
  use mod_types
  use mod_legacy_io
  implicit none

  type(CEA_Problem), allocatable:: cea(:)

  character(MAX_FILENAME):: inp_filename
  character(MAX_FILENAME):: out_filename
  character(MAX_FILENAME):: plt_filename
  logical:: legacy_mode

  call parse_arguments(inp_filename, out_filename, plt_filename, legacy_mode)

  if (legacy_mode) then
     call read_legacy_input(cea, inp_filename)

     call run_all_cases(cea, out_filename, plt_filename)

     deallocate(cea)

  else
     write(stderr, *) '[ERROR] Not implemented yet.'
     error stop

  end if

  stop
end program main


subroutine parse_arguments(inp_filename, out_filename, plt_filename, legacy_mode)
  use mod_constants, only: MAX_FILENAME, stderr
  implicit none

  character(MAX_FILENAME), intent(out):: inp_filename
  character(MAX_FILENAME), intent(out):: out_filename
  character(MAX_FILENAME), intent(out):: plt_filename
  logical, intent(out):: legacy_mode

  character(MAX_FILENAME-4):: basename
  character(MAX_FILENAME), allocatable:: args(:)
  integer:: i, nargs
  logical:: exists

  inp_filename = ''
  out_filename = ''
  plt_filename = ''
  legacy_mode = .false.

  nargs = command_argument_count()
  allocate(args(0:nargs))

  do i = 0, nargs
     call get_command_argument(i, args(i))

     if (trim(args(i)) == '-l' .or. trim(args(i)) == '--legacy-mode') legacy_mode = .true.
  end do


  if (legacy_mode) then
     if (nargs > 2) then
        write(stderr, *) '[WARNING] Command line arguments are ignored in legacy mode.'
     end if

     write(*, '(//" ENTER INPUT FILE NAME WITHOUT .inp EXTENSION."/  &
          & "   THE OUTPUT FILES FOR LISTING AND PLOTTING WILL HAVE", / &
          & " THE SAME NAME WITH EXTENSIONS .out AND .plt RESPECTIVELY" &
          & //)')

     read(*, '(a)') basename

     inp_filename = trim(basename) // '.inp'
     out_filename = trim(basename) // '.out'
     plt_filename = trim(basename) // '.plt'

     return
  end if


  i = 0
  do while (i < nargs)
     select case (args(i))
     case ('-o', '--output')
        out_filename = trim(adjustl(args(i+1)))

     case ('-p', '--plt-file')
        plt_filename = trim(adjustl(args(i+1)))

     case default
        if (len_trim(inp_filename) == 0) then
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
  end do

  return
end subroutine parse_arguments
