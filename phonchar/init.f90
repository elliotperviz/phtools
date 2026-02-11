!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! phonchar. Copyright (C) 2022 Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Program to calculate the atomic character
! of phonon eigenvectors obtained from PHONOPY
! ( https://phonopy.github.io/phonopy )
!
!    This file is part of phonchar.
!
!    phonchar is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    phonchar is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with phonchar.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init
  use functions,   only: i2a
  use pars,        only: error_string
  use io_units,    only: input, inp_yaml
  use in_user
  use in_yaml
  use refconf
  use String_mod,  only: String_type
  use omp_lib

  implicit none

  integer :: i, j, ios
  integer(8) :: l
  character(len=256) :: infile, inposcar, line, key, dum
  type(String_type) :: inputfiles
  logical :: have_files, file_exists
  logical :: have_scan, have_groups

  call show_logo

  !$omp parallel
  max_num_othreads = omp_get_num_threads()
  !$omp end parallel
  write(*,'(*(a))') ' Using maximum ', i2a(max_num_othreads), ' OpenMP threads'

  ! ------------------------------------------------------------
  ! Read input filename
  ! ------------------------------------------------------------
  call get_command_argument(1, infile)
  if (trim(infile) == '' .or. trim(infile) == '-h') then
     write(*,'(a)') ' Syntax: phonchar <input file>'
     stop
  end if

  inquire(file=infile, exist=file_exists)
  if (.not. file_exists) then
     write(0,'(*(a))') error_string, 'input file ', trim(infile), ' not found.'
     stop
  end if

  open(unit=input, file=infile, action='READ')

  have_files  = .false.
  have_scan   = .false.
  have_groups = .false.

  ! ------------------------------------------------------------
  ! Parse phonchar.inp
  ! ------------------------------------------------------------
  do
    call read_line(input, line, ios)
    if (ios < 0) exit

    read(line,*) key

    select case (trim(key))

    case ('run_flags')
      read(line,*) key, run_flag(:)
      if (all(run_flag == 0)) then
        stop 'QUIT: all run_flags are zero'
      end if

    case ('files')
      have_files = .true.
      inposcar = ''   ! default: not provided
      read(line,*,iostat=ios) key, inyaml, inposcar
      if (ios /= 0) then
        read(line,*) key, inyaml
      end if

    case ('[scan]')
      if (run_flag(2) == 0) then
        call skip_block(input)
      else
        call read_scan_block(input)
        have_scan = .true.
      end if

    case ('[groups]')
      call read_groups_block(input)
      have_groups = .true.

    case default
      stop 'ERROR: unknown keyword or block: '//trim(key)

    end select
  end do

  close(input)

  ! ------------------------------------------------------------
  ! Validation
  ! ------------------------------------------------------------

  if (run_flag(2) == 1 .and. .not. have_scan) then
    stop 'ERROR: scan block required but missing'
  end if

  if (.not. have_groups) then
    stop 'ERROR: atomic groups block is required'
  end if

  if (.not. have_files) then
    stop 'ERROR: files block is required'
  else
    call init_geometry(inyaml, inposcar)
  end if

contains

  ! ------------------------------------------------------------
  logical function ignorable(line)
    character(*), intent(in) :: line
    character(len=:), allocatable :: t
    t = adjustl(line)
    ignorable = (len_trim(t) == 0) .or. t(1:1) == '#' .or. t(1:1) == '!'
  end function ignorable
  ! ------------------------------------------------------------

  subroutine read_line(unit, line, ios)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: line
    integer, intent(out) :: ios
    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) return
      if (.not. ignorable(line)) return
    end do
  end subroutine read_line
  ! ------------------------------------------------------------

  subroutine skip_block(unit)
    integer, intent(in) :: unit
    character(len=256) :: line
    integer :: ios
    do
      call read_line(unit, line, ios)
      if (ios /= 0) exit
      if (line(1:1) == '[') then
        backspace(unit)
        exit
      end if
    end do
  end subroutine skip_block
  ! ------------------------------------------------------------

  subroutine read_scan_block(unit)
    integer, intent(in) :: unit
    character(len=256) :: line, key
    integer :: ios

    call read_line(unit, line, ios)
    read(line,*) key, dir(:)

    if (dot_product(dir,dir) < tiny(1.d0)) then
      stop 'ERROR: direction vector has zero length'
    end if

    call read_line(unit, line, ios)
    read(line,*) key, rotax(:)

    call read_line(unit, line, ios)
    read(line,*) key, scanang(:)
  end subroutine read_scan_block
  ! ------------------------------------------------------------

  subroutine read_groups_block(unit)
    integer, intent(in) :: unit
    character(len=256) :: line
    integer :: ios, i

    call read_line(unit, line, ios)
    read(line,*) nag
    if (nag < 2) stop 'ERROR: at least two atomic groups required'

    allocate(ag(nag,natmax))

    do i = 1, nag
      call read_line(unit, line, ios)
      read(line,*) natg(i), ag(i,1:natg(i))
      if (natg(i) > natmax) stop 'ERROR: natg exceeds natmax'
    end do
  end subroutine read_groups_block
  ! ------------------------------------------------------------

  subroutine init_geometry(inyaml, inposcar)
    use pars,     only : error_string
    use io_units, only : inp_yaml
    use refconf
    implicit none

    character(len=*), intent(in) :: inyaml, inposcar
    logical :: geo_from_yaml
    logical :: file_exists

    geo_from_yaml = yaml_has_points(inyaml)

    if (geo_from_yaml) then
      write(*,'(a)') ' Reading geometry from YAML'
      call read_geometry_from_yaml(inyaml)
    else
      inquire(file=inposcar, exist=file_exists)
      if (.not. file_exists) then
        write(0,'(*(a))') error_string, 'input file ', trim(inposcar), ' not found.'
        stop
      end if
      write(*,'(a)') ' Reading geometry from POSCAR'
      call read_poscar(inposcar)
    end if
  end subroutine init_geometry 
    
  logical function yaml_has_points(filename)
    use io_units, only : inp_yaml
    implicit none

    character(len=*), intent(in) :: filename
    character(len=256) :: dum
    integer :: ios

    yaml_has_points = .false.

    open(unit=inp_yaml, file=filename, action='READ')

    do
      read(inp_yaml,*,iostat=ios) dum
      if (ios < 0) exit
      if (dum == 'points:') then
        yaml_has_points = .true.
        exit
      end if
    end do

    close(inp_yaml)
  end function yaml_has_points
 
  subroutine read_geometry_from_yaml(inyaml)
    use io_units, only : inp_yaml
    use refconf
    use functions, only : i2a
    implicit none

    character(len=*), intent(in) :: inyaml
    character(len=256) :: dum
    integer :: i

    open(unit=inp_yaml, file=inyaml, action='READ')

    ! ---- natom
    do
      read(inp_yaml,*,iostat=i) dum
      if (i < 0) stop 'reached end of file, natom string not found.'
      if (dum == 'natom:') then
        backspace(inp_yaml)
        read(inp_yaml,*) dum, atoms_UC
        write(*,'(2a)') ' Number of atoms: ', i2a(atoms_UC)
        exit
      end if
    end do

    rewind(inp_yaml)

    ! ---- lattice
    do
      read(inp_yaml,*,iostat=i) dum
      if (i < 0) stop 'reached end of file, lattice parameters not found.'
      if (dum == 'lattice:') then
        do i = 1, 3
          read(inp_yaml,*) dum, dum, side_UC(i,:)
        end do
        exit
      end if
    end do

    rewind(inp_yaml)

    allocate(pos_eq_UC(atoms_UC,3))
    allocate(mass_UC(atoms_UC))

    ! ---- points
    do
      read(inp_yaml,*,iostat=i) dum
      if (i < 0) stop 'reached end of file, points: not found.'
      if (dum == 'points:') then
        do i = 1, atoms_UC
          read(inp_yaml,*)
          read(inp_yaml,*) dum, dum, pos_eq_UC(i,:)
          read(inp_yaml,*) dum, mass_UC(i)
        end do
        exit
      end if
    end do

    close(inp_yaml)
end subroutine read_geometry_from_yaml

end subroutine init

