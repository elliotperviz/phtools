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
  use functions, only: i2a
  use pars, only: error_string
  use io_units, only: input, inp_yaml
  use in_user
  use in_yaml
  use refconf
  use String_mod, only: String_type
  use omp_lib

  implicit none
  integer :: i, j
  integer(8) :: l
  character(256) :: infile, inposcar, dum
  type(String_type) :: inputfiles
  logical :: file_exists, geo_from_yaml

  call show_logo

  !$omp parallel
  max_num_othreads = omp_get_num_threads()
  !$omp end parallel
  write(*,'(*(a))') ' Using maximum ', i2a(max_num_othreads), ' OpenMP threads' 

  call get_command_argument(1,infile)
  if ( (infile == '-h') .or. (infile == '' ) ) then
     write(*,'(a)') ' Syntax: phonchar <setting file>'
     write(*,*)
     stop
  end if

  inquire(file=infile,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(0,'(*(a))') error_string,'input file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  open(unit=input,file=infile,action='READ')

  write(*,'(2a)') ' Reading settings from file: ', trim(infile)
  read(input, *) run_flag(:)
  if ( run_flag(1) == 0 .and. run_flag(2) == 0 .and. run_flag(3) == 0 ) then
    write(0,*) error_string, 'all run flags are set to 0, no action to perform.'
    write(*,*)
    stop
  end if

  read(input, '(a)') dum
  inputfiles%value = trim(dum)
  inputfiles%Parts = inputfiles%split(inputfiles%value, delim = " ")
  write(inyaml,'(a)') trim(inputfiles%Parts(1)%record)

  write(*,'(*(a))',advance='no') ' Checking if ',trim(inyaml),' contains information on atoms: '
  inquire(file=inyaml,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,*)
    write(0,'(*(a))') error_string,'input file ',trim(inyaml),' not found.'
    write(*,*)
    stop
  end if

  open(unit=inp_yaml,file=inyaml,action='READ')

  geo_from_yaml = .false.
  do
    read(inp_yaml,*,iostat=i) dum
    if ( i<0 ) then
      write(inposcar,'(a)') trim(inputfiles%Parts(2)%record)
      write(*,'(*(a))') 'no, reading input geometry from ', trim(inposcar)
      inquire(file=inposcar,exist=file_exists)
      if ( .not. (file_exists) ) then
        write(0,'(*(a))') error_string, 'input file ',trim(infile),' not found.'
        write(*,*)
        stop
      end if
      call read_poscar(inposcar)
      exit
    else if ( dum == 'points:' ) then
      write(*,'(a)') 'yes'
      geo_from_yaml = .true.
      exit
    end if
  end do

  close(inp_yaml)

  read(input,*) dir(:)
  if ( dot_product(dir(:),dir(:))<tiny(1.d0) ) then
    write(0,*) error_string, 'the specified direction has null length.'
    write(*,*)
    stop
  end if 

  read(input,*) rotax(:)
  read(input,*) scanang(:)

  read(input,*) nag  ! number of atomic groups
  if ( nag < 2 ) then
     write(*,'(a)') ' ERROR: the number of atomic groups must be greater than 1.'
     write(*,*)
     stop
  end if
  allocate ( ag(nag,natmax), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for ag'

  do i = 1, nag
    read(input,*) natg(i) ! number of atoms in i-th group
    if ( natg(i) > natmax ) then
       write(*,'(a)') ' ERROR: the number of atoms in the groups exceeds natmax.'
       write(*,*)
       stop
    end if
    do j = 1, natg(i)
      read(input,*) ag(i,j)    ! atom j in group i
    end do
  end do

  close(input)

  open(unit=inp_yaml,file=inyaml,action='READ')

  if (geo_from_yaml) then 

    do
      read(inp_yaml,*,iostat=i) dum
      if ( i<0 ) then
        write(*,'(a)') 'reached end of file, natom string not found.'
        write(*,*)
        stop 
      else if ( dum == 'natom:' ) then
        backspace(inp_yaml)
        read(inp_yaml,*) dum, atoms_UC
        write(*,'(2a)') ' Number of atoms: ', i2a(atoms_UC)
        exit
      end if
    end do
    call fseek(inp_yaml, 0, 0, i)
    l=ftell(inp_yaml)
  
    do
      read(inp_yaml,*,iostat=i) dum
      if ( i<0 ) then
        write(*,'(a)') 'reached end of file, lattice parameters not found.'
        write(*,*)
        stop 
      else if ( dum == 'lattice:' ) then
        do i = 1, 3
          read(inp_yaml,*) dum, dum, side_UC(i,:)
        end do
        exit
      end if
    end do
    call fseek(inp_yaml, 0, 0, i)
    l=ftell(inp_yaml)
  
    allocate ( pos_eq_UC(atoms_UC,3), stat = i )
    if ( i /= 0 ) stop 'Allocation failed for pos_eq_UC'
    allocate ( mass_UC(atoms_UC), stat = i )
    if ( i /= 0 ) stop 'Allocation failed for mass_UC'
    do
      read(inp_yaml,*,iostat=i) dum
      if ( i<0 ) then
        write(*,'(a)') 'reached end of file, points: not found.'
        write(*,*)
        stop 
      else if ( dum == 'points:' ) then
        do i = 1, atoms_UC
          read(inp_yaml,*) 
          read(inp_yaml,*) dum, dum, pos_eq_UC(i,:)
          read(inp_yaml,*) dum, mass_UC(i)
        end do
        exit
      end if
    end do
  
    close(inp_yaml)

  end if

  return
end subroutine init


