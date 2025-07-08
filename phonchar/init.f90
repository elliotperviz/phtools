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
  use io_units, only: inp_unit, inp_band
  use var, only: nag, ag, natmax, natg
  use refconf
  use functions, only: i2a
  use omp_lib

  implicit none
  integer :: i, j, threads
  integer(8) :: l
  character(256) :: infile, in2file, dum
  logical :: file_exists

  call show_logo

  !$omp parallel
  threads = omp_get_num_threads()
  !$omp end parallel
  write(*,'(*(a))') ' Running on ', i2a(threads), ' OpenMP threads' 

  call get_command_argument(1,infile)
  if ( (infile == '-h') .or. (infile == '' ) ) then
     write(*,'(a)') ' Syntax: phonchar <setting file> <band.yaml file>'
     write(*,*)
     stop
  end if

  inquire(file=infile,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  call get_command_argument(2,in2file)
  if ( (in2file == '-h') .or. (in2file == '' ) ) then
     write(*,'(a)') ' Syntax: phonchar <setting file> <band.yaml file>'
     write(*,*)
     stop
  end if
  inquire(file=infile,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input yaml file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  open(unit=inp_unit,file=infile,action='READ')

  write(*,'(2a)') ' Reading settings from file: ', trim(infile)
  read(inp_unit,*) nag  ! number of atomic groups
  if ( nag < 2 ) then
     write(*,'(a)') ' ERROR: the number of atomic groups must be greater than 1.'
     write(*,*)
     stop
  end if
  allocate ( ag(nag,natmax), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for ag'

  do i = 1, nag
    read(inp_unit,*) natg(i) ! number of atoms in i-th group
    if ( natg(i) > natmax ) then
       write(*,'(a)') ' ERROR: the number of atoms in the groups exceeds natmax.'
       write(*,*)
       stop
    end if

    do j = 1, natg(i)
      read(inp_unit,*) ag(i,j)    ! atom j in group i
    end do
  end do

  close(inp_unit)

  open(unit=inp_band,file=in2file,action='READ')

  do
    read(inp_band,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(a)') 'reached end of file, natom string not found.'
      write(*,*)
      stop 
    else if ( dum == 'natom:' ) then
      backspace(inp_band)
      read(inp_band,*) dum, atoms_UC
      write(*,'(2a)') ' Number of atoms: ', i2a(atoms_UC)
      exit
    end if
  end do
  call fseek(inp_band, 0, 0, i)
  l=ftell(inp_band)

  do
    read(inp_band,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(a)') 'reached end of file, lattice parameters not found.'
      write(*,*)
      stop 
    else if ( dum == 'lattice:' ) then
      do i = 1, 3
        read(inp_band,*) dum, dum, side_UC(i,:)
      end do
      exit
    end if
  end do
  call fseek(inp_band, 0, 0, i)
  l=ftell(inp_band)

  allocate ( pos_eq_UC(atoms_UC,3), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for pos_eq_UC'
  allocate ( mass_UC(atoms_UC), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for mass_UC'
  do
    read(inp_band,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(a)') 'reached end of file, points: not found.'
      write(*,*)
      stop 
    else if ( dum == 'points:' ) then
      do i = 1, atoms_UC
        read(inp_band,*) 
        read(inp_band,*) dum, dum, pos_eq_UC(i,:)
        read(inp_band,*) dum, mass_UC(i)
      end do
      exit
    end if
  end do

  close(inp_band)

  return
end subroutine init
