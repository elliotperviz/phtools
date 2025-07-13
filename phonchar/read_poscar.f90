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

subroutine read_poscar(inposcar)
  use functions, only: wordcount
  use pars, only: error_string
  use refconf
  use io_units, only: inp_poscar
  use masses
  implicit none
!  integer, external :: wordcount
  integer :: i, j, k, l
  real(8) :: scalevol, x_tmp, y_tmp, z_tmp
  character(1) :: word
  character(256), intent(in) :: inposcar
  character(256) :: dum
  character(2), allocatable :: at_pertype(:)

  open(inp_poscar,file=inposcar,action='read')
  ! comment
  read(inp_poscar,*) 

  ! scale factor
  read(inp_poscar,*) scalevol
  if ( abs(scalevol-1.d0) > tiny(1.d0) ) then
     write(0,*) error_string, 'POSCAR scale factor must be equal to 1.0'
     stop
  end if
  
  ! read UC box
  do i = 1, 3
     read(inp_poscar,*) (side_UC(i,j),j=1,3) ! Ang
  end do

  read(inp_poscar,'(a)') dum
  backspace(inp_poscar)

  atom_types = wordcount(dum)

  allocate ( at_pertype(atom_types), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array at_pertype.'
     stop
  end if

  ! read atomic types
  read(inp_poscar,*) at_pertype(:) ! atom symbols
  
  ! allocate number of atoms for each type
  allocate ( natoms_UC(atom_types), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array natoms_UC.'
     stop
  end if
  
  ! read number of atoms for each type
  read(inp_poscar,*) natoms_UC(:)
  
  ! calculate total number of atoms
  atoms_UC = 0
  do i = 1, atom_types
     atoms_UC = atoms_UC + natoms_UC(i) 
  end do

  ! set atomic masses from phonopy database
  allocate ( mass_UC(atoms_UC), stat = i )
  if ( i /= 0 ) then
    write(0,*) error_string, 'Allocation failed for array mass_UC.'
    stop
  end if
  k = 0
  do i = 1, atom_types
    ! loop over all the atomic symbols to find the corresponding mass
    do l = 1, n_phonopy
      if ( at_pertype(i) == at_phonopy(l) ) then
        do j = 1, natoms_UC(i)
           k = k + 1
           mass_UC(k) = mass_phonopy(l) ! amu
        end do
      end if
    end do
  end do

  ! allocate equilibrium position array
  allocate ( pos_eq_UC(atoms_UC,3), stat = i )
  if ( i /= 0 ) then
     write(0,*) error_string, 'Allocation failed for array pos_eq_UC.'
     stop
  end if

  ! read keyword
  read(inp_poscar,*) word
  if ( word == 'S' .or. word == 's' ) read(inp_poscar,*) word
  
  ! read atomic positions and store them in cartesian
  do i = 1, atoms_UC
     read(inp_poscar,*) (pos_eq_UC(i,j),j=1,3) ! Ang or adim
     if ( word /= 'C' .and. word /= 'c' .and. word /= 'K' .and. word /= 'k' ) then
        x_tmp = side_UC(1,1)*pos_eq_UC(i,1) + side_UC(2,1)*pos_eq_UC(i,2) + side_UC(3,1)*pos_eq_UC(i,3) 
        y_tmp = side_UC(1,2)*pos_eq_UC(i,1) + side_UC(2,2)*pos_eq_UC(i,2) + side_UC(3,2)*pos_eq_UC(i,3) 
        z_tmp = side_UC(1,3)*pos_eq_UC(i,1) + side_UC(2,3)*pos_eq_UC(i,2) + side_UC(3,3)*pos_eq_UC(i,3) 
        pos_eq_UC(i,1) = x_tmp
        pos_eq_UC(i,2) = y_tmp
        pos_eq_UC(i,3) = z_tmp
     end if
  end do
  close(inp_poscar)

  return
end subroutine read_poscar
