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
! If used for production, you should cite
! Phys. Rev. B 103, 035406 (2021)
! https://doi.org/10.1103/PhysRevB.103.035406
! where the formulation is reported in section V "Atomic character of the phonon modes" 
! of the supplemental material.
!
! Only the real part of the eigenvector is considered, as the purpose
! is to analyse the contribution to the atomic motions.
! 
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

module functions
contains

  function i2a(i) result(out)
    character(:), allocatable :: out
    integer(4), intent(in) :: i
    character(range(i)+2) :: x

    write(x,'(i0)') i

    out = trim(x)
  end function i2a

end module functions

module pars
  ! strings
  character(6), parameter :: version = 'v1.0.1', output_format = 'g22.14', screen_format = 'f12.4'
  character(7), parameter :: error_string = 'ERROR: '
  character(9), parameter :: warning_string = 'WARNING: '
end module pars

module var
  ! global parameters
  integer, parameter :: natmax = 10000
  character(3), parameter :: version = '2.3'
  character(8), parameter :: progname = 'PHONCHAR'
  ! from band file
  integer, save :: nqp, nq, nag
  integer, save :: natg(natmax)
  integer, save, allocatable :: ag(:,:)
  real(8), save, allocatable :: freq(:,:), vec(:,:)
  complex(8), save, allocatable :: eig(:,:,:)
end module var

module refconf
  integer, save :: atoms_UC
  real(8), save :: side_UC(3,3)
  real(8), save, allocatable :: mass_UC(:), pos_eq_UC(:,:)
end module refconf

module io_units
  integer, parameter :: inp_unit          = 10
  integer, parameter :: inp_band          = 11
end module io_units
