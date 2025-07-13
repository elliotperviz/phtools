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

subroutine rotvec(dir, rotax, rotang, rotdir)
  implicit none
  real(8), intent(in) :: dir(3)       ! Input vector to rotate
  real(8), intent(in) :: rotax(3)     ! Rotation axis (not necessarily normalized)
  real(8), intent(in) :: rotang       ! Rotation angle in radians
  real(8), intent(out) :: rotdir(3)   ! Output: rotated vector

  real(8) :: k(3)
  real(8) :: norm, cost, sint, dotkv
  real(8) :: crosskv(3)

  ! Normalize the rotation axis
  norm = sqrt(dot_product(rotax, rotax))
  if ( abs(norm) < tiny(1.d0) ) then
    rotdir = dir
    return
  end if
  k = rotax / norm

  ! Compute dot and cross products
  dotkv = dot_product(k, dir)
  crosskv(1) = k(2)*dir(3) - k(3)*dir(2)
  crosskv(2) = k(3)*dir(1) - k(1)*dir(3)
  crosskv(3) = k(1)*dir(2) - k(2)*dir(1)

  cost = cos(rotang)
  sint = sin(rotang)

  ! Rodrigues' rotation formula
  rotdir = dir * cost + crosskv * sint + k * dotkv * (1.0d0 - cost)

end subroutine rotvec

