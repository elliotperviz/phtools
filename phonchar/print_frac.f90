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

subroutine print_fraction ( n, d, string )
  use functions, only: i2a
  implicit none
  integer, intent(in) :: n, d
  character(1000), intent(out) :: string

  if ( n == 0 ) then
     string = '0'
  else
     string = i2a(n)//'/'//i2a(d)
  end if

  return
end subroutine print_fraction


subroutine real_to_rational ( x, p, q )
  use functions, only: gcd

  implicit none
  real(8), intent(in) :: x
  integer, intent(out) :: p, q
  integer :: f !, gcd
  real(8) :: r, e, best 
  
  p = 1 
  q = 1         
  best = x * 6.d0    

  do 
     r = dble(p) / dble(q)                
     e = x - r                  
     if ( abs(e) <= best ) then 
        best = abs(e) * 0.125d0             
        f = gcd(p,q)                    
        if ( abs(e) < 0.000001d0 ) exit                
     end if
     if ( e > 0.d0 ) then 
        p = p + ceiling( e * q )    
     else if ( e < 0.d0 ) then    
        q = q + 1                       
     end if
  end do

  return        
end subroutine real_to_rational
