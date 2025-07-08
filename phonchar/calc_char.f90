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

subroutine calc_char
  use var, only: nqp, nq, eig, freq, nag, ag, natg, version, vec
  use functions
  use refconf 

  implicit none
  integer, external :: lcm
  integer :: i, j, k, ix, l, ncells_tot, i1, i2 ,i3, m, n, sign_vec, mcm
  integer :: ncells(3), q(3), p(3)
  real(8) :: xt, yt, zt
  real(8) :: s(nqp), w(nqp,nq), s1(nag,nqp,nq), maxl(nag)
  real(8), allocatable :: pos_eq_EC(:,:,:)
  complex(8) :: eexp, q_amp
  character(1000) :: frac(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate phonon character based on contribution to unit displacement
! unit displacement is calculated by setting same mode amplitude q_amp
! for all the modes and considering the quadratic displacement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  write(*,'(a)') ' Calculating phonon character... '
  write(*,'(a)') '  # q-point, commensurate supercell'


  s1(:,:,:) = 0.d0
  s(:) = 0.d0
  do k = 1, nqp
    ! the amplitude is the same for all the modes at fixed q
    ! if q=Gamma, the amplitude is real
    if ( vec(k,1)<tiny(1.d0) .and. vec(k,2)<tiny(1.d0) .and. vec(k,3)<tiny(1.d0) ) then
      q_amp = cmplx ( 1.d0, 0.d0, KIND = 8 )
    else
      q_amp = cmplx( 1.17d0, 3.47d0, KIND=8 ) 
    end if

    ! print the q-point as a fraction
    do i = 1, 3
       if ( abs(vec(k,i)) > tiny(1.d0) ) then
          sign_vec = 1
          if ( vec(k,i) < 0.d0 ) sign_vec = -1
          call real_to_rational ( dble(sign_vec)*vec(k,i), p(i), q(i) )
          p(i) = sign_vec*p(i)
       else
          p(i) = 0
          q(i) = 1
       end if
       call print_fraction ( p(i), q(i), frac(i) )
    end do
    write(*,'(2x,8a)',advance='no') i2a(k), ' (', trim(frac(1)), ',', trim(frac(2)), ',', trim(frac(3)), '), '

    ! supercell commensurate with vec(k,:)
    do i = 1, 3
       if ( p(i) == 0  ) then
          mcm = 1
       else
          mcm = lcm(p(i), q(i))
       end if
       ncells(i) = abs(mcm)
    end do

    write(*,'(6a)') i2a(ncells(1)),' x ', i2a(ncells(2)),' x ', i2a(ncells(3))
    ncells_tot = ncells(1)*ncells(2)*ncells(3)

    allocate ( pos_eq_EC(atoms_UC,3,ncells_tot), stat = i )
    if ( i /= 0 ) stop 'Allocation failed for pos_eq_EC'

    ! generate supercell atomic positions according to the reference configuration
    do i = 1, atoms_UC
       n = 0
       do i3 = 1, ncells(3)
          do i2 = 1, ncells(2)
             do i1 = 1, ncells(1)
                n = n + 1
                xt = dble(i1-1)*side_UC(1,1) + dble(i2-1)*side_UC(2,1) + dble(i3-1)*side_UC(3,1)
                yt = dble(i1-1)*side_UC(1,2) + dble(i2-1)*side_UC(2,2) + dble(i3-1)*side_UC(3,2)
                zt = dble(i1-1)*side_UC(1,3) + dble(i2-1)*side_UC(2,3) + dble(i3-1)*side_UC(3,3)
                pos_eq_EC(i,1,n) = pos_eq_UC(i,1) + xt
                pos_eq_EC(i,2,n) = pos_eq_UC(i,2) + yt
                pos_eq_EC(i,3,n) = pos_eq_UC(i,3) + zt
             end do
          end do
       end do
    end do

    do j = 1, nq
  
      do l = 1, nag ! loop over number of atomic groups
        do m = 1, natg(l) ! loop over number of atoms in group l
          do ix = 1, 3 ! loop over Cartesian components

            !$omp parallel do schedule ( dynamic ) private (eexp) reduction(+:s1)
            do n = 1, ncells_tot 
              eexp = exp( dcmplx( 0.d0, dot_product( vec(k,:), pos_eq_EC(ag(l,m),:,n) ) ) )
              ! it is not necessary to consider 2 times the real part of the displacement
              ! to calculate the weight; also, better to consider the quadratic displacement as
              ! the direction is not relevant
              s1(l,k,j) = s1(l,k,j) + realpart((q_amp*eig(k,j,(ag(l,m)-1)*3+ix))*eexp)*realpart((q_amp*eig(k,j,(ag(l,m)-1)*3+ix))*eexp)
            end do

          end do
          s1(l,k,j) = s1(l,k,j)/sqrt(mass_UC(m)*dble(ncells_tot))
        end do
        s(k) = s(k) + s1(l,k,j)
      end do
  
    end do

    deallocate ( pos_eq_EC, stat = i )
    if ( i /= 0 ) stop 'Deallocation failed for pos_eq_EC'

  end do

  if ( nag == 2 ) then
    do k = 1, nqp
      do j = 1, nq
        w(k,j) = (-s1(1,k,j)+s1(2,k,j)) / s(k)
      end do
    end do
  else
    do k = 1, nqp
      do j = 1, nq
        maxl(:) = s1(:,k,j)
        w(k,j) = dble(maxloc(maxl,dim=1))
      end do
    end do
  end if
  
  open(unit=30,file='phchar.dat')
  write(30,'(*(a))') '# phonchar v. ', version
  write(30,'(*(a))') '# q-points: ',i2a(nqp),' ; bands: ',i2a(nq),' ; groups: ',i2a(nag)
  if ( nag > 2 ) then
    write(30,'(a)') '# the weight value is equal to the group label with greatest projection'
  else
    write(30,'(a)') '# smaller weights correspond to largest projections on group 1'
  end if

  write(30,'(a28)') '# q-point, freq[THz], weight'
  do j = 1, nq
    do k = 1, nqp
      write(30,'(a,1x,E12.6,1x,E12.6)') i2a(k), freq(k,j), w(k,j)
    end do
    write(30,*) 
    write(30,*) 
  end do
  close(30)

  write(*,'(a)') ' done.'

  write(*,'(a)') ' Output written in phchar.dat'

  return
end subroutine calc_char

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
  implicit none
  real(8), intent(in) :: x
  integer, intent(out) :: p, q
  integer :: f, gcd
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


  integer function gcd ( a, b )
    implicit none
    integer, intent(in) :: a, b
    integer :: aa, bb, t
    
    aa = a
    bb = b
    do while ( bb /= 0 )
       t = bb
       bb = mod(aa,bb)
       aa = t
    end do
    gcd = abs(aa)
    
    return
  end function gcd

  integer function lcm ( a, b )
    implicit none
    integer, intent(in) :: a, b
    integer :: gcd
    
    lcm = a * b / gcd(a,b)
    
    return
  end function lcm
