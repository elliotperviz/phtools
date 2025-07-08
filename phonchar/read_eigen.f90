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

subroutine read_eigen
  use io_units, only: inp_band
  use var, only: nqp, nq, eig, freq, vec
  use functions, only: i2a
  use refconf, only: atoms_UC
                 
  implicit none
  integer(4) :: i, j, k, ix
  integer(8) :: l
  real(8) :: tmpr, tmpi
  character(200) :: dum, in2file

  call get_command_argument(2,in2file)

  open(unit=inp_band,file=in2file, action='READ')

  write(*,'(2a)') ' Reading eigenvectors and frequencies from file: ', trim(in2file)
  write(*,'(a)', advance='no') ' Checking if input file contains eigenvectors: '
  do
    read(inp_band,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(a)') 'reached end of file, no eigenvectors found.'
      write(*,*)
      stop 
    else if ( dum == 'eigenvector:' ) then
      write(*,'(a)') 'yes'
      exit
    end if
  end do
  call fseek(inp_band, 0, 0, i)
  l=ftell(inp_band)
  do
    read(inp_band,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(a)') 'reached end of file, nqpoints string not found.'
      write(*,*)
      stop 
    else if ( dum == 'nqpoint:' ) then
      backspace(inp_band)
      read(inp_band,*) dum, nqp
      exit
    end if
  end do

  write(*,'(2a)') ' Number of q-points: ', i2a(nqp)
  call fseek(inp_band, 0, 0, i)
  l=ftell(inp_band)

  nq=atoms_UC*3

  allocate ( eig(nqp,nq,nq), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for eig'
  allocate ( freq(nqp,nq), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for freq'
  allocate ( vec(nqp,3), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for freq'

  do
    read(inp_band,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(*(a))') 'reached end of file, the input file ', trim(in2file), ' is not complete.'
      write(*,*)
      stop 
    else if ( dum == 'phonon:' ) then
      exit
    end if
  end do
  
  do i = 1, nqp
    read(inp_band,*) dum, dum, dum, vec(i,:)
    do k = 1,2
      read(inp_band,*)
    end do

    do j = 1, nq
      read(inp_band,*)
      read(inp_band,*) dum, freq(i,j)
      read(inp_band,*)

      do k = 1, atoms_UC
        read(inp_band,*)
        do ix = 1, 3
          read(inp_band,*) dum, dum, tmpr, tmpi
          eig(i,j,(k-1)*3+ix) = cmplx(tmpr, tmpi, 8) 
        end do
      end do

    end do
    read(inp_band,*)
  end do

  close(inp_band)

  write(*,*) 'Reading input file done.'

  return
end subroutine read_eigen
