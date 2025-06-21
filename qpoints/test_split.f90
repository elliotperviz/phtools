module String_mod

    use iso_fortran_env, only: IK => int32
    implicit none
    public

    type :: CharVec_type
        character(:)    , allocatable   :: record
    end type CharVec_type

    type :: String_type
        character(:)      , allocatable   :: value          !< The string value.
        type(CharVec_type), allocatable   :: Parts(:)       !< The string parts.
        integer(IK)                       :: nPart = 0_IK   !< The number of parts in the string.
    contains
        procedure, nopass :: split
    end type String_type

contains

    function split(string,delim,npart) result(Parts)

        implicit none
        character(len=*)    , intent(in)            :: string, delim

        integer(IK)         , intent(out), optional :: npart

        type(CharVec_type)  , allocatable           :: Parts(:)
        integer(IK)         , allocatable           :: PartEnd(:)
        integer(IK)         , allocatable           :: PartBegin(:)
        integer(IK)                                 :: dlmlenMinusOne
        integer(IK)                                 :: strlen, dlmlen, npartMax, ipart, ibeg, iend, i
        logical                                     :: npartIsPresent

        dlmlen = len(delim)
        strlen = len(string)
        npartIsPresent = present(npart)

        ! if dlm is empty, return the whole string split character by character

        if (dlmlen==0_IK) then
            allocate(Parts(strlen))
            do ipart = 1, strlen
                Parts(ipart)%record = string(ipart:ipart)
            end do
            if (npartIsPresent) npart = strlen
            return
        end if

        npartMax = 1_IK + strlen / dlmlen ! There can be at most strlen + 1 splits
        allocate(PartBegin(npartMax), PartEnd(npartMax)) ! This will contain the beginning and the ends of the splits.
        dlmlenMinusOne = dlmlen - 1_IK

        ibeg = 0_IK
        ipart = 1_IK
        PartBegin(ipart) = 1_IK
        loopParseString: do

            ibeg = ibeg + 1_IK
            iend = ibeg + dlmlenMinusOne

            if (strlen<iend) then ! the remaining part of the string is shorter than the delim
                PartEnd(ipart) = strlen
                exit loopParseString
            elseif ( string(ibeg:iend) == delim ) then
                PartEnd(ipart) = ibeg - 1_IK
                ipart = ipart + 1_IK
                PartBegin(ipart) = iend + 1_IK
                ibeg = iend
            end if

        end do loopParseString

        allocate(Parts(ipart))
        do i = 1, ipart
            Parts(i)%record = string(PartBegin(i):PartEnd(i))
        end do
        if (present(npart)) npart = ipart

        deallocate(PartBegin, PartEnd)

    end function split
    
end module String_mod

program test

use, intrinsic :: iso_fortran_env, only: output_unit
use String_mod, only: String_type
implicit none
type(String_type) :: string

string%value = "command file/path # prova"
string%Parts = string%split(string = string%value, delim = "#")
write(output_unit,"(*(g0,:,' '))") "fname =", string%Parts(1)%record

end program test
