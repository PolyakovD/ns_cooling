!This module implements tridiagonal matrix algorithm
module tridiagonal 

    implicit none

    type, public :: tridiagonalMatrix
        logical, private :: isInit = .false.
        integer, private :: mSize
        real(8), private, allocatable, dimension(:,:) :: table
        real(8), private, allocatable, dimension(:) :: rightPart
    end type tridiagonalMatrix
    
contains
 
    subroutine tMConstructor(tM, newTable, newRightPart, newSize)
        type(tridiagonalMatrix) :: tM
        integer :: newSize
        real(8) :: newTable(3, newSize)
        real(8) :: newRightPart(newSize)
        tM%mSize = newSize
        allocate(tM%table(3, newSize), tM%rightPart(newSize))
        tM%table = newTable
        tM%rightPart = newRightPart
        tM%isInit = .true.
    end subroutine tMConstructor
	
    function tMSolution(tM) 
        real(8), allocatable, dimension(:) :: tMSolution 
        type(tridiagonalMatrix) :: tM
    
        integer(4) i, newSize
    
        real(8) :: multy
	
        newSize = tM%mSize
        allocate(tMSolution(newSize))
        if (tM%isInit .eqv. .true.) then
            do i = 1, newSize - 1
                tM%rightPart(i) = tM%rightPart(i) / tM%table(2, i)
                tM%table(3, i) = tM%table(3, i) / tM%table(2, i)
                tM%table(2, i) = 1.0
		
                multy = tM%table(1, i + 1)
                tM%rightPart(i + 1) = tM%rightPart(i + 1) - (multy * tM%rightPart(i))
                tM%table(1, i + 1) = 0.0
                tM%table(2, i + 1) = tM%table(2, i + 1) - (multy * tM%table(3, i))
            end do
					
            tM%rightPart(newSize) = tM%rightPart(newSize) / tM%table(2, newSize)
            tM%table(3, newSize) = tM%table(3, newSize) / tM%table(2, newSize)
            tM%table(2, newSize) = 1.0

            tMSolution(newSize) = tM%rightPart(newSize)
	
            do i = newSize - 1, 1, -1
                tMSolution(i) = tM%rightPart(i) - tM%table(3, i) * tMSolution(i + 1)
            end do
            
            deallocate(tM%table)
            deallocate(tM%rightPart)
            tM%isInit = .false.
        else
            tMSolution = 0.0
            write(*,*) 'tridiagonalMatrix have not initialized'
        endif
    end function tMSolution
	
end module tridiagonal

	