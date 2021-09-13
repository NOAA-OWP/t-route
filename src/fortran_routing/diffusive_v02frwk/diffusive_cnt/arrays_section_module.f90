module arrays_section_module

    implicit none
    save

    double precision, allocatable :: elevTable(:),areaTable(:)
    double precision, allocatable :: pereTable(:),rediTable(:)
    double precision, allocatable :: convTable(:),topwTable(:)
    double precision, allocatable :: nwi1Table(:),dPdATable(:)
    double precision, allocatable :: ncompElevTable(:), ncompAreaTable(:), skkkTable(:)
! only in version 2 20191511
    double precision, allocatable :: I2Tablep(:),I2Tablec(:)
    double precision, allocatable :: upstreamI2Tablec(:), downstreamI2Tablep(:)
    double precision, allocatable :: currentSquareDepth(:), currentCubicDepth(:), downstreamSquareDepth(:), upstreamSquareDepth(:)

    integer :: maxTableLength, nel

    character*4 :: file_num
contains

    subroutine setup_arrays_section

        implicit none

        allocate(elevTable(nel))
        allocate(areaTable(nel))
        allocate(pereTable(nel))
        allocate(rediTable(nel))
        allocate(convTable(nel))
        allocate(topwTable(nel))
        allocate(nwi1Table(nel))
        allocate(dPdATable(nel))
        allocate(skkkTable(nel))
		allocate(ncompElevTable(nel))
        allocate(ncompAreaTable(nel))
! only in version 2 20151115
        allocate(I2Tablep(nel))
        allocate(I2Tablec(nel))
        allocate(upstreamI2Tablec(nel))
        allocate(downstreamI2Tablep(nel))
        allocate(currentSquareDepth(nel))
        allocate(currentCubicDepth(nel))
        allocate(downstreamSquareDepth(nel))
        allocate(upstreamSquareDepth(nel))

    end subroutine setup_arrays_section
end module
