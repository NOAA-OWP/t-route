module ini
    use var
    use subtools
    implicit none
    integer :: ncomp

contains

    subroutine iniydsbd
    !* This subroutine computes the values of y at time=1 (that is, simulation beginning time) for
    !* all the segments of a downstream terminal reach, p.F<->P_init.cond.-1~2
        implicit none
        integer :: i
        real ::dpth, Axs, vncomp
        ncomp= ncompy
        !* lateral inflow given
        !qlat(n,i,j)= 0.5*real(i) !* move to python
        !* headbasin segment discharge given
        !do n=1,ntim        !* move to python
        !j=1; q(n,1,j)=10.0 !* move to python
        !j=2; q(n,1,j)=12.0 !* move to python
        !enddo              !* move to python
        !+++---------------------------------------------------------------------------------+
        !+                      ESTIMATE of q at time=1
        !+++---------------------------------------------------------------------------------+
        !+ 1) Headbasin link: initialize q along all the nodes with the measured data at n=1
        !+ until the link meets a junction.
        !n=1 !* move to python
        !j=1 !* move to python
        !ncomp=nx1(j)    !* move to python
        !do i=1, ncomp    !* move to python
        !    q(n,i,j) = q(n,1,j) + qlat(n,i,j) !* move to python
        !end do    !* move to python
        !j=2    !* move to python
        !ncomp=nx1(j)    !* move to python
        !do i=1, ncomp    !* move to python
        !    q(n,i,j) = q(n,1,j) + qlat(n,i,j)    !* move to python
        !end do    !* move to python
        !+ 2) At a junction: Treat q for link j that receive water from more than one immediate
        !+ upstream link.
        !do j=1, nlinks    !* move to python
        !    ncomp=nx1(j)    !* move to python
        !    if (ndep(j)>1) then    !* move to python
        !        n=1    !* move to python
        !        q(n,1,j)=0    !* move to python
        !        do k=1, ndep(j)    !* move to python
        !            linknb=uslinks(k,j) !*link number of k_th link in the immediate upstream of link j    !* move to python
        !            q(n,1,j)=q(n,1,j) + q(n, nx1(linknb), linknb)    !* move to python
        !        end do    !* move to python
        !        do i=1, ncomp    !* move to python
        !            q(n,i,j) = q(n,1,j) + qlat(n,i,j)    !* move to python
        !        end do    !* move to python
        !    end if    !* move to python
        !end do    !* move to python

        !*test2. Set normal depth at ds boundary with constant q and qlat in time.
        !dsc=0.0
        !do j=1, nlinks
        !    ncomp=nx1(j)
        !    do i=1,ncomp
        !        dsc= dsc + q(n,i,j)
        !    end do
        !enddo
        !ynm0=4.4
        !j=3 !5
        !ncomp=nx1(j)
        !call nmaldepcalc(ncomp, j, ynm0, dsc, ynm)
        !* given stage data at the downstream boundary
        !do n=1, ntim    !* move to python
        !    ncomp=nx1(j)    !* move to python
        !    y(n,ncomp,j)=ynm + z(ncomp,j)    !* move to python
        !end do    !* move to python
        !+++----------------------------------------------------------------------------------------+
        !+                      ESTIMATE OF y AT TIME=1, p68-2~3,RM1    !+
        !+
        !+ Up to this point, all nodes of links that are connected at a common junction already take
        !+ initial estimates of q.  Now, initial values of y is estimated with the assumption that
        !+ velocity is consistent throughout all the nodes of a link and its value is computed by
        !+ Q=AV where A is first computed using known stage value at the last node of the most
        !+ downstream link.
        !+++----------------------------------------------------------------------------------------+
        !n=1        !* move to python
        !j=3    !* move to python
        !* First, using known stage at the last node of a downstream terminal reach, determine initial velocity
        !* with the known stage data and then compute stage of the rest nodes of the reach.
        !do j=1, nlinks
        !   if (instrdflag(j,2)==1) then
        !* when water depth data is available at the lower end of a link
            !* velocity at node ncomp, p68-2,RM1
        !ncomp=nx1(j) !* move to python
        dpth=y(ncomp) - z(ncomp) !dpth=y(n,ncomp,j) - z(ncomp,j)
        call areacalc(ncomp, dpth, Axs) !call areacalc(ncomp, j, dpth, Axs)
        vncomp=q(ncomp)/Axs
        !* estimate y at the rest nodes of the link using the estimated velocity vncomp
        do i=1, ncomp-1
            Axs=q(i)/vncomp !Axs=q(n,i,j)/vncomp
            call depthcalc(i, Axs, dpth) !call depthcalc(i, j, Axs, dpth)
            y(i)=dpth+z(i)  !y(n,i,j)=dpth+z(i,j)
        end do
        !    end if
        !end do
    end subroutine iniydsbd
    !#-----------------------------------------------------

    !#                  NEXT SUBROUTINE

    !#-----------------------------------------------------
    subroutine iniy(depth_ds)
    !* This subroutine computes the values of y at time=1 (that is, simulation beginning time) for
    !* all the segments of reaches before the downstream terminal reach, p.F<->P_init.cond.-1~2
        implicit none
        real, intent(in) :: depth_ds
        integer :: i
        real :: Axs, dpth, vncomp
        ncomp= ncompy

        !* assume that depths at a junction are all the same for any segments directly connected to it.
        y(ncomp)= depth_ds + z(ncomp)
        call areacalc(ncomp, depth_ds, Axs) !call areacalc(ncomp, j, dpth, Axs)
        vncomp=q(ncomp)/Axs ! vncomp=q(n,ncomp,j)/Axs
        do i=1, ncomp-1
            Axs=q(i)/vncomp   !Axs=q(n,i,j)/vncomp
            call depthcalc(i, Axs, dpth) !call depthcalc(i, j, Axs, dpth)
            y(i)=dpth +  z(i) !y(n,i,j)=dpth +  z(i,j)
        end do
    end subroutine iniy

end module ini
