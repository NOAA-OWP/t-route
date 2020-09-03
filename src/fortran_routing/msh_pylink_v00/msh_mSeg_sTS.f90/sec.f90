module sec_
    use var
    use subtools
    implicit none

contains
!+---------------------------------------------------------

!+                          SECTION

!+---------------------------------------------------------
    subroutine section(ncomp, depth, ci1, co, ci2, aso, gso)

    use var
    use subtools

    implicit none

    integer, intent(in) :: ncomp
    real, dimension(ncomp), intent(out) :: depth, ci1, co, ci2, aso, gso
    integer :: i
    integer :: iup, idw
    real :: dpth, hy, Axs, hydR, I1
    real :: I2im1, I2i, I2
    real :: cnvey,  sfi, soi, sfim1, soim1, asoi, asoim1
    real :: delx,  wtus, wtds, hoI2


    do i=1,ncomp
        dpth=y(i)-z(i) !dpth=y(n,i,j)-z(i,j)
        depth(i)=dpth
        call areacalc(i, dpth, Axs) !call areacalc(i, j, dpth, Axs)
        area(i)=Axs !area(i,j)=Axs
        !* I1, p112~113,RM3
        call I1calc(i, dpth, Axs, I1) !call I1calc(i, j, dpth, Axs, I1)
        ci1(i)=I1
        !* R=A/P
        call hydRcalc(i, dpth, Axs, hydR)  !call hydRcalc(i, j, dpth, Axs, hydR)
        hy=hydR
        !* K=(1/n)*A*R**(2/3)
        call Kcalc(i, dpth, Axs, hydR, cnvey) !call Kcalc(i, j, dpth, Axs, hydR, cnvey)
        co(i)=cnvey
    end do

    !* compute I2 at i=1
    i=1
    iup=1
    idw=2
    delx= dx(iup) !delx= dx(iup,j)
    dpth=depth(i)
    call I2calc(i, iup, idw, delx, dpth, I2im1) !call I2calc(i, iup, idw, j, delx, dpth, I2im1)

    do i=2,ncomp
        !if(ityp(i-1) == 1) then
            iup=i-1; idw=i
            !------------------------------------------------------------------------------------------
            !* ci2(i)=0.5*[I2(i) + I2(i-1)], p29-c,RM2; only used in RHS2, p86,RM2; I2 calc, p.114-115, RM3
            !------------------------------------------------------------------------------------------
            delx= dx(iup) !delx= dx(iup,j)
            dpth=depth(i)
            call I2calc(i, iup, idw, delx, dpth, I2i) !call I2calc(i, iup, idw, j, delx, dpth, I2i)
!            wtus= area(iup)/(area(iup)+area(idw))
!            wtds= area(idw)/(area(iup)+area(idw))
!            hoI2= wtus*depth(iup) + wtds*depth(idw)
!            call I2calc(i, iup, idw, delx, hoI2, dpth, I2i)
            ci2(i)= 0.5*(I2i + I2im1)
            I2im1= I2i
            !---------------------------------------------------------------------------------------------------
            !* aso(i)= 0.5*{A(i)*[So(i)-Sf(i)] + A(i-1)*[So(i-1)-Sf(i-1)]}, p29-c,RM2; only used in RHS2,p86,RM2
            !---------------------------------------------------------------------------------------------------
            !* 1) channel bed slope
            soim1= chbtslp(i-1) !soim1= chbtslp(i-1,j)
            soi= chbtslp(i) !soi= chbtslp(i,j)
            !* 2) Sf= Q**2/K**2
            sfim1= f*q(i-1)*abs(q(i-1))/(co(i-1)**2.0)  !sfim1= f*q(n,i-1,j)*abs(q(n,i-1,j))/(co(i-1)**2.0)
            sfi= f*q(i)*abs(q(i))/(co(i)**2.0)  !sfi= f*q(n,i,j)*abs(q(n,i,j))/(co(i)**2.0)
            aso(i) = 0.5*(area(i)*(soi-sfi) + area(i-1)*(soim1-sfim1))  !aso(i) = 0.5*(area(i,j)*(soi-sfi) + area(i-1,j)*(soim1-sfim1))
            !--------------------------------
            !* g[So(i) - Sf(i)],p29-c,RM2
            !--------------------------------
            gso(i)= grav*(soi-sfi)
            if (i==2) then
                gso(1)= grav*(soim1-sfim1)
            endif
        !end if
    end do
end subroutine section
!+---------------------------------------------------------

!+                          SECPRED

!+---------------------------------------------------------
    subroutine secpred(ncomp, depth, ci1, co, ci2, aso, gso)

        implicit none

        integer,intent(in) :: ncomp
        !real, dimension(ncomp), intent(in) :: areap, qp
        real, dimension(ncomp), intent(out) :: depth, ci1, co, ci2, aso, gso
        integer :: i, iup, idw
        real :: beds, fs, hy
        real :: dpth, Axs, I1, hydR, ybf
        real :: cnvey, delx
        real :: I2, I2i, I2ip1, asoi, asoip1, naval
        real :: sfi, soi, sfip1, soip1

        do i=1,ncomp
            Axs=areap(i)
            call depthcalc(i, Axs, dpth)
            depth(i)=dpth
            !* I1, P112~113,RM3
            call I1calc(i, dpth, Axs, I1)
            ci1(i)=I1
            !* R=A/P
            call hydRcalc(i, dpth, Axs, hydR)
            hy=hydR
            !* K
            call Kcalc(i, dpth, Axs, hydR, cnvey)
            co(i)=cnvey
        end do

        !* compute I2 at i=ncomp
        i=ncomp
        iup=ncomp-1
        idw=ncomp
        delx= dx(iup)
        dpth=depth(i)
        call I2calc(i, iup, idw, delx, dpth, I2ip1)
        I2ip1=0.0
        do i=ncomp-1, 1, -1
            !if(ityp(i) == 1) then
            iup=i; idw=i+1
            !---------------------------------------------------------------------------------------------
            !* ci2(i)=0.5*[I2(i+1) + I2(i)], p35-c,RM2; only used in RHS2, p90,RM2; I2 calc, p.114~5, RM3
            !---------------------------------------------------------------------------------------------
            delx= dx(iup)
            dpth= depth(i)
            call I2calc(i, iup, idw, delx, dpth, I2i)
            ci2(i)= 0.5*(I2ip1 + I2i)
            I2ip1= I2i
            !---------------------------------------------------------------------------------------------------
            !* aso(i)= 0.5*{A(i+1)*[So(i+1)-Sf(i+1)] + A(i)*[So(i)-Sf(i)]}, p35-c,RM2; only used in RHS2,p90,RM2
            !---------------------------------------------------------------------------------------------------
            !* 1) channel bed slope
            soip1= chbtslp(i+1)
            soi= chbtslp(i)
            !* 2) Sf= Q**2/K**2
            sfip1= f*qp(i+1)*abs(qp(i+1))/(co(i+1)**2.0)
            sfi= f*qp(i)*abs(qp(i))/(co(i)**2.0)
            aso(i) = 0.5*(areap(i+1)*(soip1-sfip1) + areap(i)*(soi-sfi))
            !--------------------------------
            !* g[So(i) - Sf(i)],p35-c,RM2
            !--------------------------------
            gso(i)= grav*(soi-sfi)
            if (i.eq.(ncomp-1)) then
                gso(ncomp)= grav*(soip1-sfip1)
            endif
            !endif
        end do
    end subroutine secpred

end module sec_
