
subroutine secpred(j)

    use constants
    use arrays
    use var
    use subtools

    implicit none

    integer,intent(in) :: j
    ! Locals
    integer :: i, iup, idw,iadj
    real :: beds, fs, hy
    real :: dpth, Axs, I1, hydR, ybf
    real :: dbdxm, dbdxf, dbdxval
    real :: cnvey
    real :: I2, I2i, I2ip1, sfi, asoi, asoip1, naval

    do i=1,ncomp
        Axs=areap(i,j)
        call depthcalc(i, j, Axs, dpth)
        depth(i)=dpth
        !* I1
        call I1calc(i, j, dpth, Axs, I1)
        ci1(i)=I1
        !* R=A/P
        call hydRcalc(i, j, dpth, Axs, hydR)
        hy=hydR
        !* K
        call Kcalc(i, j, dpth, Axs, hydR, cnvey)
        co(i)=cnvey
    end do

    naval=-1000.0

    do i=ncomp-1, 1, -1
        if(ityp(i) == 1) then
            iup=i; idw=i+1; iadj=i+1
            !* I2, p81, RM1
            dbdxval=naval; dbdxm= dbdxval; dbdxf= dbdxval
            ybf= (Tw(i,j) - bo(i,j))/(2.0*traps(i,j))
            dpth= depth(i)
            if (dpth<=ybf) then
                !* regular or trapezoidal channel, p81, p29-C, p86, RM1
                !* db/bx for main channel at current node, p81-2
                call dbdxcalc(i, iup, idw, j, dpth, dbdxval)
                I2i=0.5*dbdxval*dpth**2.0
            else
                !* compound channel db/dx, P81-, RM1
                call dbdxCCcalc(iup, idw, iadj, j, dbdxm, dbdxf)
                call I2CCcalc(i, j, dbdxm, dbdxf, dpth, I2)
                I2i=I2
            endif
            if (i==ncomp-1) I2ip1=I2i
            ci2(i)=0.5*(I2i+I2ip1)
            I2ip1=I2i
            !* So
            beds=(z(i,j)-z(i+1,j))/dx(i,j)
            !* Sf
            fs=f*0.5*qp(i,j)*abs(qp(i,j))/(co(i)**2)+f*0.5*qp(i+1,j)*abs(qp(i+1,j))/(co(i+1)**2)
            !* A(So - Sf), p35-c,p90,RM1
            sfi=fs  !*take average friction slope
            asoi=areap(i,j)*(beds-sfi)
            if (i==ncomp-1) asoip1=asoi
            aso(i)=0.5*(asoi+asoip1)
            asoip1=asoi
            !* g(So - Sf)
            gso(i)=grav*(beds-fs)
            !* db/dx
            dbdx(i)=dbdxval
            dbdxCCm(i)=dbdxm
            dbdxCCf(i)=dbdxf
        endif
    end do
    !* assume value of db/dx at i=1 zero
    dbdx(ncomp)= 0.0
    dbdxCCm(ncomp)= 0.0
    dbdxCCf(ncomp)= 0.0
    !** added to fill in values at i=ncomp
    gso(ncomp)=gso(ncomp-1)

end subroutine secpred
