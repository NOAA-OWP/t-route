subroutine section(n,j)

    use constants
    use arrays
    use var
    use subtools

    implicit none

    integer, intent(in) :: n, j
    ! Locals
    integer :: i
    integer :: iup, idw, iadj
    real :: beds, fs, hy, Axs, hydR, I1
    real :: dpth, ybf
    real :: dbdxval
    real :: naval, dbdxm, dbdxf, I2im1, I2i, I2
    real :: cnvey,  sfi, asoi, asoim1

    do i=1,ncomp
        depth(i)=y(n,i,j)-z(i,j)
        dpth=depth(i)
        call areacalc(i, j, dpth, Axs)
        area(i,j)=Axs
        !* I1=yc*A
        call I1calc(i, j, dpth, Axs, I1)
        ci1(i)=I1
        !* R=A/P
        call hydRcalc(i, j, dpth, Axs, hydR)
        hy=hydR
        !* K=(1/n)*A*R**(2/3)
        call Kcalc(i, j, dpth, Axs, hydR, cnvey)
        co(i)=cnvey
    end do

    naval=-1000.0
    do i=2,ncomp
        if(ityp(i-1) == 1) then
            iup=i-1; idw=i; iadj=i-1
            !* I2
            dbdxval=naval; dbdxm= dbdxval; dbdxf= dbdxval
            ybf= (Tw(i,j) - bo(i,j))/(2.0*traps(i,j))
            dpth= depth(i)
            if (dpth<=ybf) then
            !* inbank flow in regular or trapezoidal channel, p81, p29-C, p86, RM1
            !* db/bx for main channel at current node, p81-2
                call dbdxcalc(i, iup, idw, j, dpth, dbdxval)
                I2i=0.5*dbdxval*dpth**2.0
            else
            !* overbank flow in compound channel db/dx at current node, P81-, RM1
                call dbdxCCcalc(iup, idw, iadj, j, dbdxm, dbdxf)
                call I2CCcalc(i, j, dbdxm, dbdxf, dpth, I2)
                I2i=I2
            endif
            !* ci2(i)=0.5*(I2(i)+I2(I-1))
            if (i==2) I2im1=I2i
            ci2(i)=0.5*(I2i+I2im1)
            I2im1=I2i
            !* So
            beds=(z(i-1,j)-z(i,j))/dx(i-1,j)
            !* Sf
            fs=f*0.5*q(n,i-1,j)*abs(q(n,i-1,j))/(co(i-1)**2.0)+f*0.5*q(n,i,j)*abs(q(n,i,j))/(co(i)**2.0)
            !* A(So - Sf), p29-c,P86,RM1
            sfi=fs  !*take average friction slope
            asoi=area(i,j)*(beds-sfi)
            if (i==2) asoim1=asoi
            aso(i)=0.5*(asoi+asoim1)
            asoim1=asoi
            !* g(So - Sf)
            gso(i)=grav*(beds-fs)
            !* db/dx of either inbank or overbank flow
            dbdx(i)=dbdxval
            dbdxCCm(i)=dbdxm
            dbdxCCf(i)=dbdxf
        end if
    end do
    !* assume value of db/dx at i=1 zero
    dbdx(1)= 0.0
    dbdxCCm(1)= 0.0
    dbdxCCf(1)= 0.0
    !** added to fill in values at i=1
    gso(1)=gso(2)

end subroutine section
