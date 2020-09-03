module predictor
    use var
    use subtools
    use sec_
    use matrix_

    implicit none

contains
!+---------------------------------------------------------

!+                          MAIN

!+---------------------------------------------------------
    subroutine main

        implicit none
        ! Local storage
        integer :: i !n, i, j, k, linknb, nodenb
        real :: areanew, cour, da, dq, dnm0, dnm, dxs
        real :: t, sslp
        real :: Axs, Axs0, Axs1, dsc
        real :: d1junc, hav, h1z, dmy, dpth, areacr
        real :: areasum, hk_ncomp, areak_ncomp, qsum, qk_ncomp
        real :: rhs1, rhs2, c11, c12, c21, c22
        real :: vncomp, FrCC
        real :: aply, bply, cply, idpth_aply, qscl, hmin
        real :: vel, Cn1, Cn2, ynm0, ynm
        integer :: linknb_ds, linknb_us
        real :: qnp1_ds, qnp1_us
        integer :: ncomp
        real, dimension(ncompy) :: depth, ci1, co, ci2, aso, gso
        real,dimension(ncompy) :: b11, b12, b21, b22, g11inv, g12inv, g21inv, g22inv
        real,dimension(ncompy) :: f1, f2, d1, d2, sigma_dtdx

        !dtini=2.0 !;  dxini=20.0;  tfin=300060.0
        !ntim=100 !5001

        !phi=1.0 !* source term treatment (0:explicit, 1:implicit)
        !alfa2=0.5 !* em parameter for artificial diffusion
        !alfa4=0.1 !*maximum value for artificial diffusion

        !theta=1.0;  thetas=1.0;  thesinv=1.0;
        !f=1.0;  skk=33;  yy=2.423;  qq=100.0;  cfl=1.0;  ots=0;
        !ityp = 1


        !dongha_A=1

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

        !
        ! Loop on time
        !
        !do n=1, ntim-1

        ncomp= ncompy
    !+---------------------------------------------------------------------------------+

    !+                              PREDICTOR STEP

    !+---------------------------------------------------------------------------------+
        !do j=1, nlinks
        !ncomp=nx1(j)

        call section(ncomp, depth, ci1, co, ci2, aso, gso)

        !+---------------------------------------------------------------------------------+
        !+                           UPSTREAM BOUNDARY CONDITION

        !+ Predictor step must start first from links where measured or estimated discharge
        !+ data at times n and n+1 is available.  Assuming that discharge values in the upper
        !+ end of a river network are always known, determine whether upstream flow is
        !+ supercritical or subcritical based on directly from y an yc instead of
        !+ steep (yn<yc) or mild channel (yn>yc), p66,RM1
        !+---------------------------------------------------------------------------------+
        !* critical depth or Froud number at the first node of a link where measured discharge
        !* at current time n is available.
        if (usbdflag==1) then
            i=1  ! i_usbd
            sslp= traps(i) !sslp= traps(i,j)
            dsc= q(i)   !dsc= q(n,i,j)
            !* current depth
            dpth= depth(i) !y(i)-z(i) !y(n,i,j)-z(i,j)
            Axs=area(i)  !Axs=area(i,j)  !* result from sub section

            call FrCCcalc(i, dpth, Axs, dsc, FrCC) !call FrCCcalc(i, j, dpth, Axs, dsc, FrCC)

            !dqp(i)= qnp1(i) - q(i)  ! dqp(i,j)= q(n+1,i,j) - q(n,i,j)
            dqp(i)= msrqnp1 -  msrq
            if (FrCC.gt.1.0) then
                !* supercritical flow
                !* Knowing that y(n+1) at this location is unknown, assume the difference of area at n and n+1 is similar
                !* to n-1 and n

!* CHECK HERE if n is indeed > 1.

                !if (n>1) then
                dpth=ynm1 - z(i) !dpth=y(n-1,i,j) - z(i,j)
                !else
                !    dpth=y(n,i,j) - z(i,j)
                !endif
                call areacalc(i, dpth, Axs0)  !call areacalc(i, j, dpth, Axs0)
                dpth=y(i) - z(i) !dpth=y(n,i,j) - z(i,j)
                call areacalc(i, dpth, Axs1) !call areacalc(i, j, dpth, Axs1)
                dap(i)= Axs1-Axs0 !dap(i,j)= Axs1-Axs0
            else
                dap(i)=0.0 !dap(i,j)=0.0
            end if
        endif

        thes=thetas
        call matrixp(ncomp, depth, co, ci1, gso, b11, b12, b21, b22, &
                    g11inv, g12inv, g21inv, g22inv, f1, f2, d1, d2, sigma_dtdx)

!* START: move to Python
        !+-------------------------------------------------------------------------------
        !+                  INTERNAL JUNCTION BOUNDARY CONDITION

        !+ Hand over water from upstream to downstream properly according
        !+ to the nature of link connections, i.e., branching.
        !+ Refer to p.52,RM1_MESH
        !+-------------------------------------------------------------------------------
!        if (ndep(j).gt.1) then  !* ndep(j): the number of links that are immediately upstream of link j.
!            !*total water areas at n+1 at the end nodes of upstream links that join link j
!            areasum=0.0
!            qsum= 0.0
!            do k=1, ndep(j)
!                linknb=uslinks(k,j); nodenb=nx1(linknb)
!                areasum=areasum + area(nodenb,linknb) + dap(nodenb,linknb)
!                qsum= qsum + q(n,nodenb,linknb) + dqp(nodenb,linknb)
!            end do
!
!            dqp(1,j)=0.0;
!            hav=0.0
!            do k=1, ndep(j)
!                linknb=uslinks(k,j); nodenb=nx1(linknb)
                !**dqp(1,j)
                !* discharge at the first node of link j, i.e.,dqp(1,j), is the sum of discharges of at
                !* the last nodes of upstream links that connect to link j at their junction.
!                dqp(1,j)=dqp(1,j)+dqp(nodenb,linknb)
!                !**dap(1,j)
!                !*area at the end nod of link k at time n+1
!                areak_ncomp = area(nodenb,linknb) + dap(nodenb,linknb)  !*may cause different chshp determined in sub section
!                qk_ncomp=  q(n,nodenb,linknb) + dqp(nodenb,linknb)
!                !* depth at the end nod of link k at time n+1
!                call depthcalc(nodenb, linknb, areak_ncomp, dxs)
!                hk_ncomp= dxs
!                !* weighted average based on areas at the end nodes of upstream link ks
!                hav = hav + (areak_ncomp/areasum)*hk_ncomp
!                !* weighted average based on discharge at the end nodes of upstream link ks
!                !hav = hav + (qk_ncomp/qsum)*hk_ncomp
!            end do
!            !* Area estimated for time n+1
!            i=1
!            call areacalc(i, j, hav, Axs)
!            dap(1,j)= Axs - area(1,j)
!        end if
!* END: move to Python

        do i=2,ncomp
            !cour=dt(i)/dx(i-1,j)
            cour= sigma_dtdx(i)
            !* rhs1 with lateral inflow, p86,RM2
            rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1)) + 0.5*dt(i)*(qlat(i)+qlat(i-1))
            !* rhs2, p86,RM2
            rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1)) + dt(i)*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
            c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
            c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
            c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
            dap(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1)-c12*dqp(i-1)
            dqp(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1)-c22*dqp(i-1)
        enddo

        do i=1, ncomp
            areap(i) = area(i) + dap(i)
            qp(i) = q(i) + dqp(i)
        end do

        !+++-----------------------------------------------------------------+
        !+                  Downstream Boundary Condition

        !+ ASSUMPTION: stage data at future time n+1 is available at this node
        !+++-----------------------------------------------------------------+
        if (dsbdflag==1) then
            call dsbc(ncomp, qp, areap, dap, dqp, dac, dqc)
        endif
    end subroutine main
!+---------------------------------------------------------

!+                          DSBC

!+---------------------------------------------------------
    !subroutine dsbc(n,j)
    subroutine dsbc(ncomp, qp, areap, dap, dqp, dac, dqc)
    !subroutine dsbc(mxnd, nlk, ncomp)
    ! Downstream boundary condition. Refer to p.66~67, RM1_MESH

        implicit none

        integer, intent(in) :: ncomp
        !integer, intent(in) :: mxnd, nlk, ncomp
        real, dimension(ncomp), intent(in) :: dap, dqp, qp, areap
        real, dimension(ncomp), intent(out) :: dac, dqc
        real :: dsc, FrCC
        real :: dpcomp, dconj, dtw, Ancomp0, Ancomp1, Axs, dxs


        !* Ancomp0 and Ancomp1, area at time n and n+1 computed from known stage data at the times.
        !dxs= y(ncomp)-z(ncomp) != y(n,ncomp,j)-z(ncomp,j)
        dxs= msrh
        call areacalc(ncomp, dxs, Axs)
        Ancomp0= Axs
        !dxs= ynp1(ncomp)-z(ncomp) != y(n+1,ncomp,j)-z(ncomp,j)
        dxs= msrhnp1
        call areacalc(ncomp, dxs, Axs)
        Ancomp1= Axs

        !* inbank or overbank flow: when overbank flow, sub-/super-critical is determined by Fround number,p93,RM1
        Axs=areap(ncomp)   !* result from sub section
        call depthcalc(ncomp, Axs, dxs)
        dpcomp= dxs
        dsc=qp(ncomp)
        call FrCCcalc(ncomp, dpcomp, Axs, dsc, FrCC)

        !+++----------------------------------------------------------------------
        !+ Determine supercritical or subcritical depending on the size difference
        !+ dconj and yTW (known tailwater stage).
        !+ It is assumed that tailwater stage at the most downstream link is always
        !+ known.
        !+++----------------------------------------------------------------------
        if (FrCC<=1.0) then
        !* subcritical flow
            dac(ncomp)=Ancomp1-Ancomp0
            dqc(ncomp)=dqp(ncomp)
        else
        !* downstream conjugate depth to a computed value of qp
            dsc=qp(ncomp)
            call conjugatedep(ncomp, dpcomp, dsc, dconj)
            !dtw=y(ncomp) - z(ncomp)
            dtw= msrh
            if (dconj<dtw) then
            !*still subcritical
                dac(ncomp)=Ancomp1-Ancomp0
                dqc(ncomp)=dqp(ncomp)
            else
            !*indeed supercritical
                dac(ncomp)=dap(ncomp)
                dqc(ncomp)=dqp(ncomp)
            end if
        end if

    end subroutine dsbc
end module predictor
