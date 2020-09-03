module corrector
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
        integer :: n, i, j, k, linknb, nodenb
        integer :: dmyi
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

        ncomp= ncompy
    !+---------------------------------------------------------------------------------+
    !+                              CORRECTOR STEP
    !+---------------------------------------------------------------------------------+
        call secpred(ncomp, depth, ci1, co, ci2, aso, gso)

        thes=thesinv

        call matrixc(ncomp, co, ci1, gso, depth, b11, b12, b21, b22, &
                        g11inv, g12inv, g21inv, g22inv, f1, f2, d1, d2)

!       !test
!       do i=1,ncomp
!            cop(i)= co(i)
!            ci1p(i)= ci1(i)
!            gsop(i)= gso(i)
!            depthp(i)=depth(i)
!            b11p(i)=b11(i)
!            b12p(i)=b12(i)
!            b21p(i)=b21(i)
!            b22p(i)=b22(i)
!            g11invp(i)= g11inv(i)
!            g12invp(i)= g12inv(i)
!            g21invp(i)= g21inv(i)
!            g22invp(i)= g22inv(i)
!            f1p(i)=f1(i)
!            f2p(i)=f2(i)
!            d1p(i)=d1(i)
!            d2p(i)=d2(i)
!       end do

!* start: move to python
        !+------------------------------------------------------------------------------
        !+          Corrector: INTERNAL JUNCTION BOUNDARY CONDITION

        !+ Handle internal downstream boundary condition for a link j that has a link
        !+ immediately downstream, either in serial or branching channels.
        !+ **Note. dac and dqc at the last node of the most downstream link are computed
        !+ in sub dsbc during predictor step.
        !+ Refer to p.53-1,RM1_MESH
        !+------------------------------------------------------------------------------
!        if (dslink(j).ne.NAnum) then
!            !*when downstream link is available,
!                !* link ID in the downstream of link j
!                linknb=dslink(j)
!                !--------------
!                !*dac(ncomp,j)
!                !--------------
!                areacr=area(1,linknb)+0.5*(dap(1,linknb)+dac(1,linknb))
!                !* depth at the last node of link j for time n+1 estimated from the first node of the next downstream link
!                i=1
!                call depthcalc(i, linknb, areacr, dxs)
!                h1z= dxs
!                call areacalc(ncomp, j, h1z, Axs)
!                dac(ncomp,j)= 2.0*(Axs - area(ncomp,j)) - dap(ncomp,j)
!                !---------------
!                !* dqc(ncomp,j)
!                !---------------
!                !dqc(ncomp,j)=dqc(1,linknb)*q(n,ncomp,j)/q(n,1,linknb)
!                !dqc(ncomp,j)=dqc(1,linknb)*qp(ncomp,j)/qp(1,linknb)
!                !* p.120,RM3
!                qsum= 0.0
!                linknb_ds= linknb
!                do k=1, ndep(linknb_ds)
!                    !* uslinks(k,j): k_th link ID that is immediately upstream of link j
!                    linknb_us=uslinks(k,linknb_ds); nodenb=nx1(linknb_us)
!                    qsum= qsum + qp(nodenb,linknb_us)
!                end do
!                qnp1_ds= q(n,1,linknb_ds) +0.5*(dqp(1,linknb_ds)+dqc(1,linknb_ds))
!                !* est. q(n+1, ncomp, link j_i), p120_RM
!                qnp1_us= qnp1_ds*qp(ncomp,j)/qsum
!                dqc(ncomp,j)= 2.0*(qnp1_us - q(n,ncomp,j)) - dqp(ncomp,j)
!
!        end if
!* end: move to python

        do i=ncomp-1,1,-1
            cour=dt(i)/dx(i)
            !* rhs1 with lateral inflow
            rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+ 0.5*dt(i)*(qlat(i+1)+qlat(i))
            rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dt(i)*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i+1)+g12inv(i)*b21(i+1)
            c12=g11inv(i)*b12(i+1)+g12inv(i)*b22(i+1)
            c21=g21inv(i)*b11(i+1)+g22inv(i)*b21(i+1)
            c22=g21inv(i)*b12(i+1)+g22inv(i)*b22(i+1)
            dac(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dac(i+1)-c12*dqc(i+1)
            dqc(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dac(i+1)-c22*dqc(i+1)
        end do

        ! Final update q an area (thus y) at time n+1
        do i=1,ncomp
            !* update y(n+1,i,j)
            da=(dap(i)+dac(i))/2.0
            areanew=da+area(i)
            if(areanew <= 0.0) areanew=0.0001
            !areafnal(n+1,i,j)=areanew
            call depthcalc(i, areanew, dpth)
            ynp1(i)= dpth + z(i) !y(n+1,i,j)= dpth + z(i,j)
            !* update q(n+1,i,j)
            dq=(dqp(i)+dqc(i))/2.0
            qnp1(i)=q(i)+dq  !q(n+1,i,j)=q(n,i,j)+dq
        end do

        !do i=1, ncomp
            !* Courant number for check,p111,RM3
        !    vel= q(n+1,i,j)/areafnal(n+1,i,j)
        !    dpth= y(n+1,i,j)-z(i,j)
        !    Cn1=(abs(vel)+(grav*dpth)**0.5)/(dx(i,j)/dt(i))
        !    Cn2=(abs(vel)-(grav*dpth)**0.5)/(dx(i,j)/dt(i))
        !enddo
    end subroutine main

end module corrector
