module matrix_
    use var
    use subtools
    implicit none

contains
!+---------------------------------------------------------------------

!+                                  MATRIXP

!+----------------------------------------------------------------------
    subroutine matrixp(ncomp, depth, co, ci1, gso, b11, b12, b21, b22, &
                    g11inv, g12inv, g21inv, g22inv, f1, f2, d1, d2, sigma_dtdx)

        implicit none

        integer, intent(in) :: ncomp
        real,dimension(ncomp),intent(in) ::  depth, co, ci1, gso
        real,dimension(ncomp),intent(out) :: b11, b12, b21, b22, g11inv, g12inv, g21inv, g22inv
        real,dimension(ncomp),intent(out) :: f1, f2, d1, d2, sigma_dtdx

        integer :: i, dm1, dm2, ii
        real, allocatable, dimension(:,:) :: e, f_mat, a, st, g
        real, allocatable, dimension(:) :: d
        real, allocatable, dimension(:) ::  u, c, eps2, eps4
        real :: cour, crmax, crmin, di, dim1, dip1, dKdA, ei, ei1, eia
        real :: dpth, Axs, cA, gdI2dA, Ai, hi
        integer :: iup, idw
        real :: delx, dI1dA, dI2dA, dI1dAprev, Cn1, Cn2

        dm1=2; dm2=2
        allocate(e(dm1,dm2), f_mat(dm1,dm2), a(dm1,dm2), d(dm1), st(dm1,dm2), g(dm1,dm2))
        allocate(u(ncomp), c(ncomp), eps2(ncomp), eps4(ncomp))

        do i=1,ncomp
            u(i)=q(i)/area(i) !u(i)=q(n,i,j)/area(i,j)
            !* c in matrix A, p80,RM1_MESH
            dpth=depth(i)
            Axs=area(i) !area(i,j)
            !* for dI1/dA, dI2/dA, and dSfdA
            if (i==1) then
                iup=i
                idw=i+1
                if (ncomp==1) idw=i
            else
                iup=i-1
                idw=i
             end if
            delx=dx(iup)

            !* dI1/dA for c in matrix A based on p112~3, RM3
            call dI1dAcalc(i, dpth, dI1dA) !call dI1dAcalc(i, j, dpth, dI1dA)

            c(i)=(grav*dI1dA)**0.5
            if (eqn_A==0) then
                ! This is the matrix L (left eigenvector matrix - eq 13)
                e(1,1)=1.0
                if(abs(u(i) - c(i)) < TOLERANCE) then
                    c(i)=c(i)+0.00001
                end if
                e(1,2)=-1.0/(u(i)-c(i))
                e(2,1)=1.0
                e(2,2)=-1.0/(u(i)+c(i))
            else
                !* L changed by Dong Ha
                e(1,1)= -(u(i)+c(i))*(u(i)-c(i))/(2.0*c(i))
                e(1,2)= (u(i)+c(i))/(2.0*c(i))
                e(2,1)= -e(1,1)
                e(2,2)= -(u(i)-c(i))/(2.0*c(i))
            endif

            if (eqn_A==0) then
                ! L^{-1} (inverse of Left eigenvector matrix)
                f_mat(1,1)=-(u(i)-c(i))/(2.0*c(i))
                f_mat(1,2)=(u(i)+c(i))/(2.0*c(i))
                f_mat(2,1)=-((u(i)**2.0)-(c(i)**2.0))/(2.0*c(i))
                f_mat(2,2)=((u(i)**2.0)-(c(i)**2.0))/(2.0*c(i))
            else
                !* L^{-1} changed by Dong Ha
                f_mat(1,1)= 1.0/(u(i)+c(i))
                if(abs(u(i) - c(i)) < TOLERANCE) then
                    f_mat(1,2)= 1.0/(u(i)-c(i)-0.00001)
                else
                    f_mat(1,2)= 1.0/(u(i)-c(i))
                endif
                f_mat(2,1)= 1.0
                f_mat(2,2)= 1.0
            endif

            ! Diagonal wave matrix D (eq 12)
            d(1)=abs(u(i)+c(i))
            d(2)=abs(u(i)-c(i))

            ! Equation 11 (L^{-1} D L)
            a(1,1)=e(1,1)*f_mat(1,1)*d(1)+e(2,1)*f_mat(1,2)*d(2)
            a(1,2)=e(1,2)*f_mat(1,1)*d(1)+e(2,2)*f_mat(1,2)*d(2)
            a(2,1)=e(1,1)*f_mat(2,1)*d(1)+e(2,1)*f_mat(2,2)*d(2)
            a(2,2)=e(1,2)*f_mat(2,1)*d(1)+e(2,2)*f_mat(2,2)*d(2)

            !* ots currently set to zero.
            if(ots == 1) then
                dt(i)=cfl*dx(i-1)/max(d(1),d(2)) !dt(i)=cfl*dx(i-1,j)/max(d(1),d(2))
            else
                dt(i)=dtini
            endif

            !* Courant number for check,p111,RM3
            Cn1=(abs(u(i))+(grav*dpth)**0.5)/(delx/dt(i))
            Cn2=(abs(u(i))-(grav*dpth)**0.5)/(delx/dt(i))

            ! Matrix S (eq 14)
            st(1, 1)=0.0
            st(1, 2)=0.0
            !* gdI2/dA of st(2,1), p114~5,RM3
            call dI2dAcalc(i, iup, idw, delx, dpth, dI2dA) !call dI2dAcalc(i, iup, idw, j, delx, dpth, dI2dA)
            gdI2dA= grav*dI2dA
            !* dK/dA in dSf/dA of S(2,1),p116~8,RM3
            call dKdAcalc_p116(i, dpth, dKdA) !call dKdAcalc_p116(i, j, dpth, dKdA)
            st(2,1)= gdI2dA + gso(i) + f*2.0*grav*area(i)*q(i)*abs(q(i))/(co(i)**3.0)*dKdA
            !st(2,1)=gdI2dA + gso(i) + f*2.0*grav*area(i,j)*q(n,i,j)*abs(q(n,i,j))/(co(i)**3.0)*dKdA
            st(2,2)= -2.0*f*grav*q(i)*area(i)/(co(i)**2.0) !st(2,2)=-2.0*f*grav*q(n,i,j)*area(i,j)/(co(i)**2.0)

            !if(dx(i,j) < TOLERANCE) then
            !    cour=dt(i)
            !else
            !    cour=dt(i)/dx(i,j)
            !    crmax=max(crmax,cour*max(d(1),d(2))) !* not used
            !    crmin=min(crmin,cour*max(d(1),d(2))) !* not used
            !endif
            !* Dong Ha version, chmay
            !* cour in p19,RM1
            if(dx(i) < TOLERANCE) then !if(dx(i,j) < TOLERANCE) then
                write(*,*) "cour in matrixp faces the issue of dx(i,j)<TOLERANCE"
                !pause
            else
                if (i==1) then
                    cour= dt(i)/dx(i) !cour= dt(i)/dx(i,j)
                else
                    cour=0.5*(dt(i-1)/dx(i-1) + dt(i)/dx(i)) !cour=0.5*(dt(i-1)/dx(i-1,j) + dt(i)/dx(i,j))
                endif
            endif
            sigma_dtdx(i)= cour

            ! LHS of eq 7 - this is done after the next check in matrixc
            b11(i)=0.5-phi-theta*cour*a(1,1)-0.5*thes*st(1,1)*dt(i)
            b12(i)=-theta*cour*a(1,2)-0.5*thes*st(1,2)*dt(i)
            b21(i)=-theta*cour*a(2,1)-0.5*thes*st(2,1)*dt(i)
            b22(i)=0.5-phi-theta*cour*a(2,2)-0.5*thes*st(2,2)*dt(i)

            ! This is switched with the dx check about in matrixc...
            !if(i == 1) then
            !    cour=dt(i)
            !else if(dx(i-1,j) < TOLERANCE) then
            !    cour=dt(i)
            !else
            !    cour=dt(i)/dx(i-1,j)
            !endif

            g(1,1)=0.5+phi+theta*cour*a(1,1)-0.5*thes*st(1,1)*dt(i)
            g(1,2)=theta*cour*a(1,2)-0.5*thes*st(1,2)*dt(i)
            g(2,1)=theta*cour*a(2,1)-0.5*thes*st(2,1)*dt(i)
            g(2,2)=0.5+phi+theta*cour*a(2,2)-0.5*thes*st(2,2)*dt(i)

            g11inv(i)= g(2,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
            g12inv(i)=-g(1,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
            g21inv(i)=-g(2,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
            g22inv(i)= g(1,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))

            !* (1,1) component of F matrix, p85,RM2
            f1(i)=q(i) !f1(i)=q(n,i,j)
            !* (2,1) component of F matrix, p85,RM2
            f2(i)=(q(i)**2.0)/area(i)+grav*ci1(i) !f2(i)=(q(n,i,j)**2.0)/area(i,j)+grav*ci1(i)

            eps2(i)=-5 !* just initialize
            if(i >= 2 .and. i < ncomp) then
                ii=i+1
                Ai=area(ii) !Ai=area(ii,j)
                call depthcalc(ii, Ai, hi) !call depthcalc(ii, j, Ai, hi)
                dip1=hi

                Ai=area(i)
                call depthcalc(i, Ai, hi)
                di=2.0*hi

                ii=i-1
                Ai=area(ii)
                call depthcalc(ii, Ai, hi)
                dim1=hi

                if (abs(dip1 + di + dim1) < TOLERANCE) then
                    eps2(i) = 0.0
                else
                    eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1)
                endif
            endif
        end do

        eps2(1)=eps2(2)
        eps2(ncomp)=eps2(ncomp-1)

        do i=2,ncomp-1
            if(ityp(i).ne.1) then  !* ityp is set to 1.
                eps2(i)=eps2(i-1)
                eps2(i+1)=eps2(i+2)
            endif
        end do
        do i=2,ncomp-1
            eps2(i)=max(eps2(i+1),eps2(i))
            eps4(i)=max(0.,alfa4-eps2(i)/(u(i)+c(i)))
        end do
        d1(1)=0.0
        d2(1)=0.0
        d1(ncomp)=0.0
        d2(ncomp)=0.0

        do i=2,ncomp-1
            d(1)=abs(u(i)+c(i))
            d(2)=abs(u(i)-c(i))
            ei=max(d(1),d(2))
            d(1)=abs(u(i+1)+c(i+1))
            d(2)=abs(u(i+1)-c(i+1))
            ei1=max(d(1),d(2))
            eia=(ei+ei1)/2.0
            if(ityp(i) /= 1) then
                d1(i)=0.0
                d2(i)=0.0
            !elseif(i.eq.2.or.i.eq.(ncomp-1)) then
            elseif((i.lt.3).or.(i.gt.(ncomp-2))) then
                d1(i)=eps2(i)*eia*(area(i+1)-area(i)) !=eps2(i)*eia*(area(i+1,j)-area(i,j))
                d2(i)=eps2(i)*eia*(q(i+1)-q(i))
            else
                d1(i)=eps2(i)*eia*(area(i+1)-area(i)) - eps4(i)*(area(i+2)-3.0*area(i+1)+3.0*area(i)-area(i-1))
                d2(i)=eps2(i)*eia*(q(i+1)-q(i)) - eps4(i)*(q(i+2)-3.0*q(i+1)+3.0*q(i)-q(i-1))
            endif
        enddo

        deallocate(e, f_mat, a, d, st, g)
        deallocate(u, c, eps2, eps4)

    end subroutine matrixp
!+---------------------------------------------------------------------

!+                                  MATRIXC

!+----------------------------------------------------------------------
    subroutine matrixc(ncomp, co, ci1, gso, depth, b11, b12, b21, b22, &
                        g11inv, g12inv, g21inv, g22inv, f1, f2, d1, d2)
        implicit none

        integer, intent(in) :: ncomp
        real,dimension(ncomp),intent(in) :: depth, co, ci1, gso
        real,dimension(ncomp),intent(out) :: b11, b12, b21, b22, g11inv, g12inv, g21inv, g22inv
        real,dimension(ncomp),intent(out) :: f1, f2, d1, d2
        integer :: i, ii, dm1, dm2
        real, allocatable, dimension(:,:) :: e, f_mat, a, st, g
        real, allocatable, dimension(:) ::  d, u, c, eps2, eps4
        real :: cour, crmax, crmin, di, dim1, dip1, dKdA, ei, ei1, eia !, ter
        real :: dpth, Axs, cA, gdI2dA, Ai, hi
        integer :: iup, idw
        real :: delx, dI1dA, dI2dA, dI1dAprev, Cn1, Cn2

        dm1=2; dm2=2
        allocate(e(dm1,dm2), f_mat(dm1,dm2), a(dm1,dm2), d(dm1), st(dm1,dm2), g(dm1,dm2))
        allocate(u(ncomp), c(ncomp), eps2(ncomp), eps4(ncomp))

        !do i=1,ncomp
        do i=ncomp,1,-1
            u(i)=qp(i)/areap(i)
            !* c in matrix A, p80,RM1_MESH
            !* depth(i) from sub secpred for a given link value j
            dpth=depth(i)  !y(n,i,j)-z(i,j)
            Axs=areap(i)
            !* for dI1/dA, dI2/dA, and dSfdA
            if (i==ncomp) then
                iup=i
                idw=i
            else
                iup=i
                idw=i+1
            end if
            delx=dx(iup)
            !* dI1/dA for c in matrix A based on p112~3, RM3
            call dI1dAcalc(i, dpth, dI1dA)
            c(i)=(grav*dI1dA)**0.5 !cA

            !* L matrix
            if (eqn_A==0) then
                e(1,1)=1.0
                if(abs(u(i) - c(i)) < TOLERANCE) then
                    c(i)=c(i)+0.00001
                end if
                e(1,2)=-1.0/(u(i)-c(i))
                e(2,1)=1.0
                e(2,2)=-1.0/(u(i)+c(i))
            else
            !* L changed by Dong Ha
                e(1,1)= -(u(i)+c(i))*(u(i)-c(i))/(2.0*c(i))
                e(1,2)= (u(i)+c(i))/(2.0*c(i))
                e(2,1)= -e(1,1)
                e(2,2)= -(u(i)-c(i))/(2.0*c(i))
            endif

            !* L^{-1} matrix
            if (eqn_A==0) then
                f_mat(1,1)=-(u(i)-c(i))/(2.0*c(i))
                f_mat(1,2)=(u(i)+c(i))/(2.0*c(i))
                f_mat(2,1)=-((u(i)**2.0)-(c(i)**2.0))/(2.0*c(i))
                f_mat(2,2)=((u(i)**2.0)-(c(i)**2.0))/(2.0*c(i))
            else
            !* L^{-1} changed by Dong Ha
                f_mat(1,1)= 1.0/(u(i)+c(i))
                if(abs(u(i) - c(i)) < TOLERANCE) then
                    f_mat(1,2)= 1.0/(u(i)-c(i)-0.00001)
                else
                    f_mat(1,2)= 1.0/(u(i)-c(i))
                endif
                f_mat(2,1)= 1.0
                f_mat(2,2)= 1.0
            endif

            !* D matrix
            d(1)=abs(u(i)+c(i))
            d(2)=abs(u(i)-c(i))

            !* A matrix
            a(1,1)=e(1,1)*f_mat(1,1)*d(1)+e(2,1)*f_mat(1,2)*d(2)
            a(1,2)=e(1,2)*f_mat(1,1)*d(1)+e(2,2)*f_mat(1,2)*d(2)
            a(2,1)=e(1,1)*f_mat(2,1)*d(1)+e(2,1)*f_mat(2,2)*d(2)
            a(2,2)=e(1,2)*f_mat(2,1)*d(1)+e(2,2)*f_mat(2,2)*d(2)

            if(ots == 1) then
                dt(i)=cfl*dx(i)/max(d(1),d(2)) !dt(i)=cfl*dx(i,j)/max(d(1),d(2))
            else
                dt(i)=dtini
            endif

            !* Courant number for check,p111,RM3
            Cn1=(abs(u(i))+(grav*dpth)**0.5)/(delx/dt(i))
            Cn2=(abs(u(i))-(grav*dpth)**0.5)/(delx/dt(i))

            !* Matrix S (eq 14)
            st(1, 1)=0.0
            st(1, 2)=0.0
            !* gdI2/dA of st(2,1), p114~5,RM3
            call dI2dAcalc(i, iup, idw, delx, dpth, dI2dA)
            gdI2dA= grav*dI2dA
            !* dK/dA in dSf/dA of S(2,1),p116~8,RM3
            call dKdAcalc_p116(i, dpth, dKdA)
            st(2,1)=gdI2dA + gso(i)+ f*2.0*grav*areap(i)*qp(i)*abs(qp(i))/(co(i)**3.0)*dKdA
            st(2,2)=-2.0*f*grav*qp(i)*areap(i)/(co(i)**2.0)



            !if(i == 1) then
            !    cour=dt(i)
            !else if(dx(i-1,j) < TOLERANCE) then
            !    cour=dt(i)
            !else
            !    cour=dt(i)/dx(i-1,j)
            !endif
            !* Dong Ha version, chmay
            !* cour in p19,RM1
            if(dx(i) < TOLERANCE) then
                write(*,*) "cour in matrixp faces the issue of dx(i,j)<TOLERANCE"
                !pause
            else
                if (i==ncomp) then
                    cour= dt(i)/dx(i)
                else
                    cour=0.5*(dt(i)/dx(i) + dt(i+1)/dx(i+1))
                endif
            endif



            b11(i)=0.5-phi-theta*cour*a(1,1)+0.5*thes*st(1,1)*dt(i)
            b12(i)=-theta*cour*a(1,2)+0.5*thes*st(1,2)*dt(i)
            b21(i)=-theta*cour*a(2,1)+0.5*thes*st(2,1)*dt(i)
            b22(i)=0.5-phi-theta*cour*a(2,2)+0.5*thes*st(2,2)*dt(i)

            !if(dx(i,j) < TOLERANCE) then
            !    cour=dt(i)
            !else
            !    cour=dt(i)/dx(i,j)
            !    crmax=max(crmax,cour*max(d(1),d(2))) !*not used
            !    crmin=min(crmin,cour*max(d(1),d(2))) !*not used
            !endif
            !* cour in p19,RM1
            !if(dx(i,j) < TOLERANCE) then
            !    write(*,*) "cour in matrixp faces the issue of dx(i,j)<TOLERANCE"
            !    pause
            !else
            !    if (i==ncomp) then
            !        cour= dt(i)/dx(i,j)
            !    else
            !        cour=0.5*(dt(i)/dx(i,j) + dt(i+1)/dx(i+1,j))
            !    endif
            !endif

            g(1,1)=0.5+phi+theta*cour*a(1,1)+0.5*thes*st(1,1)*dt(i)
            g(1,2)=theta*cour*a(1,2)+0.5*thes*st(1,2)*dt(i)
            g(2,1)=theta*cour*a(2,1)+0.5*thes*st(2,1)*dt(i)
            g(2,2)=0.5+phi+theta*cour*a(2,2)+0.5*thes*st(2,2)*dt(i)

            g11inv(i)= g(2,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
            g12inv(i)=-g(1,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
            g21inv(i)=-g(2,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
            g22inv(i)= g(1,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))

            f1(i)=qp(i)
            f2(i)=(qp(i)**2.0)/areap(i)+grav*ci1(i)

            eps2(i)=-5 !* just initialize
            if(i >= 2 .and. i < ncomp) then
                ii=i+1
                Ai=areap(ii)
                call depthcalc(ii, Ai, hi)
                dip1=hi

                Ai=areap(i)
                call depthcalc(i, Ai, hi)
                di=2.0*hi

                ii=i-1
                Ai=areap(ii)
                call depthcalc(ii, Ai, hi)
                dim1=hi

                if (abs(dip1 + di + dim1) < TOLERANCE) then
                    eps2(i) = 0.0
                else
                    eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1)
                endif
            endif
        end do

        eps2(1)=eps2(2)
        eps2(ncomp)=eps2(ncomp-1)

        do i=2,ncomp-1
            if(ityp(i).ne.1) then
                eps2(i)=eps2(i-1)
                eps2(i+1)=eps2(i+2)
            endif
        end do
        do i=ncomp-1,1,-1
            eps2(i+1)=max(eps2(i+1),eps2(i))
            eps4(i+1)=max(0.,alfa4-eps2(i+1)/(u(i+1)+c(i+1)))
        end do
        d1(1)=0.0
        d2(1)=0.0
        d1(ncomp)=0.0
        d2(ncomp)=0.0

        do i=2,ncomp-1
            d(1)=abs(u(i)+c(i))
            d(2)=abs(u(i)-c(i))
            ei=max(d(1),d(2))
            d(1)=abs(u(i-1)+c(i-1))
            d(2)=abs(u(i-1)-c(i-1))
            ei1=max(d(1),d(2))
            eia=(ei+ei1)/2.0
            if(ityp(i-1) /= 1) then
                d1(i)=0.0
                d2(i)=0.0
            !elseif(i == 2 .or. i == (ncomp-1)) then
            elseif((i.lt.3).or.(i.gt.(ncomp-2))) then
                d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))
                d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))
            else
                d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))-eps4(i)*(areap(i+1)  &
                            -3.0*areap(i)+3.0*areap(i-1)-areap(i-2))
                d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))-eps4(i)*(qp(i+1)-3.0*qp(i)   &
                            +3.0*qp(i-1)-qp(i-2))
            endif
        end do

        deallocate(e, f_mat, a, d, st, g)
        deallocate(u, c, eps2, eps4)

    end subroutine matrixc
end module matrix_
