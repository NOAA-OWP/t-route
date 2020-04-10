subroutine matrixc(n,j)

    use constants
    use arrays
    use var
    use subtools

    implicit none

    integer, intent(in) :: n, j
    integer :: i, ii, dm1, dm2
    real, allocatable, dimension(:,:) :: e, f_mat, a, st, g
    real, allocatable, dimension(:) :: d
    real :: cour, crmax, crmin, di, dim1, dip1, dkda, ei, ei1, eia !, ter
    real :: dpth, Axs, cA, gdI2dA, Ai, hi

    dm1=2; dm2=2
    allocate(e(dm1,dm2), f_mat(dm1,dm2), a(dm1,dm2), d(dm1), st(dm1,dm2), g(dm1,dm2))

    do i=1,ncomp
        u(i)=qp(i,j)/areap(i,j)
        !* c in matrix A, p80,RM1_MESH
        !* depth(i) from sub secpred for a given link value j
        dpth=depth(i)  !y(n,i,j)-z(i,j)

        Axs=areap(i,j)

        call c_mtrxAcalc(i, j, dpth, Axs, cA)
        c(i)=cA

        e(1,1)=1.0
        if(abs(u(i) - c(i)) < TOLERANCE) then
            c(i)=c(i)+0.00001
        end if
        e(1,2)=-1.0/(u(i)-c(i))
        e(2,1)=1.0
        e(2,2)=-1.0/(u(i)+c(i))

        f_mat(1,1)=-(u(i)-c(i))/(2.0*c(i))
        f_mat(1,2)=(u(i)+c(i))/(2.0*c(i))
        f_mat(2,1)=-((u(i)**2.0)-(c(i)**2.0))/(2.0*c(i))
        f_mat(2,2)=((u(i)**2.0)-(c(i)**2.0))/(2.0*c(i))

        d(1)=abs(u(i)+c(i))
        d(2)=abs(u(i)-c(i))

        a(1,1)=e(1,1)*f_mat(1,1)*d(1)+e(2,1)*f_mat(1,2)*d(2)
        a(1,2)=e(1,2)*f_mat(1,1)*d(1)+e(2,2)*f_mat(1,2)*d(2)
        a(2,1)=e(1,1)*f_mat(2,1)*d(1)+e(2,1)*f_mat(2,2)*d(2)
        a(2,2)=e(1,2)*f_mat(2,1)*d(1)+e(2,2)*f_mat(2,2)*d(2)

        if(ots == 1) then
           dt(i)=cfl*dx(i,j)/max(d(1),d(2))
        else
           dt(i)=dtini
        endif

        !* Matrix S (eq 14)
        st(1, 1)=0.0
        st(1, 2)=0.0
        !* gdI1/dA of st(2,1)
        call gdI2dAcalc(i, j, dpth, Axs, gdI2dA)
        !* dK/dA of st(2,1)
        call dKdAcalc(i, j, dpth, Axs, dkda)
        st(2, 1)=gdI2dA + gso(i)+ f*2.0*grav*areap(i,j)*qp(i,j)*abs(qp(i,j))/(co(i)**3.0)*dkda
        st(2,2)=-2.0*f*grav*qp(i,j)*areap(i,j)/(co(i)**2.0)

        if(i == 1) then
            cour=dt(i)
        else if(dx(i-1,j) < TOLERANCE) then
            cour=dt(i)
        else
            cour=dt(i)/dx(i-1,j)
        endif
        b11(i)=0.5-phi-theta*cour*a(1,1)+0.5*thes*st(1,1)*dt(i)
        b12(i)=-theta*cour*a(1,2)+0.5*thes*st(1,2)*dt(i)
        b21(i)=-theta*cour*a(2,1)+0.5*thes*st(2,1)*dt(i)
        b22(i)=0.5-phi-theta*cour*a(2,2)+0.5*thes*st(2,2)*dt(i)

        if(dx(i,j) < TOLERANCE) then
            cour=dt(i)
        else
            cour=dt(i)/dx(i,j)
            crmax=max(crmax,cour*max(d(1),d(2)))
            crmin=min(crmin,cour*max(d(1),d(2)))
        endif
        g(1,1)=0.5+phi+theta*cour*a(1,1)+0.5*thes*st(1,1)*dt(i)
        g(1,2)=theta*cour*a(1,2)+0.5*thes*st(1,2)*dt(i)
        g(2,1)=theta*cour*a(2,1)+0.5*thes*st(2,1)*dt(i)
        g(2,2)=0.5+phi+theta*cour*a(2,2)+0.5*thes*st(2,2)*dt(i)

        g11inv(i)= g(2,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g12inv(i)=-g(1,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g21inv(i)=-g(2,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g22inv(i)= g(1,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))

        f1(i)=qp(i,j)
        f2(i)=(qp(i,j)**2.0)/areap(i,j)+grav*ci1(i)

        if(i >= 2 .and. i < ncomp) then
            ii=i+1
            Ai=areap(ii,j)
            call depthcalc(ii, j, Ai, hi)
            dip1=hi

            Ai=areap(i,j)
            call depthcalc(i, j, Ai, hi)
            di=2.0*hi

            ii=i-1
            Ai=areap(ii,j)
            call depthcalc(ii, j, Ai, hi)
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
        elseif(i == 2 .or. i == (ncomp-1)) then
            d1(i)=eps2(i)*eia*(areap(i,j)-areap(i-1,j))
            d2(i)=eps2(i)*eia*(qp(i,j)-qp(i-1,j))
        else
          d1(i)=eps2(i)*eia*(areap(i,j)-areap(i-1,j))-eps4(i)*(areap(i+1,j)  &
                            -3.0*areap(i,j)+3.0*areap(i-1,j)-areap(i-2,j))
          d2(i)=eps2(i)*eia*(qp(i,j)-qp(i-1,j))-eps4(i)*(qp(i+1,j)-3.0*qp(i,j)   &
                            +3.0*qp(i-1,j)-qp(i-2,j))
        endif
    end do

    deallocate(e, f_mat, a, d, st, g)

end subroutine matrixc
