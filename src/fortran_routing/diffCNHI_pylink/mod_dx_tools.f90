module mdxtools

    use var
    implicit none

contains
    !*--------------------------------------------------------------------------
    !* Compute qlatj at a given time for the changed length delxj_m
    !* Refer to p.2-5, RM5.
    !*--------------------------------------------------------------------------
    subroutine mod_qlatj(j, ncomp, ncomp_m, delxj_m)
        implicit none
        integer(KIND=i4b),intent(in) :: j, ncomp, ncomp_m
        real(KIND=dp),intent(in) :: delxj_m

        integer(KIND=i4b) :: i, i_m, i_s, i2, x_mflag
        real(KIND=dp) :: x_m, sumx1, sumx2, sumqlatdx, dx_s, sumqlatx_m, sumqlatdx_m

        x_m=0.0
        do i_m=1, ncomp_m-1
            x_m= x_m + dx_ar_m_g(i_m)
            x_mflag=0
            i=1
            sumx1=0.0
            sumx2=0.0
            do while ((x_mflag.eq.0).and.(i.le.ncomp-1))
                sumx1= sumx1 + dx_ar_g(i,j)
                sumx2= sumx1 + dx_ar_g(i+1,j)
                if (x_m.le.dx_ar_g(1,j)) then
                    i_s=1
                    sumx2= dx_ar_g(1,j)
                    x_mflag=1
                elseif ((x_m.gt.sumx1).and.(x_m.le.sumx2)) then
                    i_s= i+1
                    x_mflag=1
                endif
                i=i+1
            enddo

            sumqlatdx=0.0
            do i2= 1, i_s
                sumqlatdx= sumqlatdx + qlatj_g(i2)*dx_ar_g(i2,j)
            enddo
            dx_s= sumx2 - x_m
            !* sum of qlat up to distance x_m from node1(=distance zero)
            sumqlatx_m= sumqlatdx - qlatj_g(i_s)*dx_s
            !* sum of qlat_m up to index i_m - 1
            sumqlatdx_m=0.0
            if (i_m>1) then
                do i2=1, i_m-1
                    sumqlatdx_m = sumqlatdx_m + qlatj_m_g(i2)*delxj_m
                enddo
            endif
            qlatj_m_g(i_m)= (sumqlatx_m - sumqlatdx_m)/delxj_m
        enddo !* do i_m=1, ncomp_m-1

    end subroutine mod_qlatj

    !* ----------------------------------------------------------------------------
    !* compute celty, diffty, q, and qpx at time ts for the changed length delxj_m
    !* p.2-5, RM5
    !* ----------------------------------------------------------------------------
    subroutine mod_cdqqpx(j, ncomp, ncomp_m, delxj_m)
        implicit none
        integer(KIND=i4b),intent(in) :: j, ncomp, ncomp_m
        real(KIND=dp),intent(in) :: delxj_m

        integer(KIND=i4b) :: i, i_m, i_s, ifc
        real(KIND=dp) :: x_m, x_mflag, sumx1, sumx2, dx_s, dx_i_s, delx, u1, u2, u_m

        !* first node
        celty_m_g(1) = celty_g(1,j)
        diffty_m_g(1) = diffty_g(1,j)
        q_m_g(1) = q_g(1,j)
        qpx_m_g(1)= qpx_g(1,j)
        !* last node
        celty_m_g(ncomp_m) = celty_g(ncomp,j)
        diffty_m_g(ncomp_m) = diffty_g(ncomp,j)
        q_m_g(ncomp_m) = q_g(ncomp,j)
        qpx_m_g(ncomp_m)= qpx_g(ncomp,j)
        !* in-between nodes
        x_m=0.0
        do i_m=2, ncomp_m-1
            x_m= x_m + dx_ar_m_g(i_m-1) !delxj_m*real(i_m-1,KIND(delxj_m))
            x_mflag=0
            i=1
            sumx1=0.0
            sumx2=0.0
            do while ((x_mflag.eq.0).and.(i.le.ncomp-1))
                sumx1= sumx1 + dx_ar_g(i,j)
                sumx2= sumx1 + dx_ar_g(i+1,j)
                if (x_m.le.dx_ar_g(1,j)) then
                    i_s=1
                    sumx1=0.0
                    sumx2= dx_ar_g(1,j)
                    x_mflag=1
                elseif ((x_m.gt.sumx1).and.(x_m.le.sumx2)) then
                    i_s= i+1
                    x_mflag=1
                endif
                i=i+1
            enddo
            dx_s= x_m- sumx1
            !dx_i_s= dx_ar_g(i_s,j)
            !delx= dx_i_s - dx_s
            do ifc=1,4
                if (ifc==1) then
                    u1= celty_g(i_s,j)
                    u2= celty_g(i_s+1,j)
                elseif (ifc==2) then
                    u1= diffty_g(i_s,j)
                    u2= diffty_g(i_s+1,j)
                elseif (ifc==3) then
                    u1= q_g(i_s,j)
                    u2= q_g(i_s+1,j)
                elseif (ifc==4) then
                    u1= qpx_g(i_s,j)
                    u2= qpx_g(i_s+1,j)
                endif

                !u_m= (u2-u1)*(delx-dx_i_s)/dx_i_s + u2
                u_m= (u2-u1)*dx_s/dx_ar_g(i_s,j) + u1

                if (ifc==1) then
                    celty_m_g(i_m) = u_m
                elseif (ifc==2) then
                    diffty_m_g(i_m) = u_m
                elseif (ifc==3) then
                    q_m_g(i_m) = u_m
                elseif (ifc==4) then
                    qpx_m_g(i_m) = u_m
                endif
            enddo
        enddo !* do i_m=2, ncomp_m-1

    end subroutine mod_cdqqpx

    !* ----------------------------------------------------------------------------------
    !* Map q_m_g and qpx_m_g at ts+1 back to q_g and qpx_g at ts+1.  q_g at tst+1 becomes
    !* input to computed elv, celty, diffty at ts+1.  qpx_g at ts+1 is also used in the
    !* next run of ef_calc or ef_calc_m. p2-10,RM5
    !* ----------------------------------------------------------------------------
    subroutine mapback_qqpx(j, ncomp, ncomp_m, delxj_m)
        implicit none
        integer(KIND=i4b),intent(in) :: j, ncomp, ncomp_m
        real(KIND=dp),intent(in) :: delxj_m

        integer(KIND=i4b) :: xflag, i, i_m, i_m_s, ifc
        real(KIND=dp) :: x_o, sumx1_m, sumx2_m, dx_s, delx_m, u1, u2, u_o

        q_g(1,j)= q_m_g(1)
        qpx_g(1,j)= qpx_m_g(1)
        q_g(ncomp,j)= q_m_g(ncomp_m)
        qpx_g(ncomp,j)= qpx_m_g(ncomp_m)

        x_o= 0
        do i=2, ncomp-1
            x_o= x_o + dx_ar_g(i-1,j) !* x distance from upstream end node of original segment lengths.
            xflag=0
            i_m=0
            sumx1_m=0.0
            sumx2_m=0.0
            do while ((xflag.eq.0).and.(i_m.lt.ncomp_m))
                sumx1_m= delxj_m*real(i_m,KIND(delxj_m))
                sumx2_m= delxj_m*real(i_m+1,KIND(delxj_m))
                if ((x_o.gt.sumx1_m).and.(x_o.le.sumx2_m)) then
                    i_m_s= i_m+1
                    xflag=1
                endif
                i_m= i_m + 1
            enddo
            dx_s=  x_o - sumx1_m
            !delx_m= delxj_m - dx_s
            do ifc=1,2
                if (ifc==1) then
                    u1= q_m_g(i_m_s)
                    u2= q_m_g(i_m_s+1)
                elseif (ifc==2) then
                    u1= qpx_m_g(i_m_s)
                    u2= qpx_m_g(i_m_s+1)
                endif

                !u_o= (u2-u1)*(delx_m - delxj_m)/delxj_m + u2
                u_o= (u2-u1)*dx_s/delxj_m + u1

                if (ifc==1) then
                    q_g(i,j) = u_o
                elseif (ifc==2) then
                    qpx_g(i,j) = u_o
                endif
            enddo !* ifc=1,2
        enddo !* i=2, ncomp-1

    end subroutine mapback_qqpx

end module mdxtools
