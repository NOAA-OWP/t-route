!Written  by Md Nazmul Azim Beg*, Ehab A Meselhe**
!*:     Postdoctoral Fellow, Tulane River and Coastal Center,
!       Department of River-Coastal Science and Engineering,
!       Tulane University
!**:    Professor, Department of River-Coastal Science and Engineering,
!       Tulane University

!Modified by Dong Ha Kim, NOAA's Office of Water Prediction, National Water Center
module attrtable
    use var
    use subtools
    use nrtype

    implicit none

contains
!*--------------------------------------------------
!       intent(out) variables: xsec_tab
!
!*--------------------------------------------------
    subroutine read_xsec
        implicit none

        integer(KIND=i4b)  :: i_area, i_find, i, j,  num, i1, i2 !jj,
        real(KIND=dp) :: el_min, el_max, el_range, el_incr, el_now, x1, y1, x2, y2, x_start, x_end
        real(KIND=dp) ::  cal_area, cal_peri, cal_topW
        real(KIND=dp), dimension(nel) :: el1, a1, peri1, redi1
        real(KIND=dp), dimension(nel) ::  conv1, tpW1
        integer(KIND=i4b) , dimension(nel) :: i_start, i_end
        real(KIND=dp) :: area_tz, hydr_tz, conv_tz, topwd_tz, hbf_tz, depth_tz, mann_tz, uniflw

        num= nxsecpt
       !* compute min. and max. of ycs
        el_min=minval(ycs)
        el_max=maxval(ycs)
        !* nel: # of lines in x-sec attribute table. So, nel-1 numbers of rectangles
        !el_range=(el_max-el_min)*timesDepth !*timesDepth is a multiplier to max. height of given x-section
        el_range=(el_max-el_min)
        el_incr=el_range/real(nel-1.0)
        !open(unit=11,file="./output/readxsection_result.txt",status="unknown")
        j=1
        xsec_attr_seg(1,j) = el_min
        xsec_attr_seg(2,j) = 0.0   !*5 -> 2
        xsec_attr_seg(3,j) = 0.0    !*6 -> 3
        xsec_attr_seg(4,j) = 0.0    !* for uniform flow

        do j=2,nel
            el_now=el_min+real(j-1)*el_incr
            el1(j)=el_now

            if (tzeq_flag==0) then
            !* use original readXsection method developed by Nazmul, Tulane University.
                i_start(1)=-999
                i_end(1)=-999
                i_area=0
                i_find=0
                do i=1,num-1
                    y1=ycs(i)
                    y2=ycs(i+1)
                    if(el_now.le.y1 .and. el_now.gt.y2 .and. i_find.eq.0)then
                        i_find=1
                        i_area=i_area+1
                        i_start(i_area)=i
                    endif

                    if(el_now.gt.y1 .and. el_now.le.y2 .and. i_find.eq.1)then
                        i_find=0
                        i_end(i_area)=i
                    endif
                enddo

                cal_area=0.
                cal_peri=0.
                cal_topW=0.
                !newI1=0.0 !Hu changed
                do i=1,i_area
                    x1=xcs(i_start(i))
                    x2=xcs(i_start(i)+1)
                    y1=ycs(i_start(i))
                    y2=ycs(i_start(i)+1)
                    if(y1.eq.y2)then
                        x_start=x1
                    else
                        x_start=x1+(el_now-y1)/(y2-y1)*(x2-x1)
                    endif

                    x1=xcs(i_end(i))
                    x2=xcs(i_end(i)+1)
                    y1=ycs(i_end(i))
                    y2=ycs(i_end(i)+1)
                    if(y1.eq.y2)then
                      x_end=x1
                    else
                      x_end=x1+(el_now-y1)/(y2-y1)*(x2-x1)
                    endif

                    cal_topW=x_end-x_start+cal_topW

                    i1=i_start(i)
                    i2=i_end(i)
                    x1= xcs(i1+1); y1= ycs(i1+1)
                    x2= xcs(i2); y2= ycs(i2)
                    cal_area = cal_area    &
                             +cal_tri_area(el_now, x_start, x1, y1)    &
                             +cal_multi_area(el_now, i1+1, i2)    &
                             +cal_tri_area(el_now, x_end, x2, y2)

                    x1= xcs(i1+1); y1= ycs(i1+1)
                    x2= xcs(i2); y2= ycs(i2)
                    cal_peri = cal_peri    &
                            +cal_dist(x_start, el_now, x1, y1)    &
                            +cal_perimeter(i1+1, i2)    &
                            +cal_dist(x_end, el_now, x2, y2)

                    x1= xcs(i1+1); y1= ycs(i1+1)
                    if (i1.eq.1) then
                        cal_peri=cal_peri - cal_dist(x_start, el_now, x1, y1)
                    endif

                    x2= xcs(i2); y2= ycs(i2)
                    if (i2.eq.(num-1)) then
                        cal_peri=cal_peri-cal_dist(x_end, el_now, x2, y2)
                    endif
                enddo

                !el1(j)=el_now
                a1(j)=cal_area
                peri1(j)=cal_peri
                redi1(j)=a1(j)/peri1(j)

                !* when overbankflow, manning's N of floodplain is usually different
                !* from that of main channel, conveyance is computed separately.
                !* This computation applies only to trapezoidal main & rectangular floodplains
                hbf_tz=(tw0 - bo0)/(2.0*traps0 )
                depth_tz= el_now - z0
                if (depth_tz>hbf_tz) then
                !* overbank
                    call areacalc(depth_tz, area_tz)
                    call hydRcalc(depth_tz, area_tz, hydr_tz)
                    call Kcalc(depth_tz, area_tz, hydr_tz, conv_tz)
                    conv1(j)= conv_tz
                else
                !* inbank
                    conv1(j)=1./mann*a1(j)*(redi1(j))**(2./3.)
                endif
                tpW1(j)=cal_topW

            elseif (tzeq_flag==1) then
            !* use x-section equations for trapezoidal main and rectangular floodplain
                depth_tz= el_now - z0
                call areacalc(depth_tz, area_tz)
                call hydRcalc(depth_tz, area_tz, hydr_tz)
                call Kcalc(depth_tz, area_tz, hydr_tz, conv_tz)
                hbf_tz= (tw0 - bo0)/(2.0*traps0) !* bankfull depth
                if (depth_tz>hbf_tz) then
                    topwd_tz= twcc0
                else
                    topwd_tz= bo0 + 2.0*traps0*depth_tz
                end if
                conv1(j)= conv_tz
                tpW1(j)=topwd_tz
            endif
            !* uniform flow
            uniflw= conv1(j)*So0**0.5
            !* For this diffusive model only columns of 1, 5, and 6 are used.
            xsec_attr_seg(1,j) = el1(j)
            xsec_attr_seg(2,j) = conv1(j)   !*5 -> 2
            xsec_attr_seg(3,j) = tpW1(j)    !*6 -> 3
            xsec_attr_seg(4,j) = uniflw     !* uniform flow
        enddo
        !close(11)

    end subroutine read_xsec
!*--------------------------------------------
!           compute triangle area
!
!*--------------------------------------------
    real function cal_tri_area(el,x0,x1,y1)
        implicit none
        real(KIND=dp), intent(in) :: el, x0, x1, y1

        cal_tri_area=abs(0.5*(x1-x0)*(el-y1))

    end function cal_tri_area
!*--------------------------------------------
!          compute trapezoidal area
!
!*--------------------------------------------
    real function cal_trap_area(el,x1,y1,x2,y2)
        implicit none
        real(KIND=dp), intent(in) :: el,x1,y1,x2,y2

        cal_trap_area=abs(0.5*(x2-x1)*(el-y1+el-y2))

    end function cal_trap_area
!*--------------------------------------------
!      compute multiple trapezoidal areas
!
!*--------------------------------------------
    !real function cal_multi_area(el,xx,yy,n,i1,i2)
    real function cal_multi_area(el, i1, i2)
        implicit none
        integer(KIND=i4b) , intent(in) :: i1, i2
        real(KIND=dp), intent(in) :: el
        integer(KIND=i4b)  :: i
        real(KIND=dp) :: area, x1, y1, x2, y2

        area=0
        do i=i1,i2-1
            x1=xcs(i)
            y1=ycs(i)
            x2=xcs(i+1)
            y2=ycs(i+1)
            area=area+cal_trap_area(el,x1,y1,x2,y2)
        enddo
        cal_multi_area=area
    end function cal_multi_area
!*--------------------------------------------
!           compute linear distance
!
!*--------------------------------------------
    real function cal_dist(x1,y1,x2,y2)
        implicit none
        real(KIND=dp), intent(in) :: x1,y1,x2,y2

        cal_dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+1.e-32)

    end function cal_dist
!*--------------------------------------------
!           compute perimeter
!
!*--------------------------------------------
!    real function cal_perimeter(xx,yy,n,i1,i2)
    real function cal_perimeter(i1,i2)
        implicit none
        integer(KIND=i4b) , intent(in) :: i1, i2
        integer(KIND=i4b)  :: i
        real(KIND=dp) :: p, x1, y1, x2, y2

        p=0.0
        do i=i1,i2-1
            x1=xcs(i)
            y1=ycs(i)
            x2=xcs(i+1)
            y2=ycs(i+1)
            p=p+cal_dist(x1,y1,x2,y2)
        enddo
        cal_perimeter=p

    end function cal_perimeter
end module attrtable

