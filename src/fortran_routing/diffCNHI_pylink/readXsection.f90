module attrtable
    !use constants_module
    !use arrays_module
    !use arrays_section_module
    !use xsec_attribute_module
    use var
    use subtools
    use nrtype

    implicit none

!    save
!    real, dimension(:), allocatable :: xcs, ycs -> var_module

contains
!*--------------------------------------------------
!       intent(out) variables: xsec_tab
!
!*--------------------------------------------------
    !subroutine readXsection(k,rmanning,timesDepth):
    !* k: segment ID, rmanning: Manning's N, timesDepth: multiplier to the highest possible water depth
    subroutine read_xsec
        implicit none

        integer(KIND=i4b)  :: i_area, i_find, i, j,  num, i1, i2, il !jj,
        real(KIND=dp) :: el_min, el_max, el_range, el_incr, el_now, x1, y1, x2, y2, x_start, x_end
        real(KIND=dp) ::  cal_area, cal_peri, cal_topW
        real(KIND=dp), dimension(nel_g) :: el1, a1, peri1, redi1
        real(KIND=dp), dimension(nel_g) ::  conv1, tpW1
        integer(KIND=i4b) , dimension(nel_g) :: i_start, i_end
        real(KIND=dp) :: area_tz, hydr_tz, conv_tz, topwd_tz, hbf_tz, depth_tz, mann_tz, uniflw


        num= nxsecpt_g ! num=i

       !* compute min. and max. of ycs
        !el_min=99999.
        !el_max=-99999.
        !do i=2,num-1
        !    if(ycs(i).lt.el_min)el_min=ycs(i)
        !    if(ycs(i).gt.el_max)el_max=ycs(i)
        !enddo
        el_min=minval(ycs_g)
        el_max=maxval(ycs_g)
        !* nel_g: # of lines in x-sec attribute table. So, nel_g-1 numbers of rectangles
        !el_range=(el_max-el_min)*timesDepth !*timesDepth is a multiplier to max. height of given x-section
        el_range=(el_max-el_min)
        el_incr=el_range/real(nel_g-1.0)

!        xcs(1)=xcs(2)
!        ycs(1)=el_min+el_range+1.
!        xcs(num)=xcs(num-1)
!        ycs(num)=el_min+el_range+1.
        !*output cs data
        open(unit=11,file="./output/readxsection_result.txt",status="unknown")
        !do i=1,num
        !    write(11,*)xcs(i),ycs(i) !changing all into m unit
        !enddo
        write(11,"(A5,6A20)") "il", "el1", "conv1", "tpW1", "uniflw",  "hbf_tz", "mann"

        il=1 !j=1
        xsec_attr_seg_g(1,il) = el_min
        xsec_attr_seg_g(2,il) = 0.0   !*5 -> 2
        xsec_attr_seg_g(3,il) = 0.0    !*6 -> 3
        xsec_attr_seg_g(4,il) = 0.0    !* for uniform flow

        do il= 2,nel_g
            el_now=el_min+real(il-1)*el_incr
                !if(abs(el_now - el_min) < TOLERANCE) then
                !    el_now=el_now+0.00001
                !end if
            el1(il)=el_now

            if (tzeq_flag_g==0) then
            !* use original readXsection method developed by Nazmul, Tulane University.
                i_start(1)=-999
                i_end(1)=-999
                i_area=0
                i_find=0
                do i=1,num-1
                    y1=ycs_g(i)
                    y2=ycs_g(i+1)
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
                    x1=xcs_g(i_start(i))
                    x2=xcs_g(i_start(i)+1)
                    y1=ycs_g(i_start(i))
                    y2=ycs_g(i_start(i)+1)
                    if(y1.eq.y2)then
                        x_start=x1
                    else
                        x_start=x1+(el_now-y1)/(y2-y1)*(x2-x1)
                    endif

                    x1=xcs_g(i_end(i))
                    x2=xcs_g(i_end(i)+1)
                    y1=ycs_g(i_end(i))
                    y2=ycs_g(i_end(i)+1)
                    if(y1.eq.y2)then
                      x_end=x1
                    else
                      x_end=x1+(el_now-y1)/(y2-y1)*(x2-x1)
                    endif

                    cal_topW=x_end-x_start+cal_topW

                    i1=i_start(i)
                    i2=i_end(i)
                    x1= xcs_g(i1+1); y1= ycs_g(i1+1)
                    x2= xcs_g(i2); y2= ycs_g(i2)
                    cal_area = cal_area    &
                             +cal_tri_area(el_now, x_start, x1, y1)    &
                             +cal_multi_area(el_now, i1+1, i2)    &
                             +cal_tri_area(el_now, x_end, x2, y2)

                    x1= xcs_g(i1+1); y1= ycs_g(i1+1)
                    x2= xcs_g(i2); y2= ycs_g(i2)
                    cal_peri = cal_peri    &
                            +cal_dist(x_start, el_now, x1, y1)    &
                            +cal_perimeter(i1+1, i2)    &
                            +cal_dist(x_end, el_now, x2, y2)

                    x1= xcs_g(i1+1); y1= ycs_g(i1+1)
                    if (i1.eq.1) then
                        cal_peri=cal_peri - cal_dist(x_start, el_now, x1, y1)
                    endif

                    x2= xcs_g(i2); y2= ycs_g(i2)
                    if (i2.eq.(num-1)) then
                        cal_peri=cal_peri-cal_dist(x_end, el_now, x2, y2)
                    endif
                enddo

                !el1(j)=el_now
                a1(il)=cal_area
                peri1(il)=cal_peri
                redi1(il)=a1(il)/peri1(il)

                !* when overbankflow, manning's N of floodplain is usually different
                !* from that of main channel, conveyance is computed separately.
                !* This computation applies only to trapezoidal main & rectangular floodplains
                hbf_tz=(tw_g - bo_g)/(2.0*traps_g )
                depth_tz= el_now - z_g
                if (depth_tz>hbf_tz) then
                !* overbank
                    call areacalc(depth_tz, area_tz)
                    call hydRcalc(depth_tz, area_tz, hydr_tz)
                    call Kcalc(depth_tz, area_tz, hydr_tz, conv_tz)
                    conv1(il)= conv_tz
                else
                !* inbank
                    conv1(il)=1./mann_g*a1(il)*(redi1(il))**(2./3.)
                endif
                tpW1(il)=cal_topW

            elseif (tzeq_flag_g==1) then
            !* use x-section equations for trapezoidal main and rectangular floodplain
                depth_tz= el_now - z_g
                call areacalc(depth_tz, area_tz)
                call hydRcalc(depth_tz, area_tz, hydr_tz)
                call Kcalc(depth_tz, area_tz, hydr_tz, conv_tz)
                hbf_tz= (tw_g - bo_g )/(2.0*traps_g) !* bankfull depth
                if (depth_tz>hbf_tz) then
                    topwd_tz= twcc_g
                else
                    topwd_tz= bo_g + 2.0*traps_g*depth_tz
                end if
                conv1(il)= conv_tz
                tpW1(il)=topwd_tz

            endif
            !* uniform flow
            uniflw= conv1(il)*So_g**0.5

!            if(j.eq.1) then
!                diffArea(j)=a1(j)
!                deffPere(j)=peri1(j)
!            else
!                diffArea(j)=a1(j)-a1(j-1)
!                deffPere(j)=peri1(j)-peri1(j-1)
!            endif
!
!            newdPdA(j)=deffPere(j)/diffArea(j)
!
!            waterElev=el1(j)
!            do jj=2,j
!                diffAreaCenter=el1(jj)-el_incr*0.5
!                newI1(j)=newI1(j)+diffArea(jj)*(waterElev-diffAreaCenter)
!            enddo
!            write(22,10)el1(j),a1(j),peri1(j),redi1(j),conv1(j),    &
!                       tpW1(j),newI1(j),newdPdA(j)
!    10      format(f9.2,3f12.2,2f20.3,2f16.3)
            !* new change 20200107 to make the model faster
            !* j: index for each water elevation
            !* k: segment ID
!            xsec_tab(1,j,k) = el1(j)
!            xsec_tab(2,j,k) = a1(j)
!            xsec_tab(3,j,k) = peri1(j)
!            xsec_tab(4,j,k) = redi1(j)
!            xsec_tab(5,j,k) = conv1(j)
!            xsec_tab(6,j,k) = tpW1(j)
!            xsec_tab(7,j,k) = newI1(j)
!            xsec_tab(8,j,k) = newdPdA(j)
!python:
            !* For this diffusive model only columns of 1, 5, and 6 are used.
            xsec_attr_seg_g(1,il) = el1(il)
!            xsec_tab(2,j,k) = a1(j)
!            xsec_tab(3,j,k) = peri1(j)
!            xsec_tab(4,j,k) = redi1(j)
            xsec_attr_seg_g(2,il) = conv1(il)   !*5 -> 2
            xsec_attr_seg_g(3,il) = tpW1(il)    !*6 -> 3
            xsec_attr_seg_g(4,il) = uniflw     !* uniform flow
!            xsec_tab(7,j,k) = newI1(j)
!            xsec_tab(8,j,k) = newdPdA(j)
            !write(11,*) j, el1(j),a1(j),peri1(j),redi1(j),conv1(j),tpW1(j), conv_tz, topwd_tz

            write(11,"(I5,6F20.8)") il, el1(il) ,conv1(il),tpW1(il), uniflw,  hbf_tz, mann_g
            !write(11,*) j, el1(j) ,conv1(j),tpW1(j), conv_tz, topwd_tz
        enddo
        !close(11)

    end subroutine read_xsec
!*--------------------------------------------
!           compute triangle area
!
!*--------------------------------------------
    real(KIND=dp) function cal_tri_area(el,x0,x1,y1)
        implicit none
        real(KIND=dp), intent(in) :: el, x0, x1, y1

        cal_tri_area=abs(0.5*(x1-x0)*(el-y1))

    end function cal_tri_area
!*--------------------------------------------
!          compute trapezoidal area
!
!*--------------------------------------------
    real(KIND=dp) function cal_trap_area(el,x1,y1,x2,y2)
        implicit none
        real(KIND=dp), intent(in) :: el,x1,y1,x2,y2

        cal_trap_area=abs(0.5*(x2-x1)*(el-y1+el-y2))

    end function cal_trap_area
!*--------------------------------------------
!      compute multiple trapezoidal areas
!
!*--------------------------------------------
    !real function cal_multi_area(el,xx,yy,n,i1,i2)
    real(KIND=dp) function cal_multi_area(el, i1, i2)
        implicit none
        integer(KIND=i4b) , intent(in) :: i1, i2
        real(KIND=dp), intent(in) :: el
        integer(KIND=i4b)  :: i
        real(KIND=dp) :: area, x1, y1, x2, y2

        area=0
        do i=i1,i2-1
            x1=xcs_g(i)
            y1=ycs_g(i)
            x2=xcs_g(i+1)
            y2=ycs_g(i+1)
            area=area+cal_trap_area(el,x1,y1,x2,y2)
        enddo
        cal_multi_area=area
    end function cal_multi_area
!*--------------------------------------------
!           compute linear distance
!
!*--------------------------------------------
    real(KIND=dp) function cal_dist(x1,y1,x2,y2)
        implicit none
        real(KIND=dp), intent(in) :: x1,y1,x2,y2

        cal_dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+1.e-32)

    end function cal_dist
!*--------------------------------------------
!           compute perimeter
!
!*--------------------------------------------
!    real function cal_perimeter(xx,yy,n,i1,i2)
    real(KIND=dp) function cal_perimeter(i1,i2)
        implicit none
        integer(KIND=i4b) , intent(in) :: i1, i2
        integer(KIND=i4b)  :: i
        real(KIND=dp) :: p, x1, y1, x2, y2

        p=0.0
        do i=i1,i2-1
            x1=xcs_g(i)
            y1=ycs_g(i)
            x2=xcs_g(i+1)
            y2=ycs_g(i+1)
            p=p+cal_dist(x1,y1,x2,y2)
        enddo
        cal_perimeter=p

    end function cal_perimeter
end module attrtable

