module lookupTB
    use var
    implicit none

contains

    subroutine readXsection(k,lftBnkMann,rmanning_main,rgtBnkMann,leftBnkX_given,rghtBnkX_given,timesDepth,num_reach)
        implicit none
        save

        integer, intent(in) :: k, num_reach
        real, intent(in) :: rmanning_main,lftBnkMann,rgtBnkMann,leftBnkX_given,rghtBnkX_given, timesDepth
        real, dimension(:), allocatable :: xcs, ycs
        real, dimension(:,:), allocatable :: el1, a1, peri1, redi1
        real, dimension(:), allocatable :: redi1All
        real, dimension(:,:), allocatable :: conv1, tpW1, diffArea, newI1, diffPere
        real, dimension(:), allocatable :: newdPdA, diffAreaAll, diffPereAll, newdKdA       ! change Nazmul 20210601
        real, dimension(:), allocatable :: compoundSKK, elev
        integer, dimension(:), allocatable :: i_start, i_end, totalNodes
        real, dimension(:,:), allocatable :: allXcs, allYcs
        integer :: i_area, i_find, i, j, jj, num
        real:: el_min, el_max, el_range, el_incr, el_now, x1, y1, x2, y2, x_start, x_end, waterElev, leftBnkX,rghtBnkX
        real:: f2m, cal_area, cal_peri, cal_topW,  diffAreaCenter
        real:: compoundMann, el_min_1
        integer:: i1, i2
        real:: leftBnkY, rghtBnkY,rmanning
        integer:: mainChanStrt, mainChanEnd,  kkk, startFound, endFound
        real :: hbf

        allocate (el1(nel,3),a1(nel,3),peri1(nel,3),redi1(nel,3),redi1All(nel))
        allocate (conv1(nel,3), tpW1(nel,3), diffArea(nel,3), newI1(nel,3), diffPere(nel,3))
        allocate (newdPdA(nel), diffAreaAll(nel), diffPereAll(nel), newdKdA(nel))       ! change Nazmul 20210601
        allocate (compoundSKK(nel), elev(nel))
        allocate (i_start(nel), i_end(nel))
        allocate (totalNodes(3))


        leftBnkX=leftBnkX_given
        rghtBnkX=rghtBnkX_given
        startFound = 0
        endFound = 0
        !* channel geometry at a given segment
        z_g = z_ar_g(k, num_reach)
        bo_g=bo_ar_g(k, num_reach)
        traps_g=traps_ar_g(k, num_reach)
        tw_g= tw_ar_g(k, num_reach)
        twcc_g= twcc_ar_g(k, num_reach)
        hbf= (tw_g-bo_g)/(2.0*traps_g) !* bankfull depth
        maxTableLength=8
        f2m=1.0
        allocate (xcs(maxTableLength), ycs(maxTableLength))
        allocate (allXcs(maxTableLength,3), allYcs(maxTableLength,3))
        do i=1, maxTableLength
            !* channel x-section vertices at a given segment
            if (i==1) then
                x1=0.0; y1=z_g + timesDepth*hbf
            elseif (i==2) then
                x1=0.0; y1=z_g + hbf
            elseif (i==3) then
                x1=(twcc_g-tw_g)/2.0; y1= z_g + hbf
            elseif (i==4) then
                x1= xcs(3) + traps_g*hbf; y1= z_g
            elseif (i==5) then
                x1= xcs(4) + bo_g; y1= z_g
            elseif (i==6) then
                x1= xcs(5) + traps_g*hbf; y1= z_g + hbf
            elseif (i==7) then
                x1= twcc_g; y1= z_g + hbf
            elseif (i==8) then
                x1= xcs(7); y1= z_g + timesDepth*hbf
            endif

            xcs(i)=x1*f2m
            ycs(i)=y1*f2m
            if ((xcs(i) .ge. leftBnkX) .and. (startFound .eq. 0)) then
                mainChanStrt = i-1
                startFound = 1
            end if
            if ((xcs(i) .ge. rghtBnkX) .and. (endFound .eq. 0)) then
                mainChanEnd = i-1
                endFound = 1
            end if
        enddo
        mainChanStrt=3
        mainChanEnd=6
        num=i

        if (leftBnkX .lt. minval(xcs(2:num-1))) leftBnkX = minval(xcs(2:num-1))
        if (rghtBnkX .gt. maxval(xcs(2:num-1))) rghtBnkX = maxval(xcs(2:num-1))

        leftBnkY = ycs(mainChanStrt)+(leftBnkX-xcs(mainChanStrt))/&
          (xcs(mainChanStrt+1)-xcs(mainChanStrt))*(ycs(mainChanStrt+1)-ycs(mainChanStrt))
        rghtBnkY = ycs(mainChanEnd)+(rghtBnkX-xcs(mainChanEnd))/&
          (xcs(mainChanEnd+1)-xcs(mainChanEnd))*(ycs(mainChanEnd+1)-ycs(mainChanEnd))
        el_min=99999.
        el_max=-99999.
        do i=2,num-1
            if(ycs(i).lt.el_min)el_min=ycs(i)
            if(ycs(i).gt.el_max)el_max=ycs(i)
        enddo
        el_range=(el_max-el_min)*2.0 ! change Nazmul 20210601

        do i=1, 3
            allXcs(i+1,1)=xcs(i) !x1*f2m
            allYcs(i+1,1)=ycs(i) !y1*f2m
        enddo
        allXcs(1,1)=xcs(1)
        allYcs(1,1)=el_min+el_range+1.
        allXcs(mainChanStrt+2,1)=xcs(3)
        allYcs(mainChanStrt+2,1)=el_min+el_range+1.

        do i=3,4
            allXcs(i-1,2)=xcs(i) !x1*f2m
            allYcs(i-1,2)=ycs(i) !y1*f2m
        enddo

        do i=5,6
            allXcs(i,2)=xcs(i) !x1*f2m
            allYcs(i,2)=ycs(i) !y1*f2m
        enddo
        allXcs(1,2)=xcs(3)
        allYcs(1,2)=el_min+el_range+1.
        allXcs(7,2)=xcs(6)
        allYcs(7,2)=el_min+el_range+1.

        do i=6,8
            allXcs(i-4,3)=xcs(i) !x1*f2m
            allYcs(i-4,3)=ycs(i) !y1*f2m
        enddo
        allXcs(1,3) = allXcs(2,3)
        allYcs(1,3) = el_min+el_range+1.
        i=5
        allXcs(i,3) = allXcs(i-1,3)
        allYcs(i,3) = el_min+el_range+1.

        totalNodes(1) = 5
        totalNodes(2) = 7
        totalNodes(3) = 5

        allXcs(4,2) = (allXcs(3,2)+allXcs(5,2))/2.0
        allYcs(4,2) = allYcs(3,2) - 0.01

        el_min_1 = el_min
        el_min = allYcs(4,2)    ! change Nazmul 20210601 ! corrected

        elev(1) = el_min
        elev(2) = el_min + 0.01/4.
        elev(3) = el_min + 0.01/4.*2.
        elev(4) = el_min + 0.01/4.*3.
        elev(5) = el_min + 0.01

        el_incr=el_range/real(nel-6.0)

        do kkk = 6,nel
            elev(kkk) = elev(5)+el_incr * (kkk-5)
        end do

        xcs = 0.
        ycs = 0.
        newI1=0.0 !Hu changed
        do kkk=1,3
            num = totalNodes(kkk)
            xcs(1:num) = allXcs(1:num,kkk)
            ycs(1:num) = allYcs(1:num,kkk)
            if (kkk .eq. 1) rmanning = lftBnkMann
            if (kkk .eq. 2) rmanning = rmanning_main
            if (kkk .eq. 3) rmanning = rgtBnkMann
            do j=1,nel
                el_now = elev(j)
                if(abs(el_now - el_min) < TOLERANCE) then
                    el_now=el_now+0.00001
                end if
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
                    write(11,*)x_start,el_now

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

                    write(11,*)x_end,el_now
                    write(11,*)'NaN NaN'

                    i1=i_start(i)
                    i2=i_end(i)

                    cal_area = cal_area    &
                             +cal_tri_area(el_now,x_start,xcs(i1+1),ycs(i1+1))    &
                             +cal_multi_area(el_now,xcs,ycs,maxTableLength,i1+1,i2)    &
                             +cal_tri_area(el_now,x_end,xcs(i2),ycs(i2))
                    cal_peri = cal_peri    &
                            +cal_dist(x_start,el_now,xcs(i1+1),ycs(i1+1))    &
                            +cal_perimeter(xcs,ycs,maxTableLength,i1+1,i2)    &
                            +cal_dist(x_end,el_now,xcs(i2),ycs(i2))
                    if(i1.eq.1)cal_peri=cal_peri    &
                             -cal_dist(x_start,el_now,xcs(i1+1),ycs(i1+1))
                    if(i2.eq.(num-1))cal_peri=cal_peri    &
                             -cal_dist(x_end,el_now,xcs(i2),ycs(i2))

                enddo

                el1(j,kkk)=el_now
                a1(j,kkk)=cal_area
                peri1(j,kkk)=cal_peri
                redi1(j,kkk)=a1(j,kkk)/peri1(j,kkk)

                conv1(j,kkk)=1./rmanning*a1(j,kkk)*(redi1(j,kkk))**(2./3.)
                !print*, j, kkk, conv1(j,kkk)
                if (peri1(j,kkk) .le. TOLERANCE) then
                    redi1(j,kkk) =0.0; conv1(j,kkk)=0.0
                endif
                tpW1(j,kkk)=cal_topW

                if(j.eq.1) then !* Dongha added
                    diffArea(j,kkk)=a1(j,kkk) !* Dongha added
                    diffPere(j,kkk)=peri1(j,kkk) !* Dongha added
                else
                    if (el_now .le. minval(ycs(1:num))) then
                      diffArea(j,kkk)=a1(j,kkk)
                      diffPere(j,kkk)=peri1(j,kkk)
                    else
                      diffArea(j,kkk)=a1(j,kkk)-a1(j-1,kkk)
                      diffPere(j,kkk)=peri1(j,kkk)-peri1(j-1,kkk)
                    endif
                endif

                waterElev=el1(j,kkk)
                do jj=2,j
                  diffAreaCenter=el1(jj,kkk)-(el1(jj,kkk)-el1(jj-1,kkk))*0.5
                  newI1(j,kkk)=newI1(j,kkk)+diffArea(jj,kkk)*(waterElev-diffAreaCenter)
                enddo
            end do
        end do

        do j = 1,nel
            el_now=el1(j,1)
            if (j .eq. 1) then
                newdPdA(j) = sum(peri1(j,:)) / sum(a1(j,:))
                newdKdA(j) = sum(conv1(j,:)) / sum(a1(j,:))     ! change Nazmul 20210601
            else
                newdPdA(j)= (sum(peri1(j,:)) - sum(peri1(j-1,:))) / (sum(a1(j,:)) - sum(a1(j-1,:)))
                newdKdA(j)= (sum(conv1(j,:)) - sum(conv1(j-1,:))) / (sum(a1(j,:)) - sum(a1(j-1,:)))
            end if

            compoundMann = sqrt((abs(peri1(j,1))*lftBnkMann ** 2. + abs(peri1(j,2))*rmanning_main ** 2.+&
             abs(peri1(j,3))*rgtBnkMann ** 2.) / (abs(peri1(j,1))+abs(peri1(j,2))+abs(peri1(j,3))))
            compoundSKK(j) = 1. / compoundMann

            redi1All(j)=sum(a1(j,:)) /sum(peri1(j,:))
            xsec_tab(1,j,k,num_reach) = el1(j,1)
            xsec_tab(2,j,k,num_reach) = sum(a1(j,:))
            xsec_tab(3,j,k,num_reach) = sum(peri1(j,:))
            xsec_tab(4,j,k,num_reach) = redi1All(j)
            xsec_tab(5,j,k,num_reach) = sum(conv1(j,:))
            xsec_tab(6,j,k,num_reach) = abs(tpW1(j,1))+abs(tpW1(j,2))+abs(tpW1(j,3))
            xsec_tab(7,j,k,num_reach) = sum(newI1(j,:))
            xsec_tab(8,j,k,num_reach) = newdPdA(j)
            xsec_tab(9,j,k,num_reach) = newdKdA(j)
            xsec_tab(11,j,k,num_reach) = compoundSKK(j)
        end do
        z(k,num_reach)= el_min

        deallocate (el1, a1, peri1, redi1, redi1All)
        deallocate (conv1, tpW1, diffArea, newI1, diffPere)
        deallocate (newdPdA, diffAreaAll, diffPereAll, newdKdA)       ! change Nazmul 20210601
        deallocate (compoundSKK, elev)
        deallocate (i_start, i_end)
        deallocate (totalNodes)
        deallocate (xcs, ycs)
        deallocate (allXcs, allYcs)

    end subroutine


    real function cal_tri_area(el,x0,x1,y1)
        implicit none
          real, intent(in) :: el,x0,x1,y1

          cal_tri_area=abs(0.5*(x1-x0)*(el-y1))
          return
    end function cal_tri_area


    real function cal_trap_area(el,x1,y1,x2,y2)
        implicit none
        real, intent(in) :: el,x1,y1,x2,y2
        cal_trap_area=abs(0.5*(x2-x1)*(el-y1+el-y2))
        return
    end function cal_trap_area

    real function cal_multi_area(el,xx,yy,n,i1,i2)
        implicit none
        integer, intent(in) :: n,i1,i2
        real, intent(in) :: el
        real, intent(in) :: xx(n),yy(n)
        integer :: i
        real :: area, x1, x2, y1, y2

        area=0
        do i=i1,i2-1
            x1=xx(i)
            y1=yy(i)
            x2=xx(i+1)
            y2=yy(i+1)
            area=area+cal_trap_area(el,x1,y1,x2,y2)
        enddo
        cal_multi_area=area
        return
    endfunction cal_multi_area

    real function cal_dist(x1,y1,x2,y2)
        implicit none
        real, intent(in) :: x1,y1,x2,y2
        cal_dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+1.e-32)
        return
    end function cal_dist

    real function cal_perimeter(xx,yy,n,i1,i2)
        implicit none
        integer, intent(in) :: n,i1,i2
        real, intent(in) :: xx(n),yy(n)
        integer :: i
        real :: p, x1, x2, y1, y2

        p=0.
        do i=i1,i2-1
            x1=xx(i)
            y1=yy(i)
            x2=xx(i+1)
            y2=yy(i+1)
            p=p+cal_dist(x1,y1,x2,y2)
        enddo
        cal_perimeter=p
        return
    endfunction cal_perimeter
endmodule lookupTB
