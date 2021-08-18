module moussa_network

	use constants_module
	use arrays_module
	use var_module
	use arrays_section_module
	use xsec_attribute_module

contains

	subroutine diffusive_CNT(j,ntim,repeatInterval)

		use subtools

		implicit none

		integer, intent(in) :: j, ntim, repeatInterval
		double precision,allocatable :: E_cnt(:,:), F_cnt(:,:)
		integer :: i, n, kk
		double precision :: hi, gi, ki, pj, qj, rj, pj_p, qj_p, rj_p, sj_p, mi, ni, qp_ghost, qp_ghost_1


	!!++++++++++++++++++++ Diffusive wave Forward sweep starts +++++++++++++++++++!!

		ncomp = frnw_g(j,1)

		allocate(E_cnt(ncomp,ntim))
		allocate(F_cnt(ncomp,ntim))

		qp_ghost_1 = qp(1,ntim,j)

		E_cnt(1:ncomp,1) = ini_E(1:ncomp,j)
		F_cnt(1:ncomp,1) = ini_F(1:ncomp,j)

		do i = 2,ncomp

			do n = 2, ntim

				hi = dx(i-1,j) / (dtini * celerity(i-1,j))
				gi = diffusivity(i-1,j) * dx(i-1,j) / ((celerity(i-1,j)**3.0)*(dtini**2.0))
				ki = (gi ** 2.0)/(hi ** 2.0)
				mi = diffusivity(i-1,j) / (celerity(i-1,j)**2.0) * dx(i-1,j) / dtini
				ni = dx(i-1,j)

				pj = - hi / 4.0 - gi / 2.0 - 2.0 * ki
				qj = 1.0 + gi + 4.0 * ki
				rj = hi / 4.0 - gi / 2.0 - 2.0 * ki

				pj_p = hi / 4.0 + gi / 2.0 - 2.0 * ki
				qj_p = 1.0 - gi + 4.0 * ki
				rj_p = -hi / 4.0 + gi / 2.0 - 2.0 * ki

				if (n .eq. ntim) then
					! applying ghost node !&
					sj_p = pj_p * qp(i-1,n-1,j) + qj_p * qp(i-1,n,j) + rj_p * qp_ghost_1 &
					+ mi / 2.0 * (lateralFlow(i-1,n,j) - lateralFlow(i-1,n-1,j) ) &
					+ ni * lateralFlow(i-1,n,j)
				else
					sj_p = pj_p * qp(i-1,n-1,j) + qj_p * qp(i-1,n,j) + rj_p * qp(i-1,n+1,j) &
					+ mi / 2.0 * (lateralFlow(i-1,n+1,j) - lateralFlow(i-1,n-1,j) ) &
					+ ni * lateralFlow(i-1,n,j)
				end if

				E_cnt(i,n) = -1.0 * rj / (pj * E_cnt(i,n-1) + qj)
				F_cnt(i,n) = ( sj_p - pj * F_cnt(i,n-1) ) / ( pj * E_cnt(i,n-1) + qj )

			end do

			qp_ghost = qp(i-1,ntim,j)+lateralFlow(i-1,ntim,j)*dx(i-1,j)

			qp_ghost_1 = qp_ghost

			! boundary at t=ntim, calculated from the ghost node at t=ntim+1
			qp(i,ntim,j) = E_cnt(i,ntim) * qp_ghost + F_cnt(i,ntim)

			do n = ntim-1,1,-1
				qp(i,n,j) = E_cnt(i,n) * qp(i,n+1,j) + F_cnt(i,n)

				if (qp(i,n,j) .lt. min_Q) then
					added_Q(i,n,j) = min_Q - qp(i,n,j)
					qp(i,n,j) = max(qp(i,n,j),min_Q)
				end if

			end do

		end do

		! replacing with the initial value at all nodes
		qp(1:ncomp,1,j) = ini_q_repeat(1:ncomp,j)

		! taking the new initial value for the next cycle
		ini_E(1:ncomp,j) = E_cnt(1:ncomp,repeatInterval+1)
		ini_F(1:ncomp,j) = F_cnt(1:ncomp,repeatInterval+1)
		ini_q_repeat(1:ncomp,j) = qp(1:ncomp,repeatInterval+1,j)

		deallocate (E_cnt, F_cnt)

	end subroutine diffusive_CNT

! ###########################################################################################
! ###########################################################################################
! ###########################################################################################

	subroutine mesh_diffusive_backward(j,ntim,repeatInterval)

		use subtools

		implicit none
		
		integer, intent(in) :: j,ntim,repeatInterval
		double precision :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
		double precision :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width, slope
		double precision :: cour, cour2, q_sk_multi, sfi, temp, dkdh
		double precision :: D_lim1, D_lim2, y_norm, y_crit, area_n, area_c, chnWidth, vel,y_norm1

		double precision :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ, stg1, stg2
		integer :: tableLength, jj, iii

		double precision :: elevTable_1(nel),areaTable_1(nel),rediTable_1(nel)
		double precision ::convTable_1(nel),topwTable_1(nel),currentSquaredDepth_1(nel)
		double precision :: dKdATable_1(nel)
		double precision :: pereTable_1(nel),depthYi,tempDepthi_1,tempCo_1,temp_q_sk_multi_1,tempY_1,tempArea_1,tempRadi_1,tempbo_1,ffy
		double precision :: ffy_1, ffy_2, ffy_3, tempCo_2
		double precision :: tempCo_3, tempsfi_2, tempsfi_3
		double precision :: ffprime,tempDepthi_1_new,tempsfi_1,toll, dkda, tempPere_1
		double precision :: tempdKdA_1, tempY_2, tempY_3, temp_v_1, tempDepthi_2
		double precision :: skkkTable_1(nel), tempsk_1
		double precision :: usFroud, dsFroud, dUdt, eHds, y_alt, area_alt, y_cnj, area_cnj, y_crit_test, y_norm_test
		integer :: depthCalOk(ncomp), wlCalcMethod
		integer :: i, pp
		double precision :: unity
		unity = 1.0		

	!!++++++++++++++++++++ Diffusive wave Backward sweep starts +++++++++++++++++++!!

			D_lim1 = -10.
			D_lim2 = -15.
			wlCalcMethod = 3

			! wlCalcMethod = 1 : Mid-point bisection method
			! wlCalcMethod = 2 : simple iterative method method with contribution from dU/dX
			! wlCalcMethod = 3 : Newton Raphson method
			! wlCalcMethod = 6 : Only Normal depth

			ncomp = frnw_g(j,1)

			depthCalOk(ncomp) = 1
			q_sk_multi = 1

			do i=ncomp,1,-1

				currentQ = qp(i,repeatInterval+1,j)
				! call calc_q_sk_multi(i,j,currentQ,q_sk_multi)

				! Calculating : read all attributes from tab file
				elevTable = xsec_tab(1,:,i,j)
				convTable = xsec_tab(5,:,i,j)
				areaTable = xsec_tab(2,:,i,j)
				pereTable = xsec_tab(3,:,i,j)
				topwTable = xsec_tab(6,:,i,j)
				skkkTable = xsec_tab(11,:,i,j)
				
				! interpolate the cross section attributes based on water elevation
				xt=newY(i,j)

				currentSquareDepth=(elevTable-z(i,j))**2.

				call r_interpol(currentSquareDepth,convTable,nel,(newY(i,j)-z(i,j))**2.0,co(i))
				if (co(i) .eq. -9999) then
					print*, 'At j = ',j,', i = ',i, 'interpolation of conveyence was not possible, wl', &
					newY(i,j), 'z',z(i,j),'previous wl',newY(i+1,j), 'previous z',z(i+1,j)
					stop
				end if
				co(i) =q_sk_multi * co(i)

				call r_interpol(elevTable,areaTable,nel,xt,newArea(i,j))
				if (newArea(i,j) .eq. -9999) then
					print*, 'At j = ',j,', i = ',i, 'interpolation of newArea was not possible'
					stop
				end if
				call r_interpol(elevTable,pereTable,nel,xt,pere(i,j))
				call r_interpol(elevTable,topwTable,nel,xt,bo(i,j))
				call r_interpol(elevTable,skkkTable,nel,xt,sk(i,j))

				sfi = qp(i,repeatInterval+1,j) * abs(qp(i,repeatInterval+1,j)) / ( co(i)** 2.0 )
				if (abs(qp(i,repeatInterval+1,j)) .lt. min_Q) then
					sfi = min_Q ** 2.0 * dsign(unity,qp(i,repeatInterval+1,j)) / ( co(i)** 2.0 )
					! at some head water basin, the Q boundary becomes 0.0 m3/s at some time. If Q = 0, then
					! sfi = 0. This leads to diffusivity = NaN. To avoid NaN diffusivity, those sfi is calculated using the min_Q
					! for those cases.
				end if

				chnWidth = rightBank(i,j)-leftBank(i,j)
				chnWidth = min(chnWidth,bo(i,j))

				if (depthCalOk(i) .eq. 1) then
					celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,repeatInterval+1,j)) ** &
						0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6
					diffusivity2(i) = abs(qp(i,repeatInterval+1,j)) / 2.0 / bo(i,j) / abs(sfi)
					vel = qp(i,repeatInterval+1,j)/newArea(i,j)

					velocity(i,j) = vel

					! Applying the upper limit of celerity
					if (celerity2(i) .gt. 3.0*vel) celerity2(i) = vel*3.0
				else
					if (qp(i,repeatInterval+1,j) .lt. 1) then
						celerity2(i)=0.5
					else
						celerity2(i)=1.0
					end if

					diffusivity2(i)=diffusivity(i,j)
				end if

				! Calculating Depth at the upstream node
				if (i .gt. 1) then

					slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
					depthYi = newY(i,j) - z(i,j)

					tempDepthi_1 = oldY(i-1,j)-z(i-1,j)

					elevTable_1 = xsec_tab(1,:,i-1,j)
					areaTable_1 = xsec_tab(2,:,i-1,j)
					pereTable_1 = xsec_tab(3,:,i-1,j)
					rediTable_1 = xsec_tab(4,:,i-1,j)
					convTable_1 = xsec_tab(5,:,i-1,j)
					topwTable_1 = xsec_tab(6,:,i-1,j)
					dKdATable_1 = xsec_tab(9,:,i-1,j)

					currentSquaredDepth_1=(elevTable_1-z(i-1,j))**2.

					toll = 1.0
					iii = 0

					dsFroud = abs(qp(i,repeatInterval+1,j))/sqrt(grav*newArea(i,j)**3.0/bo(i,j))

					! Applying mid point bisection
					if (wlCalcMethod .eq. 1) then
						tempY_1 = elevTable_1(2)
						tempY_2 = depthYi * 3. + z(i-1,j)
						tempY_3 = (tempY_1 + tempY_2) / 2.
						do while ( abs(toll) .gt. 0.001)
							iii = iii +1

							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_1-z(i-1,j))**2.0,tempCo_1)
							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_2-z(i-1,j))**2.0,tempCo_2)
							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_3-z(i-1,j))**2.0,tempCo_3)

							temp_q_sk_multi_1 = 1
							! call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
							tempCo_1 = tempCo_1 * temp_q_sk_multi_1
							tempCo_2 = tempCo_2 * temp_q_sk_multi_1
							tempCo_3 = tempCo_3 * temp_q_sk_multi_1

							tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )
							tempsfi_2 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_2** 2.0 )
							tempsfi_3 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_3** 2.0 )

							ffy_1 = (tempY_1-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1)
							ffy_2 = (tempY_2-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_2)
							ffy_3 = (tempY_3-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_3)

							if ((ffy_1 * ffy_2) .gt. 0.) then
								tempY_2 = (tempY_2 - z(i-1,j)) * 2.0 + z(i-1,j)
							elseif ((ffy_1 * ffy_3) .le. 0.) then
								tempY_2 = tempY_3
							elseif ((ffy_2 * ffy_3) .le. 0.) then
								tempY_1 = tempY_3
							end if
							tempY_3 = (tempY_1 + tempY_2) / 2.0
							toll = tempY_2 - tempY_1
							tempDepthi_1 = tempY_3 - z(i-1,j)
							depthCalOk(i-1) = 1

						end do
					end if

					! Applying simple iteration method with inclusion of dU/dX
					if (wlCalcMethod .eq. 2) then
						vel = qp(i,repeatInterval+1,j)/newArea(i,j)
						toll = 1.0

						do while ( abs(toll) .gt. 0.001)
							iii = iii +1
							tempY_1 = tempDepthi_1 + z(i-1,j)

							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempDepthi_1)**2.0,tempCo_1)

							temp_q_sk_multi_1 = 1
							! call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
							tempCo_1 = tempCo_1 * temp_q_sk_multi_1

							call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)

							temp_v_1 = qp(i-1,repeatInterval+1,j) / tempArea_1

							tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )!

							tempDepthi_2 = depthYi - dx(i-1,j)*slope + 0.5*dx(i-1,j)*(sfi + tempsfi_1)&
								+0.5*(vel+temp_v_1)*(vel-temp_v_1)/grav

							toll = tempDepthi_2 - tempDepthi_1

							tempDepthi_1 = tempDepthi_2

						end do
						depthCalOk(i-1) = 1

					end if

					! Applying NewtonRaphson method corrected
					if (wlCalcMethod .eq. 3) then
						vel = qp(i,repeatInterval+1,j)/newArea(i,j)
						do while ( abs(toll) .gt. 0.001)
							iii = iii +1
							tempY_1 = tempDepthi_1 + z(i-1,j)

							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempDepthi_1)**2.0,tempCo_1)
							temp_q_sk_multi_1 = 1
							! call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
							tempCo_1 = tempCo_1 * temp_q_sk_multi_1

							call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)
							call r_interpol(elevTable_1,pereTable_1,nel,tempY_1,tempPere_1)
							call r_interpol(elevTable_1,rediTable_1,nel,tempY_1,tempRadi_1)
							call r_interpol(elevTable_1,topwTable_1,nel,tempY_1,tempbo_1)
							call r_interpol(elevTable_1,dKdATable_1,nel,tempY_1,tempdKdA_1)!    Change 20210520

							tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )!

							temp_v_1 = qp(i-1,repeatInterval+1,j)/tempArea_1

							ffy = tempDepthi_1- depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1)

							dkda = tempdKdA_1

							ffprime = 1 + dx(i-1,j) * tempbo_1 *  qp(i-1,repeatInterval+1,j) * &
								abs(qp(i-1,repeatInterval+1,j)) / (tempCo_1 ** 3.0) * dkda

							if (ffprime .eq. 0.) then
								! print*, 'ffprime = 0'
								! print*, j, i, qp(i-1,repeatInterval+1,j), tempCo_1
								! pause
							end if

							tempDepthi_1_new = tempDepthi_1 - ffy / ffprime

							tempDepthi_1_new = max(tempDepthi_1_new,0.005)

							toll = abs(tempDepthi_1_new - tempDepthi_1)

							if(iii .gt. 30)then
								print*, 'Warning: Depth iteration reached maximum trial at j=', j, 'i=', i-1 , 'and',i, &
								'depths are',tempDepthi_1, tempDepthi_1_new, 'slope=', slope, 'tempbo_1=',tempbo_1, 'dx=',&
								dx(i-1,j), 'Q=', qp(i-1,repeatInterval+1,j)
								print*, 'depth at d/s',depthYi
								depthCalOk(i-1) = 0
								EXIT
							endif
							tempDepthi_1 = tempDepthi_1_new
							depthCalOk(i-1) = 1
						end do

					end if      ! calculation of Wl based on d/s wl ends


					usFroud = abs(qp(i-1,repeatInterval+1,j))/sqrt(grav*tempArea_1**3.0/tempbo_1)

					newY(i-1,j) = tempDepthi_1 + z(i-1,j)

					! calculating the normal depth at i=i-1 using the slope between i and i-1
					! Calculating the normal depth and critical depth at the river reach upstream as an output
					! very small slope is creating an issue. Applying a min limit of bed slope.
					! for Normal depth calculation, applying the absolute value of the Q
					call normal_crit_y(i-1, j, q_sk_multi, max(slope,0.0001), max(abs(qp(i-1,repeatInterval+1,j)),min_Q),&
						y_norm, y_crit, area_norm, area_crit)

					! to avoid sudden jump in water level calculation
					! if the WL jumps abruptly, the tempsfi_1 would be very small compared to sfi.
					! we are replacing the water level using normal water level
					if (sfi / tempsfi_1 .gt. 1e5)  newY(i-1,j) = y_norm

					! limiting very high value of WL
					if (newY(i-1,j) .gt. 10.0**5.) newY(i-1,j) = 10.0**5.

					if (newY(i-1,j)-z(i-1,j) .le. 0.) then
						! print*, 'depth is negative at j= ', j,'i=',i-1,'newY=',(newY(jj,j),jj=1,ncomp)
						! print*, 'dimensionless_D',(dimensionless_D(jj,j),jj=1,ncomp)
						! print*, 'newQ',(qp(jj,repeatInterval+1,j),jj=1,ncomp)
						! print*, 'Bed',(z(jj,j),jj=1,ncomp)
						! print*, 'dx',(dx(jj,j),jj=1,ncomp-1)
						! pause 777
					end if

				end if      ! end of if (i .gt. 1) || end of WL calculation at j reach

				! Book-keeping: Counting the number as for how many time steps the routing method is unchanged
				if (i .gt. 1) then
					routingNotChanged(i-1,j) = routingNotChanged(i-1,j) + 1
				end if 

			end do

			celerity(1:ncomp,j) =  sum(celerity2(1:ncomp)) / ncomp  ! change 20210524

			do i = 1, ncomp
				celerity(i,j) = max(celerity(i,j),0.5)   !!! Applying Celerity lower limit
			end do

			diffusivity(1:ncomp,j)=sum(diffusivity2(1:ncomp)) / ncomp

			do i = 1, ncomp
				if (diffusivity(i,j) .gt. maxDiffuLm) diffusivity(i,j) = maxDiffuLm !!! Applying diffusivity upper limit
				if (diffusivity(i,j) .lt. minDiffuLm) diffusivity(i,j) = minDiffuLm !!! Applying diffusivity lower limit
			end do

	!!++++++++++++++++++++ Diffusive wave Backward sweep ends +++++++++++++++++++!!

	end subroutine mesh_diffusive_backward

end module