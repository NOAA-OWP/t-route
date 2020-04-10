!            FINITE DIFFERENCE METHOD
!
!  A program for one dimensional flow in open channel
!
program mesh

    use constants
    use arrays
    use var
    use subtools

    implicit none

    integer :: n, i, j, k, linknb, nodenb
    real :: areanew, cour, da, dq, dnm0, dnm, dxs, x
    real :: t, sslp
    real :: Axs, Axs0, Axs1, dsc
    real :: d1junc, hav, h1z, dmy, dpth, areacr
    real :: areasum, hk_ncomp, areak_ncomp, qsum, qk_ncomp
    real :: rhs1, rhs2, c11, c12, c21, c22
    real :: vncomp, FrCC

    !*temporarily
    mxncomp=6;   !** Make sure that ncomp value should be the maximum value among
                   !** the numbers of nodes of links in consideration
    dtini=60.0;  dxini=20.0;  tfin=300060.0;  ntim=5001
    phi=0.0;  theta=1.0;  thetas=1.0;  thesinv=1.0;  alfa2=0.5;  alfa4=0.1
    f=1.0;  skk=33;  yy=2.423;  qq=100.0;  cfl=1.0;  ots=0;

    !*three links; make sure that the value of nlinks is also equivalent to link number in the most downstream.
    nlinks=3
    NAnum=-100  !* proxy for not avail
    mxnbrch=10  !* maximum number of branches that a lower link can hold
    nends=2     !* it counts upper and lower ends of each link

    call setup_arrays(ntim, mxncomp, nlinks, mxnbrch, nends)

    !+++---------------------------------------------------------------------------------+
    !+ Input channel geometries & lateral inflows; Initialize y & q;
    !+ ; Set up boundary conditions with known data, if exists.
    !+++---------------------------------------------------------------------------------+
    !* inputs from hydrofabric. For simplicity the link numbers for a branching channel,
    !* ,as shown in p.47, MESH of RM1, are 1, 2, and 3 for upper left, upper right, and
    !* lower links, respectively.
    !* the number of links that are immediately upstream of link j
    ndep(1)=0; ndep(2)=0; ndep(3)=2
    !* uslinks(k,j): k_th link ID that is immediately upstream of link j
    uslinks(1,1)=NAnum  !not available
    uslinks(1,2)=NAnum  !not available
    uslinks(1,3)=1;  uslinks(2,3)=2
    !* link ID that is immediately downstream of link j
    dslink(1)=3; dslink(2)=3; dslink(3)=NAnum
    !* the number of nodes in link j
    nx1(1)=5;  nx1(2)=5;  nx1(3)=5

    !* when data at either upper or lower end of link j is available,
    !* instrdflag(j,1)=1 when water level data is known at the upper end of link j
    !* instrdflag(j,1)=2 when discharge data is known
    !* instrdflag(j,1)=3 when rating curve is known
    !* instrdflag(j,2)=1 when water level data is known at the lower end of link j
    !* instrdflag(j,2)=2 when discharge data is known
    !* instrdflag(j,2)=3 when rating curve is known
    !* Otherwise, instrdflag(j,1/2)=0
    instrdflag(1,1)=2; instrdflag(1,2)=0    !*discharge data is known at the upper end of link 1
    instrdflag(2,1)=2; instrdflag(2,2)=0    !*discharge is known at the upper end of link 2
    instrdflag(3,1)=0; instrdflag(3,2)=1    !*stage data is known at the lower end of link 3

    ityp = 1

    do j=1, nlinks
        ncomp=nx1(j)
        !* channel width
        open(unit=100, file="./input/channel_width.txt")
        do i=1,ncomp
            read(100,*) x, bo(i,j)
        end do
        close(100)
        !* channel inbank and overbank flow top width
        do i=1, ncomp
            Tw(i,j)= 100.0 + 4.0*4.50 !* the last number is assumed bankfull depth
            TwCC(i,j)= Tw(i,j)+100.0
        end do
        !* bed elevation
        open(unit=97, file="./input/bed_elevation.txt")
        do i=1,ncomp
            read(97, *) x, z(i,j)
        end do
        close(97)

        do i=1, ncomp
            !* Manning's n: skk=1/n
            sk(i,j) = skk
            !* spatial interval
            dx(i,j) = dxini
            !* side slope of trapezoidal channel
            traps(i,j)=2.0
            !* lateral inflows if exist
            do n=1, ntim
                qlat(n,i,j)=0.0
            end do
        end do
        !++--------------------------------------------------------------+
        !+ Assume that hydrographs are always available at the upper ends
        !+ of the most upstream links.
        !++--------------------------------------------------------------+
        i=1 !*first node of a link
        if (instrdflag(j,1)==1) then
        !* stage data at the upper end of link j
        elseif (instrdflag(j,1)==2) then
        !* discharge data at the upper end of link j
            open(unit=981, file="./input/upstream_hydrograph.txt")
            do n=1,ntim
                read(981,*) t, q(n, 1, j)
            !* double the size of discharge for link j=2 for a test
            if (j==2) then
                q(n,i,j)=2.0*q(n,i,j-1)
            end if

            end do
            close(981)
            !* initialize q along all the nodes of link j with the measured data at n=1
            !* until the link meets a junction.
            n=1
            do i=1, ncomp
                q(n,i,j) = q(n,1,j)
            end do
        elseif (instrdflag(j,1)==3) then
            !* Add lines here for rating curve boundary condition!!
        end if
        !++--------------------------------------------------------------+
        !+ when stage data of lower end of link j are known, plug them in.
        !+ Assume that stage time-series data is always available at the
        !+ lower end of the most downstream link.
        !++--------------------------------------------------------------+
        i=ncomp
        if (instrdflag(j,2)==1) then
        !* stage data at the lower end of link j. Assume that stage data is always measured from channel bottom.
            open(unit=99, file="./input/downstream_stage.txt")
            do n=1,ntim
                read(99,*) t, dmy
                y(n, ncomp, j)=dmy+z(i,j)  !*y is stage from a datum not channel bottom
            end do
            close(99)
        elseif (instrdflag(j,2)==2) then
        !* discharge data at the lower end of link j
        elseif (instrdflag(j,2)==3) then
        !* Add lines here for rating curve boundary condition!!
        end if
    end do  !*do j=1, nlinks


    !+++---------------------------------------------------------------------------------+
    !+ Treat q for link j that receive water from more than one immediate upstream link
    !+++---------------------------------------------------------------------------------+
    do j=1, nlinks
        ncomp=nx1(j)
        if (ndep(j)>=1) then    !*here >=1 mean that this code treat serial reaches let alone branching reaches
        n=1
            q(n,1,j)=0
            do k=1, ndep(j)
                linknb=uslinks(k,j) !*link number of k_th link in the immediate upstream of link j
                q(n,1,j)=q(n,1,j) + q(n, nx1(linknb), linknb)
            end do

            do i=2, ncomp
                q(n,i,j) = q(n,1,j)
            end do
        end if
    end do

    !+++----------------------------------------------------------------------------------------+
    !+                   INITIAL ESTIMATE OF y AT TIME=1, p68-2~3,RM1
    !+ Up to this point, all nodes of links that are connected at a common junction already take
    !+ initial estimates of q.  Now, initial values of y is estimated with the assumption that
    !+ velocity is consistent throughout all the nodes of a link and its value is computed by
    !+ Q=AV where A is computed using known stage value at the last node of the most
    !+ downstream link.
    !+++----------------------------------------------------------------------------------------+
    n=1
    !* First, determine initial velocity with known stage data at the end node of the most downstream link.
    do j=1, nlinks
        if (instrdflag(j,2)==1) then
        !* when water depth data is available at the lower end of a link
            !* velocity at node ncomp, p68-2,RM1
            ncomp=nx1(j)
            dpth=y(n,ncomp,j) - z(ncomp,j)
            call areacalc(ncomp, j, dpth, Axs)
            vncomp=q(n,ncomp,j)/Axs
            !* estimate initial velocity at the rest nodes of the link.
            do i=1, ncomp-1
                Axs=q(n,i,j)/vncomp
                call depthcalc(i, j, Axs, dpth)
                y(n,i,j)=dpth+z(i,j)
            end do
            d1junc=y(n,1,j)-z(1,j)
        end if
    end do

    !* Second, estimate initial velocity through adjacent links connected to a common junction, p68-3,RM1
    do j=1,nlinks
        if (instrdflag(j,2)/=1) then
            !* velocity at node ncomp
            ncomp=nx1(j)
            dpth=d1junc
            call areacalc(ncomp, j, dpth, Axs)
            vncomp=q(n,ncomp,j)/Axs

            do i=1, ncomp
                Axs=q(n,i,j)/vncomp
                call depthcalc(i, j, Axs, dpth)
                y(n,i,j)=dpth +  z(i,j)
            end do
        end if
    end do

    open(unit=85, file="./output/areafinal.txt", status="unknown")
    open(unit=101, file="./output/f_depth0.txt", status="unknown")
    open(unit=102, file="./output/f_q0.txt", status="unknown")
    !
    ! Loop on time
    !
    do n=1, ntim-1
    !+++---------------------------------------------------------------------------------+
    !+                                     PREDICTOR
    !+
    !+
    !+                  Run predictor though all the links for a given time n
    !+                  with previously computed y(->area) and q variables
    !+++---------------------------------------------------------------------------------+
        do j=1, nlinks
            ncomp=nx1(j)

            call section(n,j)

            !++--------------------------------------------------------------------------------+
            !+                     UPSTREAM BOUNDARY CONDTION
            !+
            !+ Assuming that discharge values in the upper end of a river network
            !+ are always known, determine whether upstream flow is supercritical
            !+ or subcritical based on directly Froud number of inbank or overbank flow.
            !+ Predictor step must start first from links where measured or estimated discharge
            !+ data at time n and n+1 are available, for example, headbasin streams.
            !++--------------------------------------------------------------------------------+
            if (instrdflag(j,1)==2) then
            !* Froud number at the first node of a link where measured discharge
            !* at current time n is available.
            i=1
                sslp= traps(i,j)
                dsc= q(n,i,j)
                !* current depth
                dpth= y(n,i,j)-z(i,j)
                Axs=area(i,j)   !* result from sub section
                call FrCCcalc(i, j, dpth, Axs, dsc, FrCC)

                if (FrCC.gt.1.0) then
                !* supercritical flow
                    !** dqp(1,j)
                    dqp(i,j)=q(n+1,i,j)-q(n,i,j)
                    !** dap(1,j)
                    !* Area(1,j) at time n
                    Axs0=area(i,j)
                    !* Area(1,j) at time n+1.
                    !* Normal depth is estimated from future discharge at time n+1, which may be a proper candidate of y(n+1).
                    dnm0=y(n,i,j)-z(i,j)  !*only for Newton-Rapson or Bisection method
                    dsc=q(n+1,i,j)
                    call nmaldepcalc(i, j, dnm0, dsc, dnm)
                    call areacalc(i, j, dnm, Axs)
                    Axs1=Axs
                    dap(i,j)=Axs1-Axs0
                else
                !* subcritical flow
                    dqp(i,j)=q(n+1,i,j)-q(n,i,j)
                    dap(i,j)=0.0
                end if
            end if

            thes=thetas

            call matrixp(n,j)

            !+++-----------------------------------------------------------------------------
            !+             INTERNAL JUNCTION BOUNDARY CONDITION in Predictor step
            !+
            !+ Hand over water from upstream to downstream properly according
            !+ to the nature of link connections, i.e., branching.
            !+ Refer to p.52,RM1_MESH
            !+++-----------------------------------------------------------------------------
            if (ndep(j).gt.1) then  !* ndep(j): the number of links that are immediately upstream of link j.
                !*total water areas at n+1 at the end nodes of upstream links that join link j
                areasum=0.0
                qsum= 0.0
                do k=1, ndep(j)
                    linknb=uslinks(k,j); nodenb=nx1(linknb)
                    areasum=areasum + area(nodenb,linknb) + dap(nodenb,linknb)
                    qsum= qsum + q(n,nodenb,linknb) + dqp(nodenb,linknb)
                end do

                dqp(1,j)=0.0;
                hav=0.0
                do k=1, ndep(j)
                    linknb=uslinks(k,j); nodenb=nx1(linknb)
                    !**dqp(1,j)
                    !* discharge at the first node of link j, i.e.,dqp(1,j), is the sum of discharges of at
                    !* the last nodes of upstream links that connect to link j at their junction.
                    dqp(1,j)=dqp(1,j)+dqp(nodenb,linknb)
                    !**dap(1,j)
                    !*area at the end nod of link k at time n+1
                    areak_ncomp = area(nodenb,linknb) + dap(nodenb,linknb)  !*may cause different chshp determined in sub section
                    qk_ncomp=  q(n,nodenb,linknb) + dqp(nodenb,linknb)
                    !* depth at the end nod of link k at time n+1
                    call depthcalc(nodenb, linknb, areak_ncomp, dxs)
                    hk_ncomp= dxs
                    !* weighted average based on areas at the end nodes of upstream link ks
                    hav = hav + (areak_ncomp/areasum)*hk_ncomp
                    !* weighted average based on discharge at the end nodes of upstream link ks
                    !hav = hav + (qk_ncomp/qsum)*hk_ncomp
                end do
                !* Area estimated for time n+1
                i=1
                call areacalc(i, j, hav, Axs)
                !* update dap at the first node of lower link of a junction, i.e., dap(1,j)
                !* **Note:  dap at the first node of lower link of a junction is updated before
                !*          the update of dap and dqp in the next do loop, while dap at the last
                !*          nodes of upper links at the junction are updated after the do loop
                dap(1,j)= Axs - area(1,j)
            end if

            do i=2,ncomp
                cour=dt(i)/dx(i-1,j)
                !* rhs1 with lateral inflow
                rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1)) + 0.5*dt(i)*(qlat(n,i,j)+qlat(n,i-1,j))
                rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1))+dt(i)*grav*(ci2(i)+aso(i))
                c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
                c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
                c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
                c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
                dap(i,j)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1,j)-c12*dqp(i-1,j)
                dqp(i,j)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1,j)-c22*dqp(i-1,j)
            end do

            !+++----------------------------------------------------------------------------
            !+          INTERNAL JUNCTION BOUNDARY CONDITION in Predictor step
            !+
            !+ Hand over water from upstream to downstream properly according
            !+ to the nature of link connections, i.e., branching.
            !+ Refer to p.52,RM1_MESH
            !+++-----------------------------------------------------------------------------
            if (ndep(j).gt.1) then  !* ndep(j): the number of links that are immediately upstream of link j.
                !* **Note:  dap at the first node of lower link of a junction is updated before
                !*          the update of dap and dqp in the previous do loop, while dap at the last
                !*          nodes of upper links at the junction are updated after the do loop.
                do k=1, ndep(j)
                    linknb=uslinks(k,j); nodenb=nx1(linknb)
                    call areacalc(nodenb, linknb, hav, Axs)
                    dap(nodenb,linknb)= Axs - area(nodenb,linknb)
                    !* update areap as well
                    areap(nodenb,linknb) = area(nodenb,linknb) + dap(nodenb,linknb)
                end do
            end if

            !+++------------------------------------------------------------+
            !+ Update via predictor for each link
            !+++------------------------------------------------------------+
            do i=1, ncomp
                areap(i,j) = area(i,j) + dap(i,j)
                qp(i,j) = q(n,i,j) + dqp(i,j)
            end do

            !+++-----------------------------------------------------------------+
            !+         Downstream Boundary Condition in Predictor step
            !+
            !+ ASSUMPTION: stage data at future time n+1 is available at this node
            !+
            !+++-----------------------------------------------------------------+
            if (instrdflag(j,2)==1) then            !
                dxs= y(n+1,ncomp,j)-z(ncomp,j)
                call areacalc(ncomp, j, dxs, Axs)
                areap(ncomp,j)= Axs
                dap(ncomp,j)= areap(ncomp,j) - area(ncomp,j)
                !* downstream boundary conditions of the most downstream link
                !* eventually for corrector step
                call dsbc(n, j)
            end if
        end do  !*do j=1, nlinks


    !+++---------------------------------------------------------------------------------+
    !+                                    CORRECTOR
    !+
    !+
    !+                  Run predictor though all the links for a given time n
    !+                  with PREDICTOR-computed areap and qp
    !+++---------------------------------------------------------------------------------+
        do j=nlinks,1,-1
            ncomp=nx1(j)

            call secpred(j)

            thes=thesinv

            call matrixc(n,j)

            !+++----------------------------------------------------------------------------
            !+          Corrector: INTERNAL JUNCTION BOUNDARY CONDITION
            !+
            !+ Handle internal downstream boundary condition for a link j that has a link
            !+ immediately downstream, either in serial or branching channels.
            !+ **Note. dac and dqc at the last node of the most downstream link are computed
            !+ in sub dsbc during predictor step.
            !+ Refer to p.53-1,RM1_MESH
            !+++----------------------------------------------------------------------------
            if (dslink(j).ne.NAnum) then
            !*when downstream link is available,
                !* link ID in the downstream of link j
                linknb=dslink(j)
                !*dac(ncomp,j)
                areacr=area(1,linknb)+0.5*(dap(1,linknb)+dac(1,linknb))
                !* depth at the last node of link j for time n+1 estimated from the first node of the next downstream link
                i=1
                call depthcalc(i, linknb, areacr, dxs)
                h1z= dxs
                call areacalc(ncomp, j, h1z, Axs)
                dac(ncomp,j)= 2.0*(Axs - area(ncomp,j)) - dap(ncomp,j)
                dqc(ncomp,j)=dqc(1,linknb)*q(n,ncomp,j)/q(n,1,linknb)
            end if

            do i=ncomp-1,1,-1
                cour=dt(i)/dx(i,j)
                !* rhs1 with lateral inflow
                rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+ 0.5*dt(i)*(qlat(n,i+1,j)+qlat(n,i,j))
                rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dt(i)*grav*(ci2(i)+aso(i))
                c11=g11inv(i)*b11(i+1)+g12inv(i)*b21(i+1)
                c12=g11inv(i)*b12(i+1)+g12inv(i)*b22(i+1)
                c21=g21inv(i)*b11(i+1)+g22inv(i)*b21(i+1)
                c22=g21inv(i)*b12(i+1)+g22inv(i)*b22(i+1)
                dac(i,j)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dac(i+1,j)-c12*dqc(i+1,j)
                dqc(i,j)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dac(i+1,j)-c22*dqc(i+1,j)
            end do

            !* when discharge data is known at upper ends of links, use the info
            if (instrdflag(j,1)==2) then
                dqc(1,j)=dqp(1,j)
            end if

            ! Final update
            do i=1,ncomp
                !* update y(n+1,i,j)
                da=(dap(i,j)+dac(i,j))/2.0
                areanew=da+area(i,j)
                if(areanew <= 0.0) areanew=0.0001
                areafnal(n,i,j)=areanew
                call depthcalc(i, j, areanew, dpth)
                y(n+1,i,j)= dpth + z(i,j)
                !* update q(n+1,i,j)
                dq=(dqp(i,j)+dqc(i,j))/2.0
                q(n+1,i,j)=q(n,i,j)+dq
            end do
            t = t + dtini
        end do  !*do j=nlinks,1,-1
    end do  !*do n=1, ntim-1
                        do n=1,ntim
                        do j=1,nlinks
                            write(101,"(2I6,5F15.3)") n, j, (y(n, i, j)-z(i,j),i=1,ncomp)
                            write(102,"(2I6,5F15.3)") n, j, (q(n, i, j),i=1,ncomp)
                            write(85,"(2I6,5F15.3)") n, j, (areafnal(n,i, j),i=1,ncomp)
                        end do
                        end do
                        close(85); close(101); close(102)

end program mesh
