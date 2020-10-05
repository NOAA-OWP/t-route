program test_MC_program

use precis
use muskingcunge_module

    implicit none
    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios = 0
    integer :: line = 0

    ! variables in the .txt file to be read-in
    real(prec) :: dt_in, dx_in, bw_in, tw_in, twcc_in, n_in, ncc_in, cs_in, s0_in, ql_in, &
      qup_in, quc_in, qdp_in, velp_in, depthp_in    
    
    ! variables output by muskingcunge subroutine
    real(prec) :: qdc_o, velc_o, depthc_o, ck_o, cn_o, X_o
    
    ! ###############################################################################
    ! read parameters from control file
    ! ###############################################################################
    ! open file containing MC parameter variables
    open(fh, file='control.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    
    !print*, "--------------------------------------------------------------"
    !print*, " "
    !print*, "MC model inputs"
    !print*, " "
    
    do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, ' 	')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)
        case ('dt')
           read(buffer, *, iostat=ios) dt_in
           !print *, 'dt: ', dt_in
        case ('dx')
           read(buffer, *, iostat=ios) dx_in
           !print *, 'dx: ', dx_in
        case ('bw')
           read(buffer, *, iostat=ios) bw_in
           !print *, 'bw: ', bw_in
        case ('tw')
           read(buffer, *, iostat=ios) tw_in
           !print *, 'tw: ', tw_in          
        case ('twcc')
           read(buffer, *, iostat=ios) twcc_in
           !print *, 'twcc: ', twcc_in           
        case ('n')
           read(buffer, *, iostat=ios) n_in
           !print *, 'n: ', n_in           
        case ('ncc')
           read(buffer, *, iostat=ios) ncc_in
           !print *, 'ncc: ', ncc_in           
        case ('cs')
           read(buffer, *, iostat=ios) cs_in
           !print *, 'cs: ', cs_in           
        case ('s0')
           read(buffer, *, iostat=ios) s0_in
           !print *, 's0: ', s0_in             
        case ('qlat')
           read(buffer, *, iostat=ios) ql_in
           !print *, 'qlat: ', ql_in  
        case ('qup')
           read(buffer, *, iostat=ios) qup_in
           !print *, 'qup: ', qup_in             
        case ('quc')
           read(buffer, *, iostat=ios) quc_in
           !print *, 'quc: ', quc_in             
        case ('qdp')
           read(buffer, *, iostat=ios) qdp_in
           !print *, 'qdp: ', qdp_in     
        case ('depthp')
           read(buffer, *, iostat=ios) depthp_in
           !print *, 'depthp: ', depthp_in     
        end select
     end if
    end do

    ! ###############################################################################
    ! run MC model and print output
    ! ###############################################################################
    
     call muskingcungenwm(dt_in, qup_in, quc_in, qdp_in, ql_in, dx_in, bw_in, tw_in, twcc_in,&
         n_in, ncc_in, cs_in, s0_in, velp_in, depthp_in, qdc_o, velc_o, depthc_o, ck_o, cn_o, X_o)

!     call muskingcungenwm(dt_in, qup_in, quc_in, qdp_in, ql_in, dx_in, bw_in, tw_in, twcc_in,&
!         n_in, ncc_in, cs_in, s0_in, velp_in, depthp_in, qdc_o, velc_o, depthc_o)
        
        
    print*, "--------------------------------------------------------------"
    print*, " "
    print*, "MC model outputs"
    print*, " "
    print*, "flow [m^3/s] = ", qdc_o
    print*, "velocity [m/s] = ", velc_o
    print*, "depth [m] = ", depthc_o
    print*, "kinematic celerity [m/s] = ", ck_o
    print*, "Courant number [-] = ", cn_o
    print*, "X parameter [-]", X_o
    
!    open(1, file = 'output.txt', status='replace')
!    write(1,"(a,f16.8)") "qdc ",qdc_o 
!    write(1,"(a,f16.8)") "velc ",velc_o 
!    write(1,"(a,f16.8)") "depthc ",depthc_o 
!    write(1,"(a,f16.8)") "ck ", ck_o 
!    write(1,"(a,f16.8)") "cn ",cn_o 
       
end program test_MC_program