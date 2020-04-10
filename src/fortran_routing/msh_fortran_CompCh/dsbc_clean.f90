! Downstream boundary condition. Refer to p.66~67, RM1_MESH
subroutine dsbc(n,j)

    use constants
    use arrays
    use var
    use subtools

    implicit none

    integer, intent(in) :: n,j
    real :: dsc, FrCC
    real :: dpcomp, dconj, dtw, Ancomp0, Ancomp1, Axs, dxs

    Axs= areap(ncomp,j)
    call depthcalc(ncomp, j, Axs, dxs)
    dpcomp= dxs
    !* Ancomp0 and Ancomp1, area at time n and n+1 computed from known stage data at the times.
    dxs= y(n,ncomp,j)-z(ncomp,j)
    call areacalc(ncomp, j, dxs, Axs)
    Ancomp0= Axs
    dxs= y(n+1,ncomp,j)-z(ncomp,j)
    call areacalc(ncomp, j, dxs, Axs)
    Ancomp1= Axs

    !* inbank or overbank flow: when overbank flow, sub-/super-critical is determined by Fround number,p93,RM1
    Axs=areap(ncomp,j)   !* result from sub section
    dsc=qp(ncomp,j)
    call FrCCcalc(ncomp, j, dpcomp, Axs, dsc, FrCC)
    !endif

    !+++----------------------------------------------------------------------
    !+ Determine supercritical or subcritical depending on the size difference
    !+ dconj and yTW (known tailwater stage).
    !+ It is assumed that tailwater stage at the most downstream link is always
    !+ known.
    !+++----------------------------------------------------------------------
    if (FrCC<=1.0) then
    !* subcritical flow
        dac(ncomp,j)=Ancomp1-Ancomp0
        dqc(ncomp,j)=dqp(ncomp,j)
    else
    !* downstream conjugate depth to a computed value of qp
        dsc=qp(ncomp,j)
        call conjugatedep(ncomp, j, dpcomp, dsc, dconj)
        dtw=y(n,ncomp,j) - z(ncomp,j)
        if (dconj<dtw) then
        !*still subcritical
            dac(ncomp,j)=Ancomp1-Ancomp0
            dqc(ncomp,j)=dqp(ncomp,j)
        else
        !*indeed supercritical
            dac(ncomp,j)=dap(ncomp,j)
            dqc(ncomp,j)=dqp(ncomp,j)
        end if
    end if

end subroutine dsbc
